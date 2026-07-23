#!/usr/bin/env python3
"""
Validate subgenome-filtered SNP-subsampling stability, diploid-reference
sensitivity, and the genome-wide cultivated-versus-introgressed scan.

Both public VCFs were written against a combined Arabica reference and contain
sgC-labelled pseudomolecules, sgE-labelled pseudomolecules, and unassigned
scaffolds.  The workflow conservatively retains only the 11
pseudomolecules explicitly labelled for the intended VCF subgenome before any
PCA, accession-ranking, permutation, or interval analysis.  It then streams each
large VCF once, creates a reproducible filtered 12,000-SNP baseline, constructs
independent deterministic priority samples for multiple panel sizes/seeds, and
accumulates full-eligible-site PCA Gram matrices and genome-wide interval
summaries without loading complete VCFs into memory. The public workflow
writes machine-readable validation outputs and does not render manuscript
figures.
"""
from __future__ import annotations

import argparse
import collections
import gzip
import heapq
import json
import math
import random
import re
import sys
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Mapping, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

from analysis_common import (
    MASK64,
    build_metadata_lookup,
    find_file,
    interval_jaccard,
    load_metadata,
    pairwise_distances,
    parse_gt_dosage,
    parse_int_list,
    pca_from_dosage,
    pca_from_gram,
    procrustes_similarity,
    read_vcf_header,
    reconcile_sample_id,
    safe_pearson,
    spearman_correlation,
    stable_text_hash,
    standardize_variant_rows,
    upper_triangle_values,
    variant_priority_key,
    write_excel,
    write_json,
)


@dataclass(slots=True)
class VariantLite:
    contig: str
    pos: int
    ref: str
    alt: str
    call_rate: float
    maf: float
    abs_delta_af: float
    dosage: np.ndarray  # int8, -1 means missing


def classify_contig(contig: str) -> str:
    """Classify combined-reference contigs by explicit subgenome label.

    The public VCF headers contain both sgC- and sgE-labelled pseudomolecules
    plus thousands of unassigned scaffolds.  For subgenome-resolved analyses,
    only pseudomolecules explicitly labelled for the intended subgenome are
    retained.
    """
    text = str(contig)
    if re.search(r"_sg_C_", text) or re.search(r"_sgCC(?:_|$)", text):
        return "sgC"
    if re.search(r"_sg_E_", text) or re.search(r"_sgEE(?:_|$)", text):
        return "sgE"
    return "unassigned"


def expected_chr4_contig(contig: object, subgenome: str) -> bool:
    if contig is None or (isinstance(contig, float) and np.isnan(contig)):
        return False
    text = str(contig)
    expected_token = "_sg_C_" if subgenome == "sgC" else "_sg_E_"
    return text.startswith("chr_D") and expected_token in text


def chromosome_display_label(contig: str) -> str:
    match = re.match(r"chr_([A-K])_sg_[CE]_", str(contig))
    if not match:
        return str(contig)
    return f"Chr{ord(match.group(1)) - ord('A') + 1}"


class ExactReservoir:
    """Exact reservoir sampler matching the submitted Part 1 algorithm."""

    def __init__(self, size: int, seed: int) -> None:
        self.size = int(size)
        self.rng = random.Random(int(seed))
        self.records: List[VariantLite] = []
        self.seen = 0

    def decision(self) -> Optional[int]:
        self.seen += 1
        if len(self.records) < self.size:
            return len(self.records)
        replacement = self.rng.randint(0, self.seen - 1)
        return replacement if replacement < self.size else None

    def apply(self, index: Optional[int], record: VariantLite) -> None:
        if index is None:
            return
        if index == len(self.records):
            self.records.append(record)
        else:
            self.records[index] = record


class PrioritySampler:
    """Deterministic bottom-k sampler. Samples are nested across panel sizes."""

    def __init__(self, max_size: int, seed: int) -> None:
        self.max_size = int(max_size)
        self.seed = int(seed)
        self.heap: List[Tuple[int, int, VariantLite]] = []  # (-key, serial, record)
        self.serial = 0

    def qualifies(self, key: int) -> bool:
        if len(self.heap) < self.max_size:
            return True
        largest_retained_key = -self.heap[0][0]
        return key < largest_retained_key

    def add(self, key: int, record: VariantLite) -> None:
        item = (-int(key), self.serial, record)
        self.serial += 1
        if len(self.heap) < self.max_size:
            heapq.heappush(self.heap, item)
        elif key < -self.heap[0][0]:
            heapq.heapreplace(self.heap, item)

    def sorted_records(self) -> List[VariantLite]:
        ordered = sorted(((-neg_key, serial, rec) for neg_key, serial, rec in self.heap), key=lambda x: (x[0], x[1]))
        return [record for _, _, record in ordered]


class BatchGramAccumulator:
    def __init__(self, n_samples: int, subset_indices: Optional[np.ndarray] = None, batch_size: int = 4096) -> None:
        self.n_samples = int(n_samples if subset_indices is None else len(subset_indices))
        self.subset_indices = subset_indices
        self.batch_size = int(batch_size)
        self.rows: List[np.ndarray] = []
        self.gram = np.zeros((self.n_samples, self.n_samples), dtype=np.float64)
        self.n_variable_variants = 0

    def add(self, dosage: np.ndarray) -> None:
        self.rows.append(dosage)
        if len(self.rows) >= self.batch_size:
            self.flush()

    def flush(self) -> None:
        if not self.rows:
            return
        matrix = np.vstack(self.rows)
        self.rows.clear()
        if self.subset_indices is not None:
            matrix = matrix[:, self.subset_indices]
        standardized, _ = standardize_variant_rows(matrix)
        if standardized.shape[0]:
            self.gram += standardized.T @ standardized
            self.n_variable_variants += int(standardized.shape[0])


class GenomeWideAccumulator:
    def __init__(self, window_size_bp: int, top_n: int = 250) -> None:
        self.window_size_bp = int(window_size_bp)
        self.top_n = int(top_n)
        self.contigs: Dict[str, Dict[str, float]] = collections.defaultdict(
            lambda: {"n_markers": 0, "sum_abs_delta_af": 0.0, "max_abs_delta_af": 0.0, "n_ge_075": 0, "n_ge_090": 0}
        )
        self.windows: Dict[Tuple[str, int], Dict[str, float]] = collections.defaultdict(
            lambda: {"n_markers": 0, "sum_abs_delta_af": 0.0, "max_abs_delta_af": 0.0, "n_ge_075": 0, "n_ge_090": 0}
        )
        self.top_heap: List[Tuple[float, int, str, int, float]] = []
        self.serial = 0

    def add(self, contig: str, pos: int, abs_delta_af: float) -> None:
        if not np.isfinite(abs_delta_af):
            return
        stats = self.contigs[contig]
        stats["n_markers"] += 1
        stats["sum_abs_delta_af"] += float(abs_delta_af)
        stats["max_abs_delta_af"] = max(stats["max_abs_delta_af"], float(abs_delta_af))
        stats["n_ge_075"] += int(abs_delta_af >= 0.75)
        stats["n_ge_090"] += int(abs_delta_af >= 0.90)
        start = ((int(pos) - 1) // self.window_size_bp) * self.window_size_bp + 1
        wstats = self.windows[(contig, start)]
        wstats["n_markers"] += 1
        wstats["sum_abs_delta_af"] += float(abs_delta_af)
        wstats["max_abs_delta_af"] = max(wstats["max_abs_delta_af"], float(abs_delta_af))
        wstats["n_ge_075"] += int(abs_delta_af >= 0.75)
        wstats["n_ge_090"] += int(abs_delta_af >= 0.90)
        item = (float(abs_delta_af), self.serial, contig, int(pos), float(abs_delta_af))
        self.serial += 1
        if len(self.top_heap) < self.top_n:
            heapq.heappush(self.top_heap, item)
        elif abs_delta_af > self.top_heap[0][0]:
            heapq.heapreplace(self.top_heap, item)

    def contig_frame(self, subgenome: str) -> pd.DataFrame:
        rows = []
        for contig, stats in self.contigs.items():
            n = int(stats["n_markers"])
            rows.append(
                {
                    "subgenome": subgenome,
                    "contig": contig,
                    "n_markers": n,
                    "mean_abs_delta_af": stats["sum_abs_delta_af"] / n if n else np.nan,
                    "max_abs_delta_af": stats["max_abs_delta_af"],
                    "n_markers_delta_ge_0_75": int(stats["n_ge_075"]),
                    "n_markers_delta_ge_0_90": int(stats["n_ge_090"]),
                }
            )
        df = pd.DataFrame(rows)
        if not df.empty:
            df = df.sort_values(
                ["n_markers_delta_ge_0_90", "n_markers_delta_ge_0_75", "mean_abs_delta_af", "max_abs_delta_af"],
                ascending=[False, False, False, False],
            ).reset_index(drop=True)
            df["genomewide_contig_rank"] = np.arange(1, len(df) + 1)
        return df

    def window_frame(self, subgenome: str) -> pd.DataFrame:
        rows = []
        for (contig, start), stats in self.windows.items():
            n = int(stats["n_markers"])
            rows.append(
                {
                    "subgenome": subgenome,
                    "contig": contig,
                    "window_start": int(start),
                    "window_end": int(start + self.window_size_bp - 1),
                    "n_markers": n,
                    "mean_abs_delta_af": stats["sum_abs_delta_af"] / n if n else np.nan,
                    "max_abs_delta_af": stats["max_abs_delta_af"],
                    "n_markers_delta_ge_0_75": int(stats["n_ge_075"]),
                    "n_markers_delta_ge_0_90": int(stats["n_ge_090"]),
                }
            )
        df = pd.DataFrame(rows)
        if not df.empty:
            df = df.sort_values(
                ["n_markers_delta_ge_0_90", "n_markers_delta_ge_0_75", "mean_abs_delta_af"],
                ascending=[False, False, False],
            ).reset_index(drop=True)
            df["genomewide_window_rank"] = np.arange(1, len(df) + 1)
        return df

    def top_marker_records(self) -> List[Tuple[str, int, float]]:
        return [(contig, pos, delta) for _, _, contig, pos, delta in sorted(self.top_heap, reverse=True)]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Validate subgenome-filtered SNP sampling and population-structure stability.")
    parser.add_argument("--project-dir", required=True, type=Path, help="Directory containing Accession_info.xlsx and both VCFs")
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--panel-sizes", default="6000,12000,24000,48000")
    parser.add_argument("--seeds", default="20250416,20250417,20250418,20250419,20250420")
    parser.add_argument("--submitted-panel-size", type=int, default=12000)
    parser.add_argument("--submitted-seed", type=int, default=20250416)
    parser.add_argument("--min-site-call-rate", type=float, default=0.80)
    parser.add_argument("--maf-threshold", type=float, default=0.05)
    parser.add_argument("--window-size-bp", type=int, default=1_000_000)
    parser.add_argument("--top-markers", type=int, default=250)
    parser.add_argument("--n-permutations", type=int, default=500)
    parser.add_argument("--min-window-markers", type=int, default=5)
    parser.add_argument("--n-components", type=int, default=5)
    parser.add_argument("--progress-every", type=int, default=250_000)
    parser.add_argument("--batch-size", type=int, default=4096)
    parser.add_argument("--max-variant-lines", type=int, default=0, help="For environment tests only; 0 scans the complete VCF")
    parser.add_argument(
        "--max-retained-variant-lines",
        type=int,
        default=0,
        help="For filtered quick tests only; stop after this many records on matching labelled pseudomolecules (0 scans all).",
    )
    parser.add_argument(
        "--contig-filter-mode",
        choices=["matching_pseudomolecules_only", "matching_plus_unassigned", "none"],
        default="matching_pseudomolecules_only",
        help="Restrict each VCF to contigs explicitly assigned to its intended subgenome. The conservative default excludes opposite-labelled and unassigned contigs.",
    )
    return parser.parse_args()


def map_vcf_samples(sample_names: Sequence[str], metadata: pd.DataFrame, common_ids: Sequence[str]) -> Tuple[List[int], pd.DataFrame]:
    lookup = build_metadata_lookup(metadata)
    canonical_to_raw: Dict[str, Tuple[int, str, str]] = {}
    audit_rows = []
    for idx, sample in enumerate(sample_names):
        canonical, method = reconcile_sample_id(sample, lookup)
        audit_rows.append({"vcf_sample": sample, "canonical_seq_id": canonical, "reconciliation_method": method})
        if canonical is not None and canonical not in canonical_to_raw:
            canonical_to_raw[canonical] = (idx, sample, method)
    missing = [sid for sid in common_ids if sid not in canonical_to_raw]
    if missing:
        raise ValueError(f"VCF is missing common canonical IDs after reconciliation: {missing}")
    indices = [canonical_to_raw[sid][0] for sid in common_ids]
    return indices, pd.DataFrame(audit_rows)


def determine_common_panel(metadata: pd.DataFrame, sgc_samples: Sequence[str], sge_samples: Sequence[str]) -> Tuple[pd.DataFrame, List[str], pd.DataFrame]:
    lookup = build_metadata_lookup(metadata)
    rows = []
    sets = []
    for label, samples in (("sgC", sgc_samples), ("sgE", sge_samples)):
        canonical = []
        for sample in samples:
            sid, method = reconcile_sample_id(sample, lookup)
            rows.append({"subgenome": label, "vcf_sample": sample, "canonical_seq_id": sid, "reconciliation_method": method})
            if sid is not None:
                canonical.append(sid)
        sets.append(set(canonical))
    common = sets[0] & sets[1]
    ordered = [str(x) for x in metadata["seq_id"].astype(str).tolist() if str(x) in common]
    panel = metadata[metadata["seq_id"].astype(str).isin(ordered)].copy()
    order_map = {sid: i for i, sid in enumerate(ordered)}
    panel["_order"] = panel["seq_id"].astype(str).map(order_map)
    panel = panel.sort_values("_order").drop(columns="_order").reset_index(drop=True)
    return panel, ordered, pd.DataFrame(rows)


def make_record(contig: str, pos: int, ref: str, alt: str, call_rate: float, maf: float, abs_delta_af: float, dosage: np.ndarray) -> VariantLite:
    return VariantLite(contig, int(pos), ref, alt, float(call_rate), float(maf), float(abs_delta_af), dosage.copy())


def update_full_gram(accumulators: Sequence[BatchGramAccumulator], dosage: np.ndarray) -> None:
    for acc in accumulators:
        acc.add(dosage)


def interval_summary(records: Sequence[VariantLite], top_n: int, condition: str, subgenome: str) -> Tuple[Dict[str, object], pd.DataFrame]:
    if not records:
        return {
            "condition": condition,
            "subgenome": subgenome,
            "n_sampled_variants": 0,
            "top_contig": None,
            "top_contig_marker_count": 0,
            "interval_start": np.nan,
            "interval_end": np.nan,
            "interval_span_bp": np.nan,
            "mean_abs_delta_af": np.nan,
            "max_abs_delta_af": np.nan,
            "expected_chromosome4_match": False,
        }, pd.DataFrame()
    ordered = sorted(records, key=lambda r: (r.abs_delta_af, r.maf), reverse=True)[:top_n]
    marker_df = pd.DataFrame(
        [{"condition": condition, "subgenome": subgenome, "contig": r.contig, "pos": r.pos, "abs_delta_af": r.abs_delta_af, "maf": r.maf} for r in ordered]
    )
    grouped = (
        marker_df.groupby("contig")
        .agg(
            n_top_markers=("pos", "size"),
            interval_start=("pos", "min"),
            interval_end=("pos", "max"),
            mean_abs_delta_af=("abs_delta_af", "mean"),
            max_abs_delta_af=("abs_delta_af", "max"),
        )
        .reset_index()
        .sort_values(["n_top_markers", "max_abs_delta_af", "mean_abs_delta_af", "contig"], ascending=[False, False, False, True])
        .reset_index(drop=True)
    )
    row = grouped.iloc[0]
    summary = {
        "condition": condition,
        "subgenome": subgenome,
        "n_sampled_variants": len(records),
        "top_contig": row["contig"],
        "top_contig_marker_count": int(row["n_top_markers"]),
        "interval_start": int(row["interval_start"]),
        "interval_end": int(row["interval_end"]),
        "interval_span_bp": int(row["interval_end"] - row["interval_start"]),
        "mean_abs_delta_af": float(row["mean_abs_delta_af"]),
        "max_abs_delta_af": float(row["max_abs_delta_af"]),
        "expected_chromosome4_match": expected_chr4_contig(row["contig"], subgenome),
    }
    return summary, marker_df


def interval_summary_from_top_markers(markers: Sequence[Tuple[str, int, float]], condition: str, subgenome: str) -> Tuple[Dict[str, object], pd.DataFrame]:
    records = [VariantLite(c, p, "", "", np.nan, np.nan, d, np.empty(0, dtype=np.int8)) for c, p, d in markers]
    return interval_summary(records, len(records), condition, subgenome)


def pca_introgression_scores(scores: np.ndarray, groups: Sequence[str], sample_ids: Sequence[str]) -> pd.DataFrame:
    groups_array = np.asarray(groups)
    cultivated = np.where(groups_array == "arabica_cultivated")[0]
    introgressed = np.where(groups_array == "arabica_introgressed")[0]
    if cultivated.size == 0 or introgressed.size == 0:
        return pd.DataFrame()
    c_centroid = scores[cultivated].mean(axis=0)
    i_centroid = scores[introgressed].mean(axis=0)
    d_c = np.linalg.norm(scores - c_centroid, axis=1)
    d_i = np.linalg.norm(scores - i_centroid, axis=1)
    affinity = d_c - d_i
    return pd.DataFrame({"sample_id": sample_ids, "pca_introgression_affinity": affinity})


def centroid_distance(scores: np.ndarray, groups: Sequence[str]) -> float:
    groups_array = np.asarray(groups)
    c = scores[groups_array == "arabica_cultivated"]
    i = scores[groups_array == "arabica_introgressed"]
    if len(c) == 0 or len(i) == 0:
        return float("nan")
    return float(np.linalg.norm(c.mean(axis=0) - i.mean(axis=0)))


def run_permutation_test(
    records: Sequence[VariantLite],
    groups: Sequence[str],
    window_size_bp: int,
    n_permutations: int,
    min_window_markers: int,
    seed: int,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    if not records:
        return pd.DataFrame(), pd.DataFrame()
    matrix = np.vstack([r.dosage for r in records]).astype(float)
    matrix[matrix < 0] = np.nan
    groups_array = np.asarray(groups)
    c_idx = np.where(groups_array == "arabica_cultivated")[0]
    i_idx = np.where(groups_array == "arabica_introgressed")[0]
    eligible = np.concatenate([c_idx, i_idx])
    if c_idx.size < 2 or i_idx.size < 2:
        return pd.DataFrame(), pd.DataFrame()
    keys = [(r.contig, ((r.pos - 1) // window_size_bp) * window_size_bp + 1) for r in records]
    unique_keys: List[Tuple[str, int]] = []
    code_lookup: Dict[Tuple[str, int], int] = {}
    codes = np.empty(len(keys), dtype=int)
    for idx, key in enumerate(keys):
        if key not in code_lookup:
            code_lookup[key] = len(unique_keys)
            unique_keys.append(key)
        codes[idx] = code_lookup[key]
    counts = np.bincount(codes, minlength=len(unique_keys))
    valid_windows = counts >= min_window_markers

    def group_af(columns: np.ndarray) -> np.ndarray:
        values = matrix[:, columns] / 2.0
        finite = np.isfinite(values)
        counts = finite.sum(axis=1)
        sums = np.nansum(values, axis=1)
        return np.divide(sums, counts, out=np.full(matrix.shape[0], np.nan, dtype=float), where=counts > 0)

    def deltas_for_groups(c_columns: np.ndarray, i_columns: np.ndarray) -> np.ndarray:
        return np.abs(group_af(c_columns) - group_af(i_columns))

    observed_delta = deltas_for_groups(c_idx, i_idx)
    sums = np.bincount(codes, weights=np.nan_to_num(observed_delta, nan=0.0), minlength=len(unique_keys))
    valid_counts = np.bincount(codes, weights=np.isfinite(observed_delta).astype(float), minlength=len(unique_keys))
    observed_means = np.divide(sums, valid_counts, out=np.full_like(sums, np.nan, dtype=float), where=valid_counts > 0)
    observed_mask = valid_windows & np.isfinite(observed_means)
    observed_max = float(np.nanmax(observed_means[observed_mask])) if observed_mask.any() else float("nan")
    observed_argmax = int(np.nanargmax(np.where(observed_mask, observed_means, np.nan))) if observed_mask.any() else -1

    rng = np.random.default_rng(seed)
    null_maxima = np.full(n_permutations, np.nan, dtype=float)
    for perm in range(n_permutations):
        permuted = rng.permutation(eligible)
        c_perm = permuted[: len(c_idx)]
        i_perm = permuted[len(c_idx) :]
        delta = deltas_for_groups(c_perm, i_perm)
        psums = np.bincount(codes, weights=np.nan_to_num(delta, nan=0.0), minlength=len(unique_keys))
        pcounts = np.bincount(codes, weights=np.isfinite(delta).astype(float), minlength=len(unique_keys))
        means = np.divide(psums, pcounts, out=np.full_like(psums, np.nan, dtype=float), where=pcounts > 0)
        mask = valid_windows & np.isfinite(means)
        if mask.any():
            null_maxima[perm] = np.nanmax(means[mask])
    finite_null = null_maxima[np.isfinite(null_maxima)]
    p_empirical = (1 + np.sum(finite_null >= observed_max)) / (1 + finite_null.size) if finite_null.size and np.isfinite(observed_max) else np.nan
    top_contig, top_start = unique_keys[observed_argmax] if observed_argmax >= 0 else (None, None)
    summary = pd.DataFrame(
        [
            {
                "n_permutations": int(n_permutations),
                "n_valid_permutations": int(finite_null.size),
                "n_windows": int(len(unique_keys)),
                "n_windows_meeting_min_marker_count": int(valid_windows.sum()),
                "min_window_markers": int(min_window_markers),
                "observed_max_window_mean_abs_delta_af": observed_max,
                "observed_top_contig": top_contig,
                "observed_top_window_start": top_start,
                "observed_top_window_end": (top_start + window_size_bp - 1) if top_start is not None else None,
                "null_95th_percentile": float(np.quantile(finite_null, 0.95)) if finite_null.size else np.nan,
                "null_99th_percentile": float(np.quantile(finite_null, 0.99)) if finite_null.size else np.nan,
                "empirical_p_value": float(p_empirical),
            }
        ]
    )
    null_df = pd.DataFrame({"permutation": np.arange(1, n_permutations + 1), "max_window_mean_abs_delta_af": null_maxima})
    return summary, null_df


def scan_vcf(
    vcf_path: Path,
    subgenome: str,
    raw_indices: Sequence[int],
    sample_ids: Sequence[str],
    groups: Sequence[str],
    arabica_indices: np.ndarray,
    panel_sizes: Sequence[int],
    seeds: Sequence[int],
    submitted_panel_size: int,
    submitted_seed: int,
    min_site_call_rate: float,
    maf_threshold: float,
    window_size_bp: int,
    top_markers: int,
    progress_every: int,
    batch_size: int,
    max_variant_lines: int = 0,
    max_retained_variant_lines: int = 0,
    contig_filter_mode: str = "matching_pseudomolecules_only",
) -> Dict[str, object]:
    start_time = time.time()
    effective_submitted_seed = submitted_seed + (1 if subgenome == "sgE" else 0)
    submitted_reservoir = ExactReservoir(submitted_panel_size, effective_submitted_seed)
    max_size = max(panel_sizes)
    priority_samplers = {seed: PrioritySampler(max_size, seed + (1 if subgenome == "sgE" else 0)) for seed in seeds}
    full_acc = BatchGramAccumulator(len(sample_ids), subset_indices=None, batch_size=batch_size)
    arabica_acc = BatchGramAccumulator(len(sample_ids), subset_indices=arabica_indices, batch_size=batch_size)
    genomewide = GenomeWideAccumulator(window_size_bp=window_size_bp, top_n=top_markers)
    group_array = np.asarray(groups)
    cultivated_idx = np.where(group_array == "arabica_cultivated")[0]
    introgressed_idx = np.where(group_array == "arabica_introgressed")[0]
    if cultivated_idx.size < 2 or introgressed_idx.size < 2:
        raise ValueError("Cultivated or introgressed group has fewer than two common samples")

    counters = collections.Counter()
    contig_lengths: Dict[str, Optional[int]] = {}
    format_cache: Dict[str, int] = {}
    contig_hash_cache: Dict[str, int] = {}

    with gzip.open(vcf_path, "rt", encoding="utf-8", errors="replace") as handle:
        for raw_line in handle:
            if raw_line.startswith("##contig="):
                import re
                id_match = re.search(r"ID=([^,>]+)", raw_line)
                length_match = re.search(r"length=(\d+)", raw_line)
                if id_match:
                    contig_lengths[id_match.group(1)] = int(length_match.group(1)) if length_match else None
                continue
            if raw_line.startswith("#"):
                continue
            counters["variant_lines_total"] += 1
            if max_variant_lines > 0 and counters["variant_lines_total"] > max_variant_lines:
                break
            if counters["variant_lines_total"] % progress_every == 0:
                elapsed = (time.time() - start_time) / 60.0
                print(f"[{subgenome}] parsed {counters['variant_lines_total']:,} variants ({elapsed:.1f} min)", flush=True)
            fields = raw_line.rstrip("\n").split("\t")
            if len(fields) < 10:
                continue
            contig, pos_text, _, ref, alt, _, _, _, fmt = fields[:9]
            contig_class = classify_contig(contig)
            counters[f"records_on_{contig_class}_contigs"] += 1
            retain_contig = (
                contig_filter_mode == "none"
                or contig_class == subgenome
                or (contig_filter_mode == "matching_plus_unassigned" and contig_class == "unassigned")
            )
            if not retain_contig:
                if contig_class in {"sgC", "sgE"}:
                    counters["opposite_subgenome_records_excluded"] += 1
                else:
                    counters["unassigned_records_excluded"] += 1
                continue
            counters["retained_contig_records"] += 1
            if max_retained_variant_lines > 0 and counters["retained_contig_records"] > max_retained_variant_lines:
                break
            if "," in alt or len(ref) != 1 or len(alt) != 1 or alt == "*":
                continue
            counters["biallelic_snp_lines"] += 1
            if fmt not in format_cache:
                fmt_fields = fmt.split(":")
                if "GT" not in fmt_fields:
                    raise RuntimeError(f"FORMAT field without GT in {vcf_path}: {fmt}")
                format_cache[fmt] = fmt_fields.index("GT")
            gt_index = format_cache[fmt]
            sample_values = fields[9:]
            dosage = np.empty(len(raw_indices), dtype=np.int8)
            called = 0
            alt_sum = 0
            for out_idx, raw_idx in enumerate(raw_indices):
                value = parse_gt_dosage(sample_values[raw_idx], gt_index)
                dosage[out_idx] = value
                if value >= 0:
                    called += 1
                    alt_sum += int(value)
            call_rate = called / len(raw_indices)
            if called == 0:
                continue
            alt_af = alt_sum / (2.0 * called)
            maf = min(alt_af, 1.0 - alt_af)
            if call_rate < min_site_call_rate or maf < maf_threshold:
                continue
            counters["eligible_variants"] += 1
            with np.errstate(invalid="ignore"):
                c_values = dosage[cultivated_idx].astype(float)
                i_values = dosage[introgressed_idx].astype(float)
                c_values[c_values < 0] = np.nan
                i_values[i_values < 0] = np.nan
                c_af = np.nanmean(c_values / 2.0)
                i_af = np.nanmean(i_values / 2.0)
            abs_delta = abs(float(c_af) - float(i_af)) if np.isfinite(c_af) and np.isfinite(i_af) else np.nan
            update_full_gram((full_acc, arabica_acc), dosage)
            genomewide.add(contig, int(pos_text), abs_delta)

            submitted_index = submitted_reservoir.decision()
            contig_hash = contig_hash_cache.setdefault(contig, stable_text_hash(contig))
            priority_keys = {seed: variant_priority_key(contig_hash, int(pos_text), ref, alt, sampler.seed) for seed, sampler in priority_samplers.items()}
            needs_record = submitted_index is not None or any(priority_samplers[seed].qualifies(key) for seed, key in priority_keys.items())
            if not needs_record:
                continue
            record = make_record(contig, int(pos_text), ref, alt, call_rate, maf, abs_delta, dosage)
            submitted_reservoir.apply(submitted_index, record)
            for seed, key in priority_keys.items():
                sampler = priority_samplers[seed]
                if sampler.qualifies(key):
                    sampler.add(key, record)

    full_acc.flush()
    arabica_acc.flush()
    priority_records = {seed: sampler.sorted_records() for seed, sampler in priority_samplers.items()}
    return {
        "subgenome": subgenome,
        "sample_ids": list(sample_ids),
        "groups": list(groups),
        "arabica_indices": arabica_indices,
        "submitted_records": submitted_reservoir.records,
        "priority_records": priority_records,
        "full_gram": full_acc.gram,
        "full_n_variants": full_acc.n_variable_variants,
        "arabica_gram": arabica_acc.gram,
        "arabica_n_variants": arabica_acc.n_variable_variants,
        "genomewide": genomewide,
        "contig_lengths": contig_lengths,
        "counters": dict(counters),
        "contig_filter_mode": contig_filter_mode,
        "elapsed_minutes": (time.time() - start_time) / 60.0,
    }


def analyze_subgenome(
    result: Dict[str, object],
    panel_sizes: Sequence[int],
    seeds: Sequence[int],
    top_markers: int,
    n_components: int,
    n_permutations: int,
    window_size_bp: int,
    min_window_markers: int,
) -> Dict[str, pd.DataFrame]:
    subgenome = str(result["subgenome"])
    sample_ids = list(result["sample_ids"])
    groups = list(result["groups"])
    arabica_indices = np.asarray(result["arabica_indices"], dtype=int)
    arabica_sample_ids = [sample_ids[i] for i in arabica_indices]
    arabica_groups = [groups[i] for i in arabica_indices]

    condition_records: Dict[str, Sequence[VariantLite]] = {"filtered_reservoir_12000": result["submitted_records"]}
    for seed in seeds:
        records = result["priority_records"][seed]
        for size in panel_sizes:
            condition_records[f"priority_seed_{seed}_{size}"] = records[: min(size, len(records))]

    pca_rows = []
    variance_rows = []
    rank_rows = []
    interval_rows = []
    top_marker_frames = []
    score_store: Dict[Tuple[str, str], np.ndarray] = {}
    rank_store: Dict[Tuple[str, str], pd.DataFrame] = {}
    pca_score_rows = []

    for condition, records in condition_records.items():
        matrix = np.vstack([r.dosage for r in records]) if records else np.empty((0, len(sample_ids)), dtype=np.int8)
        for panel_type, indices, ids, panel_groups in (
            ("full_panel", np.arange(len(sample_ids)), sample_ids, groups),
            ("arabica_only", arabica_indices, arabica_sample_ids, arabica_groups),
        ):
            if matrix.shape[0] == 0:
                continue
            scores, ratio = pca_from_dosage(matrix[:, indices], n_components=n_components)
            score_store[(condition, panel_type)] = scores
            for sample_id, group, score_row in zip(ids, panel_groups, scores):
                record = {
                    "subgenome": subgenome,
                    "condition": condition,
                    "panel_type": panel_type,
                    "sample_id": sample_id,
                    "analysis_group": group,
                }
                for pc_idx, value in enumerate(score_row, start=1):
                    record[f"PC{pc_idx}"] = float(value)
                pca_score_rows.append(record)
            pca_rows.append(
                {
                    "subgenome": subgenome,
                    "condition": condition,
                    "panel_type": panel_type,
                    "n_variants": int(matrix.shape[0]),
                    "n_samples": int(len(ids)),
                    "cultivated_introgressed_centroid_distance": centroid_distance(scores, panel_groups),
                }
            )
            for pc_idx, value in enumerate(ratio, start=1):
                variance_rows.append(
                    {
                        "subgenome": subgenome,
                        "condition": condition,
                        "panel_type": panel_type,
                        "component": f"PC{pc_idx}",
                        "explained_variance_ratio": float(value),
                    }
                )
            ranking = pca_introgression_scores(scores, panel_groups, ids)
            if not ranking.empty:
                ranking["subgenome"] = subgenome
                ranking["condition"] = condition
                ranking["panel_type"] = panel_type
                ranking = ranking.merge(pd.DataFrame({"sample_id": ids, "analysis_group": panel_groups}), on="sample_id", how="left")
                ranking["within_introgressed_rank"] = np.nan
                mask = ranking["analysis_group"] == "arabica_introgressed"
                ranking.loc[mask, "within_introgressed_rank"] = ranking.loc[mask, "pca_introgression_affinity"].rank(method="min", ascending=False)
                rank_store[(condition, panel_type)] = ranking
                rank_rows.append(ranking)
        summary, markers = interval_summary(records, top_markers, condition, subgenome)
        interval_rows.append(summary)
        if not markers.empty:
            top_marker_frames.append(markers)

    full_scores, full_ratio = pca_from_gram(result["full_gram"], result["full_n_variants"], n_components=n_components)
    arabica_scores, arabica_ratio = pca_from_gram(result["arabica_gram"], result["arabica_n_variants"], n_components=n_components)
    score_store[("all_eligible", "full_panel")] = full_scores
    score_store[("all_eligible", "arabica_only")] = arabica_scores
    for panel_type, scores, ratio, ids, panel_groups, nvar in (
        ("full_panel", full_scores, full_ratio, sample_ids, groups, result["full_n_variants"]),
        ("arabica_only", arabica_scores, arabica_ratio, arabica_sample_ids, arabica_groups, result["arabica_n_variants"]),
    ):
        for sample_id, group, score_row in zip(ids, panel_groups, scores):
            record = {
                "subgenome": subgenome,
                "condition": "all_eligible",
                "panel_type": panel_type,
                "sample_id": sample_id,
                "analysis_group": group,
            }
            for pc_idx, value in enumerate(score_row, start=1):
                record[f"PC{pc_idx}"] = float(value)
            pca_score_rows.append(record)
        pca_rows.append(
            {
                "subgenome": subgenome,
                "condition": "all_eligible",
                "panel_type": panel_type,
                "n_variants": int(nvar),
                "n_samples": int(len(ids)),
                "cultivated_introgressed_centroid_distance": centroid_distance(scores, panel_groups),
            }
        )
        for pc_idx, value in enumerate(ratio, start=1):
            variance_rows.append({"subgenome": subgenome, "condition": "all_eligible", "panel_type": panel_type, "component": f"PC{pc_idx}", "explained_variance_ratio": float(value)})
        ranking = pca_introgression_scores(scores, panel_groups, ids)
        if not ranking.empty:
            ranking["subgenome"] = subgenome
            ranking["condition"] = "all_eligible"
            ranking["panel_type"] = panel_type
            ranking = ranking.merge(pd.DataFrame({"sample_id": ids, "analysis_group": panel_groups}), on="sample_id", how="left")
            ranking["within_introgressed_rank"] = np.nan
            mask = ranking["analysis_group"] == "arabica_introgressed"
            ranking.loc[mask, "within_introgressed_rank"] = ranking.loc[mask, "pca_introgression_affinity"].rank(method="min", ascending=False)
            rank_store[("all_eligible", panel_type)] = ranking
            rank_rows.append(ranking)
    full_interval, full_markers = interval_summary_from_top_markers(result["genomewide"].top_marker_records(), "all_eligible", subgenome)
    full_interval["n_sampled_variants"] = int(result["counters"].get("eligible_variants", 0))
    interval_rows.append(full_interval)
    if not full_markers.empty:
        top_marker_frames.append(full_markers)

    baseline_condition = "filtered_reservoir_12000"
    baseline_interval = next(row for row in interval_rows if row["condition"] == baseline_condition)
    pca_df = pd.DataFrame(pca_rows)
    stability_rows = []
    diploid_rows = []
    rank_stability_rows = []
    for (condition, panel_type), scores in score_store.items():
        baseline_scores = score_store.get((baseline_condition, panel_type))
        if baseline_scores is not None and condition != baseline_condition:
            distance_corr = safe_pearson(upper_triangle_values(pairwise_distances(baseline_scores)), upper_triangle_values(pairwise_distances(scores)))
            proc, rmse = procrustes_similarity(baseline_scores, scores)
        else:
            distance_corr = 1.0 if condition == baseline_condition else np.nan
            proc, rmse = (1.0, 0.0) if condition == baseline_condition else (np.nan, np.nan)
        stability_rows.append(
            {
                "subgenome": subgenome,
                "condition": condition,
                "panel_type": panel_type,
                "pairwise_distance_correlation_vs_filtered_12k": distance_corr,
                "procrustes_similarity_vs_filtered_12k": proc,
                "procrustes_rmse_vs_filtered_12k": rmse,
            }
        )
        if panel_type == "full_panel" and (condition, "arabica_only") in score_store:
            restricted = scores[arabica_indices]
            arabica_only = score_store[(condition, "arabica_only")]
            diploid_corr = safe_pearson(upper_triangle_values(pairwise_distances(restricted)), upper_triangle_values(pairwise_distances(arabica_only)))
            proc2, rmse2 = procrustes_similarity(restricted, arabica_only)
            full_sep = centroid_distance(restricted, arabica_groups)
            arabica_sep = centroid_distance(arabica_only, arabica_groups)
            diploid_rows.append(
                {
                    "subgenome": subgenome,
                    "condition": condition,
                    "n_arabica_samples": len(arabica_indices),
                    "arabica_pairwise_distance_correlation_full_vs_arabica_only": diploid_corr,
                    "arabica_procrustes_similarity_full_vs_arabica_only": proc2,
                    "arabica_procrustes_rmse_full_vs_arabica_only": rmse2,
                    "cultivated_introgressed_distance_full_panel_restricted": full_sep,
                    "cultivated_introgressed_distance_arabica_only": arabica_sep,
                    "arabica_only_to_full_separation_ratio": arabica_sep / full_sep if np.isfinite(full_sep) and full_sep > 0 else np.nan,
                }
            )

    # Rank stability is evaluated separately for the full and Arabica-only panels.
    for (condition, panel_type), ranking in rank_store.items():
        baseline_ranking = rank_store.get((baseline_condition, panel_type))
        if baseline_ranking is None:
            continue
        merged = baseline_ranking[["sample_id", "pca_introgression_affinity"]].merge(
            ranking[["sample_id", "pca_introgression_affinity"]], on="sample_id", suffixes=("_baseline", "_condition")
        )
        intro_base = baseline_ranking[baseline_ranking["analysis_group"] == "arabica_introgressed"].nlargest(6, "pca_introgression_affinity")["sample_id"].tolist()
        intro_cond = ranking[ranking["analysis_group"] == "arabica_introgressed"].nlargest(6, "pca_introgression_affinity")["sample_id"].tolist()
        rank_stability_rows.append(
            {
                "subgenome": subgenome,
                "condition": condition,
                "panel_type": panel_type,
                "spearman_accession_affinity_vs_filtered_12k": spearman_correlation(merged["pca_introgression_affinity_baseline"], merged["pca_introgression_affinity_condition"]),
                "top6_introgressed_overlap_count": len(set(intro_base) & set(intro_cond)),
                "baseline_top6_introgressed_ids": ";".join(map(str, intro_base)),
                "condition_top6_introgressed_ids": ";".join(map(str, intro_cond)),
            }
        )

    interval_df = pd.DataFrame(interval_rows)
    interval_df["same_top_contig_as_filtered_12k"] = interval_df["top_contig"] == baseline_interval["top_contig"]
    interval_df["interval_jaccard_vs_filtered_12k"] = [
        interval_jaccard(baseline_interval["interval_start"], baseline_interval["interval_end"], row["interval_start"], row["interval_end"])
        if row["top_contig"] == baseline_interval["top_contig"] else 0.0
        for _, row in interval_df.iterrows()
    ]
    interval_df["interval_start_shift_vs_filtered_12k_bp"] = interval_df["interval_start"] - baseline_interval["interval_start"]
    interval_df["interval_end_shift_vs_filtered_12k_bp"] = interval_df["interval_end"] - baseline_interval["interval_end"]

    perm_summary, perm_null = run_permutation_test(
        result["submitted_records"], groups, window_size_bp, n_permutations, min_window_markers, seed=91731 + (1 if subgenome == "sgE" else 0)
    )
    if not perm_summary.empty:
        perm_summary.insert(0, "subgenome", subgenome)
    if not perm_null.empty:
        perm_null.insert(0, "subgenome", subgenome)

    return {
        "pca_conditions": pca_df,
        "pca_scores": pd.DataFrame(pca_score_rows),
        "pca_variance": pd.DataFrame(variance_rows),
        "pca_stability": pd.DataFrame(stability_rows),
        "diploid_sensitivity": pd.DataFrame(diploid_rows),
        "accession_rankings": pd.concat(rank_rows, ignore_index=True) if rank_rows else pd.DataFrame(),
        "accession_rank_stability": pd.DataFrame(rank_stability_rows),
        "interval_stability": interval_df,
        "top_markers": pd.concat(top_marker_frames, ignore_index=True) if top_marker_frames else pd.DataFrame(),
        "genomewide_contigs": result["genomewide"].contig_frame(subgenome),
        "genomewide_windows": result["genomewide"].window_frame(subgenome),
        "permutation_summary": perm_summary,
        "permutation_null": perm_null,
    }


def build_status(all_tables: Mapping[str, pd.DataFrame]) -> Dict[str, object]:
    stability = all_tables["pca_stability"]
    diploid = all_tables["diploid_sensitivity"]
    intervals = all_tables["interval_stability"]
    ranks = all_tables["accession_rank_stability"]
    nonbaseline_stability = stability[(stability["condition"] != "filtered_reservoir_12000") & (stability["panel_type"] == "arabica_only")]
    priority_intervals = intervals[intervals["condition"].str.startswith("priority_seed_")]
    arabica_ranks = ranks[ranks["panel_type"] == "arabica_only"] if "panel_type" in ranks.columns else ranks
    baseline_rows = intervals[intervals["condition"] == "filtered_reservoir_12000"]
    metrics = {
        "median_pca_distance_correlation": float(nonbaseline_stability["pairwise_distance_correlation_vs_filtered_12k"].median()),
        "minimum_pca_distance_correlation": float(nonbaseline_stability["pairwise_distance_correlation_vs_filtered_12k"].min()),
        "median_full_vs_arabica_only_distance_correlation": float(diploid["arabica_pairwise_distance_correlation_full_vs_arabica_only"].median()),
        "minimum_full_vs_arabica_only_distance_correlation": float(diploid["arabica_pairwise_distance_correlation_full_vs_arabica_only"].min()),
        "top_contig_stability_fraction": float(priority_intervals["same_top_contig_as_filtered_12k"].mean()),
        "expected_chr4_recovery_fraction": float(priority_intervals["expected_chromosome4_match"].mean()),
        "baseline_expected_chr4_fraction": float(baseline_rows["expected_chromosome4_match"].mean()),
        "median_interval_jaccard": float(priority_intervals["interval_jaccard_vs_filtered_12k"].median()),
        "median_accession_rank_spearman": float(arabica_ranks["spearman_accession_affinity_vs_filtered_12k"].median()),
        "median_top6_overlap": float(arabica_ranks["top6_introgressed_overlap_count"].median()),
    }
    pass_flags = {
        "pca_sampling_stable": metrics["median_pca_distance_correlation"] >= 0.95,
        "diploid_reference_geometry_change_detected": metrics["median_full_vs_arabica_only_distance_correlation"] < 0.95,
        "chromosome4_interval_recovered": metrics["baseline_expected_chr4_fraction"] == 1.0 and metrics["expected_chr4_recovery_fraction"] >= 0.80,
        "interval_boundaries_reasonably_stable": metrics["median_interval_jaccard"] >= 0.60,
        "accession_ranking_stable": metrics["median_accession_rank_spearman"] >= 0.90 and metrics["median_top6_overlap"] >= 5,
    }
    scientific_pass_keys = ["pca_sampling_stable", "chromosome4_interval_recovered", "interval_boundaries_reasonably_stable", "accession_ranking_stable"]
    n_scientific_pass = sum(bool(pass_flags[k]) for k in scientific_pass_keys)
    overall = "PASS" if n_scientific_pass == len(scientific_pass_keys) else ("WARN" if n_scientific_pass >= 3 else "FAIL")
    return {"overall_status": overall, "metrics": metrics, "checks": pass_flags}


def main() -> None:
    args = parse_args()
    project_dir = args.project_dir.expanduser().resolve()
    output_dir = args.output_dir.expanduser().resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    panel_sizes = sorted(set(parse_int_list(args.panel_sizes)))
    seeds = sorted(set(parse_int_list(args.seeds)))
    if args.submitted_panel_size > max(panel_sizes):
        panel_sizes.append(args.submitted_panel_size)
        panel_sizes.sort()

    metadata_path = find_file(project_dir, ["Accession_info.xlsx"])
    sgc_vcf = find_file(project_dir, ["Arabica_sgC.TIP.BB.vcf.gz"])
    sge_vcf = find_file(project_dir, ["Arabica_sgE.TIP.BB.vcf.gz"])
    metadata = load_metadata(metadata_path)
    sgc_samples, _, _ = read_vcf_header(sgc_vcf)
    sge_samples, _, _ = read_vcf_header(sge_vcf)
    panel, common_ids, reconciliation_audit = determine_common_panel(metadata, sgc_samples, sge_samples)
    if len(common_ids) < 10:
        raise RuntimeError(f"Only {len(common_ids)} samples were reconciled across subgenomes; expected approximately 44")
    sgc_raw_indices, sgc_audit = map_vcf_samples(sgc_samples, metadata, common_ids)
    sge_raw_indices, sge_audit = map_vcf_samples(sge_samples, metadata, common_ids)
    sample_ids = [str(x) for x in panel["seq_id"]]
    groups = panel["analysis_group"].astype(str).tolist()
    arabica_indices = np.where(panel["species_name"].astype(str).to_numpy() == "Coffea arabica")[0]

    run_parameters = pd.DataFrame(
        {
            "parameter": [
                "project_dir", "metadata_path", "sgC_vcf", "sgE_vcf", "panel_sizes", "seeds",
                "submitted_panel_size", "submitted_seed", "min_site_call_rate", "maf_threshold",
                "window_size_bp", "top_markers", "n_permutations", "min_window_markers",
                "common_panel_size", "arabica_only_size", "max_variant_lines", "max_retained_variant_lines",
                "contig_filter_mode", "unassigned_contig_policy",
            ],
            "value": [
                str(project_dir), str(metadata_path), str(sgc_vcf), str(sge_vcf), ",".join(map(str, panel_sizes)),
                ",".join(map(str, seeds)), args.submitted_panel_size, args.submitted_seed,
                args.min_site_call_rate, args.maf_threshold, args.window_size_bp, args.top_markers,
                args.n_permutations, args.min_window_markers, len(common_ids), len(arabica_indices), args.max_variant_lines, args.max_retained_variant_lines,
                args.contig_filter_mode, "excluded from primary subgenome-resolved analysis",
            ],
        }
    )

    results = []
    scan_results_by_subgenome: Dict[str, Dict[str, object]] = {}
    for subgenome, vcf, raw_indices in (("sgC", sgc_vcf, sgc_raw_indices), ("sgE", sge_vcf, sge_raw_indices)):
        print(f"\n=== Streaming {subgenome}: {vcf} ===", flush=True)
        scan_result = scan_vcf(
            vcf_path=vcf,
            subgenome=subgenome,
            raw_indices=raw_indices,
            sample_ids=sample_ids,
            groups=groups,
            arabica_indices=arabica_indices,
            panel_sizes=panel_sizes,
            seeds=seeds,
            submitted_panel_size=args.submitted_panel_size,
            submitted_seed=args.submitted_seed,
            min_site_call_rate=args.min_site_call_rate,
            maf_threshold=args.maf_threshold,
            window_size_bp=args.window_size_bp,
            top_markers=args.top_markers,
            progress_every=args.progress_every,
            batch_size=args.batch_size,
            max_variant_lines=args.max_variant_lines,
            max_retained_variant_lines=args.max_retained_variant_lines,
            contig_filter_mode=args.contig_filter_mode,
        )
        scan_results_by_subgenome[subgenome] = scan_result
        results.append(analyze_subgenome(scan_result, panel_sizes, seeds, args.top_markers, args.n_components, args.n_permutations, args.window_size_bp, args.min_window_markers))
        write_json(output_dir / f"{subgenome}_scan_counters.json", {"counters": scan_result["counters"], "contig_filter_mode": scan_result["contig_filter_mode"], "elapsed_minutes": scan_result["elapsed_minutes"], "full_pca_variable_variants": scan_result["full_n_variants"]})

    keys = results[0].keys()
    combined = {key: pd.concat([r[key] for r in results if not r[key].empty], ignore_index=True) if any(not r[key].empty for r in results) else pd.DataFrame() for key in keys}
    status = build_status(combined)
    status_rows = []
    for name, value in status["metrics"].items():
        status_rows.append({"category": "metric", "item": name, "value": value})
    for name, value in status["checks"].items():
        status_rows.append({"category": "check", "item": name, "value": value})
    status_rows.insert(0, {"category": "overall", "item": "overall_status", "value": status["overall_status"]})

    panel_out = panel.copy()
    panel_out.insert(0, "panel_order", np.arange(1, len(panel_out) + 1))
    audit = pd.concat([reconciliation_audit, sgc_audit.assign(subgenome="sgC"), sge_audit.assign(subgenome="sgE")], ignore_index=True, sort=False)
    contig_audit_rows = []
    # Use the raw scan results retained above through scan_results_by_subgenome.
    for subgenome, scan_result in scan_results_by_subgenome.items():
        counters = scan_result["counters"]
        contig_audit_rows.append({
            "subgenome": subgenome,
            "filter_mode": scan_result["contig_filter_mode"],
            "raw_variant_records": counters.get("variant_lines_total", 0),
            "records_on_sgC_contigs": counters.get("records_on_sgC_contigs", 0),
            "records_on_sgE_contigs": counters.get("records_on_sgE_contigs", 0),
            "records_on_unassigned_contigs": counters.get("records_on_unassigned_contigs", 0),
            "retained_contig_records": counters.get("retained_contig_records", 0),
            "opposite_subgenome_records_excluded": counters.get("opposite_subgenome_records_excluded", 0),
            "unassigned_records_excluded": counters.get("unassigned_records_excluded", 0),
            "eligible_variants_after_filtering": counters.get("eligible_variants", 0),
        })
    contig_filter_audit = pd.DataFrame(contig_audit_rows)
    workbook = output_dir / "sampling_population_structure_validation.xlsx"
    write_excel(
        workbook,
        {
            "Run_parameters": run_parameters,
            "Validation_status": pd.DataFrame(status_rows),
            "Common_panel": panel_out,
            "Sample_reconciliation": audit,
            "Contig_filter_audit": contig_filter_audit,
            "PCA_conditions": combined["pca_conditions"],
            "PCA_scores": combined["pca_scores"],
            "PCA_explained_variance": combined["pca_variance"],
            "PCA_sampling_stability": combined["pca_stability"],
            "Diploid_reference_sensitivity": combined["diploid_sensitivity"],
            "Accession_rank_stability": combined["accession_rank_stability"],
            "Accession_rankings": combined["accession_rankings"],
            "Interval_stability": combined["interval_stability"],
            "Top_marker_support": combined["top_markers"],
            "Genomewide_contig_scan": combined["genomewide_contigs"],
            "Genomewide_window_scan": combined["genomewide_windows"],
            "Permutation_summary": combined["permutation_summary"],
            "Permutation_null": combined["permutation_null"],
        },
    )
    write_json(output_dir / "sampling_validation_status.json", status)

    print("\nSampling/population-structure validation completed.")
    print(f"Workbook: {workbook}")
    print(f"Overall status: {status['overall_status']}")
    print(json.dumps(status, indent=2))


if __name__ == "__main__":
    main()
