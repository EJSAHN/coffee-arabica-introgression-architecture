#!/usr/bin/env python3
"""
Integrate matched subgenome sample sets and summarize population structure.

The script reconciles sample identifiers, rebuilds the common accession panel,
computes subgenome-resolved genetic similarity summaries, derives PCA centroid
summaries, estimates group affinities, and ranks accessions for introgression
and subgenome-asymmetry signals.
"""

from __future__ import annotations

import argparse
import math
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd


REQUIRED_ANALYSIS_SHEETS = [
    "13_sgC_sample_match",
    "14_sgE_sample_match",
    "16_sgC_sample_qc",
    "17_sgE_sample_qc",
    "20_subgenome_delta",
    "23_sgC_pca_scores",
    "24_sgE_pca_scores",
    "27_sgC_ibs",
    "28_sgE_ibs",
    "31_sgC_top_contrasts",
    "32_sgE_top_contrasts",
]

REQUIRED_SUPPLEMENTARY_SHEETS = [
    "S1_accession_master",
    "S1_column_guide",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Create Part 2 coffee population analysis workbooks from Part 1 outputs."
    )
    parser.add_argument(
        "--analysis-workbook",
        required=True,
        help="Path to the Part 1 manuscript workbook."
    )
    parser.add_argument(
        "--supplementary-workbook",
        required=True,
        help="Path to the Part 1 supplementary accession workbook."
    )
    parser.add_argument(
        "--output-workbook",
        required=True,
        help="Path to the Part 2 manuscript workbook to write."
    )
    parser.add_argument(
        "--supplementary-output-workbook",
        required=True,
        help="Path to the Part 2 supplementary workbook to write."
    )
    parser.add_argument(
        "--min-group-size",
        type=int,
        default=2,
        help="Minimum number of samples required to retain a group in the rescued common analysis set."
    )
    parser.add_argument(
        "--top-markers-per-contrast",
        type=int,
        default=50,
        help="Number of top markers per contrast to retain in manuscript-focused summaries."
    )
    return parser.parse_args()


def ensure_workbook_has_sheets(path: Path, required_sheets: Sequence[str]) -> None:
    xl = pd.ExcelFile(path)
    missing = [sheet for sheet in required_sheets if sheet not in xl.sheet_names]
    if missing:
        raise ValueError(
            f"Workbook '{path}' is missing required sheets: {', '.join(missing)}"
        )


def normalize_token(value: object) -> Optional[str]:
    if value is None or (isinstance(value, float) and np.isnan(value)):
        return None
    token = str(value).strip().lower()
    if token == "":
        return None
    token = re.sub(r"[^a-z0-9]+", "", token)
    return token or None


def build_metadata_lookup(accession_master: pd.DataFrame) -> Dict[str, str]:
    lookup: Dict[str, str] = {}
    for _, row in accession_master.iterrows():
        canonical_id = str(row["seq_id"])
        for source_value in (row.get("seq_id"), row.get("accession_name")):
            token = normalize_token(source_value)
            if token:
                lookup[token] = canonical_id
    return lookup


def candidate_sample_names(sample_id: str) -> List[str]:
    value = str(sample_id).strip()
    candidates = [value]

    # Common suffix patterns observed in Part 1 outputs
    candidates.append(re.sub(r"_(sgc|sge)$", "", value, flags=re.IGNORECASE))
    candidates.append(re.sub(r"_eugenioides$", "", value, flags=re.IGNORECASE))
    candidates.append(re.sub(r"_canephora$", "", value, flags=re.IGNORECASE))

    # Remove trailing subgenome labels even when additional tokens exist
    candidates.append(value.replace("_sgC", "").replace("_sgE", ""))
    candidates.append(value.replace("_Eugenioides", "").replace("_Canephora", ""))

    # Remove repeated separators
    cleaned = []
    seen = set()
    for candidate in candidates:
        candidate = candidate.strip("_ ")
        if candidate and candidate not in seen:
            cleaned.append(candidate)
            seen.add(candidate)
    return cleaned


@dataclass
class RescueResult:
    canonical_seq_id: Optional[str]
    rescue_status: str
    rescue_method: Optional[str]


def rescue_sample_id(sample_id: object, metadata_lookup: Dict[str, str]) -> RescueResult:
    if sample_id is None or (isinstance(sample_id, float) and np.isnan(sample_id)):
        return RescueResult(None, "unresolved", None)

    raw = str(sample_id)

    # Exact normalized match
    token = normalize_token(raw)
    if token in metadata_lookup:
        return RescueResult(metadata_lookup[token], "matched", "exact")

    for candidate in candidate_sample_names(raw):
        token = normalize_token(candidate)
        if token in metadata_lookup:
            return RescueResult(metadata_lookup[token], "rescued", f"normalized:{candidate}")

    return RescueResult(None, "unresolved", None)


def apply_rescue(
    df: pd.DataFrame,
    sample_col: str,
    metadata_lookup: Dict[str, str]
) -> pd.DataFrame:
    out = df.copy()
    rescued = out[sample_col].apply(lambda x: rescue_sample_id(x, metadata_lookup))
    out["canonical_seq_id"] = [r.canonical_seq_id for r in rescued]
    out["rescue_status"] = [r.rescue_status for r in rescued]
    out["rescue_method"] = [r.rescue_method for r in rescued]
    return out


def build_identity_rescue_audit(
    sgc_match: pd.DataFrame,
    sge_match: pd.DataFrame,
    accession_master: pd.DataFrame
) -> pd.DataFrame:
    metadata_cols = [
        "seq_id",
        "accession_name",
        "species_name",
        "variety",
        "analysis_group",
        "country_of_origin",
        "genome_structure",
    ]
    meta = accession_master[metadata_cols].copy()

    sgc = sgc_match[["sample_id", "canonical_seq_id", "rescue_status", "rescue_method"]].copy()
    sgc["subgenome_source"] = "sgC"

    sge = sge_match[["sample_id", "canonical_seq_id", "rescue_status", "rescue_method"]].copy()
    sge["subgenome_source"] = "sgE"

    audit = pd.concat([sgc, sge], ignore_index=True)
    audit = audit.merge(meta, left_on="canonical_seq_id", right_on="seq_id", how="left")
    audit["matched_metadata"] = audit["seq_id"].notna()
    audit = audit.sort_values(
        by=["matched_metadata", "subgenome_source", "sample_id"],
        ascending=[False, True, True]
    ).reset_index(drop=True)
    return audit


def load_ibs_matrix(workbook_path: Path, sheet_name: str) -> pd.DataFrame:
    raw = pd.read_excel(workbook_path, sheet_name=sheet_name)
    first_col = raw.columns[0]
    raw = raw.rename(columns={first_col: "sample_id"})
    raw["sample_id"] = raw["sample_id"].astype(str)
    column_map = {c: str(c) for c in raw.columns if c != "sample_id"}
    raw = raw.rename(columns=column_map)
    return raw


def collapse_duplicate_matrix_entries(df: pd.DataFrame) -> pd.DataFrame:
    """
    Collapse duplicate row and column labels by arithmetic mean.
    """
    numeric = df.copy()
    numeric = numeric.set_index("sample_id")
    numeric.index = numeric.index.astype(str)
    numeric.columns = numeric.columns.astype(str)
    numeric = numeric.apply(pd.to_numeric, errors="coerce")

    if numeric.index.duplicated().any():
        numeric = numeric.groupby(level=0, sort=False).mean()

    if pd.Index(numeric.columns).duplicated().any():
        numeric = numeric.T.groupby(level=0, sort=False).mean().T

    numeric = numeric.loc[~numeric.index.duplicated(keep="first")]
    numeric = numeric.loc[:, ~pd.Index(numeric.columns).duplicated(keep="first")]

    keep = [idx for idx in numeric.index if idx in numeric.columns]
    numeric = numeric.loc[keep, keep]
    numeric.index.name = "sample_id"
    return numeric.reset_index()


def canonicalize_ibs_matrix(
    df: pd.DataFrame,
    metadata_lookup: Dict[str, str],
    allowed_ids: Optional[Sequence[str]] = None
) -> pd.DataFrame:
    row_map = {str(v): rescue_sample_id(v, metadata_lookup).canonical_seq_id for v in df["sample_id"].tolist()}
    col_map = {str(c): rescue_sample_id(c, metadata_lookup).canonical_seq_id for c in df.columns if c != "sample_id"}

    renamed = df.copy()
    renamed["sample_id"] = renamed["sample_id"].map(row_map)
    renamed = renamed.rename(columns=col_map)

    renamed = renamed[renamed["sample_id"].notna()].copy()
    keep_columns = ["sample_id"] + [c for c in renamed.columns if c != "sample_id" and c is not None]
    renamed = renamed[keep_columns]

    renamed = collapse_duplicate_matrix_entries(renamed)

    if allowed_ids is not None:
        allowed = [str(x) for x in allowed_ids]
        numeric = renamed.set_index("sample_id")
        available = [x for x in allowed if x in numeric.index and x in numeric.columns]
        numeric = numeric.loc[available, available]
        numeric.index.name = "sample_id"
        return numeric.reset_index()

    return renamed


def matrix_to_long(
    matrix_df: pd.DataFrame,
    group_lookup: pd.Series,
    remove_self: bool = True
) -> pd.DataFrame:
    numeric = matrix_df.set_index("sample_id")
    numeric.index = numeric.index.astype(str)
    numeric.columns = numeric.columns.astype(str)

    rows = []
    index_ids = list(numeric.index)
    for i, sample_a in enumerate(index_ids):
        for j, sample_b in enumerate(index_ids):
            if remove_self and sample_a == sample_b:
                continue
            if sample_a > sample_b:
                continue
            value = numeric.loc[sample_a, sample_b]
            rows.append(
                {
                    "sample_a": sample_a,
                    "sample_b": sample_b,
                    "ibs": value,
                    "group_a": group_lookup.get(sample_a),
                    "group_b": group_lookup.get(sample_b),
                }
            )
    return pd.DataFrame(rows)


def summarize_group_ibs(
    matrix_df: pd.DataFrame,
    group_lookup: pd.Series
) -> pd.DataFrame:
    long_df = matrix_to_long(matrix_df, group_lookup=group_lookup, remove_self=True)
    records = []

    groups = sorted(pd.Series(group_lookup.dropna().unique()).astype(str).tolist())
    for group_a in groups:
        for group_b in groups:
            subset = long_df[
                (
                    ((long_df["group_a"] == group_a) & (long_df["group_b"] == group_b)) |
                    ((long_df["group_a"] == group_b) & (long_df["group_b"] == group_a))
                )
            ].copy()

            if group_a == group_b:
                subset = subset[subset["sample_a"] != subset["sample_b"]]

            if subset.empty:
                continue

            records.append(
                {
                    "group_a": group_a,
                    "group_b": group_b,
                    "n_pairs": int(subset.shape[0]),
                    "mean_ibs": subset["ibs"].mean(),
                    "median_ibs": subset["ibs"].median(),
                    "min_ibs": subset["ibs"].min(),
                    "max_ibs": subset["ibs"].max(),
                    "mean_distance_1_minus_ibs": 1.0 - subset["ibs"].mean(),
                }
            )

    out = pd.DataFrame(records)
    out = out.sort_values(
        by=["group_a", "group_b", "mean_ibs"],
        ascending=[True, True, False]
    ).reset_index(drop=True)
    return out


def compute_group_affinity(
    matrix_df: pd.DataFrame,
    sample_metadata: pd.DataFrame,
    subgenome_label: str
) -> pd.DataFrame:
    numeric = matrix_df.set_index("sample_id")
    numeric.index = numeric.index.astype(str)
    numeric.columns = numeric.columns.astype(str)

    sample_meta = sample_metadata.copy()
    sample_meta["canonical_seq_id"] = sample_meta["canonical_seq_id"].astype(str)

    groups = sorted(sample_meta["analysis_group"].dropna().astype(str).unique().tolist())
    records = []

    for _, row in sample_meta.iterrows():
        sample_id = str(row["canonical_seq_id"])
        if sample_id not in numeric.index:
            continue

        result = {
            "canonical_seq_id": sample_id,
            "accession_name": row.get("accession_name"),
            "species_name": row.get("species_name"),
            "variety": row.get("variety"),
            "analysis_group": row.get("analysis_group"),
            "country_of_origin": row.get("country_of_origin"),
            "genome_structure": row.get("genome_structure"),
        }

        for group_name in groups:
            member_ids = (
                sample_meta.loc[sample_meta["analysis_group"] == group_name, "canonical_seq_id"]
                .astype(str)
                .tolist()
            )
            if group_name == row.get("analysis_group"):
                member_ids = [x for x in member_ids if x != sample_id]

            member_ids = [x for x in member_ids if x in numeric.columns]
            column_name = f"{subgenome_label}_mean_ibs_to_{group_name}"

            if len(member_ids) == 0:
                result[column_name] = np.nan
            else:
                result[column_name] = numeric.loc[sample_id, member_ids].astype(float).mean()

        records.append(result)

    return pd.DataFrame(records)


def compute_pca_centroids(
    pca_df: pd.DataFrame,
    sample_metadata: pd.DataFrame,
    subgenome_label: str,
    n_components: int = 5
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    merged = pca_df.merge(
        sample_metadata,
        on="canonical_seq_id",
        how="inner",
        suffixes=("", "_meta")
    ).copy()

    component_cols = [c for c in merged.columns if re.fullmatch(r"PC\d+", str(c))]
    component_cols = component_cols[:n_components]

    centroid_records = []
    for group_name, sub in merged.groupby("analysis_group", dropna=True):
        centroid = {
            "analysis_group": group_name,
            "n_samples": int(sub.shape[0]),
        }
        for col in component_cols:
            centroid[f"{col}_centroid"] = sub[col].mean()
        centroid_records.append(centroid)

    centroids = pd.DataFrame(centroid_records).sort_values("analysis_group").reset_index(drop=True)

    distance_records = []
    if not centroids.empty:
        for i, row_a in centroids.iterrows():
            for j, row_b in centroids.iterrows():
                if j < i:
                    continue
                vec_a = np.array([row_a[f"{col}_centroid"] for col in component_cols], dtype=float)
                vec_b = np.array([row_b[f"{col}_centroid"] for col in component_cols], dtype=float)
                dist = float(np.linalg.norm(vec_a - vec_b))
                distance_records.append(
                    {
                        "group_a": row_a["analysis_group"],
                        "group_b": row_b["analysis_group"],
                        f"{subgenome_label}_centroid_distance_pc1_to_pc{len(component_cols)}": dist,
                    }
                )

    distances = pd.DataFrame(distance_records)
    return centroids, distances


def prepare_sample_qc_rescued(
    sgc_qc: pd.DataFrame,
    sge_qc: pd.DataFrame,
    accession_master: pd.DataFrame,
    metadata_lookup: Dict[str, str],
    min_group_size: int
) -> pd.DataFrame:
    sgc = apply_rescue(sgc_qc, "sample_id", metadata_lookup)
    sge = apply_rescue(sge_qc, "sample_id", metadata_lookup)

    sgc = sgc.rename(columns={c: f"sgC_{c}" for c in sgc.columns if c != "canonical_seq_id"})
    sge = sge.rename(columns={c: f"sgE_{c}" for c in sge.columns if c != "canonical_seq_id"})

    merged = sgc.merge(
        sge,
        on="canonical_seq_id",
        how="inner",
        suffixes=("", "")
    )

    metadata_cols = [
        "seq_id", "accession_name", "species_name", "variety", "location_code",
        "country_of_origin", "Latitude", "Longitude", "altitude", "ploidy_level",
        "genome_size_gb", "genome_structure", "donor_institute", "collection_location",
        "notes", "ref", "full_reference", "analysis_group", "has_latitude",
        "has_longitude", "has_altitude"
    ]
    metadata_subset = accession_master.loc[:, [c for c in metadata_cols if c in accession_master.columns]].copy()

    merged = merged.merge(
        metadata_subset,
        left_on="canonical_seq_id",
        right_on="seq_id",
        how="left"
    )

    group_counts = (
        merged["analysis_group"]
        .dropna()
        .astype(str)
        .value_counts()
        .to_dict()
    )
    merged["retained_for_common_analysis"] = merged["analysis_group"].astype(str).map(group_counts).fillna(0) >= min_group_size
    merged["present_in_rescued_common_set"] = True

    merged["delta_heterozygosity_sgC_minus_sgE"] = (
        merged["sgC_heterozygosity_rate_biallelic_snps"] - merged["sgE_heterozygosity_rate_biallelic_snps"]
    )
    merged["delta_call_rate_sgC_minus_sgE"] = (
        merged["sgC_call_rate_biallelic_snps"] - merged["sgE_call_rate_biallelic_snps"]
    )
    merged["delta_non_reference_rate_sgC_minus_sgE"] = (
        merged["sgC_non_reference_rate_biallelic_snps"] - merged["sgE_non_reference_rate_biallelic_snps"]
    )

    return merged


def summarize_group_counts(analysis_set: pd.DataFrame) -> pd.DataFrame:
    out = (
        analysis_set.groupby(["species_name", "variety", "analysis_group"], dropna=False)
        .agg(
            n_samples=("canonical_seq_id", "count"),
            n_countries=("country_of_origin", pd.Series.nunique),
        )
        .reset_index()
        .sort_values(["species_name", "analysis_group"])
        .reset_index(drop=True)
    )
    return out


def build_group_lookup(analysis_set: pd.DataFrame) -> pd.Series:
    lookup = analysis_set.set_index("canonical_seq_id")["analysis_group"].astype(str)
    return lookup


def build_pca_rescued_table(
    pca_df: pd.DataFrame,
    accession_master: pd.DataFrame,
    metadata_lookup: Dict[str, str],
    allowed_ids: Sequence[str]
) -> pd.DataFrame:
    rescued = apply_rescue(pca_df, "sample_id", metadata_lookup)
    metadata_cols = [
        "seq_id", "accession_name", "species_name", "variety", "analysis_group",
        "country_of_origin", "genome_structure"
    ]
    metadata_subset = accession_master.loc[:, [c for c in metadata_cols if c in accession_master.columns]].copy()
    out = rescued.merge(metadata_subset, left_on="canonical_seq_id", right_on="seq_id", how="left")
    out = out[out["canonical_seq_id"].astype(str).isin([str(x) for x in allowed_ids])].copy()
    return out


def summarize_hotspots(
    top_contrasts_df: pd.DataFrame,
    subgenome_label: str
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    by_contrast = (
        top_contrasts_df.groupby("contrast", dropna=False)
        .agg(
            n_markers=("contig", "count"),
            n_contigs=("contig", pd.Series.nunique),
            mean_abs_delta_af=("abs_delta_af", "mean"),
            median_abs_delta_af=("abs_delta_af", "median"),
            max_abs_delta_af=("abs_delta_af", "max"),
            n_markers_delta_ge_0_90=("abs_delta_af", lambda x: int((x >= 0.90).sum())),
            n_markers_delta_ge_0_75=("abs_delta_af", lambda x: int((x >= 0.75).sum())),
        )
        .reset_index()
        .sort_values(["max_abs_delta_af", "n_markers"], ascending=[False, False])
        .reset_index(drop=True)
    )
    by_contrast["subgenome"] = subgenome_label

    by_contrast_contig = (
        top_contrasts_df.groupby(["contrast", "contig"], dropna=False)
        .agg(
            n_markers=("pos", "count"),
            min_pos=("pos", "min"),
            max_pos=("pos", "max"),
            span_bp=("pos", lambda x: int(x.max() - x.min()) if len(x) > 1 else 0),
            mean_abs_delta_af=("abs_delta_af", "mean"),
            max_abs_delta_af=("abs_delta_af", "max"),
        )
        .reset_index()
        .sort_values(
            ["contrast", "n_markers", "max_abs_delta_af", "mean_abs_delta_af"],
            ascending=[True, False, False, False]
        )
        .reset_index(drop=True)
    )
    by_contrast_contig["subgenome"] = subgenome_label
    return by_contrast, by_contrast_contig


def merge_distance_summaries(
    sgc_df: pd.DataFrame,
    sge_df: pd.DataFrame,
    value_column_sgc: str,
    value_column_sge: str,
    merged_value_name: str
) -> pd.DataFrame:
    merged = sgc_df.merge(
        sge_df,
        on=["group_a", "group_b"],
        how="outer",
        suffixes=("_sgC", "_sgE")
    )
    left_col = f"{value_column_sgc}_sgC" if f"{value_column_sgc}_sgC" in merged.columns else value_column_sgc
    right_col = f"{value_column_sge}_sgE" if f"{value_column_sge}_sgE" in merged.columns else value_column_sge
    merged[f"delta_{merged_value_name}_sgC_minus_sgE"] = merged[left_col] - merged[right_col]
    return merged


def compute_priority_tables(
    analysis_set: pd.DataFrame,
    sgc_affinity: pd.DataFrame,
    sge_affinity: pd.DataFrame
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    merged = sgc_affinity.merge(
        sge_affinity,
        on=[
            "canonical_seq_id", "accession_name", "species_name", "variety",
            "analysis_group", "country_of_origin", "genome_structure"
        ],
        how="inner",
        suffixes=("", "")
    ).merge(
        analysis_set[
            [
                "canonical_seq_id",
                "delta_heterozygosity_sgC_minus_sgE",
                "delta_call_rate_sgC_minus_sgE",
                "delta_non_reference_rate_sgC_minus_sgE",
                "sgC_mean_depth",
                "sgE_mean_depth",
            ]
        ],
        on="canonical_seq_id",
        how="left"
    )

    arabica_only = merged[merged["species_name"] == "Coffea arabica"].copy()

    def safe_col(df: pd.DataFrame, name: str) -> pd.Series:
        return df[name] if name in df.columns else pd.Series(np.nan, index=df.index)

    arabica_only["sgC_introgressed_shift"] = (
        safe_col(arabica_only, "sgC_mean_ibs_to_arabica_introgressed") -
        safe_col(arabica_only, "sgC_mean_ibs_to_arabica_cultivated")
    )
    arabica_only["sgE_introgressed_shift"] = (
        safe_col(arabica_only, "sgE_mean_ibs_to_arabica_introgressed") -
        safe_col(arabica_only, "sgE_mean_ibs_to_arabica_cultivated")
    )
    arabica_only["mean_introgressed_shift"] = arabica_only[
        ["sgC_introgressed_shift", "sgE_introgressed_shift"]
    ].mean(axis=1)

    arabica_only["sgC_wild_shift"] = (
        safe_col(arabica_only, "sgC_mean_ibs_to_arabica_wild") -
        safe_col(arabica_only, "sgC_mean_ibs_to_arabica_cultivated")
    )
    arabica_only["sgE_wild_shift"] = (
        safe_col(arabica_only, "sgE_mean_ibs_to_arabica_wild") -
        safe_col(arabica_only, "sgE_mean_ibs_to_arabica_cultivated")
    )

    arabica_only["mean_canephora_affinity"] = arabica_only[
        [c for c in ["sgC_mean_ibs_to_canephora", "sgE_mean_ibs_to_canephora"] if c in arabica_only.columns]
    ].mean(axis=1)
    arabica_only["mean_eugenioides_affinity"] = arabica_only[
        [c for c in ["sgC_mean_ibs_to_eugenioides", "sgE_mean_ibs_to_eugenioides"] if c in arabica_only.columns]
    ].mean(axis=1)

    arabica_only["abs_delta_heterozygosity"] = arabica_only["delta_heterozygosity_sgC_minus_sgE"].abs()
    arabica_only["abs_delta_non_reference_rate"] = arabica_only["delta_non_reference_rate_sgC_minus_sgE"].abs()
    arabica_only["abs_delta_introgressed_shift"] = (
        arabica_only["sgC_introgressed_shift"] - arabica_only["sgE_introgressed_shift"]
    ).abs()

    rank_specs = {
        "rank_mean_introgressed_shift_desc": ("mean_introgressed_shift", False),
        "rank_mean_canephora_affinity_desc": ("mean_canephora_affinity", False),
        "rank_abs_delta_heterozygosity_desc": ("abs_delta_heterozygosity", False),
        "rank_abs_delta_non_reference_rate_desc": ("abs_delta_non_reference_rate", False),
        "rank_abs_delta_introgressed_shift_desc": ("abs_delta_introgressed_shift", False),
    }

    for rank_col, (value_col, ascending) in rank_specs.items():
        arabica_only[rank_col] = arabica_only[value_col].rank(
            method="average",
            ascending=ascending,
            na_option="bottom"
        )

    rank_cols = list(rank_specs.keys())
    arabica_only["mean_rank_unweighted"] = arabica_only[rank_cols].mean(axis=1)
    arabica_only["median_rank_unweighted"] = arabica_only[rank_cols].median(axis=1)
    arabica_only["n_rank_metrics_available"] = arabica_only[rank_cols].notna().sum(axis=1)
    arabica_only = arabica_only.sort_values(
        ["mean_rank_unweighted", "median_rank_unweighted", "analysis_group", "accession_name"],
        ascending=[True, True, True, True]
    ).reset_index(drop=True)

    asymmetry = arabica_only[
        [
            "canonical_seq_id", "accession_name", "analysis_group", "country_of_origin", "genome_structure",
            "delta_heterozygosity_sgC_minus_sgE", "delta_non_reference_rate_sgC_minus_sgE",
            "delta_call_rate_sgC_minus_sgE", "sgC_introgressed_shift", "sgE_introgressed_shift",
            "abs_delta_introgressed_shift", "mean_canephora_affinity", "mean_eugenioides_affinity"
        ]
    ].copy()

    asymmetry["rank_abs_delta_introgressed_shift_desc"] = asymmetry["abs_delta_introgressed_shift"].rank(
        method="average", ascending=False, na_option="bottom"
    )
    asymmetry["rank_abs_delta_heterozygosity_desc"] = asymmetry["delta_heterozygosity_sgC_minus_sgE"].abs().rank(
        method="average", ascending=False, na_option="bottom"
    )
    asymmetry["rank_abs_delta_non_reference_rate_desc"] = asymmetry["delta_non_reference_rate_sgC_minus_sgE"].abs().rank(
        method="average", ascending=False, na_option="bottom"
    )
    asymmetry["mean_rank_unweighted"] = asymmetry[
        [
            "rank_abs_delta_introgressed_shift_desc",
            "rank_abs_delta_heterozygosity_desc",
            "rank_abs_delta_non_reference_rate_desc",
        ]
    ].mean(axis=1)
    asymmetry["median_rank_unweighted"] = asymmetry[
        [
            "rank_abs_delta_introgressed_shift_desc",
            "rank_abs_delta_heterozygosity_desc",
            "rank_abs_delta_non_reference_rate_desc",
        ]
    ].median(axis=1)
    asymmetry = asymmetry.sort_values(
        ["mean_rank_unweighted", "median_rank_unweighted", "analysis_group", "accession_name"],
        ascending=[True, True, True, True]
    ).reset_index(drop=True)

    return arabica_only, asymmetry


def build_run_summary_sheet(args: argparse.Namespace) -> pd.DataFrame:
    return pd.DataFrame(
        {
            "parameter": [
                "analysis_workbook",
                "supplementary_workbook",
                "output_workbook",
                "supplementary_output_workbook",
                "min_group_size",
                "top_markers_per_contrast",
            ],
            "value": [
                str(args.analysis_workbook),
                str(args.supplementary_workbook),
                str(args.output_workbook),
                str(args.supplementary_output_workbook),
                args.min_group_size,
                args.top_markers_per_contrast,
            ]
        }
    )


def build_input_inventory(analysis_workbook: Path, supplementary_workbook: Path) -> pd.DataFrame:
    return pd.DataFrame(
        {
            "input_name": ["analysis_workbook", "supplementary_workbook"],
            "path": [str(analysis_workbook), str(supplementary_workbook)],
            "file_size_bytes": [analysis_workbook.stat().st_size, supplementary_workbook.stat().st_size],
        }
    )


def build_analysis_scope(analysis_set: pd.DataFrame) -> pd.DataFrame:
    records = [
        {"metric": "n_common_samples_rescued", "value": int(analysis_set.shape[0])},
        {"metric": "n_species", "value": int(analysis_set["species_name"].nunique())},
        {"metric": "n_groups", "value": int(analysis_set["analysis_group"].nunique())},
        {"metric": "n_arabica_wild", "value": int((analysis_set["analysis_group"] == "arabica_wild").sum())},
        {"metric": "n_arabica_cultivated", "value": int((analysis_set["analysis_group"] == "arabica_cultivated").sum())},
        {"metric": "n_arabica_introgressed", "value": int((analysis_set["analysis_group"] == "arabica_introgressed").sum())},
        {"metric": "n_canephora", "value": int((analysis_set["analysis_group"] == "canephora").sum())},
        {"metric": "n_eugenioides", "value": int((analysis_set["analysis_group"] == "eugenioides").sum())},
    ]
    return pd.DataFrame(records)


def build_column_guide() -> pd.DataFrame:
    records = [
        ("canonical_seq_id", "Canonical sequencing identifier after metadata rescue."),
        ("rescue_status", "Matching status of each sample identifier against the accession metadata."),
        ("mean_ibs", "Mean identity-by-state value for a sample or group comparison."),
        ("mean_distance_1_minus_ibs", "Mean genetic distance estimated as one minus mean IBS."),
        ("centroid_distance_pc1_to_pc5", "Euclidean distance between group centroids in PCA space using PC1-PC5."),
        ("mean_introgressed_shift", "Mean difference between affinity to introgressed arabica and cultivated arabica across subgenomes."),
        ("median_rank_unweighted", "Unweighted median rank across transparent ranking metrics; lower values indicate stronger signal."),
        ("span_bp", "Genomic span between minimum and maximum positions within a hotspot summary."),
    ]
    return pd.DataFrame(records, columns=["column_name", "description"])


def select_columns(df: pd.DataFrame, preferred_order: Sequence[str]) -> pd.DataFrame:
    columns = [c for c in preferred_order if c in df.columns] + [c for c in df.columns if c not in preferred_order]
    return df.loc[:, columns]


def main() -> None:
    args = parse_args()

    analysis_workbook = Path(args.analysis_workbook)
    supplementary_workbook = Path(args.supplementary_workbook)
    output_workbook = Path(args.output_workbook)
    supplementary_output_workbook = Path(args.supplementary_output_workbook)

    print("Validating workbook structure")
    ensure_workbook_has_sheets(analysis_workbook, REQUIRED_ANALYSIS_SHEETS)
    ensure_workbook_has_sheets(supplementary_workbook, REQUIRED_SUPPLEMENTARY_SHEETS)

    print("Loading Part 1 workbooks")
    accession_master = pd.read_excel(supplementary_workbook, sheet_name="S1_accession_master")
    accession_master["seq_id"] = accession_master["seq_id"].astype(str)
    metadata_lookup = build_metadata_lookup(accession_master)

    sgc_match = pd.read_excel(analysis_workbook, sheet_name="13_sgC_sample_match")
    sge_match = pd.read_excel(analysis_workbook, sheet_name="14_sgE_sample_match")
    sgc_qc = pd.read_excel(analysis_workbook, sheet_name="16_sgC_sample_qc")
    sge_qc = pd.read_excel(analysis_workbook, sheet_name="17_sgE_sample_qc")
    sgc_pca = pd.read_excel(analysis_workbook, sheet_name="23_sgC_pca_scores")
    sge_pca = pd.read_excel(analysis_workbook, sheet_name="24_sgE_pca_scores")
    sgc_ibs = load_ibs_matrix(analysis_workbook, "27_sgC_ibs")
    sge_ibs = load_ibs_matrix(analysis_workbook, "28_sgE_ibs")
    sgc_top = pd.read_excel(analysis_workbook, sheet_name="31_sgC_top_contrasts")
    sge_top = pd.read_excel(analysis_workbook, sheet_name="32_sgE_top_contrasts")

    print("Rescuing sample identifiers and building the common analysis set")
    sgc_match_rescued = apply_rescue(sgc_match, "sample_id", metadata_lookup)
    sge_match_rescued = apply_rescue(sge_match, "sample_id", metadata_lookup)

    identity_rescue = build_identity_rescue_audit(
        sgc_match=sgc_match_rescued,
        sge_match=sge_match_rescued,
        accession_master=accession_master,
    )

    canonical_master = prepare_sample_qc_rescued(
        sgc_qc=sgc_qc,
        sge_qc=sge_qc,
        accession_master=accession_master,
        metadata_lookup=metadata_lookup,
        min_group_size=args.min_group_size,
    )
    analysis_set = canonical_master[
        canonical_master["retained_for_common_analysis"] &
        canonical_master["analysis_group"].notna()
    ].copy()

    allowed_ids = analysis_set["canonical_seq_id"].astype(str).tolist()

    print("Computing rescued distance and affinity summaries")
    sgc_ibs_rescued = canonicalize_ibs_matrix(sgc_ibs, metadata_lookup, allowed_ids=allowed_ids)
    sge_ibs_rescued = canonicalize_ibs_matrix(sge_ibs, metadata_lookup, allowed_ids=allowed_ids)

    group_lookup = build_group_lookup(analysis_set)

    sgc_group_ibs = summarize_group_ibs(sgc_ibs_rescued, group_lookup=group_lookup)
    sge_group_ibs = summarize_group_ibs(sge_ibs_rescued, group_lookup=group_lookup)
    group_ibs_delta = merge_distance_summaries(
        sgc_group_ibs,
        sge_group_ibs,
        value_column_sgc="mean_ibs",
        value_column_sge="mean_ibs",
        merged_value_name="mean_ibs",
    )

    sgc_affinity = compute_group_affinity(sgc_ibs_rescued, analysis_set, subgenome_label="sgC")
    sge_affinity = compute_group_affinity(sge_ibs_rescued, analysis_set, subgenome_label="sgE")

    sgc_pca_rescued = build_pca_rescued_table(
        pca_df=sgc_pca,
        accession_master=accession_master,
        metadata_lookup=metadata_lookup,
        allowed_ids=allowed_ids,
    )
    sge_pca_rescued = build_pca_rescued_table(
        pca_df=sge_pca,
        accession_master=accession_master,
        metadata_lookup=metadata_lookup,
        allowed_ids=allowed_ids,
    )

    sgc_centroids, sgc_centroid_distances = compute_pca_centroids(
        sgc_pca_rescued, analysis_set, subgenome_label="sgC", n_components=5
    )
    sge_centroids, sge_centroid_distances = compute_pca_centroids(
        sge_pca_rescued, analysis_set, subgenome_label="sgE", n_components=5
    )
    centroid_distance_delta = merge_distance_summaries(
        sgc_centroid_distances,
        sge_centroid_distances,
        value_column_sgc="sgC_centroid_distance_pc1_to_pc5",
        value_column_sge="sgE_centroid_distance_pc1_to_pc5",
        merged_value_name="centroid_distance_pc1_to_pc5",
    )

    arabica_priority, asymmetry_ranking = compute_priority_tables(
        analysis_set=analysis_set,
        sgc_affinity=sgc_affinity,
        sge_affinity=sge_affinity,
    )

    sgc_contrast_summary, sgc_hotspots = summarize_hotspots(sgc_top, subgenome_label="sgC")
    sge_contrast_summary, sge_hotspots = summarize_hotspots(sge_top, subgenome_label="sgE")

    signal_overview = pd.concat([sgc_contrast_summary, sge_contrast_summary], ignore_index=True)
    hotspot_manuscript_sgc = (
        sgc_hotspots.groupby("contrast", group_keys=False)
        .head(args.top_markers_per_contrast)
        .reset_index(drop=True)
    )
    hotspot_manuscript_sge = (
        sge_hotspots.groupby("contrast", group_keys=False)
        .head(args.top_markers_per_contrast)
        .reset_index(drop=True)
    )

    run_info = build_run_summary_sheet(args)
    input_inventory = build_input_inventory(analysis_workbook, supplementary_workbook)
    analysis_scope = build_analysis_scope(analysis_set)
    group_counts = summarize_group_counts(analysis_set)
    column_guide = build_column_guide()

    # Manuscript workbook
    print(f"Writing workbook: {output_workbook}")
    output_workbook.parent.mkdir(parents=True, exist_ok=True)
    with pd.ExcelWriter(output_workbook, engine="openpyxl") as writer:
        run_info.to_excel(writer, sheet_name="01_run_parameters", index=False)
        input_inventory.to_excel(writer, sheet_name="02_input_inventory", index=False)
        analysis_scope.to_excel(writer, sheet_name="03_analysis_scope", index=False)
        group_counts.to_excel(writer, sheet_name="04_group_composition", index=False)
        identity_rescue.to_excel(writer, sheet_name="05_identity_rescue_audit", index=False)
        select_columns(
            analysis_set,
            [
                "canonical_seq_id", "accession_name", "species_name", "variety", "analysis_group",
                "country_of_origin", "genome_structure", "sgC_mean_depth", "sgE_mean_depth",
                "delta_heterozygosity_sgC_minus_sgE", "delta_call_rate_sgC_minus_sgE",
                "delta_non_reference_rate_sgC_minus_sgE"
            ]
        ).to_excel(writer, sheet_name="06_common_sample_master", index=False)
        sgc_group_ibs.to_excel(writer, sheet_name="07_group_ibs_sgC", index=False)
        sge_group_ibs.to_excel(writer, sheet_name="08_group_ibs_sgE", index=False)
        group_ibs_delta.to_excel(writer, sheet_name="09_group_ibs_delta", index=False)
        sgc_centroids.to_excel(writer, sheet_name="10_pca_centroids_sgC", index=False)
        sge_centroids.to_excel(writer, sheet_name="11_pca_centroids_sgE", index=False)
        centroid_distance_delta.to_excel(writer, sheet_name="12_centroid_distance_delta", index=False)
        arabica_priority.to_excel(writer, sheet_name="13_accession_priority", index=False)
        asymmetry_ranking.to_excel(writer, sheet_name="14_asymmetry_ranking", index=False)
        hotspot_manuscript_sgc.to_excel(writer, sheet_name="15_sgC_hotspots", index=False)
        hotspot_manuscript_sge.to_excel(writer, sheet_name="16_sgE_hotspots", index=False)
        signal_overview.to_excel(writer, sheet_name="17_contrast_signal_overview", index=False)
        column_guide.to_excel(writer, sheet_name="18_column_guide", index=False)

    # Supplementary workbook
    print(f"Writing workbook: {supplementary_output_workbook}")
    supplementary_output_workbook.parent.mkdir(parents=True, exist_ok=True)
    with pd.ExcelWriter(supplementary_output_workbook, engine="openpyxl") as writer:
        identity_rescue.to_excel(writer, sheet_name="S2_identity_rescue_audit", index=False)
        canonical_master.to_excel(writer, sheet_name="S2_rescued_sample_master", index=False)
        analysis_set.to_excel(writer, sheet_name="S2_common_analysis_set", index=False)
        sgc_pca_rescued.to_excel(writer, sheet_name="S2_sgC_pca_scores", index=False)
        sge_pca_rescued.to_excel(writer, sheet_name="S2_sgE_pca_scores", index=False)
        sgc_affinity.to_excel(writer, sheet_name="S2_sgC_group_affinity", index=False)
        sge_affinity.to_excel(writer, sheet_name="S2_sgE_group_affinity", index=False)
        sgc_group_ibs.to_excel(writer, sheet_name="S2_sgC_group_ibs", index=False)
        sge_group_ibs.to_excel(writer, sheet_name="S2_sgE_group_ibs", index=False)
        sgc_centroids.to_excel(writer, sheet_name="S2_sgC_pca_centroids", index=False)
        sge_centroids.to_excel(writer, sheet_name="S2_sgE_pca_centroids", index=False)
        sgc_centroid_distances.to_excel(writer, sheet_name="S2_sgC_centroid_dist", index=False)
        sge_centroid_distances.to_excel(writer, sheet_name="S2_sgE_centroid_dist", index=False)
        sgc_top.to_excel(writer, sheet_name="S2_sgC_top_markers", index=False)
        sge_top.to_excel(writer, sheet_name="S2_sgE_top_markers", index=False)
        sgc_hotspots.to_excel(writer, sheet_name="S2_sgC_hotspot_summary", index=False)
        sge_hotspots.to_excel(writer, sheet_name="S2_sgE_hotspot_summary", index=False)
        column_guide.to_excel(writer, sheet_name="S2_column_guide", index=False)

    print("Part 2 analysis finished successfully.")


if __name__ == "__main__":
    main()
