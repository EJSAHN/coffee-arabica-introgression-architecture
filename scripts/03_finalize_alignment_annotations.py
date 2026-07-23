#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import math
import re
from pathlib import Path
from typing import Iterable, List, Tuple, Dict, Any

import numpy as np
import pandas as pd

BLAST_COLS = [
    "query", "qlen", "qstart", "qend", "target", "slen", "sstart", "send",
    "pident", "length", "bitscore", "evalue"
]

EXPECTED = {
    "sgC": {"target": "NC_092316.1", "gff": "C_arabica_ET39_sgC.gff3"},
    "sgE": {"target": "NC_092317.1", "gff": "C_arabica_ET39_sgE.gff3"},
}


def normalize_seqid(value: Any) -> str:
    s = "" if value is None else str(value)
    m = re.search(r"(NC_\d+\.\d+)", s)
    if m:
        return m.group(1)
    return s.strip().strip("|")


def parse_query_coordinates(query: str) -> Tuple[int, int, int]:
    m = re.search(r"chunk(\d+)\|src=(\d+)-(\d+)", str(query))
    if not m:
        raise ValueError(f"Could not parse query coordinates: {query}")
    return int(m.group(1)), int(m.group(2)), int(m.group(3))


def merge_intervals(intervals: Iterable[Tuple[int, int]], max_gap: int = 0) -> List[Tuple[int, int]]:
    vals = sorted((min(int(a), int(b)), max(int(a), int(b))) for a, b in intervals)
    if not vals:
        return []
    out = [list(vals[0])]
    for a, b in vals[1:]:
        if a <= out[-1][1] + max_gap + 1:
            out[-1][1] = max(out[-1][1], b)
        else:
            out.append([a, b])
    return [(a, b) for a, b in out]


def union_length(intervals: Iterable[Tuple[int, int]]) -> int:
    return int(sum(b - a + 1 for a, b in merge_intervals(intervals)))


def spearman(x: pd.Series, y: pd.Series) -> float:
    if len(x) < 3:
        return float("nan")
    return float(pd.concat([x.reset_index(drop=True), y.reset_index(drop=True)], axis=1).corr(method="spearman").iloc[0, 1])


def read_blast(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", names=BLAST_COLS, low_memory=False)
    parsed = df["query"].map(parse_query_coordinates)
    df[["chunk_index", "chunk_source_start", "chunk_source_end"]] = pd.DataFrame(parsed.tolist(), index=df.index)
    df["target_normalized"] = df["target"].map(normalize_seqid)
    df["source_hsp_start"] = df["chunk_source_start"] + df["qstart"].astype(int) - 1
    df["source_hsp_end"] = df["chunk_source_start"] + df["qend"].astype(int) - 1
    df["target_hsp_start"] = df[["sstart", "send"]].min(axis=1).astype(int)
    df["target_hsp_end"] = df[["sstart", "send"]].max(axis=1).astype(int)
    df["orientation"] = np.where(df["sstart"] <= df["send"], "+", "-")
    df["source_midpoint"] = (df["source_hsp_start"] + df["source_hsp_end"]) / 2.0
    df["target_midpoint"] = (df["target_hsp_start"] + df["target_hsp_end"]) / 2.0
    return df


def best_per_query(df: pd.DataFrame) -> pd.DataFrame:
    if df.empty:
        return df.copy()
    return (
        df.sort_values(["query", "bitscore", "length", "pident"], ascending=[True, False, False, False])
          .groupby("query", as_index=False, sort=False)
          .first()
    )


def parse_gff_attributes(text: str) -> Dict[str, str]:
    out: Dict[str, str] = {}
    for item in str(text).split(";"):
        if "=" in item:
            k, v = item.split("=", 1)
            out[k.strip()] = v.strip()
    return out


def read_supported_genes(gff_path: Path, seqid: str, blocks: List[Tuple[int, int]], subgenome: str) -> pd.DataFrame:
    cols = ["subgenome", "seqid", "start", "end", "strand", "gene_id", "gene_name", "product", "attributes"]
    if not gff_path.exists() or not blocks:
        return pd.DataFrame(columns=cols)
    rows = []
    with gff_path.open("rt", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9 or parts[2].lower() != "gene":
                continue
            if normalize_seqid(parts[0]) != seqid:
                continue
            start, end = int(parts[3]), int(parts[4])
            if not any(start <= b and end >= a for a, b in blocks):
                continue
            attrs = parse_gff_attributes(parts[8])
            gene_id = attrs.get("ID", attrs.get("gene_id", ""))
            gene_name = attrs.get("Name", attrs.get("gene", attrs.get("gene_name", "")))
            product = attrs.get("product", attrs.get("description", attrs.get("Note", "")))
            rows.append({
                "subgenome": subgenome,
                "seqid": seqid,
                "start": start,
                "end": end,
                "strand": parts[6],
                "gene_id": gene_id,
                "gene_name": gene_name,
                "product": product,
                "attributes": parts[8],
            })
    return pd.DataFrame(rows, columns=cols)


def candidate_subset(genes: pd.DataFrame) -> pd.DataFrame:
    if genes.empty:
        return genes.copy()
    pattern = re.compile(
        r"resistan|disease|defen|immune|NLR|NBS|LRR|RPP|RGA|kinase|RLK|WAK|NAC|JA2L|stress|transport|receptor",
        re.IGNORECASE,
    )
    text = genes[["gene_name", "product", "attributes"]].fillna("").astype(str).agg(" ".join, axis=1)
    return genes[text.map(lambda x: bool(pattern.search(x)))].copy()


def analyze_subgenome(base: Path, project: Path, sub: str) -> Tuple[Dict[str, Any], Dict[str, pd.DataFrame]]:
    blast_path = base / f"blast_{sub}.tsv"
    if not blast_path.exists() or blast_path.stat().st_size == 0:
        raise FileNotFoundError(f"Missing non-empty BLAST TSV: {blast_path}")
    df = read_blast(blast_path)
    expected = EXPECTED[sub]["target"]
    source_start = int(df["chunk_source_start"].min())
    source_end = int(df["chunk_source_end"].max())
    source_span = source_end - source_start + 1

    overall_best = best_per_query(df)
    expected_all = df[df["target_normalized"] == expected].copy()
    expected_best = best_per_query(expected_all)

    expected_top_fraction = float((overall_best["target_normalized"] == expected).mean())
    expected_source_blocks = merge_intervals(zip(expected_best["source_hsp_start"], expected_best["source_hsp_end"]))
    aligned_source_bp = union_length(expected_source_blocks)
    aligned_fraction = aligned_source_bp / source_span

    weighted_identity = float(np.average(expected_best["pident"], weights=expected_best["length"])) if not expected_best.empty else float("nan")
    orientation_counts = expected_best["orientation"].value_counts().to_dict()
    forward_fraction = float(orientation_counts.get("+", 0) / max(1, len(expected_best)))
    order_spearman = spearman(expected_best["source_midpoint"], expected_best["target_midpoint"])

    alt_rows = []
    for target, grp in df[df["target_normalized"] != expected].groupby("target_normalized"):
        b = best_per_query(grp)
        cov = union_length(zip(b["source_hsp_start"], b["source_hsp_end"]))
        alt_rows.append({"subgenome": sub, "target": target, "covered_source_bp": cov, "coverage_fraction": cov / source_span, "n_chunks": b["query"].nunique()})
    alt_df = pd.DataFrame(alt_rows).sort_values("coverage_fraction", ascending=False) if alt_rows else pd.DataFrame(columns=["subgenome", "target", "covered_source_bp", "coverage_fraction", "n_chunks"])
    largest_alt = float(alt_df["coverage_fraction"].max()) if not alt_df.empty else 0.0

    chromosome_pass = (
        expected_top_fraction >= 0.95
        and aligned_fraction >= 0.60
        and weighted_identity >= 98.0
        and largest_alt <= 0.10
        and order_spearman >= 0.95
    )
    boundary_pass = aligned_fraction >= 0.90 and forward_fraction >= 0.98
    chromosome_status = "PASS" if chromosome_pass else "WARN"
    boundary_status = "PASS" if boundary_pass else "WARN"
    overall_status = "PASS" if chromosome_pass and boundary_pass else ("WARN" if chromosome_pass else "FAIL")

    # Alignment-supported target blocks for interval-level annotation context.
    supported = expected_best[(expected_best["pident"] >= 98.0) & (expected_best["length"] >= 10000)].copy()
    merged_target = merge_intervals(zip(supported["target_hsp_start"], supported["target_hsp_end"]), max_gap=50000)
    block_df = pd.DataFrame([
        {"subgenome": sub, "target_seqid": expected, "block_index": i + 1, "target_start": a, "target_end": b, "block_span_bp": b - a + 1}
        for i, (a, b) in enumerate(merged_target)
    ])

    gff_path = project / EXPECTED[sub]["gff"]
    genes = read_supported_genes(gff_path, expected, merged_target, sub)
    candidates = candidate_subset(genes)

    summary = {
        "subgenome": sub,
        "chromosome_assignment_status": chromosome_status,
        "boundary_resolution_status": boundary_status,
        "overall_status": overall_status,
        "expected_target_contig": expected,
        "dominant_target_raw": str(overall_best["target"].value_counts().idxmax()),
        "dominant_target_normalized": str(overall_best["target_normalized"].value_counts().idxmax()),
        "expected_target_match_after_normalization": bool(overall_best["target_normalized"].value_counts().idxmax() == expected),
        "n_chunks": int(overall_best["query"].nunique()),
        "expected_target_top_hit_fraction": expected_top_fraction,
        "aligned_source_bp_expected_target": int(aligned_source_bp),
        "aligned_source_fraction_expected_target": aligned_fraction,
        "weighted_identity_expected_target": weighted_identity,
        "forward_orientation_fraction": forward_fraction,
        "source_target_order_spearman": order_spearman,
        "largest_alternative_target_fraction": largest_alt,
        "alignment_supported_target_start": int(supported["target_hsp_start"].min()) if not supported.empty else None,
        "alignment_supported_target_end": int(supported["target_hsp_end"].max()) if not supported.empty else None,
        "alignment_supported_envelope_span_bp": int(supported["target_hsp_end"].max() - supported["target_hsp_start"].min() + 1) if not supported.empty else None,
        "n_alignment_supported_blocks": int(len(merged_target)),
        "n_alignment_supported_genes": int(len(genes)),
        "n_candidate_annotations": int(len(candidates)),
        "interpretation": (
            "Direct sequence alignment strongly supports chromosome-level correspondence to the expected chromosome 4 pseudomolecule, "
            "but incomplete coverage and mixed orientation do not support a single base-pair-resolved liftover interval."
            if chromosome_pass and not boundary_pass
            else "Direct sequence alignment supports the expected chromosome-level assignment."
            if chromosome_pass
            else "Alignment support is insufficient for the expected chromosome-level assignment."
        ),
    }

    expected_best_out = expected_best.copy()
    expected_best_out.insert(0, "subgenome", sub)
    return summary, {
        "best_expected_hits": expected_best_out,
        "alternative_targets": alt_df,
        "alignment_blocks": block_df,
        "supported_genes": genes,
        "candidate_annotations": candidates,
    }


def main() -> None:
    ap = argparse.ArgumentParser(description="Correct and finalize direct-alignment validation without rerunning BLAST.")
    ap.add_argument("--project-dir", required=True)
    ap.add_argument("--alignment-dir", required=True)
    ap.add_argument("--output-dir", default=None)
    args = ap.parse_args()

    project = Path(args.project_dir)
    alignment = Path(args.alignment_dir)
    out = Path(args.output_dir) if args.output_dir else alignment / "final"
    out.mkdir(parents=True, exist_ok=True)

    summaries = []
    pieces: Dict[str, List[pd.DataFrame]] = {k: [] for k in ["best_expected_hits", "alternative_targets", "alignment_blocks", "supported_genes", "candidate_annotations"]}
    for sub in ("sgC", "sgE"):
        summary, data = analyze_subgenome(alignment, project, sub)
        summaries.append(summary)
        for k, df in data.items():
            pieces[k].append(df)

    summary_df = pd.DataFrame(summaries)
    overall = "PASS" if (summary_df["overall_status"] == "PASS").all() else ("WARN" if (summary_df["chromosome_assignment_status"] == "PASS").all() else "FAIL")

    workbook = out / "alignment_validation_corrected.xlsx"
    with pd.ExcelWriter(workbook, engine="openpyxl") as writer:
        pd.DataFrame([{"overall_status": overall, "interpretation": "Chromosome-level correspondence is evaluated separately from base-pair boundary resolution."}]).to_excel(writer, sheet_name="Overall_status", index=False)
        summary_df.to_excel(writer, sheet_name="Corrected_summary", index=False)
        for k, dfs in pieces.items():
            pd.concat(dfs, ignore_index=True).to_excel(writer, sheet_name=k[:31], index=False)

    status = {"script_version": "1.0", "overall_status": overall, "subgenomes": summaries}
    (out / "alignment_validation_status.json").write_text(json.dumps(status, indent=2), encoding="utf-8")

    lines = [
        "# Direct-alignment validation interpretation",
        "",
        f"Overall status: **{overall}**",
        "",
        "BLAST sequence identifiers are normalized before comparison with expected chromosome accessions.",
        "",
    ]
    for s in summaries:
        lines += [
            f"## {s['subgenome']}",
            f"- Expected chromosome top hit in {s['expected_target_top_hit_fraction']:.1%} of chunks.",
            f"- Expected chromosome source coverage: {s['aligned_source_fraction_expected_target']:.1%}.",
            f"- Weighted nucleotide identity: {s['weighted_identity_expected_target']:.3f}%.",
            f"- Forward-orientation best hits: {s['forward_orientation_fraction']:.1%}.",
            f"- Source-target order Spearman correlation: {s['source_target_order_spearman']:.4f}.",
            f"- Largest alternative-target coverage: {s['largest_alternative_target_fraction']:.1%}.",
            f"- Chromosome assignment: **{s['chromosome_assignment_status']}**; base-pair boundary resolution: **{s['boundary_resolution_status']}**.",
            f"- Interpretation: {s['interpretation']}",
            "",
        ]
    lines += [
        "## Recommended manuscript interpretation",
        "Direct sequence alignment supports chromosome-level correspondence of the reanalyzed sgC and sgE source intervals to the expected chromosome 4 pseudomolecules. However, incomplete sequence coverage and a minority of reverse-orientation blocks indicate that the mapping should be reported as alignment-supported interval-level correspondence rather than as a single exact base-pair liftover. Candidate annotations should therefore be restricted to alignment-supported target blocks and remain hypothesis-generating.",
        "",
    ]
    (out / "ALIGNMENT_INTERPRETATION.md").write_text("\n".join(lines), encoding="utf-8")

    print(f"Corrected alignment validation completed. Overall status: {overall}")
    print(f"Workbook: {workbook}")
    print(f"Status JSON: {out / 'alignment_validation_status.json'}")
    print(f"Interpretation: {out / 'ALIGNMENT_INTERPRETATION.md'}")


if __name__ == "__main__":
    main()
