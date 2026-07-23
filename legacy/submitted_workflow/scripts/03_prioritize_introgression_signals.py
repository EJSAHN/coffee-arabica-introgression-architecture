#!/usr/bin/env python3
"""
Aggregate accession-level and hotspot-level evidence for introgression analysis.

The script combines population-structure summaries, accession affinity metrics,
subgenome-asymmetry rankings, and high-difference marker intervals to identify
shared, same-contig, and subgenome-specific hotspot patterns.
"""

from __future__ import annotations

import argparse
import datetime as dt
from pathlib import Path
from typing import Dict, Iterable, List, Tuple

import numpy as np
import pandas as pd


REQUIRED_MANUSCRIPT_SHEETS = [
    "04_group_composition",
    "06_common_sample_master",
    "07_group_ibs_sgC",
    "08_group_ibs_sgE",
    "09_group_ibs_delta",
    "12_centroid_distance_delta",
    "13_accession_priority",
    "14_asymmetry_ranking",
    "15_sgC_hotspots",
    "16_sgE_hotspots",
    "17_contrast_signal_overview",
]

REQUIRED_SUPPLEMENTARY_SHEETS = [
    "S2_common_analysis_set",
    "S2_sgC_group_affinity",
    "S2_sgE_group_affinity",
    "S2_sgC_top_markers",
    "S2_sgE_top_markers",
    "S2_sgC_hotspot_summary",
    "S2_sgE_hotspot_summary",
    "S2_sgC_pca_scores",
    "S2_sgE_pca_scores",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Build evidence-focused manuscript tables from coffee introgression Part 2 outputs."
    )
    parser.add_argument(
        "--analysis-workbook",
        required=True,
        help="Path to the Part 2 manuscript workbook.",
    )
    parser.add_argument(
        "--supplementary-workbook",
        required=True,
        help="Path to the Part 2 supplementary workbook.",
    )
    parser.add_argument(
        "--output-workbook",
        required=True,
        help="Path to the manuscript-oriented Part 3 workbook to be written.",
    )
    parser.add_argument(
        "--supplementary-output-workbook",
        required=True,
        help="Path to the supplementary Part 3 workbook to be written.",
    )
    return parser.parse_args()


def require_file(path_string: str) -> Path:
    path = Path(path_string)
    if not path.exists():
        raise FileNotFoundError(f"Input file not found: {path}")
    return path


def list_excel_sheets(path: Path) -> List[str]:
    return pd.ExcelFile(path).sheet_names


def validate_workbook_structure(path: Path, required_sheets: List[str]) -> None:
    sheet_names = list_excel_sheets(path)
    missing = [sheet for sheet in required_sheets if sheet not in sheet_names]
    if missing:
        missing_string = ", ".join(missing)
        raise ValueError(f"Workbook is missing required sheets: {missing_string}")


def read_required_sheets(path: Path, sheets: Iterable[str]) -> Dict[str, pd.DataFrame]:
    excel_file = pd.ExcelFile(path)
    return {sheet: pd.read_excel(excel_file, sheet_name=sheet) for sheet in sheets}


def safe_mean(values: Iterable[float]) -> float:
    values = [value for value in values if pd.notna(value)]
    if not values:
        return np.nan
    return float(np.mean(values))


def safe_percentile_rank(series: pd.Series, higher_is_better: bool = True) -> pd.Series:
    numeric = pd.to_numeric(series, errors="coerce")
    if higher_is_better:
        return numeric.rank(method="average", pct=True, ascending=True)
    return numeric.rank(method="average", pct=True, ascending=False)


def add_rank(series: pd.Series, higher_is_better: bool = True) -> pd.Series:
    numeric = pd.to_numeric(series, errors="coerce")
    return numeric.rank(method="min", ascending=not higher_is_better)


def build_panel_overview(common_panel: pd.DataFrame) -> pd.DataFrame:
    rows = []
    rows.append(
        {
            "summary_type": "overall",
            "category": "all_accessions",
            "n_accessions": int(len(common_panel)),
            "n_unique_accessions": int(common_panel["accession_name"].nunique()),
            "n_countries": int(common_panel["country_of_origin"].dropna().nunique()),
        }
    )
    for column in ["analysis_group", "species_name", "variety", "country_of_origin", "genome_structure"]:
        counts = (
            common_panel.groupby(column, dropna=False)
            .size()
            .reset_index(name="n_accessions")
            .rename(columns={column: "category"})
        )
        counts.insert(0, "summary_type", column)
        counts["n_unique_accessions"] = counts["n_accessions"]
        counts["n_countries"] = np.nan
        rows.extend(counts.to_dict(orient="records"))
    output = pd.DataFrame(rows)
    return output


def prepare_group_pair_key(df: pd.DataFrame, group_a_col: str = "group_a", group_b_col: str = "group_b") -> pd.DataFrame:
    out = df.copy()
    pair_keys = out.apply(
        lambda row: tuple(sorted([str(row[group_a_col]), str(row[group_b_col])])), axis=1
    )
    out["pair_key"] = pair_keys
    out["pair_group_1"] = out["pair_key"].apply(lambda value: value[0])
    out["pair_group_2"] = out["pair_key"].apply(lambda value: value[1])
    return out


def collapse_group_relationships(group_ibs_delta: pd.DataFrame, centroid_delta: pd.DataFrame) -> pd.DataFrame:
    ibs = prepare_group_pair_key(group_ibs_delta)
    numeric_ibs_cols = [column for column in ibs.columns if column not in {"group_a", "group_b", "pair_key", "pair_group_1", "pair_group_2"}]
    ibs_unique = (
        ibs.groupby(["pair_group_1", "pair_group_2"], as_index=False)[numeric_ibs_cols]
        .mean(numeric_only=True)
    )

    centroid = prepare_group_pair_key(centroid_delta)
    numeric_centroid_cols = [
        column for column in centroid.columns
        if column not in {"group_a", "group_b", "pair_key", "pair_group_1", "pair_group_2"}
    ]
    centroid_unique = (
        centroid.groupby(["pair_group_1", "pair_group_2"], as_index=False)[numeric_centroid_cols]
        .mean(numeric_only=True)
    )

    merged = ibs_unique.merge(
        centroid_unique,
        on=["pair_group_1", "pair_group_2"],
        how="outer",
    )

    merged["subgenome_with_higher_similarity_by_mean_ibs"] = np.where(
        merged["delta_mean_ibs_sgC_minus_sgE"] > 0,
        "sgC",
        np.where(merged["delta_mean_ibs_sgC_minus_sgE"] < 0, "sgE", "tie"),
    )
    merged["subgenome_with_greater_centroid_separation"] = np.where(
        merged["delta_centroid_distance_pc1_to_pc5_sgC_minus_sgE"] > 0,
        "sgC",
        np.where(merged["delta_centroid_distance_pc1_to_pc5_sgC_minus_sgE"] < 0, "sgE", "tie"),
    )
    merged["absolute_delta_mean_ibs"] = merged["delta_mean_ibs_sgC_minus_sgE"].abs()
    merged["absolute_delta_centroid_distance"] = merged[
        "delta_centroid_distance_pc1_to_pc5_sgC_minus_sgE"
    ].abs()
    merged = merged.sort_values(
        ["absolute_delta_mean_ibs", "absolute_delta_centroid_distance", "pair_group_1", "pair_group_2"],
        ascending=[False, False, True, True],
    ).reset_index(drop=True)
    return merged


def build_accession_evidence(
    common_panel: pd.DataFrame,
    sgc_affinity: pd.DataFrame,
    sge_affinity: pd.DataFrame,
) -> pd.DataFrame:
    common = common_panel.copy()
    sgc = sgc_affinity.copy()
    sge = sge_affinity.copy()

    merged = common.merge(
        sgc,
        on=[
            "canonical_seq_id",
            "accession_name",
            "species_name",
            "variety",
            "analysis_group",
            "country_of_origin",
            "genome_structure",
        ],
        how="left",
    ).merge(
        sge,
        on=[
            "canonical_seq_id",
            "accession_name",
            "species_name",
            "variety",
            "analysis_group",
            "country_of_origin",
            "genome_structure",
        ],
        how="left",
    )

    merged["sgC_introgressed_shift"] = (
        merged["sgC_mean_ibs_to_arabica_introgressed"] - merged["sgC_mean_ibs_to_arabica_cultivated"]
    )
    merged["sgE_introgressed_shift"] = (
        merged["sgE_mean_ibs_to_arabica_introgressed"] - merged["sgE_mean_ibs_to_arabica_cultivated"]
    )
    merged["mean_introgressed_shift"] = merged[
        ["sgC_introgressed_shift", "sgE_introgressed_shift"]
    ].mean(axis=1)

    merged["sgC_wild_shift"] = (
        merged["sgC_mean_ibs_to_arabica_wild"] - merged["sgC_mean_ibs_to_arabica_cultivated"]
    )
    merged["sgE_wild_shift"] = (
        merged["sgE_mean_ibs_to_arabica_wild"] - merged["sgE_mean_ibs_to_arabica_cultivated"]
    )
    merged["mean_wild_shift"] = merged[["sgC_wild_shift", "sgE_wild_shift"]].mean(axis=1)

    merged["mean_canephora_affinity"] = merged[
        ["sgC_mean_ibs_to_canephora", "sgE_mean_ibs_to_canephora"]
    ].mean(axis=1)
    merged["mean_eugenioides_affinity"] = merged[
        ["sgC_mean_ibs_to_eugenioides", "sgE_mean_ibs_to_eugenioides"]
    ].mean(axis=1)

    merged["sgC_parental_advantage_canephora_minus_eugenioides"] = (
        merged["sgC_mean_ibs_to_canephora"] - merged["sgC_mean_ibs_to_eugenioides"]
    )
    merged["sgE_parental_advantage_eugenioides_minus_canephora"] = (
        merged["sgE_mean_ibs_to_eugenioides"] - merged["sgE_mean_ibs_to_canephora"]
    )
    merged["parental_consistency_score"] = merged[
        [
            "sgC_parental_advantage_canephora_minus_eugenioides",
            "sgE_parental_advantage_eugenioides_minus_canephora",
        ]
    ].mean(axis=1)
    merged["parental_distortion_score"] = (
        merged["sgC_parental_advantage_canephora_minus_eugenioides"]
        - merged["sgE_parental_advantage_eugenioides_minus_canephora"]
    ).abs()

    merged["abs_delta_heterozygosity"] = merged["delta_heterozygosity_sgC_minus_sgE"].abs()
    merged["abs_delta_non_reference_rate"] = merged["delta_non_reference_rate_sgC_minus_sgE"].abs()
    merged["abs_delta_call_rate"] = merged["delta_call_rate_sgC_minus_sgE"].abs()
    merged["abs_delta_introgressed_shift"] = (
        merged["sgC_introgressed_shift"] - merged["sgE_introgressed_shift"]
    ).abs()

    arabica_mask = merged["analysis_group"].isin(
        ["arabica_cultivated", "arabica_introgressed", "arabica_wild"]
    )
    arabica = merged.loc[arabica_mask].copy()

    arabica["introgression_shift_percentile"] = safe_percentile_rank(arabica["mean_introgressed_shift"])
    arabica["canephora_affinity_percentile"] = safe_percentile_rank(arabica["mean_canephora_affinity"])
    arabica["parental_distortion_percentile"] = safe_percentile_rank(arabica["parental_distortion_score"])
    arabica["asymmetry_introgressed_shift_percentile"] = safe_percentile_rank(
        arabica["abs_delta_introgressed_shift"]
    )
    arabica["heterozygosity_delta_percentile"] = safe_percentile_rank(
        arabica["abs_delta_heterozygosity"]
    )
    arabica["non_reference_delta_percentile"] = safe_percentile_rank(
        arabica["abs_delta_non_reference_rate"]
    )

    arabica["introgression_affinity_percentile"] = arabica[
        ["introgression_shift_percentile", "canephora_affinity_percentile"]
    ].mean(axis=1)
    arabica["subgenome_asymmetry_percentile"] = arabica[
        [
            "parental_distortion_percentile",
            "asymmetry_introgressed_shift_percentile",
            "heterozygosity_delta_percentile",
            "non_reference_delta_percentile",
        ]
    ].mean(axis=1)
    arabica["combined_evidence_percentile"] = arabica[
        ["introgression_affinity_percentile", "subgenome_asymmetry_percentile"]
    ].mean(axis=1)

    arabica["overall_combined_evidence_rank"] = add_rank(
        arabica["combined_evidence_percentile"], higher_is_better=True
    )
    arabica["overall_introgression_affinity_rank"] = add_rank(
        arabica["introgression_affinity_percentile"], higher_is_better=True
    )
    arabica["overall_subgenome_asymmetry_rank"] = add_rank(
        arabica["subgenome_asymmetry_percentile"], higher_is_better=True
    )
    arabica["within_group_combined_evidence_rank"] = (
        arabica.groupby("analysis_group")["combined_evidence_percentile"]
        .rank(method="min", ascending=False)
    )
    arabica["within_group_introgression_affinity_rank"] = (
        arabica.groupby("analysis_group")["introgression_affinity_percentile"]
        .rank(method="min", ascending=False)
    )
    arabica["within_group_subgenome_asymmetry_rank"] = (
        arabica.groupby("analysis_group")["subgenome_asymmetry_percentile"]
        .rank(method="min", ascending=False)
    )

    merged = merged.merge(
        arabica[
            [
                "canonical_seq_id",
                "introgression_shift_percentile",
                "canephora_affinity_percentile",
                "parental_distortion_percentile",
                "asymmetry_introgressed_shift_percentile",
                "heterozygosity_delta_percentile",
                "non_reference_delta_percentile",
                "introgression_affinity_percentile",
                "subgenome_asymmetry_percentile",
                "combined_evidence_percentile",
                "overall_combined_evidence_rank",
                "overall_introgression_affinity_rank",
                "overall_subgenome_asymmetry_rank",
                "within_group_combined_evidence_rank",
                "within_group_introgression_affinity_rank",
                "within_group_subgenome_asymmetry_rank",
            ]
        ],
        on="canonical_seq_id",
        how="left",
    )

    merged = merged.sort_values(
        ["combined_evidence_percentile", "mean_introgressed_shift", "parental_distortion_score", "accession_name"],
        ascending=[False, False, False, True],
    ).reset_index(drop=True)
    return merged


def count_markers_in_interval(marker_df: pd.DataFrame, contrast: str, contig: str, start: float, end: float) -> int:
    subset = marker_df[
        (marker_df["contrast"] == contrast)
        & (marker_df["contig"] == contig)
        & (marker_df["pos"] >= start)
        & (marker_df["pos"] <= end)
    ]
    return int(len(subset))


def build_hotspot_tables(
    sgc_hotspots: pd.DataFrame,
    sge_hotspots: pd.DataFrame,
    sgc_markers: pd.DataFrame,
    sge_markers: pd.DataFrame,
) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    overlap_rows: List[Dict[str, object]] = []
    same_contig_nonoverlap_rows: List[Dict[str, object]] = []
    subgenome_specific_rows: List[Dict[str, object]] = []

    sge_lookup = sge_hotspots.groupby(["contrast", "contig"])
    sgc_keys = set(zip(sgc_hotspots["contrast"], sgc_hotspots["contig"]))
    sge_keys = set(zip(sge_hotspots["contrast"], sge_hotspots["contig"]))

    for _, sgc_row in sgc_hotspots.iterrows():
        key = (sgc_row["contrast"], sgc_row["contig"])
        sgc_marker_count = count_markers_in_interval(
            sgc_markers, sgc_row["contrast"], sgc_row["contig"], sgc_row["min_pos"], sgc_row["max_pos"]
        )
        if key not in sge_keys:
            subgenome_specific_rows.append(
                {
                    "contrast": sgc_row["contrast"],
                    "contig": sgc_row["contig"],
                    "subgenome": "sgC",
                    "n_markers_hotspot": sgc_row["n_markers"],
                    "n_markers_top_marker_table": sgc_marker_count,
                    "min_pos": sgc_row["min_pos"],
                    "max_pos": sgc_row["max_pos"],
                    "span_bp": sgc_row["span_bp"],
                    "mean_abs_delta_af": sgc_row["mean_abs_delta_af"],
                    "max_abs_delta_af": sgc_row["max_abs_delta_af"],
                }
            )
            continue

        matches = sge_lookup.get_group(key)
        for _, sge_row in matches.iterrows():
            start = max(float(sgc_row["min_pos"]), float(sge_row["min_pos"]))
            end = min(float(sgc_row["max_pos"]), float(sge_row["max_pos"]))
            overlap_bp = max(0.0, end - start)
            union_start = min(float(sgc_row["min_pos"]), float(sge_row["min_pos"]))
            union_end = max(float(sgc_row["max_pos"]), float(sge_row["max_pos"]))
            union_bp = max(0.0, union_end - union_start)
            jaccard_bp = overlap_bp / union_bp if union_bp > 0 else np.nan

            sge_marker_count = count_markers_in_interval(
                sge_markers, sge_row["contrast"], sge_row["contig"], sge_row["min_pos"], sge_row["max_pos"]
            )
            sgc_overlap_marker_count = (
                count_markers_in_interval(sgc_markers, sgc_row["contrast"], sgc_row["contig"], start, end)
                if overlap_bp > 0
                else 0
            )
            sge_overlap_marker_count = (
                count_markers_in_interval(sge_markers, sge_row["contrast"], sge_row["contig"], start, end)
                if overlap_bp > 0
                else 0
            )

            row = {
                "contrast": sgc_row["contrast"],
                "contig": sgc_row["contig"],
                "sgC_n_markers_hotspot": sgc_row["n_markers"],
                "sgE_n_markers_hotspot": sge_row["n_markers"],
                "combined_n_markers_hotspot": int(sgc_row["n_markers"] + sge_row["n_markers"]),
                "sgC_top_marker_count_in_hotspot": sgc_marker_count,
                "sgE_top_marker_count_in_hotspot": sge_marker_count,
                "combined_top_marker_count_in_hotspot": int(sgc_marker_count + sge_marker_count),
                "sgC_min_pos": sgc_row["min_pos"],
                "sgC_max_pos": sgc_row["max_pos"],
                "sgE_min_pos": sge_row["min_pos"],
                "sgE_max_pos": sge_row["max_pos"],
                "union_min_pos": int(union_start),
                "union_max_pos": int(union_end),
                "union_span_bp": int(union_bp),
                "overlap_min_pos": int(start) if overlap_bp > 0 else np.nan,
                "overlap_max_pos": int(end) if overlap_bp > 0 else np.nan,
                "overlap_bp": int(overlap_bp),
                "jaccard_bp": jaccard_bp,
                "sgC_mean_abs_delta_af": sgc_row["mean_abs_delta_af"],
                "sgE_mean_abs_delta_af": sge_row["mean_abs_delta_af"],
                "combined_mean_abs_delta_af": safe_mean(
                    [sgc_row["mean_abs_delta_af"], sge_row["mean_abs_delta_af"]]
                ),
                "combined_max_abs_delta_af": max(
                    float(sgc_row["max_abs_delta_af"]), float(sge_row["max_abs_delta_af"])
                ),
                "sgC_overlap_fraction": overlap_bp / float(sgc_row["span_bp"]) if sgc_row["span_bp"] else np.nan,
                "sgE_overlap_fraction": overlap_bp / float(sge_row["span_bp"]) if sge_row["span_bp"] else np.nan,
                "hotspot_support_balance": (
                    min(float(sgc_row["n_markers"]), float(sge_row["n_markers"]))
                    / max(float(sgc_row["n_markers"]), float(sge_row["n_markers"]))
                    if max(float(sgc_row["n_markers"]), float(sge_row["n_markers"])) > 0
                    else np.nan
                ),
                "sgC_top_marker_count_in_overlap": sgc_overlap_marker_count,
                "sgE_top_marker_count_in_overlap": sge_overlap_marker_count,
                "combined_top_marker_count_in_overlap": int(sgc_overlap_marker_count + sge_overlap_marker_count),
                "overlap_status": "overlap" if overlap_bp > 0 else "same_contig_no_overlap",
            }

            if overlap_bp > 0:
                overlap_rows.append(row)
            else:
                same_contig_nonoverlap_rows.append(row)

    for _, sge_row in sge_hotspots.iterrows():
        key = (sge_row["contrast"], sge_row["contig"])
        sge_marker_count = count_markers_in_interval(
            sge_markers, sge_row["contrast"], sge_row["contig"], sge_row["min_pos"], sge_row["max_pos"]
        )
        if key not in sgc_keys:
            subgenome_specific_rows.append(
                {
                    "contrast": sge_row["contrast"],
                    "contig": sge_row["contig"],
                    "subgenome": "sgE",
                    "n_markers_hotspot": sge_row["n_markers"],
                    "n_markers_top_marker_table": sge_marker_count,
                    "min_pos": sge_row["min_pos"],
                    "max_pos": sge_row["max_pos"],
                    "span_bp": sge_row["span_bp"],
                    "mean_abs_delta_af": sge_row["mean_abs_delta_af"],
                    "max_abs_delta_af": sge_row["max_abs_delta_af"],
                }
            )

    overlap_df = pd.DataFrame(overlap_rows)
    same_contig_nonoverlap_df = pd.DataFrame(same_contig_nonoverlap_rows)
    subgenome_specific_df = pd.DataFrame(subgenome_specific_rows)

    if not overlap_df.empty:
        overlap_df["combined_marker_support_percentile"] = safe_percentile_rank(
            overlap_df["combined_n_markers_hotspot"]
        )
        overlap_df["combined_effect_percentile"] = safe_percentile_rank(
            overlap_df["combined_mean_abs_delta_af"]
        )
        overlap_df["interval_overlap_percentile"] = safe_percentile_rank(overlap_df["jaccard_bp"])
        overlap_df["support_balance_percentile"] = safe_percentile_rank(
            overlap_df["hotspot_support_balance"]
        )
        overlap_df["consensus_hotspot_percentile"] = overlap_df[
            [
                "combined_marker_support_percentile",
                "combined_effect_percentile",
                "interval_overlap_percentile",
                "support_balance_percentile",
            ]
        ].mean(axis=1)
        overlap_df["consensus_hotspot_rank"] = add_rank(
            overlap_df["consensus_hotspot_percentile"], higher_is_better=True
        )
        overlap_df = overlap_df.sort_values(
            [
                "consensus_hotspot_percentile",
                "combined_n_markers_hotspot",
                "combined_mean_abs_delta_af",
                "overlap_bp",
                "contig",
            ],
            ascending=[False, False, False, False, True],
        ).reset_index(drop=True)

    if not same_contig_nonoverlap_df.empty:
        same_contig_nonoverlap_df["combined_marker_support_percentile"] = safe_percentile_rank(
            same_contig_nonoverlap_df["combined_n_markers_hotspot"]
        )
        same_contig_nonoverlap_df["combined_effect_percentile"] = safe_percentile_rank(
            same_contig_nonoverlap_df["combined_mean_abs_delta_af"]
        )
        same_contig_nonoverlap_df["same_contig_support_percentile"] = same_contig_nonoverlap_df[
            ["combined_marker_support_percentile", "combined_effect_percentile"]
        ].mean(axis=1)
        same_contig_nonoverlap_df = same_contig_nonoverlap_df.sort_values(
            ["same_contig_support_percentile", "combined_n_markers_hotspot", "contig"],
            ascending=[False, False, True],
        ).reset_index(drop=True)

    if not subgenome_specific_df.empty:
        subgenome_specific_df["subgenome_specific_support_percentile"] = safe_percentile_rank(
            subgenome_specific_df["n_markers_hotspot"]
        )
        subgenome_specific_df = subgenome_specific_df.sort_values(
            ["subgenome_specific_support_percentile", "mean_abs_delta_af", "contig"],
            ascending=[False, False, True],
        ).reset_index(drop=True)

    return overlap_df, same_contig_nonoverlap_df, subgenome_specific_df


def build_contrast_consensus(
    contrast_overview: pd.DataFrame,
    overlap_df: pd.DataFrame,
    same_contig_nonoverlap_df: pd.DataFrame,
    subgenome_specific_df: pd.DataFrame,
) -> pd.DataFrame:
    pivot = contrast_overview.pivot(index="contrast", columns="subgenome")
    pivot.columns = [f"{metric}_{subgenome}" for metric, subgenome in pivot.columns]
    pivot = pivot.reset_index()

    if overlap_df.empty:
        overlap_summary = pd.DataFrame(columns=["contrast", "n_overlapping_hotspot_pairs", "n_overlapping_contigs", "total_overlap_bp", "mean_consensus_hotspot_percentile"])
    else:
        overlap_summary = overlap_df.groupby("contrast", as_index=False).agg(
            n_overlapping_hotspot_pairs=("contig", "size"),
            n_overlapping_contigs=("contig", pd.Series.nunique),
            total_overlap_bp=("overlap_bp", "sum"),
            mean_consensus_hotspot_percentile=("consensus_hotspot_percentile", "mean"),
        )

    if same_contig_nonoverlap_df.empty:
        nonoverlap_summary = pd.DataFrame(columns=["contrast", "n_same_contig_nonoverlap_pairs"])
    else:
        nonoverlap_summary = same_contig_nonoverlap_df.groupby("contrast", as_index=False).agg(
            n_same_contig_nonoverlap_pairs=("contig", "size")
        )

    if subgenome_specific_df.empty:
        specific_summary = pd.DataFrame(columns=["contrast", "n_subgenome_specific_hotspots"])
    else:
        specific_summary = subgenome_specific_df.groupby("contrast", as_index=False).agg(
            n_subgenome_specific_hotspots=("contig", "size")
        )

    merged = pivot.merge(overlap_summary, on="contrast", how="left")
    merged = merged.merge(nonoverlap_summary, on="contrast", how="left")
    merged = merged.merge(specific_summary, on="contrast", how="left")

    fill_zero_cols = [
        "n_overlapping_hotspot_pairs",
        "n_overlapping_contigs",
        "total_overlap_bp",
        "n_same_contig_nonoverlap_pairs",
        "n_subgenome_specific_hotspots",
    ]
    for column in fill_zero_cols:
        if column in merged.columns:
            merged[column] = merged[column].fillna(0)

    merged["mean_abs_delta_af_average"] = merged[
        ["mean_abs_delta_af_sgC", "mean_abs_delta_af_sgE"]
    ].mean(axis=1)
    merged["n_markers_total"] = merged[["n_markers_sgC", "n_markers_sgE"]].sum(axis=1)
    merged["subgenome_with_higher_mean_abs_delta_af"] = np.where(
        merged["mean_abs_delta_af_sgC"] > merged["mean_abs_delta_af_sgE"],
        "sgC",
        np.where(merged["mean_abs_delta_af_sgC"] < merged["mean_abs_delta_af_sgE"], "sgE", "tie"),
    )
    merged["subgenome_with_more_top_markers"] = np.where(
        merged["n_markers_sgC"] > merged["n_markers_sgE"],
        "sgC",
        np.where(merged["n_markers_sgC"] < merged["n_markers_sgE"], "sgE", "tie"),
    )
    merged["consensus_strength_percentile"] = safe_percentile_rank(merged["n_markers_total"]).add(
        safe_percentile_rank(merged["mean_abs_delta_af_average"]), fill_value=0
    )
    merged["consensus_strength_percentile"] = merged["consensus_strength_percentile"] / 2.0
    merged["consensus_strength_rank"] = add_rank(
        merged["consensus_strength_percentile"], higher_is_better=True
    )

    merged = merged.sort_values(
        ["consensus_strength_percentile", "n_overlapping_hotspot_pairs", "mean_abs_delta_af_average", "contrast"],
        ascending=[False, False, False, True],
    ).reset_index(drop=True)
    return merged


def build_group_parental_summary(all_accession_evidence: pd.DataFrame) -> pd.DataFrame:
    summary = (
        all_accession_evidence.groupby("analysis_group", as_index=False)
        .agg(
            n_accessions=("canonical_seq_id", "count"),
            mean_sgC_parental_advantage=("sgC_parental_advantage_canephora_minus_eugenioides", "mean"),
            mean_sgE_parental_advantage=("sgE_parental_advantage_eugenioides_minus_canephora", "mean"),
            mean_parental_consistency_score=("parental_consistency_score", "mean"),
            mean_parental_distortion_score=("parental_distortion_score", "mean"),
            mean_introgression_shift=("mean_introgressed_shift", "mean"),
            mean_canephora_affinity=("mean_canephora_affinity", "mean"),
            mean_eugenioides_affinity=("mean_eugenioides_affinity", "mean"),
            mean_abs_delta_introgressed_shift=("abs_delta_introgressed_shift", "mean"),
        )
    )
    summary = summary.sort_values("mean_parental_distortion_score", ascending=False).reset_index(drop=True)
    return summary


def select_top_rows(df: pd.DataFrame, sort_columns: List[str], ascending: List[bool], n: int) -> pd.DataFrame:
    if df.empty:
        return df.copy()
    return df.sort_values(sort_columns, ascending=ascending).head(n).reset_index(drop=True)


def build_manuscript_table_candidates(
    panel_overview: pd.DataFrame,
    group_relationships: pd.DataFrame,
    all_accession_evidence: pd.DataFrame,
    overlap_df: pd.DataFrame,
    contrast_consensus: pd.DataFrame,
) -> Dict[str, pd.DataFrame]:
    arabica = all_accession_evidence[
        all_accession_evidence["analysis_group"].isin(
            ["arabica_cultivated", "arabica_introgressed", "arabica_wild"]
        )
    ].copy()

    cultivated_candidates = select_top_rows(
        arabica[arabica["analysis_group"] == "arabica_cultivated"],
        ["combined_evidence_percentile", "introgression_affinity_percentile", "subgenome_asymmetry_percentile", "accession_name"],
        [False, False, False, True],
        10,
    )

    introgressed_reference = select_top_rows(
        arabica[arabica["analysis_group"] == "arabica_introgressed"],
        ["combined_evidence_percentile", "introgression_affinity_percentile", "accession_name"],
        [False, False, True],
        10,
    )

    asymmetry_candidates = select_top_rows(
        arabica,
        ["subgenome_asymmetry_percentile", "combined_evidence_percentile", "accession_name"],
        [False, False, True],
        12,
    )

    hotspot_candidates = select_top_rows(
        overlap_df,
        ["consensus_hotspot_percentile", "combined_n_markers_hotspot", "combined_mean_abs_delta_af", "overlap_bp", "contig"],
        [False, False, False, False, True],
        20,
    )

    relationship_focus = group_relationships[
        group_relationships["pair_group_1"].isin(["arabica_cultivated", "arabica_introgressed", "arabica_wild"])
        | group_relationships["pair_group_2"].isin(["arabica_cultivated", "arabica_introgressed", "arabica_wild"])
    ].copy()
    relationship_focus = relationship_focus.sort_values(
        ["absolute_delta_mean_ibs", "absolute_delta_centroid_distance", "pair_group_1", "pair_group_2"],
        ascending=[False, False, True, True],
    ).reset_index(drop=True)

    manuscript_tables = {
        "Table_01_panel_overview": panel_overview,
        "Table_02_group_relationships": relationship_focus,
        "Table_03_cultivated_candidates": cultivated_candidates,
        "Table_04_introgressed_reference": introgressed_reference,
        "Table_05_asymmetry_candidates": asymmetry_candidates,
        "Table_06_shared_hotspots": hotspot_candidates,
        "Table_07_contrast_consensus": contrast_consensus,
    }
    return manuscript_tables


def build_column_guide() -> pd.DataFrame:
    rows = [
        {
            "sheet_name": "03_panel_overview",
            "description": "Counts of accessions by group, species, variety, country and genome structure in the rescued common panel.",
        },
        {
            "sheet_name": "04_group_relationship_summary",
            "description": "Pairwise group relationships summarised from IBS and PCA centroid distance comparisons between sgC and sgE.",
        },
        {
            "sheet_name": "05_group_parental_summary",
            "description": "Mean parental-affinity and asymmetry statistics summarised by analysis group.",
        },
        {
            "sheet_name": "06_accession_evidence_all",
            "description": "Accession-level evidence matrix with affinity, asymmetry and percentile-based prioritisation metrics.",
        },
        {
            "sheet_name": "07_cultivated_priority_candidates",
            "description": "Arabica cultivated accessions ranked by combined evidence for introgression-like and asymmetric signals.",
        },
        {
            "sheet_name": "08_introgressed_reference_accessions",
            "description": "Arabica introgressed accessions ranked as a reference panel.",
        },
        {
            "sheet_name": "09_wild_transition_candidates",
            "description": "Arabica wild accessions ranked by combined evidence metrics.",
        },
        {
            "sheet_name": "10_asymmetry_priority_candidates",
            "description": "Arabica accessions ranked specifically by subgenome asymmetry metrics.",
        },
        {
            "sheet_name": "11_shared_hotspot_evidence",
            "description": "Shared hotspot intervals between sgC and sgE on the same contrast and contig with interval-overlap metrics.",
        },
        {
            "sheet_name": "12_same_contig_nonoverlap_evidence",
            "description": "Hotspots found on the same contrast and contig across subgenomes but without interval overlap.",
        },
        {
            "sheet_name": "13_subgenome_specific_hotspots",
            "description": "Hotspots without a same-contig counterpart in the other subgenome.",
        },
        {
            "sheet_name": "14_contrast_consensus",
            "description": "Contrast-level support summary integrating sgC and sgE signal strength and hotspot concordance.",
        },
        {
            "sheet_name": "15_manuscript_table_candidates",
            "description": "Concise tables intended for direct manuscript drafting or adaptation.",
        },
    ]
    return pd.DataFrame(rows)


def write_workbook(path: Path, sheet_map: Dict[str, pd.DataFrame]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with pd.ExcelWriter(path, engine="openpyxl") as writer:
        for sheet_name, df in sheet_map.items():
            safe_name = sheet_name[:31]
            df.to_excel(writer, sheet_name=safe_name, index=False)


def main() -> None:
    args = parse_args()

    analysis_path = require_file(args.analysis_workbook)
    supplementary_path = require_file(args.supplementary_workbook)
    output_path = Path(args.output_workbook)
    supplementary_output_path = Path(args.supplementary_output_workbook)

    print("Validating workbook structure")
    validate_workbook_structure(analysis_path, REQUIRED_MANUSCRIPT_SHEETS)
    validate_workbook_structure(supplementary_path, REQUIRED_SUPPLEMENTARY_SHEETS)

    print("Loading Part 2 workbooks")
    manuscript = read_required_sheets(analysis_path, REQUIRED_MANUSCRIPT_SHEETS)
    supplementary = read_required_sheets(supplementary_path, REQUIRED_SUPPLEMENTARY_SHEETS)

    common_panel = supplementary["S2_common_analysis_set"].copy()
    sgc_affinity = supplementary["S2_sgC_group_affinity"].copy()
    sge_affinity = supplementary["S2_sgE_group_affinity"].copy()
    sgc_hotspots = supplementary["S2_sgC_hotspot_summary"].copy()
    sge_hotspots = supplementary["S2_sgE_hotspot_summary"].copy()
    sgc_markers = supplementary["S2_sgC_top_markers"].copy()
    sge_markers = supplementary["S2_sgE_top_markers"].copy()

    print("Building group relationship summaries")
    group_relationships = collapse_group_relationships(
        manuscript["09_group_ibs_delta"], manuscript["12_centroid_distance_delta"]
    )

    print("Building accession evidence matrices")
    all_accession_evidence = build_accession_evidence(common_panel, sgc_affinity, sge_affinity)

    print("Deriving hotspot consensus tables")
    overlap_df, same_contig_nonoverlap_df, subgenome_specific_df = build_hotspot_tables(
        sgc_hotspots, sge_hotspots, sgc_markers, sge_markers
    )

    print("Summarising contrast-level evidence")
    contrast_consensus = build_contrast_consensus(
        manuscript["17_contrast_signal_overview"],
        overlap_df,
        same_contig_nonoverlap_df,
        subgenome_specific_df,
    )

    print("Preparing manuscript-facing tables")
    panel_overview = build_panel_overview(common_panel)
    group_parental_summary = build_group_parental_summary(all_accession_evidence)

    arabica_mask = all_accession_evidence["analysis_group"].isin(
        ["arabica_cultivated", "arabica_introgressed", "arabica_wild"]
    )
    cultivated_candidates = all_accession_evidence[
        all_accession_evidence["analysis_group"] == "arabica_cultivated"
    ].sort_values(
        ["combined_evidence_percentile", "introgression_affinity_percentile", "subgenome_asymmetry_percentile", "accession_name"],
        ascending=[False, False, False, True],
    ).reset_index(drop=True)

    introgressed_reference = all_accession_evidence[
        all_accession_evidence["analysis_group"] == "arabica_introgressed"
    ].sort_values(
        ["combined_evidence_percentile", "introgression_affinity_percentile", "accession_name"],
        ascending=[False, False, True],
    ).reset_index(drop=True)

    wild_candidates = all_accession_evidence[
        all_accession_evidence["analysis_group"] == "arabica_wild"
    ].sort_values(
        ["combined_evidence_percentile", "introgression_affinity_percentile", "accession_name"],
        ascending=[False, False, True],
    ).reset_index(drop=True)

    asymmetry_candidates = all_accession_evidence[arabica_mask].sort_values(
        ["subgenome_asymmetry_percentile", "combined_evidence_percentile", "accession_name"],
        ascending=[False, False, True],
    ).reset_index(drop=True)

    manuscript_tables = build_manuscript_table_candidates(
        panel_overview,
        group_relationships,
        all_accession_evidence,
        overlap_df,
        contrast_consensus,
    )

    run_parameters = pd.DataFrame(
        {
            "parameter": [
                "analysis_workbook",
                "supplementary_workbook",
                "execution_time_utc",
            ],
            "value": [
                str(analysis_path),
                str(supplementary_path),
                dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds").replace("+00:00", "Z"),
            ],
        }
    )

    input_inventory = pd.DataFrame(
        {
            "input_type": ["analysis_workbook", "supplementary_workbook"],
            "path": [str(analysis_path), str(supplementary_path)],
            "exists": [analysis_path.exists(), supplementary_path.exists()],
        }
    )

    manuscript_sheet_map = {
        "01_run_parameters": run_parameters,
        "02_input_inventory": input_inventory,
        "03_panel_overview": panel_overview,
        "04_group_relationship_summary": group_relationships,
        "05_group_parental_summary": group_parental_summary,
        "06_accession_evidence_all": all_accession_evidence,
        "07_cultivated_priority_candidates": cultivated_candidates,
        "08_introgressed_reference_accessions": introgressed_reference,
        "09_wild_transition_candidates": wild_candidates,
        "10_asymmetry_priority_candidates": asymmetry_candidates,
        "11_shared_hotspot_evidence": overlap_df,
        "12_same_contig_nonoverlap_evidence": same_contig_nonoverlap_df,
        "13_subgenome_specific_hotspots": subgenome_specific_df,
        "14_contrast_consensus": contrast_consensus,
        "15_column_guide": build_column_guide(),
    }

    supplementary_sheet_map = {
        "S3_common_panel": common_panel,
        "S3_accession_evidence_all": all_accession_evidence,
        "S3_group_relationship_summary": group_relationships,
        "S3_group_parental_summary": group_parental_summary,
        "S3_cultivated_candidates": cultivated_candidates,
        "S3_introgressed_reference": introgressed_reference,
        "S3_wild_candidates": wild_candidates,
        "S3_asymmetry_candidates": asymmetry_candidates,
        "S3_shared_hotspot_evidence": overlap_df,
        "S3_same_contig_nonoverlap": same_contig_nonoverlap_df,
        "S3_subgenome_specific_hotspots": subgenome_specific_df,
        "S3_contrast_consensus": contrast_consensus,
        "S3_sgC_top_markers": sgc_markers,
        "S3_sgE_top_markers": sge_markers,
        "S3_pca_scores_sgC": supplementary["S2_sgC_pca_scores"],
        "S3_pca_scores_sgE": supplementary["S2_sgE_pca_scores"],
        "S3_column_guide": build_column_guide(),
    }

    for sheet_name, table in manuscript_tables.items():
        manuscript_sheet_map[sheet_name] = table

    print(f"Writing workbook: {output_path}")
    write_workbook(output_path, manuscript_sheet_map)

    print(f"Writing workbook: {supplementary_output_path}")
    write_workbook(supplementary_output_path, supplementary_sheet_map)

    print("Part 3 evidence analysis finished successfully.")


if __name__ == "__main__":
    main()
