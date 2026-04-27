#!/usr/bin/env python3
"""
Rank candidate accessions, shared hotspots, and subgenome-specific intervals.

The script integrates accession-priority metrics and hotspot-support tables to
produce compact outputs used for downstream interpretation and annotation.
"""

import argparse
import json
import math
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np
import pandas as pd


REQUIRED_ANALYSIS_SHEETS = [
    "03_panel_overview",
    "04_group_relationship_summary",
    "05_group_parental_summary",
    "06_accession_evidence_all",
    "07_cultivated_priority_candidat",
    "08_introgressed_reference_acces",
    "09_wild_transition_candidates",
    "10_asymmetry_priority_candidate",
    "11_shared_hotspot_evidence",
    "12_same_contig_nonoverlap_evide",
    "13_subgenome_specific_hotspots",
    "14_contrast_consensus",
]

REQUIRED_SUPPLEMENTARY_SHEETS = [
    "S3_common_panel",
    "S3_accession_evidence_all",
    "S3_shared_hotspot_evidence",
    "S3_same_contig_nonoverlap",
    "S3_subgenome_specific_hotspots",
    "S3_contrast_consensus",
    "S3_sgC_top_markers",
    "S3_sgE_top_markers",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Create final manuscript-facing evidence tables for the coffee introgression analysis."
    )
    parser.add_argument("--analysis-workbook", required=True, help="Part 3 manuscript workbook (.xlsx)")
    parser.add_argument("--supplementary-workbook", required=True, help="Part 3 supplementary workbook (.xlsx)")
    parser.add_argument("--output-workbook", required=True, help="Final manuscript-facing workbook (.xlsx)")
    parser.add_argument(
        "--supplementary-output-workbook",
        required=True,
        help="Final supplementary workbook (.xlsx)",
    )
    parser.add_argument("--cultivated-top-n", type=int, default=8, help="Number of cultivated candidates to retain")
    parser.add_argument("--introgressed-top-n", type=int, default=6, help="Number of introgressed candidates to retain")
    parser.add_argument("--wild-top-n", type=int, default=8, help="Number of wild candidates to retain")
    parser.add_argument("--asymmetry-top-n", type=int, default=10, help="Number of asymmetry candidates to retain")
    parser.add_argument("--shared-hotspots-top-n", type=int, default=10, help="Number of shared hotspots to retain")
    parser.add_argument("--same-contig-top-n", type=int, default=8, help="Number of same-contig non-overlap events to retain")
    parser.add_argument("--subgenome-specific-top-n", type=int, default=15, help="Number of subgenome-specific hotspots to retain")
    parser.add_argument("--marker-support-top-n", type=int, default=10, help="Number of marker support rows to retain per selected hotspot")
    return parser.parse_args()


def validate_workbook(path: Path, required_sheets: List[str]) -> None:
    workbook = pd.ExcelFile(path)
    missing = [sheet for sheet in required_sheets if sheet not in workbook.sheet_names]
    if missing:
        raise ValueError(f"Workbook is missing required sheets: {missing}")


def load_sheets(path: Path, sheet_names: List[str]) -> Dict[str, pd.DataFrame]:
    return {sheet: pd.read_excel(path, sheet_name=sheet) for sheet in sheet_names}


def to_numeric_if_possible(df: pd.DataFrame) -> pd.DataFrame:
    result = df.copy()
    for col in result.columns:
        try:
            result[col] = pd.to_numeric(result[col])
        except Exception:
            pass
    return result


def safe_percentile_rank(series: pd.Series, ascending: bool = True) -> pd.Series:
    numeric = pd.to_numeric(series, errors="coerce")
    if numeric.notna().sum() == 0:
        return pd.Series(np.nan, index=series.index)
    if ascending:
        return numeric.rank(pct=True, method="average")
    return numeric.rank(pct=True, method="average", ascending=False)


def normalize_contrast_label(label: str) -> str:
    return str(label).replace("__vs__", " vs ").replace("_", " ")


def sheet_name_safe(name: str) -> str:
    cleaned = "".join(ch for ch in str(name) if ch not in r'[]:*?/\\')
    return cleaned[:31]


def build_input_inventory(args: argparse.Namespace) -> pd.DataFrame:
    rows = [
        {"item": "analysis_workbook", "path": str(Path(args.analysis_workbook).resolve())},
        {"item": "supplementary_workbook", "path": str(Path(args.supplementary_workbook).resolve())},
        {"item": "output_workbook", "path": str(Path(args.output_workbook).resolve())},
        {"item": "supplementary_output_workbook", "path": str(Path(args.supplementary_output_workbook).resolve())},
    ]
    return pd.DataFrame(rows)


def build_run_parameters(args: argparse.Namespace) -> pd.DataFrame:
    rows = []
    for key, value in vars(args).items():
        rows.append({"parameter": key, "value": value})
    return pd.DataFrame(rows)


def top_n(df: pd.DataFrame, n: int, sort_columns: List[Tuple[str, bool]]) -> pd.DataFrame:
    out = df.copy()
    by = [col for col, _ in sort_columns]
    ascending = [asc for _, asc in sort_columns]
    for col in by:
        if col in out.columns:
            try:
                out[col] = pd.to_numeric(out[col])
            except Exception:
                pass
    return out.sort_values(by=by, ascending=ascending).head(n).reset_index(drop=True)


def build_narrative_scoreboard(
    contrast_consensus: pd.DataFrame,
    shared_hotspots: pd.DataFrame,
    same_contig: pd.DataFrame,
    subgenome_specific: pd.DataFrame,
    cultivated: pd.DataFrame,
    introgressed: pd.DataFrame,
    wild: pd.DataFrame,
    asymmetry: pd.DataFrame,
) -> pd.DataFrame:
    rows = []

    contrast_df = contrast_consensus.copy()
    contrast_df["high_effect_markers_total"] = (
        pd.to_numeric(contrast_df["n_markers_delta_ge_0_90_sgC"], errors="coerce").fillna(0)
        + pd.to_numeric(contrast_df["n_markers_delta_ge_0_90_sgE"], errors="coerce").fillna(0)
    )
    contrast_df["overlap_hotspots_total"] = pd.to_numeric(contrast_df["n_overlapping_hotspot_pairs"], errors="coerce").fillna(0)
    contrast_df["effect_size_score"] = safe_percentile_rank(contrast_df["mean_abs_delta_af_average"])
    contrast_df["extreme_marker_score"] = safe_percentile_rank(contrast_df["high_effect_markers_total"])
    contrast_df["overlap_score"] = safe_percentile_rank(contrast_df["overlap_hotspots_total"])
    contrast_df["shared_hotspot_score"] = pd.to_numeric(contrast_df["mean_consensus_hotspot_percentile"], errors="coerce").fillna(0)
    contrast_df["base_direction_score"] = (
        pd.to_numeric(contrast_df["consensus_strength_percentile"], errors="coerce").fillna(0)
        + contrast_df["effect_size_score"].fillna(0)
        + contrast_df["extreme_marker_score"].fillna(0)
        + contrast_df["overlap_score"].fillna(0)
        + contrast_df["shared_hotspot_score"].fillna(0)
    ) / 5.0

    for _, row in contrast_df.iterrows():
        rows.append(
            {
                "direction_type": "contrast",
                "direction_label": normalize_contrast_label(row["contrast"]),
                "candidate_focus": "contrast-first manuscript",
                "support_metric_1": row["mean_abs_delta_af_average"],
                "support_metric_1_name": "mean_abs_delta_af_average",
                "support_metric_2": row["high_effect_markers_total"],
                "support_metric_2_name": "n_markers_delta_ge_0_90_total",
                "support_metric_3": row["n_overlapping_hotspot_pairs"],
                "support_metric_3_name": "n_overlapping_hotspot_pairs",
                "support_metric_4": row["mean_consensus_hotspot_percentile"],
                "support_metric_4_name": "mean_consensus_hotspot_percentile",
                "composite_support_score": row["base_direction_score"],
                "support_scope": "panel-wide contrast architecture",
                "source_sheet": "14_contrast_consensus",
                "source_key": row["contrast"],
            }
        )

    shared_agg = shared_hotspots.copy()
    if not shared_agg.empty:
        grouped = (
            shared_agg.groupby("contrast", dropna=False)
            .agg(
                n_shared_hotspots=("contig", "count"),
                max_consensus_hotspot_percentile=("consensus_hotspot_percentile", "max"),
                mean_consensus_hotspot_percentile=("consensus_hotspot_percentile", "mean"),
                max_combined_markers=("combined_n_markers_hotspot", "max"),
                mean_combined_effect=("combined_mean_abs_delta_af", "mean"),
            )
            .reset_index()
        )
        grouped["n_shared_hotspots_score"] = safe_percentile_rank(grouped["n_shared_hotspots"])
        grouped["max_hotspot_score"] = safe_percentile_rank(grouped["max_consensus_hotspot_percentile"])
        grouped["marker_support_score"] = safe_percentile_rank(grouped["max_combined_markers"])
        grouped["effect_score"] = safe_percentile_rank(grouped["mean_combined_effect"])
        grouped["composite"] = (
            grouped["n_shared_hotspots_score"].fillna(0)
            + grouped["max_hotspot_score"].fillna(0)
            + grouped["marker_support_score"].fillna(0)
            + grouped["effect_score"].fillna(0)
        ) / 4.0

        for _, row in grouped.iterrows():
            rows.append(
                {
                    "direction_type": "shared_hotspot",
                    "direction_label": f"Shared hotspot architecture: {normalize_contrast_label(row['contrast'])}",
                    "candidate_focus": "cross-subgenome hotspot architecture",
                    "support_metric_1": row["n_shared_hotspots"],
                    "support_metric_1_name": "n_shared_hotspots",
                    "support_metric_2": row["max_consensus_hotspot_percentile"],
                    "support_metric_2_name": "max_consensus_hotspot_percentile",
                    "support_metric_3": row["max_combined_markers"],
                    "support_metric_3_name": "max_combined_n_markers_hotspot",
                    "support_metric_4": row["mean_combined_effect"],
                    "support_metric_4_name": "mean_combined_mean_abs_delta_af",
                    "composite_support_score": row["composite"],
                    "support_scope": "shared hotspot architecture",
                    "source_sheet": "11_shared_hotspot_evidence",
                    "source_key": row["contrast"],
                }
            )

    asym = asymmetry.copy()
    if not asym.empty:
        top_asym = top_n(
            asym,
            min(10, len(asym)),
            [("subgenome_asymmetry_percentile", False), ("combined_evidence_percentile", False), ("overall_subgenome_asymmetry_rank", True)],
        )
        group_diversity = top_asym["analysis_group"].nunique()
        role_score = (
            top_asym["subgenome_asymmetry_percentile"].mean()
            + top_asym["combined_evidence_percentile"].mean()
            + (group_diversity / max(1, asym["analysis_group"].nunique()))
        ) / 3.0
        rows.append(
            {
                "direction_type": "accession_asymmetry",
                "direction_label": "Subgenome asymmetry across the accession panel",
                "candidate_focus": "accession-level asymmetry",
                "support_metric_1": top_asym["subgenome_asymmetry_percentile"].mean(),
                "support_metric_1_name": "mean_top_asymmetry_percentile",
                "support_metric_2": top_asym["combined_evidence_percentile"].mean(),
                "support_metric_2_name": "mean_top_combined_evidence_percentile",
                "support_metric_3": group_diversity,
                "support_metric_3_name": "n_analysis_groups_in_top_asymmetry_set",
                "support_metric_4": top_asym["overall_subgenome_asymmetry_rank"].min(),
                "support_metric_4_name": "best_overall_subgenome_asymmetry_rank",
                "composite_support_score": role_score,
                "support_scope": "accession-level asymmetric structure",
                "source_sheet": "10_asymmetry_priority_candidate",
                "source_key": "top_asymmetry_candidates",
            }
        )

    if not introgressed.empty:
        top_intro = top_n(
            introgressed,
            min(6, len(introgressed)),
            [("combined_evidence_percentile", False), ("overall_combined_evidence_rank", True)],
        )
        intro_score = (
            top_intro["combined_evidence_percentile"].mean()
            + safe_percentile_rank(top_intro["mean_introgressed_shift"]).mean()
            + safe_percentile_rank(top_intro["mean_canephora_affinity"]).mean()
        ) / 3.0
        rows.append(
            {
                "direction_type": "introgressed_reference",
                "direction_label": "Introgressed reference accessions as breeding anchors",
                "candidate_focus": "introgressed accession shortlist",
                "support_metric_1": top_intro["combined_evidence_percentile"].mean(),
                "support_metric_1_name": "mean_top_combined_evidence_percentile",
                "support_metric_2": top_intro["mean_introgressed_shift"].mean(),
                "support_metric_2_name": "mean_top_introgressed_shift",
                "support_metric_3": top_intro["mean_canephora_affinity"].mean(),
                "support_metric_3_name": "mean_top_canephora_affinity",
                "support_metric_4": top_intro["overall_combined_evidence_rank"].min(),
                "support_metric_4_name": "best_overall_combined_evidence_rank",
                "composite_support_score": intro_score,
                "support_scope": "breeding reference accessions",
                "source_sheet": "08_introgressed_reference_acces",
                "source_key": "top_introgressed_candidates",
            }
        )

    if not wild.empty:
        top_wild = top_n(
            wild,
            min(8, len(wild)),
            [("combined_evidence_percentile", False), ("overall_combined_evidence_rank", True)],
        )
        wild_score = (
            top_wild["combined_evidence_percentile"].mean()
            + safe_percentile_rank(top_wild["mean_wild_shift"]).mean()
            + safe_percentile_rank(top_wild["subgenome_asymmetry_percentile"]).mean()
        ) / 3.0
        rows.append(
            {
                "direction_type": "wild_transition",
                "direction_label": "Wild transition accessions as diversification anchors",
                "candidate_focus": "wild accession shortlist",
                "support_metric_1": top_wild["combined_evidence_percentile"].mean(),
                "support_metric_1_name": "mean_top_combined_evidence_percentile",
                "support_metric_2": top_wild["mean_wild_shift"].mean(),
                "support_metric_2_name": "mean_top_wild_shift",
                "support_metric_3": top_wild["subgenome_asymmetry_percentile"].mean(),
                "support_metric_3_name": "mean_top_subgenome_asymmetry_percentile",
                "support_metric_4": top_wild["overall_combined_evidence_rank"].min(),
                "support_metric_4_name": "best_overall_combined_evidence_rank",
                "composite_support_score": wild_score,
                "support_scope": "wild accession shortlist",
                "source_sheet": "09_wild_transition_candidates",
                "source_key": "top_wild_candidates",
            }
        )

    if not subgenome_specific.empty:
        top_specific = top_n(
            subgenome_specific,
            min(15, len(subgenome_specific)),
            [("subgenome_specific_support_percentile", False), ("mean_abs_delta_af", False), ("n_markers_hotspot", False)],
        )
        specific_score = (
            top_specific["subgenome_specific_support_percentile"].mean()
            + safe_percentile_rank(top_specific["mean_abs_delta_af"]).mean()
            + safe_percentile_rank(top_specific["n_markers_hotspot"]).mean()
        ) / 3.0
        rows.append(
            {
                "direction_type": "subgenome_specific_hotspots",
                "direction_label": "Subgenome-specific hotspot architecture",
                "candidate_focus": "subgenome-specific hotspot shortlist",
                "support_metric_1": top_specific["subgenome_specific_support_percentile"].mean(),
                "support_metric_1_name": "mean_top_subgenome_specific_support_percentile",
                "support_metric_2": top_specific["mean_abs_delta_af"].mean(),
                "support_metric_2_name": "mean_top_mean_abs_delta_af",
                "support_metric_3": top_specific["n_markers_hotspot"].sum(),
                "support_metric_3_name": "sum_top_n_markers_hotspot",
                "support_metric_4": top_specific["contrast"].nunique(),
                "support_metric_4_name": "n_contrasts_represented",
                "composite_support_score": specific_score,
                "support_scope": "subgenome-specific architecture",
                "source_sheet": "13_subgenome_specific_hotspots",
                "source_key": "top_subgenome_specific_hotspots",
            }
        )

    if not same_contig.empty:
        top_nonoverlap = top_n(
            same_contig,
            min(8, len(same_contig)),
            [("same_contig_support_percentile", False), ("combined_mean_abs_delta_af", False), ("combined_n_markers_hotspot", False)],
        )
        nonoverlap_score = (
            top_nonoverlap["same_contig_support_percentile"].mean()
            + safe_percentile_rank(top_nonoverlap["combined_mean_abs_delta_af"]).mean()
            + safe_percentile_rank(top_nonoverlap["combined_n_markers_hotspot"]).mean()
        ) / 3.0
        rows.append(
            {
                "direction_type": "same_contig_nonoverlap",
                "direction_label": "Same-contig but non-overlapping hotspot architecture",
                "candidate_focus": "offset hotspot support",
                "support_metric_1": top_nonoverlap["same_contig_support_percentile"].mean(),
                "support_metric_1_name": "mean_top_same_contig_support_percentile",
                "support_metric_2": top_nonoverlap["combined_mean_abs_delta_af"].mean(),
                "support_metric_2_name": "mean_top_combined_mean_abs_delta_af",
                "support_metric_3": top_nonoverlap["combined_n_markers_hotspot"].sum(),
                "support_metric_3_name": "sum_top_combined_n_markers_hotspot",
                "support_metric_4": top_nonoverlap["contrast"].nunique(),
                "support_metric_4_name": "n_contrasts_represented",
                "composite_support_score": nonoverlap_score,
                "support_scope": "same-contig offset hotspot architecture",
                "source_sheet": "12_same_contig_nonoverlap_evide",
                "source_key": "top_same_contig_nonoverlap",
            }
        )

    scoreboard = pd.DataFrame(rows)
    scoreboard["composite_support_score"] = pd.to_numeric(scoreboard["composite_support_score"], errors="coerce")
    scoreboard["recommended_rank"] = scoreboard["composite_support_score"].rank(ascending=False, method="dense").astype(int)
    scoreboard = scoreboard.sort_values(["recommended_rank", "direction_type", "direction_label"]).reset_index(drop=True)
    scoreboard["is_primary_recommendation"] = scoreboard["recommended_rank"] == 1
    return scoreboard


def integrate_accession_shortlist(
    all_accessions: pd.DataFrame,
    cultivated: pd.DataFrame,
    introgressed: pd.DataFrame,
    wild: pd.DataFrame,
    asymmetry: pd.DataFrame,
    top_counts: Dict[str, int],
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    cultivated_top = top_n(cultivated, top_counts["cultivated"], [("combined_evidence_percentile", False), ("overall_combined_evidence_rank", True)])
    introgressed_top = top_n(introgressed, top_counts["introgressed"], [("combined_evidence_percentile", False), ("overall_combined_evidence_rank", True)])
    wild_top = top_n(wild, top_counts["wild"], [("combined_evidence_percentile", False), ("overall_combined_evidence_rank", True)])
    asymmetry_top = top_n(asymmetry, top_counts["asymmetry"], [("subgenome_asymmetry_percentile", False), ("combined_evidence_percentile", False), ("overall_subgenome_asymmetry_rank", True)])

    role_frames = {
        "cultivated_priority": cultivated_top[["seq_id", "accession_name"]].copy(),
        "introgressed_reference": introgressed_top[["seq_id", "accession_name"]].copy(),
        "wild_transition": wild_top[["seq_id", "accession_name"]].copy(),
        "asymmetry_priority": asymmetry_top[["seq_id", "accession_name"]].copy(),
    }
    for role_name, frame in role_frames.items():
        frame["role_label"] = role_name

    role_long = pd.concat(role_frames.values(), ignore_index=True)
    role_long["role_flag"] = True

    accession_role_summary = (
        role_long.groupby(["seq_id", "accession_name"], dropna=False)
        .agg(
            role_count=("role_label", "nunique"),
            role_labels=("role_label", lambda x: "; ".join(sorted(set(map(str, x))))),
        )
        .reset_index()
    )

    for role_name, frame in role_frames.items():
        accession_role_summary[role_name] = accession_role_summary["seq_id"].isin(frame["seq_id"]).astype(int)

    rank_sources = {
        "cultivated_priority": cultivated_top[["seq_id", "overall_combined_evidence_rank"]].rename(columns={"overall_combined_evidence_rank": "cultivated_priority_rank"}),
        "introgressed_reference": introgressed_top[["seq_id", "overall_combined_evidence_rank"]].rename(columns={"overall_combined_evidence_rank": "introgressed_reference_rank"}),
        "wild_transition": wild_top[["seq_id", "overall_combined_evidence_rank"]].rename(columns={"overall_combined_evidence_rank": "wild_transition_rank"}),
        "asymmetry_priority": asymmetry_top[["seq_id", "overall_subgenome_asymmetry_rank"]].rename(columns={"overall_subgenome_asymmetry_rank": "asymmetry_priority_rank"}),
    }

    integrated = accession_role_summary.merge(all_accessions, on=["seq_id", "accession_name"], how="left")
    for frame in rank_sources.values():
        integrated = integrated.merge(frame, on="seq_id", how="left")

    rank_cols = [
        "cultivated_priority_rank",
        "introgressed_reference_rank",
        "wild_transition_rank",
        "asymmetry_priority_rank",
        "overall_combined_evidence_rank",
        "overall_introgression_affinity_rank",
        "overall_subgenome_asymmetry_rank",
    ]
    for col in rank_cols:
        if col not in integrated.columns:
            integrated[col] = np.nan

    integrated["best_role_rank"] = integrated[
        ["cultivated_priority_rank", "introgressed_reference_rank", "wild_transition_rank", "asymmetry_priority_rank"]
    ].min(axis=1, skipna=True)
    integrated["priority_tier"] = np.select(
        [
            integrated["role_count"] >= 3,
            (integrated["role_count"] >= 2) | (integrated["best_role_rank"] <= 3),
        ],
        ["core", "high"],
        default="supporting",
    )
    integrated["shortlist_sort_score"] = (
        integrated["role_count"].fillna(0) * 10
        + pd.to_numeric(integrated["combined_evidence_percentile"], errors="coerce").fillna(0) * 5
        + pd.to_numeric(integrated["subgenome_asymmetry_percentile"], errors="coerce").fillna(0) * 2
        + pd.to_numeric(integrated["introgression_affinity_percentile"], errors="coerce").fillna(0) * 3
    )

    integrated = integrated.sort_values(
        ["role_count", "priority_tier", "shortlist_sort_score", "overall_combined_evidence_rank", "overall_subgenome_asymmetry_rank"],
        ascending=[False, True, False, True, True],
    ).reset_index(drop=True)

    display_cols = [
        "seq_id",
        "accession_name",
        "analysis_group",
        "species_name",
        "variety",
        "role_count",
        "role_labels",
        "priority_tier",
        "best_role_rank",
        "overall_combined_evidence_rank",
        "overall_introgression_affinity_rank",
        "overall_subgenome_asymmetry_rank",
        "combined_evidence_percentile",
        "introgression_affinity_percentile",
        "subgenome_asymmetry_percentile",
        "mean_introgressed_shift",
        "mean_wild_shift",
        "mean_canephora_affinity",
        "mean_eugenioides_affinity",
        "parental_distortion_score",
        "cultivated_priority_rank",
        "introgressed_reference_rank",
        "wild_transition_rank",
        "asymmetry_priority_rank",
        "cultivated_priority",
        "introgressed_reference",
        "wild_transition",
        "asymmetry_priority",
        "country_of_origin",
        "donor_institute",
        "notes",
    ]
    display_cols = [col for col in display_cols if col in integrated.columns]
    integrated_display = integrated[display_cols].copy()
    return integrated_display, integrated


def integrate_hotspot_shortlists(
    shared_hotspots: pd.DataFrame,
    same_contig: pd.DataFrame,
    subgenome_specific: pd.DataFrame,
    top_counts: Dict[str, int],
) -> Tuple[pd.DataFrame, Dict[str, pd.DataFrame]]:
    shared_top = top_n(
        shared_hotspots,
        top_counts["shared"],
        [("consensus_hotspot_rank", True), ("consensus_hotspot_percentile", False), ("combined_n_markers_hotspot", False)],
    ).copy()
    shared_top["hotspot_category"] = "shared"
    shared_top["priority_score"] = pd.to_numeric(shared_top["consensus_hotspot_percentile"], errors="coerce")

    same_top = top_n(
        same_contig,
        top_counts["same_contig"],
        [("same_contig_support_percentile", False), ("combined_mean_abs_delta_af", False), ("combined_n_markers_hotspot", False)],
    ).copy()
    same_top["hotspot_category"] = "same_contig_nonoverlap"
    same_top["priority_score"] = pd.to_numeric(same_top["same_contig_support_percentile"], errors="coerce")
    if "consensus_hotspot_rank" not in same_top.columns:
        same_top["consensus_hotspot_rank"] = np.nan

    sub_top = top_n(
        subgenome_specific,
        top_counts["subgenome_specific"],
        [("subgenome_specific_support_percentile", False), ("mean_abs_delta_af", False), ("n_markers_hotspot", False)],
    ).copy()
    sub_top["hotspot_category"] = "subgenome_specific"
    sub_top["priority_score"] = pd.to_numeric(sub_top["subgenome_specific_support_percentile"], errors="coerce")
    sub_top["combined_n_markers_hotspot"] = sub_top.get("n_markers_hotspot")
    sub_top["consensus_hotspot_rank"] = np.nan

    common_cols = [
        "hotspot_category",
        "contrast",
        "contig",
        "subgenome",
        "combined_n_markers_hotspot",
        "n_markers_hotspot",
        "combined_mean_abs_delta_af",
        "mean_abs_delta_af",
        "combined_max_abs_delta_af",
        "max_abs_delta_af",
        "union_min_pos",
        "union_max_pos",
        "union_span_bp",
        "min_pos",
        "max_pos",
        "span_bp",
        "jaccard_bp",
        "overlap_status",
        "consensus_hotspot_rank",
        "priority_score",
    ]

    def ensure_cols(df: pd.DataFrame) -> pd.DataFrame:
        out = df.copy()
        for col in common_cols:
            if col not in out.columns:
                out[col] = np.nan
        return out[common_cols]

    integrated = pd.concat(
        [ensure_cols(shared_top), ensure_cols(same_top), ensure_cols(sub_top)],
        ignore_index=True,
    )
    integrated["integrated_hotspot_rank"] = integrated["priority_score"].rank(ascending=False, method="dense").astype(int)
    integrated = integrated.sort_values(
        ["priority_score", "combined_n_markers_hotspot", "n_markers_hotspot", "contrast"],
        ascending=[False, False, False, True],
    ).reset_index(drop=True)

    return integrated, {
        "shared": shared_top,
        "same_contig": same_top,
        "subgenome_specific": sub_top,
    }




def extract_primary_contrast_from_scoreboard(scoreboard: pd.DataFrame) -> str:
    primary = scoreboard.sort_values(["recommended_rank", "composite_support_score"], ascending=[True, False]).iloc[0]
    source_key = str(primary.get("source_key", ""))
    if "__vs__" in source_key:
        return source_key
    contrast_rows = scoreboard[scoreboard["direction_type"] == "contrast"].copy()
    if contrast_rows.empty:
        return ""
    return str(contrast_rows.sort_values(["recommended_rank", "composite_support_score"], ascending=[True, False]).iloc[0]["source_key"])


def build_primary_accession_shortlist(
    primary_contrast: str,
    cultivated: pd.DataFrame,
    introgressed: pd.DataFrame,
    wild: pd.DataFrame,
    asymmetry: pd.DataFrame,
    cultivated_top_n: int,
    introgressed_top_n: int,
    wild_top_n: int,
    asymmetry_top_n: int,
) -> pd.DataFrame:
    def _sorted(df: pd.DataFrame, top_n_value: int) -> pd.DataFrame:
        if df.empty:
            return df.copy()
        return top_n(df, top_n_value, [("combined_evidence_percentile", False), ("overall_combined_evidence_rank", True)])

    if "arabica_cultivated__vs__arabica_introgressed" == primary_contrast:
        role_priority_map = {"introgressed_reference": 1, "cultivated_comparator": 2, "asymmetry_support": 3}
        frames = [
            _sorted(introgressed, introgressed_top_n).assign(primary_role="introgressed_reference"),
            _sorted(cultivated, cultivated_top_n).assign(primary_role="cultivated_comparator"),
            top_n(asymmetry, min(asymmetry_top_n, len(asymmetry)), [("subgenome_asymmetry_percentile", False), ("combined_evidence_percentile", False), ("overall_subgenome_asymmetry_rank", True)]).assign(primary_role="asymmetry_support"),
        ]
    elif "arabica_cultivated__vs__arabica_wild" == primary_contrast:
        role_priority_map = {"wild_reference": 1, "cultivated_comparator": 2, "asymmetry_support": 3}
        frames = [
            _sorted(wild, wild_top_n).assign(primary_role="wild_reference"),
            _sorted(cultivated, cultivated_top_n).assign(primary_role="cultivated_comparator"),
            top_n(asymmetry, min(asymmetry_top_n, len(asymmetry)), [("subgenome_asymmetry_percentile", False), ("combined_evidence_percentile", False), ("overall_subgenome_asymmetry_rank", True)]).assign(primary_role="asymmetry_support"),
        ]
    else:
        role_priority_map = {"introgressed_reference": 1, "wild_reference": 1, "cultivated_comparator": 2, "asymmetry_support": 3}
        frames = [
            _sorted(introgressed, introgressed_top_n).assign(primary_role="introgressed_reference"),
            _sorted(wild, wild_top_n).assign(primary_role="wild_reference"),
            _sorted(cultivated, cultivated_top_n).assign(primary_role="cultivated_comparator"),
            top_n(asymmetry, min(asymmetry_top_n, len(asymmetry)), [("subgenome_asymmetry_percentile", False), ("combined_evidence_percentile", False), ("overall_subgenome_asymmetry_rank", True)]).assign(primary_role="asymmetry_support"),
        ]

    merged = pd.concat(frames, ignore_index=True).drop_duplicates(subset=["seq_id", "accession_name", "primary_role"])
    merged["primary_role_priority"] = merged["primary_role"].map(role_priority_map).fillna(99)
    role_summary = (
        merged.groupby(["seq_id", "accession_name"], dropna=False)
        .agg(
            primary_role_labels=("primary_role", lambda x: "; ".join(sorted(set(map(str, x))))),
            primary_role_count=("primary_role", "nunique"),
            best_primary_role_priority=("primary_role_priority", "min"),
        )
        .reset_index()
    )
    best_rows = (
        merged.sort_values(
            ["combined_evidence_percentile", "subgenome_asymmetry_percentile", "overall_combined_evidence_rank"],
            ascending=[False, False, True],
        )
        .drop_duplicates(subset=["seq_id", "accession_name"])
    )
    out = role_summary.merge(best_rows, on=["seq_id", "accession_name"], how="left")
    out = out.sort_values(
        ["best_primary_role_priority", "primary_role_count", "combined_evidence_percentile", "subgenome_asymmetry_percentile", "overall_combined_evidence_rank"],
        ascending=[True, False, False, False, True],
    ).reset_index(drop=True)
    preferred_cols = [
        "seq_id", "accession_name", "analysis_group", "variety", "primary_role_labels", "primary_role_count", "best_primary_role_priority",
        "combined_evidence_percentile", "overall_combined_evidence_rank", "introgression_affinity_percentile",
        "overall_introgression_affinity_rank", "subgenome_asymmetry_percentile", "overall_subgenome_asymmetry_rank",
        "mean_introgressed_shift", "mean_wild_shift", "mean_canephora_affinity", "mean_eugenioides_affinity",
        "country_of_origin", "donor_institute", "notes"
    ]
    preferred_cols = [c for c in preferred_cols if c in out.columns]
    return out[preferred_cols].copy()

def build_primary_direction_summary(scoreboard: pd.DataFrame, shared_hotspots: pd.DataFrame, accession_shortlist: pd.DataFrame) -> pd.DataFrame:
    primary = scoreboard.sort_values(["recommended_rank", "composite_support_score"], ascending=[True, False]).iloc[0]
    rows = [
        {
            "field": "recommended_direction_type",
            "value": primary["direction_type"],
            "note": "Highest composite support score in the narrative scoreboard.",
        },
        {
            "field": "recommended_direction_label",
            "value": primary["direction_label"],
            "note": "Recommended primary manuscript direction.",
        },
        {
            "field": "support_scope",
            "value": primary["support_scope"],
            "note": "Panel scope represented by the primary direction.",
        },
        {
            "field": "composite_support_score",
            "value": primary["composite_support_score"],
            "note": "Average of direction-specific support components.",
        },
        {
            "field": "top_accession_shortlist_size",
            "value": len(accession_shortlist),
            "note": "Union of top role-specific accession candidates.",
        },
        {
            "field": "top_shared_hotspot_count",
            "value": len(shared_hotspots),
            "note": "Number of retained shared hotspot rows.",
        },
    ]
    if primary["direction_type"] in {"contrast", "shared_hotspot"}:
        rows.append(
            {
                "field": "primary_source_key",
                "value": primary["source_key"],
                "note": "Primary contrast or contrast-linked source key.",
            }
        )
    return pd.DataFrame(rows)


def build_data_driven_statements(
    scoreboard: pd.DataFrame,
    contrast_consensus: pd.DataFrame,
    shared_hotspots: pd.DataFrame,
    accession_shortlist: pd.DataFrame,
    cultivated_top: pd.DataFrame,
    introgressed_top: pd.DataFrame,
    wild_top: pd.DataFrame,
    asymmetry_top: pd.DataFrame,
) -> pd.DataFrame:
    statements = []
    scoreboard = scoreboard.sort_values(["recommended_rank", "composite_support_score"], ascending=[True, False]).reset_index(drop=True)
    primary = scoreboard.iloc[0]
    top_contrast = contrast_consensus.sort_values(["consensus_strength_rank", "mean_abs_delta_af_average"], ascending=[True, False]).iloc[0]
    top_hotspot = shared_hotspots.sort_values(["consensus_hotspot_rank", "combined_n_markers_hotspot"], ascending=[True, False]).iloc[0]
    top_accession = accession_shortlist.sort_values(["role_count", "combined_evidence_percentile", "subgenome_asymmetry_percentile"], ascending=[False, False, False]).iloc[0]
    top_asymmetry = asymmetry_top.sort_values(["subgenome_asymmetry_percentile", "combined_evidence_percentile"], ascending=[False, False]).iloc[0]
    top_introgressed = introgressed_top.sort_values(["combined_evidence_percentile", "overall_combined_evidence_rank"], ascending=[False, True]).iloc[0]
    top_wild = wild_top.sort_values(["combined_evidence_percentile", "overall_combined_evidence_rank"], ascending=[False, True]).iloc[0]

    statements.append(
        {
            "statement_id": "S1",
            "statement_text": f"The highest-supported manuscript direction is '{primary['direction_label']}' with a composite support score of {primary['composite_support_score']:.3f}.",
            "source_sheet": primary["source_sheet"],
            "source_key": primary["source_key"],
        }
    )
    statements.append(
        {
            "statement_id": "S2",
            "statement_text": f"The strongest contrast in the panel is {normalize_contrast_label(top_contrast['contrast'])}, with mean_abs_delta_af_average = {top_contrast['mean_abs_delta_af_average']:.3f} and consensus strength rank = {int(top_contrast['consensus_strength_rank'])}.",
            "source_sheet": "14_contrast_consensus",
            "source_key": top_contrast["contrast"],
        }
    )
    statements.append(
        {
            "statement_id": "S3",
            "statement_text": f"The top shared hotspot is on {top_hotspot['contig']} for {normalize_contrast_label(top_hotspot['contrast'])}, supported by {int(top_hotspot['combined_n_markers_hotspot'])} combined markers and combined_mean_abs_delta_af = {top_hotspot['combined_mean_abs_delta_af']:.3f}.",
            "source_sheet": "11_shared_hotspot_evidence",
            "source_key": f"{top_hotspot['contrast']}::{top_hotspot['contig']}",
        }
    )
    statements.append(
        {
            "statement_id": "S4",
            "statement_text": f"The accession shortlist is led by {top_accession['accession_name']} ({top_accession['analysis_group']}), which appears in {int(top_accession['role_count'])} role-defined priority lists.",
            "source_sheet": "final_accession_shortlist",
            "source_key": top_accession["seq_id"],
        }
    )
    statements.append(
        {
            "statement_id": "S5",
            "statement_text": f"The highest-ranked introgressed reference accession is {top_introgressed['accession_name']}, while the top wild transition accession is {top_wild['accession_name']}.",
            "source_sheet": "08_introgressed_reference_acces;09_wild_transition_candidates",
            "source_key": f"{top_introgressed['seq_id']};{top_wild['seq_id']}",
        }
    )
    statements.append(
        {
            "statement_id": "S6",
            "statement_text": f"The strongest accession-level asymmetry signal is observed in {top_asymmetry['accession_name']} ({top_asymmetry['analysis_group']}), with subgenome_asymmetry_percentile = {top_asymmetry['subgenome_asymmetry_percentile']:.3f}.",
            "source_sheet": "10_asymmetry_priority_candidate",
            "source_key": top_asymmetry["seq_id"],
        }
    )
    return pd.DataFrame(statements)


def build_hotspot_marker_support(
    selected_shared_hotspots: pd.DataFrame,
    sgc_markers: pd.DataFrame,
    sge_markers: pd.DataFrame,
    marker_support_top_n: int,
) -> pd.DataFrame:
    frames = []
    if selected_shared_hotspots.empty:
        return pd.DataFrame()

    def subset_markers(marker_df: pd.DataFrame, contrast: str, contig: str, min_pos: float, max_pos: float, category: str) -> pd.DataFrame:
        df = marker_df.copy()
        position_col = "position" if "position" in df.columns else "pos"
        df[position_col] = pd.to_numeric(df[position_col], errors="coerce")
        mask = (
            df["contrast"].astype(str).eq(str(contrast))
            & df["contig"].astype(str).eq(str(contig))
            & df[position_col].between(min_pos, max_pos, inclusive="both")
        )
        df = df.loc[mask].copy()
        if df.empty:
            return df
        df["position"] = df[position_col]
        df["support_category"] = category
        return df.sort_values(["abs_delta_af", "position"], ascending=[False, True]).head(marker_support_top_n)

    for _, row in selected_shared_hotspots.iterrows():
        contrast = row["contrast"]
        contig = row["contig"]
        sgc_min = row.get("sgC_min_pos", np.nan)
        sgc_max = row.get("sgC_max_pos", np.nan)
        sge_min = row.get("sgE_min_pos", np.nan)
        sge_max = row.get("sgE_max_pos", np.nan)
        if pd.notna(sgc_min) and pd.notna(sgc_max):
            out = subset_markers(sgc_markers, contrast, contig, float(sgc_min), float(sgc_max), "sgC")
            if not out.empty:
                out["hotspot_key"] = f"{contrast}::{contig}"
                frames.append(out)
        if pd.notna(sge_min) and pd.notna(sge_max):
            out = subset_markers(sge_markers, contrast, contig, float(sge_min), float(sge_max), "sgE")
            if not out.empty:
                out["hotspot_key"] = f"{contrast}::{contig}"
                frames.append(out)
    if not frames:
        return pd.DataFrame()
    return pd.concat(frames, ignore_index=True)


def build_column_guide() -> pd.DataFrame:
    rows = [
        {"column_or_sheet": "03_narrative_scoreboard", "description": "Ranked manuscript directions derived from contrast, hotspot, and accession evidence."},
        {"column_or_sheet": "04_primary_direction_summary", "description": "Compact summary of the highest-supported manuscript direction."},
        {"column_or_sheet": "05_data_driven_statements", "description": "Direct manuscript-facing statements with source-traceable support."},
        {"column_or_sheet": "06_final_accession_shortlist", "description": "Union of role-based accession lists with integrated role counts and ranks."},
        {"column_or_sheet": "07_accession_role_matrix", "description": "Role flags and ranking metadata for the integrated accession shortlist."},
        {"column_or_sheet": "08_final_hotspot_shortlist", "description": "Integrated hotspot shortlist across shared, same-contig, and subgenome-specific categories."},
        {"column_or_sheet": "09_final_shared_hotspots", "description": "Top shared hotspots retained for manuscript interpretation."},
        {"column_or_sheet": "10_final_subgenome_specific", "description": "Top subgenome-specific hotspots retained for manuscript interpretation."},
        {"column_or_sheet": "11_final_same_contig_support", "description": "Top same-contig non-overlap hotspot events retained for manuscript interpretation."},
        {"column_or_sheet": "12_marker_support_for_shared_hot", "description": "Top marker rows supporting retained shared hotspots."},
        {"column_or_sheet": "13_primary_contrast_support", "description": "The highest-ranked contrast row retained from the contrast consensus table."},
        {"column_or_sheet": "14_alternative_direction_summary", "description": "Competing manuscript directions ranked below the primary recommendation."},
        {"column_or_sheet": "Table_*", "description": "Manuscript-facing tables derived from the detailed evidence sheets."},
    ]
    return pd.DataFrame(rows)


def write_workbook(path: Path, sheets: Dict[str, pd.DataFrame]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with pd.ExcelWriter(path, engine="openpyxl") as writer:
        for sheet_name, df in sheets.items():
            name = sheet_name_safe(sheet_name)
            df.to_excel(writer, sheet_name=name, index=False)


def main() -> None:
    args = parse_args()

    analysis_path = Path(args.analysis_workbook)
    supplementary_path = Path(args.supplementary_workbook)
    output_path = Path(args.output_workbook)
    supplementary_output_path = Path(args.supplementary_output_workbook)

    print("Validating workbook structure")
    validate_workbook(analysis_path, REQUIRED_ANALYSIS_SHEETS)
    validate_workbook(supplementary_path, REQUIRED_SUPPLEMENTARY_SHEETS)

    print("Loading Part 3 workbooks")
    analysis = load_sheets(analysis_path, REQUIRED_ANALYSIS_SHEETS)
    supplementary = load_sheets(supplementary_path, REQUIRED_SUPPLEMENTARY_SHEETS)

    panel_overview = to_numeric_if_possible(analysis["03_panel_overview"])
    group_relationships = to_numeric_if_possible(analysis["04_group_relationship_summary"])
    group_parental = to_numeric_if_possible(analysis["05_group_parental_summary"])
    accession_all = to_numeric_if_possible(analysis["06_accession_evidence_all"])
    cultivated = to_numeric_if_possible(analysis["07_cultivated_priority_candidat"])
    introgressed = to_numeric_if_possible(analysis["08_introgressed_reference_acces"])
    wild = to_numeric_if_possible(analysis["09_wild_transition_candidates"])
    asymmetry = to_numeric_if_possible(analysis["10_asymmetry_priority_candidate"])
    shared_hotspots = to_numeric_if_possible(analysis["11_shared_hotspot_evidence"])
    same_contig = to_numeric_if_possible(analysis["12_same_contig_nonoverlap_evide"])
    subgenome_specific = to_numeric_if_possible(analysis["13_subgenome_specific_hotspots"])
    contrast_consensus = to_numeric_if_possible(analysis["14_contrast_consensus"])

    sgc_markers = to_numeric_if_possible(supplementary["S3_sgC_top_markers"])
    sge_markers = to_numeric_if_possible(supplementary["S3_sgE_top_markers"])

    print("Scoring manuscript directions")
    scoreboard = build_narrative_scoreboard(
        contrast_consensus,
        shared_hotspots,
        same_contig,
        subgenome_specific,
        cultivated,
        introgressed,
        wild,
        asymmetry,
    )

    print("Integrating accession shortlist evidence")
    accession_display, accession_full = integrate_accession_shortlist(
        accession_all,
        cultivated,
        introgressed,
        wild,
        asymmetry,
        {
            "cultivated": args.cultivated_top_n,
            "introgressed": args.introgressed_top_n,
            "wild": args.wild_top_n,
            "asymmetry": args.asymmetry_top_n,
        },
    )

    print("Integrating hotspot shortlist evidence")
    hotspot_integrated, hotspot_parts = integrate_hotspot_shortlists(
        shared_hotspots,
        same_contig,
        subgenome_specific,
        {
            "shared": args.shared_hotspots_top_n,
            "same_contig": args.same_contig_top_n,
            "subgenome_specific": args.subgenome_specific_top_n,
        },
    )

    print("Generating manuscript-facing summaries")
    primary_contrast_key = extract_primary_contrast_from_scoreboard(scoreboard)
    primary_summary = build_primary_direction_summary(scoreboard, hotspot_parts["shared"], accession_display)
    statements = build_data_driven_statements(
        scoreboard,
        contrast_consensus,
        hotspot_parts["shared"],
        accession_display,
        top_n(cultivated, args.cultivated_top_n, [("combined_evidence_percentile", False), ("overall_combined_evidence_rank", True)]),
        top_n(introgressed, args.introgressed_top_n, [("combined_evidence_percentile", False), ("overall_combined_evidence_rank", True)]),
        top_n(wild, args.wild_top_n, [("combined_evidence_percentile", False), ("overall_combined_evidence_rank", True)]),
        top_n(asymmetry, args.asymmetry_top_n, [("subgenome_asymmetry_percentile", False), ("combined_evidence_percentile", False), ("overall_subgenome_asymmetry_rank", True)]),
    )
    marker_support = build_hotspot_marker_support(
        hotspot_parts["shared"],
        sgc_markers,
        sge_markers,
        args.marker_support_top_n,
    )


    primary_accession_shortlist = build_primary_accession_shortlist(
        primary_contrast_key,
        cultivated,
        introgressed,
        wild,
        asymmetry,
        args.cultivated_top_n,
        args.introgressed_top_n,
        args.wild_top_n,
        args.asymmetry_top_n,
    )

    primary_shared_hotspots = hotspot_parts["shared"].copy()
    if primary_contrast_key:
        primary_shared_hotspots = primary_shared_hotspots[
            primary_shared_hotspots["contrast"].astype(str).eq(primary_contrast_key)
        ].reset_index(drop=True)

    primary_subgenome_specific = hotspot_parts["subgenome_specific"].copy()
    if primary_contrast_key:
        primary_subgenome_specific = primary_subgenome_specific[
            primary_subgenome_specific["contrast"].astype(str).eq(primary_contrast_key)
        ].reset_index(drop=True)

    primary_same_contig = hotspot_parts["same_contig"].copy()
    if primary_contrast_key:
        primary_same_contig = primary_same_contig[
            primary_same_contig["contrast"].astype(str).eq(primary_contrast_key)
        ].reset_index(drop=True)


    primary_contrast = contrast_consensus.sort_values(["consensus_strength_rank", "mean_abs_delta_af_average"], ascending=[True, False]).head(1).reset_index(drop=True)
    alternative_directions = scoreboard[scoreboard["recommended_rank"] > 1].copy().reset_index(drop=True)

    manuscript_sheets = {
        "01_run_parameters": build_run_parameters(args),
        "02_input_inventory": build_input_inventory(args),
        "03_narrative_scoreboard": scoreboard,
        "04_primary_direction_summary": primary_summary,
        "05_data_driven_statements": statements,
        "06_final_accession_shortlist": accession_display,
        "06A_primary_accession_shortlist": primary_accession_shortlist,
        "07_accession_role_matrix": accession_full,
        "08_final_hotspot_shortlist": hotspot_integrated,
        "09_final_shared_hotspots": hotspot_parts["shared"],
        "09A_primary_shared_hotspots": primary_shared_hotspots,
        "10_final_subgenome_specific": hotspot_parts["subgenome_specific"],
        "10A_primary_subgenome_specific": primary_subgenome_specific,
        "11_final_same_contig_support": hotspot_parts["same_contig"],
        "11A_primary_same_contig_support": primary_same_contig,
        "12_marker_support_for_shared_hot": marker_support,
        "13_primary_contrast_support": primary_contrast,
        "14_alternative_direction_summary": alternative_directions,
        "15_panel_overview_reference": panel_overview,
        "16_group_relationship_reference": group_relationships,
        "17_group_parental_reference": group_parental,
        "18_column_guide": build_column_guide(),
        "Table_01_panel_overview": panel_overview,
        "Table_02_primary_direction": primary_summary,
        "Table_03_accession_shortlist": primary_accession_shortlist,
        "Table_04_shared_hotspots": primary_shared_hotspots,
        "Table_05_subgenome_specific": primary_subgenome_specific,
        "Table_06_direction_scoreboard": scoreboard,
        "Table_07_data_driven_statements": statements,
    }

    supplementary_sheets = {
        "S4_narrative_scoreboard": scoreboard,
        "S4_primary_direction_summary": primary_summary,
        "S4_data_driven_statements": statements,
        "S4_accession_shortlist_display": accession_display,
        "S4_primary_accession_shortlist": primary_accession_shortlist,
        "S4_accession_shortlist_full": accession_full,
        "S4_cultivated_candidates_top": top_n(cultivated, args.cultivated_top_n, [("combined_evidence_percentile", False), ("overall_combined_evidence_rank", True)]),
        "S4_introgressed_candidates_top": top_n(introgressed, args.introgressed_top_n, [("combined_evidence_percentile", False), ("overall_combined_evidence_rank", True)]),
        "S4_wild_candidates_top": top_n(wild, args.wild_top_n, [("combined_evidence_percentile", False), ("overall_combined_evidence_rank", True)]),
        "S4_asymmetry_candidates_top": top_n(asymmetry, args.asymmetry_top_n, [("subgenome_asymmetry_percentile", False), ("combined_evidence_percentile", False), ("overall_subgenome_asymmetry_rank", True)]),
        "S4_hotspot_shortlist_integrated": hotspot_integrated,
        "S4_shared_hotspots_top": hotspot_parts["shared"],
        "S4_primary_shared_hotspots": primary_shared_hotspots,
        "S4_same_contig_support_top": hotspot_parts["same_contig"],
        "S4_primary_same_contig": primary_same_contig,
        "S4_subgenome_specific_top": hotspot_parts["subgenome_specific"],
        "S4_primary_subgenome_specific": primary_subgenome_specific,
        "S4_marker_support_shared": marker_support,
        "S4_primary_contrast_support": primary_contrast,
        "S4_contrast_consensus_full": contrast_consensus,
        "S4_group_relationships_full": group_relationships,
        "S4_group_parental_full": group_parental,
        "S4_column_guide": build_column_guide(),
    }

    print(f"Writing workbook: {output_path}")
    write_workbook(output_path, manuscript_sheets)
    print(f"Writing workbook: {supplementary_output_path}")
    write_workbook(supplementary_output_path, supplementary_sheets)
    print("Part 4 final table analysis finished successfully.")


if __name__ == "__main__":
    main()
