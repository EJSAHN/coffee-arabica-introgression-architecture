#!/usr/bin/env python3
"""
Project prioritized hotspot intervals into ET-39 annotation space.

The script reads contig inventories, prioritized hotspot tables, and sgC/sgE
GFF3 annotations. It builds conservative contig-to-annotation mappings, transfers
hotspot intervals into the annotation coordinate system when possible, and
recovers overlapping gene models for downstream candidate-gene interpretation.
"""

from __future__ import annotations

import argparse
import math
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple

import pandas as pd


MAJOR_CONTIG_PATTERN = re.compile(r"^chr_([A-Z])_sg_([CE])_\(")


@dataclass
class MappingDecision:
    source_contig: str
    source_subgenome: str
    source_length: int
    target_seqid: Optional[str]
    target_length: Optional[int]
    target_chromosome_label: Optional[str]
    relative_length_difference: Optional[float]
    second_best_relative_difference: Optional[float]
    mapping_method: str
    coordinate_method: str
    mapping_confidence: str
    note: str


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Annotate hotspot intervals with gene models using GFF3 files."
    )
    parser.add_argument("--part1-workbook", required=True, help="Path to the Part 1 manuscript workbook.")
    parser.add_argument("--part4-workbook", required=True, help="Path to the Part 4 manuscript workbook.")
    parser.add_argument("--sgc-gff3", required=True, help="Path to the sgC GFF3 annotation file.")
    parser.add_argument("--sge-gff3", required=True, help="Path to the sgE GFF3 annotation file.")
    parser.add_argument("--output-workbook", required=True, help="Path to the manuscript-facing output workbook.")
    parser.add_argument(
        "--supplementary-output-workbook",
        required=True,
        help="Path to the supplementary output workbook.",
    )
    parser.add_argument(
        "--max-shared-hotspots",
        type=int,
        default=25,
        help="Maximum number of shared hotspot rows to annotate from the final hotspot table.",
    )
    parser.add_argument(
        "--max-subgenome-specific-hotspots",
        type=int,
        default=50,
        help="Maximum number of subgenome-specific hotspot rows to annotate.",
    )
    parser.add_argument(
        "--max-same-contig-hotspots",
        type=int,
        default=25,
        help="Maximum number of same-contig support rows to annotate.",
    )
    parser.add_argument(
        "--maximum-scaled-relative-length-difference",
        type=float,
        default=0.08,
        help="Maximum relative length difference allowed for scaled interval mapping.",
    )
    parser.add_argument(
        "--minimum-high-confidence-gap",
        type=float,
        default=0.01,
        help="Minimum gap between best and second-best relative length matches to call a mapping high confidence.",
    )
    parser.add_argument(
        "--keyword-priority-only",
        action="store_true",
        help="If set, primary candidate gene shortlist will only include genes with keyword hits.",
    )
    return parser.parse_args()


def validate_inputs(paths: Iterable[Path]) -> None:
    for path in paths:
        if not path.exists():
            raise FileNotFoundError(f"Missing required file: {path}")


def parse_attributes(attribute_text: str) -> Dict[str, str]:
    attributes: Dict[str, str] = {}
    if not attribute_text or attribute_text == ".":
        return attributes
    for item in attribute_text.split(";"):
        item = item.strip()
        if not item:
            continue
        if "=" in item:
            key, value = item.split("=", 1)
            attributes[key] = value
    return attributes


def parse_gff3(path: Path, file_subgenome: str) -> Tuple[pd.DataFrame, pd.DataFrame]:
    region_rows: List[Dict[str, object]] = []
    gene_rows: List[Dict[str, object]] = []

    with open(path, "r", encoding="utf-8", errors="ignore") as handle:
        for raw_line in handle:
            if not raw_line.strip() or raw_line.startswith("#"):
                continue
            parts = raw_line.rstrip("\n").split("\t")
            if len(parts) != 9:
                continue

            seqid, source, feature_type, start, end, score, strand, phase, attributes_text = parts
            attributes = parse_attributes(attributes_text)

            start_i = int(start)
            end_i = int(end)

            if feature_type == "region":
                chromosome_label = attributes.get("chromosome")
                region_rows.append(
                    {
                        "seqid": seqid,
                        "source": source,
                        "feature_type": feature_type,
                        "start": start_i,
                        "end": end_i,
                        "length": end_i - start_i + 1,
                        "subgenome": file_subgenome,
                        "chromosome_label": chromosome_label,
                        "region_name": attributes.get("Name"),
                        "attributes": attributes_text,
                    }
                )

            elif feature_type == "gene":
                description = attributes.get("description", "")
                gene_name = attributes.get("gene") or attributes.get("Name") or attributes.get("ID")
                gene_rows.append(
                    {
                        "seqid": seqid,
                        "source": source,
                        "feature_type": feature_type,
                        "start": start_i,
                        "end": end_i,
                        "strand": strand,
                        "gene_identifier": attributes.get("ID"),
                        "gene_id_numeric": attributes.get("Dbxref", ""),
                        "gene_name": gene_name,
                        "description": description,
                        "gene_biotype": attributes.get("gene_biotype"),
                        "attributes": attributes_text,
                        "subgenome": file_subgenome,
                    }
                )

    region_df = pd.DataFrame(region_rows)
    gene_df = pd.DataFrame(gene_rows)

    if region_df.empty:
        raise ValueError(f"No region rows were parsed from {path}")
    if gene_df.empty:
        raise ValueError(f"No gene rows were parsed from {path}")

    return region_df, gene_df


def load_part1_contigs(part1_workbook: Path) -> pd.DataFrame:
    sgc_df = pd.read_excel(part1_workbook, sheet_name="07_sgC_contig_header")
    sge_df = pd.read_excel(part1_workbook, sheet_name="08_sgE_contig_header")
    combined = pd.concat([sgc_df, sge_df], ignore_index=True).drop_duplicates()

    combined["vcf_subgenome"] = combined["contig_id"].str.extract(MAJOR_CONTIG_PATTERN)[1]
    combined["contig_letter"] = combined["contig_id"].str.extract(MAJOR_CONTIG_PATTERN)[0]
    combined["is_major_contig"] = combined["contig_id"].str.match(MAJOR_CONTIG_PATTERN, na=False)

    def infer_scaffold_token(contig_id: str) -> Optional[str]:
        match = re.search(r"(Scaffold_\d+;HRSCAF_\d+)", str(contig_id))
        return match.group(1) if match else None

    combined["scaffold_token"] = combined["contig_id"].apply(infer_scaffold_token)
    return combined


def classify_mapping_confidence(
    relative_difference: float,
    second_best_difference: Optional[float],
    maximum_scaled_relative_length_difference: float,
    minimum_high_confidence_gap: float,
    exact_match: bool,
) -> Tuple[str, str]:
    if exact_match:
        return "exact", "direct"

    if math.isnan(relative_difference):
        return "unresolved", "unresolved"

    gap = None
    if second_best_difference is not None and not math.isnan(second_best_difference):
        gap = second_best_difference - relative_difference

    if relative_difference <= maximum_scaled_relative_length_difference:
        if gap is not None and gap >= minimum_high_confidence_gap:
            if relative_difference <= 0.03:
                return "high", "length_scaled"
            return "medium", "length_scaled"
        if relative_difference <= 0.03:
            return "medium", "length_scaled"
        return "low", "length_scaled"

    return "unresolved", "unresolved"


def build_contig_crosswalk(
    vcf_contigs: pd.DataFrame,
    region_df: pd.DataFrame,
    subgenome: str,
    maximum_scaled_relative_length_difference: float,
    minimum_high_confidence_gap: float,
) -> pd.DataFrame:
    source_df = (
        vcf_contigs.loc[(vcf_contigs["is_major_contig"]) & (vcf_contigs["vcf_subgenome"] == subgenome), ["contig_id", "contig_length", "contig_letter", "scaffold_token"]]
        .drop_duplicates()
        .copy()
    )
    target_df = region_df.loc[region_df["subgenome"] == subgenome, ["seqid", "length", "chromosome_label"]].drop_duplicates().copy()

    records: List[Dict[str, object]] = []

    for _, source_row in source_df.iterrows():
        source_length = int(source_row["contig_length"])
        candidate_table = target_df.copy()
        candidate_table["relative_difference"] = (candidate_table["length"] - source_length).abs() / candidate_table["length"].clip(lower=1)
        candidate_table = candidate_table.sort_values(["relative_difference", "length", "seqid"]).reset_index(drop=True)

        if candidate_table.empty:
            records.append(
                MappingDecision(
                    source_contig=str(source_row["contig_id"]),
                    source_subgenome=subgenome,
                    source_length=source_length,
                    target_seqid=None,
                    target_length=None,
                    target_chromosome_label=None,
                    relative_length_difference=None,
                    second_best_relative_difference=None,
                    mapping_method="no_candidate_region",
                    coordinate_method="unresolved",
                    mapping_confidence="unresolved",
                    note="No target regions were available in the annotation file.",
                ).__dict__
            )
            continue

        best = candidate_table.iloc[0]
        second_best_difference = float(candidate_table.iloc[1]["relative_difference"]) if len(candidate_table) > 1 else math.nan
        exact_match = int(best["length"]) == source_length
        confidence, coordinate_method = classify_mapping_confidence(
            relative_difference=float(best["relative_difference"]),
            second_best_difference=second_best_difference,
            maximum_scaled_relative_length_difference=maximum_scaled_relative_length_difference,
            minimum_high_confidence_gap=minimum_high_confidence_gap,
            exact_match=exact_match,
        )

        mapping_method = "exact_length_match" if exact_match else "nearest_length_match"

        note = (
            "Exact length match between the VCF contig and the annotation region."
            if exact_match
            else "Approximate mapping based on the nearest annotation-region length. Coordinates require length scaling and should be interpreted cautiously."
        )

        records.append(
            MappingDecision(
                source_contig=str(source_row["contig_id"]),
                source_subgenome=subgenome,
                source_length=source_length,
                target_seqid=str(best["seqid"]) if pd.notna(best["seqid"]) else None,
                target_length=int(best["length"]) if pd.notna(best["length"]) else None,
                target_chromosome_label=str(best["chromosome_label"]) if pd.notna(best["chromosome_label"]) else None,
                relative_length_difference=float(best["relative_difference"]),
                second_best_relative_difference=second_best_difference,
                mapping_method=mapping_method,
                coordinate_method=coordinate_method,
                mapping_confidence=confidence,
                note=note,
            ).__dict__
        )

    crosswalk_df = pd.DataFrame(records)
    crosswalk_df["mapping_ready"] = crosswalk_df["coordinate_method"].isin(["direct", "length_scaled"])
    crosswalk_df["exact_length_match"] = crosswalk_df["mapping_method"].eq("exact_length_match")
    return crosswalk_df.sort_values(["mapping_ready", "mapping_confidence", "relative_length_difference"], ascending=[False, True, True])


def derive_paired_contig_name(source_contig: str, target_subgenome: str, contig_inventory: pd.DataFrame) -> Optional[str]:
    match = MAJOR_CONTIG_PATTERN.match(str(source_contig))
    if not match:
        return None
    letter = match.group(1)
    hits = contig_inventory.loc[
        (contig_inventory["is_major_contig"]) &
        (contig_inventory["contig_letter"] == letter) &
        (contig_inventory["vcf_subgenome"] == target_subgenome),
        "contig_id"
    ].drop_duplicates().tolist()
    if len(hits) == 1:
        return hits[0]
    return None


def build_crosswalk_lookup(crosswalk_df: pd.DataFrame) -> Dict[str, Dict[str, object]]:
    lookup: Dict[str, Dict[str, object]] = {}
    for _, row in crosswalk_df.iterrows():
        lookup[str(row["source_contig"])] = row.to_dict()
    return lookup


def transfer_interval(
    source_contig: str,
    source_start: int,
    source_end: int,
    crosswalk_lookup: Dict[str, Dict[str, object]],
) -> Dict[str, object]:
    row = crosswalk_lookup.get(source_contig)
    if row is None:
        return {
            "source_contig": source_contig,
            "target_seqid": None,
            "target_chromosome_label": None,
            "mapping_method": "missing_crosswalk",
            "coordinate_method": "unresolved",
            "mapping_confidence": "unresolved",
            "source_start": source_start,
            "source_end": source_end,
            "source_span_bp": source_end - source_start + 1,
            "target_start": None,
            "target_end": None,
            "target_span_bp": None,
            "length_scale_factor": None,
            "relative_length_difference": None,
            "annotation_ready": False,
            "mapping_note": "No contig crosswalk entry was available.",
        }

    source_length = int(row["source_length"])
    target_length = row["target_length"]
    coordinate_method = row["coordinate_method"]

    if coordinate_method == "direct":
        target_start = int(source_start)
        target_end = int(source_end)
        scale_factor = 1.0
        ready = True
    elif coordinate_method == "length_scaled" and pd.notna(target_length):
        target_start = max(1, int(round(source_start * float(target_length) / source_length)))
        target_end = max(target_start, int(round(source_end * float(target_length) / source_length)))
        scale_factor = float(target_length) / source_length
        ready = True
    else:
        target_start = None
        target_end = None
        scale_factor = None
        ready = False

    if ready and pd.notna(target_length):
        target_start = max(1, min(int(target_length), int(target_start)))
        target_end = max(1, min(int(target_length), int(target_end)))
        if target_end < target_start:
            target_end = target_start
        target_span = target_end - target_start + 1
    else:
        target_span = None

    return {
        "source_contig": source_contig,
        "target_seqid": row.get("target_seqid"),
        "target_chromosome_label": row.get("target_chromosome_label"),
        "mapping_method": row.get("mapping_method"),
        "coordinate_method": coordinate_method,
        "mapping_confidence": row.get("mapping_confidence"),
        "source_start": int(source_start),
        "source_end": int(source_end),
        "source_span_bp": int(source_end) - int(source_start) + 1,
        "target_start": target_start,
        "target_end": target_end,
        "target_span_bp": target_span,
        "length_scale_factor": scale_factor,
        "relative_length_difference": row.get("relative_length_difference"),
        "annotation_ready": ready,
        "mapping_note": row.get("note"),
    }


KEYWORD_LIBRARY = {
    "receptor_kinase": [
        "receptor kinase",
        "receptor-like kinase",
        "wall-associated receptor kinase",
        "serine/threonine-protein kinase",
        "protein kinase",
        "kinase-like",
    ],
    "resistance_like": [
        "disease resistance",
        "resistance protein",
        "nb-arc",
        "nbs-lrr",
        "leucine-rich repeat",
        "lrr",
        "rpp",
        "rpm",
        "rps",
    ],
    "defense_or_stress": [
        "defense",
        "stress",
        "pathogenesis-related",
        "peroxidase",
        "chitinase",
        "thaumatin",
        "wound",
        "immune",
    ],
    "transcription_factor": [
        "transcription factor",
        "wrky",
        "myb",
        "bhlh",
        "bzip",
        "nac",
        "ap2",
        "erf",
    ],
    "transporter": [
        "transporter",
        "abc transporter",
        "mfs transporter",
        "sugar transporter",
        "ion channel",
    ],
}


def classify_gene_keywords(gene_name: str, description: str) -> Tuple[str, str, int]:
    text = f"{gene_name} {description}".lower()
    matched_categories: List[str] = []
    matched_terms: List[str] = []

    for category, terms in KEYWORD_LIBRARY.items():
        category_hit = False
        for term in terms:
            if term in text:
                matched_terms.append(term)
                category_hit = True
        if category_hit:
            matched_categories.append(category)

    matched_categories = sorted(set(matched_categories))
    matched_terms = sorted(set(matched_terms))
    return "; ".join(matched_categories), "; ".join(matched_terms), len(matched_terms)


def intersect_genes(
    interval_table: pd.DataFrame,
    gene_df: pd.DataFrame,
    hotspot_category: str,
    hotspot_subgenome: str,
) -> pd.DataFrame:
    records: List[Dict[str, object]] = []

    for _, interval in interval_table.iterrows():
        if not bool(interval.get("annotation_ready", False)):
            continue

        seqid = interval["target_seqid"]
        start = int(interval["target_start"])
        end = int(interval["target_end"])
        overlap_genes = gene_df.loc[
            (gene_df["seqid"] == seqid) &
            (gene_df["start"] <= end) &
            (gene_df["end"] >= start)
        ].copy()

        for _, gene in overlap_genes.iterrows():
            overlap_start = max(start, int(gene["start"]))
            overlap_end = min(end, int(gene["end"]))
            overlap_bp = max(0, overlap_end - overlap_start + 1)
            category_hits, term_hits, hit_count = classify_gene_keywords(
                gene_name=str(gene.get("gene_name", "")),
                description=str(gene.get("description", "")),
            )

            row = {
                "hotspot_category": hotspot_category,
                "hotspot_subgenome": hotspot_subgenome,
                "hotspot_identifier": interval["hotspot_identifier"],
                "contrast": interval["contrast"],
                "source_contig": interval["source_contig"],
                "target_seqid": seqid,
                "target_chromosome_label": interval.get("target_chromosome_label"),
                "mapping_method": interval["mapping_method"],
                "coordinate_method": interval["coordinate_method"],
                "mapping_confidence": interval["mapping_confidence"],
                "source_start": interval["source_start"],
                "source_end": interval["source_end"],
                "target_start": interval["target_start"],
                "target_end": interval["target_end"],
                "gene_seqid": gene["seqid"],
                "gene_start": gene["start"],
                "gene_end": gene["end"],
                "gene_identifier": gene.get("gene_identifier"),
                "gene_name": gene.get("gene_name"),
                "description": gene.get("description"),
                "gene_biotype": gene.get("gene_biotype"),
                "strand": gene.get("strand"),
                "overlap_start": overlap_start,
                "overlap_end": overlap_end,
                "overlap_bp": overlap_bp,
                "overlap_fraction_of_interval": overlap_bp / max(1, int(interval["target_span_bp"])) if pd.notna(interval.get("target_span_bp")) else math.nan,
                "keyword_categories": category_hits,
                "keyword_terms": term_hits,
                "keyword_hit_count": hit_count,
            }
            records.append(row)

    df = pd.DataFrame(records)
    if not df.empty:
        df["is_protein_coding"] = df["gene_biotype"].eq("protein_coding")
        df["keyword_priority_score"] = (
            df["keyword_hit_count"].fillna(0).astype(float)
            + df["is_protein_coding"].astype(int) * 0.25
            + df["overlap_bp"].fillna(0).astype(float) / max(1.0, df["overlap_bp"].max())
        )
        df = df.sort_values(
            ["hotspot_identifier", "keyword_hit_count", "is_protein_coding", "overlap_bp", "gene_start"],
            ascending=[True, False, False, False, True],
        ).reset_index(drop=True)
    return df


def summarize_interval_table(interval_table: pd.DataFrame, gene_table: pd.DataFrame) -> pd.DataFrame:
    records: List[Dict[str, object]] = []
    if interval_table.empty:
        return pd.DataFrame(records)

    grouped_genes = gene_table.groupby("hotspot_identifier") if not gene_table.empty else None

    for _, interval in interval_table.iterrows():
        hotspot_identifier = interval["hotspot_identifier"]
        if grouped_genes is not None and hotspot_identifier in grouped_genes.groups:
            genes = grouped_genes.get_group(hotspot_identifier)
            n_genes = int(genes.shape[0])
            n_protein_coding = int(genes["is_protein_coding"].sum()) if "is_protein_coding" in genes.columns else 0
            n_keyword_hits = int((genes["keyword_hit_count"] > 0).sum()) if "keyword_hit_count" in genes.columns else 0
            unique_keyword_categories = "; ".join(
                sorted(
                    {
                        entry
                        for value in genes["keyword_categories"].fillna("")
                        for entry in str(value).split("; ")
                        if entry
                    }
                )
            )
        else:
            n_genes = 0
            n_protein_coding = 0
            n_keyword_hits = 0
            unique_keyword_categories = ""

        row = interval.to_dict()
        row["n_overlapping_genes"] = n_genes
        row["n_protein_coding_genes"] = n_protein_coding
        row["n_keyword_supported_genes"] = n_keyword_hits
        row["keyword_categories_detected"] = unique_keyword_categories
        records.append(row)

    return pd.DataFrame(records)


def build_primary_shared_interval_table(
    primary_shared: pd.DataFrame,
    contig_inventory: pd.DataFrame,
    sgc_lookup: Dict[str, Dict[str, object]],
    sge_lookup: Dict[str, Dict[str, object]],
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    sgc_records: List[Dict[str, object]] = []
    sge_records: List[Dict[str, object]] = []

    for _, row in primary_shared.iterrows():
        contrast = row["contrast"]
        sgc_contig = row["contig"]
        sgc_interval = transfer_interval(
            source_contig=str(sgc_contig),
            source_start=int(row["sgC_min_pos"]),
            source_end=int(row["sgC_max_pos"]),
            crosswalk_lookup=sgc_lookup,
        )
        sgc_interval["hotspot_identifier"] = f"{contrast}::{sgc_contig}::sgC"
        sgc_interval["contrast"] = contrast
        sgc_interval["priority_score"] = row.get("priority_score")
        sgc_records.append(sgc_interval)

        sge_contig = derive_paired_contig_name(str(sgc_contig), "E", contig_inventory)
        if sge_contig is None:
            sge_interval = {
                "source_contig": None,
                "target_seqid": None,
                "target_chromosome_label": None,
                "mapping_method": "paired_contig_not_found",
                "coordinate_method": "unresolved",
                "mapping_confidence": "unresolved",
                "source_start": int(row["sgE_min_pos"]),
                "source_end": int(row["sgE_max_pos"]),
                "source_span_bp": int(row["sgE_max_pos"]) - int(row["sgE_min_pos"]) + 1,
                "target_start": None,
                "target_end": None,
                "target_span_bp": None,
                "length_scale_factor": None,
                "relative_length_difference": None,
                "annotation_ready": False,
                "mapping_note": "The paired sgE contig could not be derived from the sgC hotspot contig.",
            }
        else:
            sge_interval = transfer_interval(
                source_contig=str(sge_contig),
                source_start=int(row["sgE_min_pos"]),
                source_end=int(row["sgE_max_pos"]),
                crosswalk_lookup=sge_lookup,
            )
        sge_interval["hotspot_identifier"] = f"{contrast}::{sgc_contig}::sgE"
        sge_interval["contrast"] = contrast
        sge_interval["paired_source_contig"] = sge_contig
        sge_interval["priority_score"] = row.get("priority_score")
        sge_records.append(sge_interval)

    return pd.DataFrame(sgc_records), pd.DataFrame(sge_records)


def build_subgenome_specific_interval_table(
    sub_df: pd.DataFrame,
    part1_contigs: pd.DataFrame,
    sgc_lookup: Dict[str, Dict[str, object]],
    sge_lookup: Dict[str, Dict[str, object]],
) -> pd.DataFrame:
    records: List[Dict[str, object]] = []

    length_lookup = part1_contigs.drop_duplicates(subset=["contig_id"]).set_index("contig_id")["contig_length"].to_dict()

    for _, row in sub_df.iterrows():
        source_contig = str(row["contig"])
        hotspot_subgenome = str(row["subgenome"])
        crosswalk_lookup = sgc_lookup if hotspot_subgenome == "sgC" else sge_lookup

        if source_contig in crosswalk_lookup:
            interval = transfer_interval(
                source_contig=source_contig,
                source_start=int(row["min_pos"]),
                source_end=int(row["max_pos"]),
                crosswalk_lookup=crosswalk_lookup,
            )
        else:
            source_length = length_lookup.get(source_contig)
            interval = {
                "source_contig": source_contig,
                "target_seqid": None,
                "target_chromosome_label": None,
                "mapping_method": "unplaced_or_unmapped_scaffold",
                "coordinate_method": "unresolved",
                "mapping_confidence": "unresolved",
                "source_start": int(row["min_pos"]),
                "source_end": int(row["max_pos"]),
                "source_span_bp": int(row["max_pos"]) - int(row["min_pos"]) + 1,
                "target_start": None,
                "target_end": None,
                "target_span_bp": None,
                "length_scale_factor": None,
                "relative_length_difference": None,
                "annotation_ready": False,
                "mapping_note": "No chromosome-level crosswalk was available for this scaffold-like hotspot.",
            }
            if source_length is not None:
                interval["source_length"] = int(source_length)

        interval["hotspot_identifier"] = f"{row['contrast']}::{source_contig}::{hotspot_subgenome}"
        interval["contrast"] = row["contrast"]
        interval["hotspot_subgenome"] = hotspot_subgenome
        interval["priority_score"] = row.get("priority_score")
        records.append(interval)

    return pd.DataFrame(records)


def build_same_contig_interval_table(
    same_df: pd.DataFrame,
    part1_contigs: pd.DataFrame,
    sgc_lookup: Dict[str, Dict[str, object]],
    sge_lookup: Dict[str, Dict[str, object]],
) -> pd.DataFrame:
    records: List[Dict[str, object]] = []
    length_lookup = part1_contigs.drop_duplicates(subset=["contig_id"]).set_index("contig_id")["contig_length"].to_dict()

    for _, row in same_df.iterrows():
        source_contig = str(row["contig"])
        entry: Dict[str, object] = {
            "hotspot_identifier": f"{row['contrast']}::{source_contig}::same_contig_support",
            "contrast": row["contrast"],
            "source_contig": source_contig,
            "sgC_min_pos": row.get("sgC_min_pos"),
            "sgC_max_pos": row.get("sgC_max_pos"),
            "sgE_min_pos": row.get("sgE_min_pos"),
            "sgE_max_pos": row.get("sgE_max_pos"),
            "priority_score": row.get("priority_score"),
            "hotspot_category": row.get("hotspot_category"),
            "mapping_note": "Same-contig support hotspots are retained for reference but are not directly annotated unless a chromosome-level crosswalk is available.",
        }
        if source_contig in length_lookup:
            entry["source_length"] = int(length_lookup[source_contig])
        records.append(entry)

    return pd.DataFrame(records)


def load_part4_tables(part4_workbook: Path) -> Dict[str, pd.DataFrame]:
    sheets = {
        "primary_direction_summary": "04_primary_direction_summary",
        "primary_shared_hotspots": "09A_primary_shared_hotspots",
        "primary_subgenome_specific": "10A_primary_subgenome_specific",
        "primary_same_contig_support": "11A_primary_same_contig_support",
        "final_shared_hotspots": "09_final_shared_hotspots",
        "final_subgenome_specific": "10_final_subgenome_specific",
        "final_same_contig_support": "11_final_same_contig_support",
    }
    return {key: pd.read_excel(part4_workbook, sheet_name=value) for key, value in sheets.items()}


def make_keyword_summary(gene_df: pd.DataFrame) -> pd.DataFrame:
    if gene_df.empty:
        return pd.DataFrame(columns=["keyword_category", "gene_count", "unique_gene_names"])
    records: List[Dict[str, object]] = []
    expanded: List[Tuple[str, str]] = []
    for _, row in gene_df.iterrows():
        for category in str(row.get("keyword_categories", "")).split("; "):
            if category:
                expanded.append((category, str(row.get("gene_name", ""))))
    if not expanded:
        return pd.DataFrame(columns=["keyword_category", "gene_count", "unique_gene_names"])
    expanded_df = pd.DataFrame(expanded, columns=["keyword_category", "gene_name"])
    for category, sub in expanded_df.groupby("keyword_category"):
        records.append(
            {
                "keyword_category": category,
                "gene_count": int(sub.shape[0]),
                "unique_gene_names": "; ".join(sorted({x for x in sub["gene_name"] if x and x != "nan"})),
            }
        )
    return pd.DataFrame(records).sort_values(["gene_count", "keyword_category"], ascending=[False, True]).reset_index(drop=True)


def build_primary_gene_shortlist(
    sgc_genes: pd.DataFrame,
    sge_genes: pd.DataFrame,
    keyword_priority_only: bool,
) -> pd.DataFrame:
    merged = pd.concat([sgc_genes, sge_genes], ignore_index=True)
    if merged.empty:
        return merged

    if keyword_priority_only:
        merged = merged.loc[merged["keyword_hit_count"] > 0].copy()

    if merged.empty:
        return merged

    keep_cols = [
        "hotspot_subgenome",
        "hotspot_identifier",
        "contrast",
        "source_contig",
        "target_seqid",
        "target_chromosome_label",
        "mapping_method",
        "coordinate_method",
        "mapping_confidence",
        "gene_identifier",
        "gene_name",
        "description",
        "gene_biotype",
        "keyword_categories",
        "keyword_terms",
        "keyword_hit_count",
        "overlap_bp",
        "keyword_priority_score",
    ]
    merged = merged[keep_cols].copy()
    merged["within_subgenome_rank"] = (
        merged.groupby("hotspot_subgenome")["keyword_priority_score"].rank(method="first", ascending=False)
    )
    merged = merged.sort_values(
        ["keyword_hit_count", "keyword_priority_score", "overlap_bp", "gene_name"],
        ascending=[False, False, False, True],
    ).reset_index(drop=True)
    merged["global_candidate_rank"] = range(1, len(merged) + 1)
    return merged


def safe_sheet_name(name: str) -> str:
    return name[:31]


def main() -> None:
    args = parse_arguments()

    part1_workbook = Path(args.part1_workbook)
    part4_workbook = Path(args.part4_workbook)
    sgc_gff3 = Path(args.sgc_gff3)
    sge_gff3 = Path(args.sge_gff3)
    output_workbook = Path(args.output_workbook)
    supplementary_output_workbook = Path(args.supplementary_output_workbook)

    validate_inputs([part1_workbook, part4_workbook, sgc_gff3, sge_gff3])

    print("Loading Part 1 contig inventories")
    part1_contigs = load_part1_contigs(part1_workbook)

    print("Loading Part 4 hotspot tables")
    part4 = load_part4_tables(part4_workbook)

    print("Parsing GFF3 annotation files")
    sgc_regions, sgc_genes = parse_gff3(sgc_gff3, "C")
    sge_regions, sge_genes = parse_gff3(sge_gff3, "E")

    print("Building contig-to-annotation crosswalks")
    sgc_crosswalk = build_contig_crosswalk(
        vcf_contigs=part1_contigs,
        region_df=sgc_regions,
        subgenome="C",
        maximum_scaled_relative_length_difference=args.maximum_scaled_relative_length_difference,
        minimum_high_confidence_gap=args.minimum_high_confidence_gap,
    )
    sge_crosswalk = build_contig_crosswalk(
        vcf_contigs=part1_contigs,
        region_df=sge_regions,
        subgenome="E",
        maximum_scaled_relative_length_difference=args.maximum_scaled_relative_length_difference,
        minimum_high_confidence_gap=args.minimum_high_confidence_gap,
    )

    sgc_lookup = build_crosswalk_lookup(sgc_crosswalk)
    sge_lookup = build_crosswalk_lookup(sge_crosswalk)

    print("Transferring primary hotspot intervals into the annotation coordinate system")
    primary_shared = part4["primary_shared_hotspots"].head(args.max_shared_hotspots).copy()
    primary_subspecific = part4["primary_subgenome_specific"].head(args.max_subgenome_specific_hotspots).copy()
    primary_same_contig = part4["primary_same_contig_support"].head(args.max_same_contig_hotspots).copy()

    primary_shared_sgc_intervals, primary_shared_sge_intervals = build_primary_shared_interval_table(
        primary_shared=primary_shared,
        contig_inventory=part1_contigs,
        sgc_lookup=sgc_lookup,
        sge_lookup=sge_lookup,
    )
    primary_subspecific_intervals = build_subgenome_specific_interval_table(
        sub_df=primary_subspecific,
        part1_contigs=part1_contigs,
        sgc_lookup=sgc_lookup,
        sge_lookup=sge_lookup,
    )
    primary_same_contig_intervals = build_same_contig_interval_table(
        same_df=primary_same_contig,
        part1_contigs=part1_contigs,
        sgc_lookup=sgc_lookup,
        sge_lookup=sge_lookup,
    )

    print("Annotating primary hotspots with overlapping genes")
    primary_shared_sgc_genes = intersect_genes(primary_shared_sgc_intervals, sgc_genes, "shared", "sgC")
    primary_shared_sge_genes = intersect_genes(primary_shared_sge_intervals, sge_genes, "shared", "sgE")
    primary_subspecific_genes_sgc = intersect_genes(
        primary_subspecific_intervals.loc[primary_subspecific_intervals["hotspot_subgenome"].eq("sgC")],
        sgc_genes,
        "subgenome_specific",
        "sgC",
    )
    primary_subspecific_genes_sge = intersect_genes(
        primary_subspecific_intervals.loc[primary_subspecific_intervals["hotspot_subgenome"].eq("sgE")],
        sge_genes,
        "subgenome_specific",
        "sgE",
    )
    primary_subspecific_genes = pd.concat(
        [primary_subspecific_genes_sgc, primary_subspecific_genes_sge],
        ignore_index=True,
    )

    print("Annotating full hotspot collections for supplementary output")
    full_shared = part4["final_shared_hotspots"].head(args.max_shared_hotspots).copy()
    full_subspecific = part4["final_subgenome_specific"].head(args.max_subgenome_specific_hotspots).copy()
    full_same_contig = part4["final_same_contig_support"].head(args.max_same_contig_hotspots).copy()

    full_shared_sgc_intervals, full_shared_sge_intervals = build_primary_shared_interval_table(
        primary_shared=full_shared,
        contig_inventory=part1_contigs,
        sgc_lookup=sgc_lookup,
        sge_lookup=sge_lookup,
    )
    full_subspecific_intervals = build_subgenome_specific_interval_table(
        sub_df=full_subspecific,
        part1_contigs=part1_contigs,
        sgc_lookup=sgc_lookup,
        sge_lookup=sge_lookup,
    )
    full_same_contig_intervals = build_same_contig_interval_table(
        same_df=full_same_contig,
        part1_contigs=part1_contigs,
        sgc_lookup=sgc_lookup,
        sge_lookup=sge_lookup,
    )

    full_shared_sgc_genes = intersect_genes(full_shared_sgc_intervals, sgc_genes, "shared", "sgC")
    full_shared_sge_genes = intersect_genes(full_shared_sge_intervals, sge_genes, "shared", "sgE")
    full_subspecific_genes_sgc = intersect_genes(
        full_subspecific_intervals.loc[full_subspecific_intervals["hotspot_subgenome"].eq("sgC")],
        sgc_genes,
        "subgenome_specific",
        "sgC",
    )
    full_subspecific_genes_sge = intersect_genes(
        full_subspecific_intervals.loc[full_subspecific_intervals["hotspot_subgenome"].eq("sgE")],
        sge_genes,
        "subgenome_specific",
        "sgE",
    )
    full_subspecific_genes = pd.concat(
        [full_subspecific_genes_sgc, full_subspecific_genes_sge],
        ignore_index=True,
    )

    print("Building summary tables")
    primary_shared_sgc_summary = summarize_interval_table(primary_shared_sgc_intervals, primary_shared_sgc_genes)
    primary_shared_sge_summary = summarize_interval_table(primary_shared_sge_intervals, primary_shared_sge_genes)
    primary_subspecific_summary = summarize_interval_table(primary_subspecific_intervals, primary_subspecific_genes)

    keyword_summary = make_keyword_summary(pd.concat([primary_shared_sgc_genes, primary_shared_sge_genes, primary_subspecific_genes], ignore_index=True))
    primary_gene_shortlist = build_primary_gene_shortlist(
        sgc_genes=primary_shared_sgc_genes,
        sge_genes=primary_shared_sge_genes,
        keyword_priority_only=args.keyword_priority_only,
    )

    readiness_summary = pd.DataFrame(
        [
            {
                "section": "sgC_crosswalk",
                "n_source_contigs": int(sgc_crosswalk.shape[0]),
                "n_mapping_ready": int(sgc_crosswalk["mapping_ready"].sum()),
                "n_exact_length_matches": int(sgc_crosswalk["exact_length_match"].sum()),
                "n_high_confidence": int(sgc_crosswalk["mapping_confidence"].eq("high").sum()),
                "n_medium_confidence": int(sgc_crosswalk["mapping_confidence"].eq("medium").sum()),
                "n_low_confidence": int(sgc_crosswalk["mapping_confidence"].eq("low").sum()),
                "n_unresolved": int(sgc_crosswalk["mapping_confidence"].eq("unresolved").sum()),
            },
            {
                "section": "sgE_crosswalk",
                "n_source_contigs": int(sge_crosswalk.shape[0]),
                "n_mapping_ready": int(sge_crosswalk["mapping_ready"].sum()),
                "n_exact_length_matches": int(sge_crosswalk["exact_length_match"].sum()),
                "n_high_confidence": int(sge_crosswalk["mapping_confidence"].eq("high").sum()),
                "n_medium_confidence": int(sge_crosswalk["mapping_confidence"].eq("medium").sum()),
                "n_low_confidence": int(sge_crosswalk["mapping_confidence"].eq("low").sum()),
                "n_unresolved": int(sge_crosswalk["mapping_confidence"].eq("unresolved").sum()),
            },
            {
                "section": "primary_shared_hotspots_sgC",
                "n_source_contigs": int(primary_shared_sgc_intervals.shape[0]),
                "n_mapping_ready": int(primary_shared_sgc_intervals["annotation_ready"].sum()),
                "n_exact_length_matches": int(primary_shared_sgc_intervals["mapping_method"].eq("exact_length_match").sum()),
                "n_high_confidence": int(primary_shared_sgc_intervals["mapping_confidence"].eq("high").sum()),
                "n_medium_confidence": int(primary_shared_sgc_intervals["mapping_confidence"].eq("medium").sum()),
                "n_low_confidence": int(primary_shared_sgc_intervals["mapping_confidence"].eq("low").sum()),
                "n_unresolved": int(primary_shared_sgc_intervals["mapping_confidence"].eq("unresolved").sum()),
            },
            {
                "section": "primary_shared_hotspots_sgE",
                "n_source_contigs": int(primary_shared_sge_intervals.shape[0]),
                "n_mapping_ready": int(primary_shared_sge_intervals["annotation_ready"].sum()),
                "n_exact_length_matches": int(primary_shared_sge_intervals["mapping_method"].eq("exact_length_match").sum()),
                "n_high_confidence": int(primary_shared_sge_intervals["mapping_confidence"].eq("high").sum()),
                "n_medium_confidence": int(primary_shared_sge_intervals["mapping_confidence"].eq("medium").sum()),
                "n_low_confidence": int(primary_shared_sge_intervals["mapping_confidence"].eq("low").sum()),
                "n_unresolved": int(primary_shared_sge_intervals["mapping_confidence"].eq("unresolved").sum()),
            },
        ]
    )

    run_parameters = pd.DataFrame(
        [
            {"parameter": "part1_workbook", "value": str(part1_workbook)},
            {"parameter": "part4_workbook", "value": str(part4_workbook)},
            {"parameter": "sgc_gff3", "value": str(sgc_gff3)},
            {"parameter": "sge_gff3", "value": str(sge_gff3)},
            {"parameter": "maximum_scaled_relative_length_difference", "value": args.maximum_scaled_relative_length_difference},
            {"parameter": "minimum_high_confidence_gap", "value": args.minimum_high_confidence_gap},
            {"parameter": "max_shared_hotspots", "value": args.max_shared_hotspots},
            {"parameter": "max_subgenome_specific_hotspots", "value": args.max_subgenome_specific_hotspots},
            {"parameter": "max_same_contig_hotspots", "value": args.max_same_contig_hotspots},
            {
                "parameter": "important_note",
                "value": "The GFF3 annotation is from the ET-39 HiFi RefSeq assembly. If VCF hotspot coordinates were derived from a different assembly version, coordinate transfer is approximate unless an exact contig-length match was found.",
            },
        ]
    )

    input_inventory = pd.DataFrame(
        [
            {"input_label": "part1_workbook", "path": str(part1_workbook)},
            {"input_label": "part4_workbook", "path": str(part4_workbook)},
            {"input_label": "sgC_gff3", "path": str(sgc_gff3)},
            {"input_label": "sgE_gff3", "path": str(sge_gff3)},
        ]
    )

    manuscript_outputs = {
        "01_run_parameters": run_parameters,
        "02_input_inventory": input_inventory,
        "03_annotation_readiness": readiness_summary,
        "04_sgC_contig_crosswalk": sgc_crosswalk,
        "05_sgE_contig_crosswalk": sge_crosswalk,
        "06_primary_shared_hotspot_sgC_mapping": primary_shared_sgc_summary,
        "07_primary_shared_hotspot_sgE_mapping": primary_shared_sge_summary,
        "08_primary_shared_hotspot_sgC_genes": primary_shared_sgc_genes,
        "09_primary_shared_hotspot_sgE_genes": primary_shared_sge_genes,
        "10_primary_subgenome_specific_mapping": primary_subspecific_summary,
        "11_primary_subgenome_specific_genes": primary_subspecific_genes,
        "12_primary_same_contig_reference": primary_same_contig_intervals,
        "13_keyword_summary": keyword_summary,
        "14_primary_candidate_gene_shortlist": primary_gene_shortlist,
        "Table_08_primary_shared_hotspot_genes": pd.concat([primary_shared_sgc_genes, primary_shared_sge_genes], ignore_index=True),
        "Table_09_primary_candidate_gene_shortlist": primary_gene_shortlist,
        "Table_10_crosswalk_diagnostics": pd.concat([sgc_crosswalk.assign(file="sgC"), sge_crosswalk.assign(file="sgE")], ignore_index=True),
    }

    supplementary_outputs = {
        "S5_run_parameters": run_parameters,
        "S5_input_inventory": input_inventory,
        "S5_sgC_region_features": sgc_regions,
        "S5_sgE_region_features": sge_regions,
        "S5_sgC_gene_features": sgc_genes,
        "S5_sgE_gene_features": sge_genes,
        "S5_sgC_crosswalk": sgc_crosswalk,
        "S5_sgE_crosswalk": sge_crosswalk,
        "S5_all_shared_hotspots_sgC_mapping": summarize_interval_table(full_shared_sgc_intervals, full_shared_sgc_genes),
        "S5_all_shared_hotspots_sgE_mapping": summarize_interval_table(full_shared_sge_intervals, full_shared_sge_genes),
        "S5_all_shared_hotspots_sgC_genes": full_shared_sgc_genes,
        "S5_all_shared_hotspots_sgE_genes": full_shared_sge_genes,
        "S5_all_subspecific_mapping": summarize_interval_table(full_subspecific_intervals, full_subspecific_genes),
        "S5_all_subspecific_genes": full_subspecific_genes,
        "S5_same_contig_reference": full_same_contig_intervals,
        "S5_keyword_summary": keyword_summary,
    }

    output_workbook.parent.mkdir(parents=True, exist_ok=True)
    supplementary_output_workbook.parent.mkdir(parents=True, exist_ok=True)

    print(f"Writing workbook: {output_workbook}")
    with pd.ExcelWriter(output_workbook, engine="openpyxl") as writer:
        for sheet_name, df in manuscript_outputs.items():
            df.to_excel(writer, sheet_name=safe_sheet_name(sheet_name), index=False)

    print(f"Writing workbook: {supplementary_output_workbook}")
    with pd.ExcelWriter(supplementary_output_workbook, engine="openpyxl") as writer:
        for sheet_name, df in supplementary_outputs.items():
            df.to_excel(writer, sheet_name=safe_sheet_name(sheet_name), index=False)

    print("Part 5 hotspot annotation analysis finished successfully.")


if __name__ == "__main__":
    main()
