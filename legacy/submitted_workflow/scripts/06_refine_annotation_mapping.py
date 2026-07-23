#!/usr/bin/env python3
"""
Refine hotspot annotation mapping with chromosome-label reconciliation.

The script uses GFF3 chromosome labels and assembly metadata to reconcile source
contig labels with RefSeq annotation identifiers. It transfers hotspot intervals
by length scaling when direct sequence aliases are unavailable and summarizes the
overlapping gene space with transparent mapping diagnostics.
"""

from __future__ import annotations

import argparse
import json
import math
import re
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple

import numpy as np
import pandas as pd


PRIMARY_SHARED_SGC_SHEET = "06_primary_shared_hotspot_sgC_m"
PRIMARY_SHARED_SGE_SHEET = "07_primary_shared_hotspot_sgE_m"
PRIMARY_SHARED_SGC_GENES_SHEET = "08_primary_shared_hotspot_sgC_g"
PRIMARY_SHARED_SGE_GENES_SHEET = "09_primary_shared_hotspot_sgE_g"
PRIMARY_SUBGENOME_SHEET = "10_primary_subgenome_specific_m"
PRIMARY_SUBGENOME_GENES_SHEET = "11_primary_subgenome_specific_g"

OUTPUT_SHEETS = {
    "run_parameters": "01_run_parameters",
    "input_inventory": "02_input_inventory",
    "jsonl_inventory": "03_jsonl_inventory",
    "annotation_readiness": "04_label_rescue_readiness",
    "sgc_labels": "05_sgC_label_dictionary",
    "sge_labels": "06_sgE_label_dictionary",
    "shared_sgc_map": "07_rescued_shared_sgC_map",
    "shared_sge_map": "08_rescued_shared_sgE_map",
    "shared_sgc_genes": "09_rescued_shared_sgC_genes",
    "shared_sge_genes": "10_rescued_shared_sgE_genes",
    "subgenome_map": "11_rescued_subgenome_specific",
    "subgenome_genes": "12_rescued_subgenome_genes",
    "candidate_shortlist": "13_rescued_candidate_shortlist",
    "shared_table": "Table_08_rescued_shared_genes",
    "candidate_table": "Table_09_rescued_candidates",
    "diagnostics": "Table_10_label_rescue_diagnostics",
}

KEYWORD_PATTERNS = {
    "defense_receptor": [
        r"\breceptor\b",
        r"\breceptor-like\b",
        r"\bwall-associated receptor kinase\b",
        r"\bleucine-rich\b",
        r"\bLRR\b",
        r"\bRLK\b",
        r"\bRLP\b",
    ],
    "kinase_signaling": [
        r"\bkinase\b",
        r"\bMAPK\b",
        r"\bsignaling\b",
        r"\bcalmodulin\b",
        r"\bphosphatase\b",
    ],
    "stress_response": [
        r"\bstress\b",
        r"\bpathogen\b",
        r"\bdisease\b",
        r"\bdefense\b",
        r"\bautophagy\b",
        r"\boxidative\b",
        r"\bwound\b",
    ],
    "transcriptional_regulation": [
        r"\btranscription factor\b",
        r"\bMYB\b",
        r"\bWRKY\b",
        r"\bNAC\b",
        r"\bERF\b",
        r"\bbHLH\b",
    ],
    "transport_membrane": [
        r"\btransporter\b",
        r"\bchannel\b",
        r"\bATPase\b",
        r"\bmembrane\b",
        r"\bABC\b",
    ],
}

LETTER_TO_NUMBER = {chr(ord("A") + i): i + 1 for i in range(12)}


def utc_now_string() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def clean_sheet_name(name: str) -> str:
    invalid = set(r'[]:*?/\\')
    cleaned = "".join("_" if c in invalid else c for c in name)
    return cleaned[:31]


def read_excel(path: Path, sheet_name: str) -> pd.DataFrame:
    return pd.read_excel(path, sheet_name=sheet_name)


def required_sheet(path: Path, sheet_name: str) -> pd.DataFrame:
    workbook = pd.ExcelFile(path)
    if sheet_name not in workbook.sheet_names:
        raise ValueError(f"Required sheet not found in {path.name}: {sheet_name}")
    return pd.read_excel(path, sheet_name=sheet_name)


def parse_jsonl_inventory(path: Path) -> pd.DataFrame:
    records = []
    if not path.exists():
        return pd.DataFrame([{
            "jsonl_present": False,
            "jsonl_path": str(path),
            "n_records": 0,
            "assembly_name": None,
            "assembly_accession": None,
            "source_database": None,
            "sequence_alias_fields_present": False,
            "note": "JSONL file was not found.",
        }])

    n_records = 0
    first_record = None
    sequence_alias_fields_present = False

    with open(path, "r", encoding="utf-8", errors="ignore") as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            n_records += 1
            record = json.loads(line)
            if first_record is None:
                first_record = record
            if contains_sequence_alias_fields(record):
                sequence_alias_fields_present = True

    if first_record is None:
        return pd.DataFrame([{
            "jsonl_present": True,
            "jsonl_path": str(path),
            "n_records": 0,
            "assembly_name": None,
            "assembly_accession": None,
            "source_database": None,
            "sequence_alias_fields_present": False,
            "note": "JSONL file was empty.",
        }])

    records.append({
        "jsonl_present": True,
        "jsonl_path": str(path),
        "n_records": n_records,
        "assembly_name": deep_get(first_record, ["assemblyInfo", "assemblyName"]),
        "assembly_accession": first_record.get("accession"),
        "source_database": first_record.get("sourceDatabase"),
        "sequence_alias_fields_present": sequence_alias_fields_present,
        "note": "This JSONL file is used as official assembly metadata context. If sequence-level aliases are absent, label-based rescue is applied.",
    })
    return pd.DataFrame(records)


def deep_get(obj, path: List[str]):
    cur = obj
    for key in path:
        if not isinstance(cur, dict) or key not in cur:
            return None
        cur = cur[key]
    return cur


def contains_sequence_alias_fields(obj) -> bool:
    if isinstance(obj, dict):
        for key, value in obj.items():
            key_lower = str(key).lower()
            if "sequence" in key_lower and ("name" in key_lower or "alias" in key_lower or "synonym" in key_lower):
                return True
            if contains_sequence_alias_fields(value):
                return True
    elif isinstance(obj, list):
        for value in obj:
            if contains_sequence_alias_fields(value):
                return True
    return False


def parse_gff3_regions_and_genes(path: Path) -> Tuple[pd.DataFrame, pd.DataFrame]:
    region_rows = []
    gene_rows = []

    with open(path, "r", encoding="utf-8", errors="ignore") as handle:
        for raw in handle:
            line = raw.rstrip("\n")
            if not line or line.startswith("#"):
                continue

            parts = line.split("\t")
            if len(parts) != 9:
                continue

            seqid, source, feature_type, start, end, score, strand, phase, attrs = parts
            attr_map = parse_attributes(attrs)

            if feature_type == "region":
                region_rows.append({
                    "seqid": seqid,
                    "start": int(start),
                    "end": int(end),
                    "length": int(end) - int(start) + 1,
                    "chromosome_label": attr_map.get("chromosome"),
                    "source": source,
                    "attributes": attrs,
                })
            elif feature_type == "gene":
                gene_rows.append({
                    "seqid": seqid,
                    "feature_type": feature_type,
                    "start": int(start),
                    "end": int(end),
                    "strand": strand,
                    "gene_id": attr_map.get("ID"),
                    "gene_symbol": attr_map.get("gene") or attr_map.get("Name"),
                    "name": attr_map.get("Name"),
                    "description": attr_map.get("description"),
                    "gene_biotype": attr_map.get("gene_biotype"),
                    "raw_attributes": attrs,
                })

    region_df = pd.DataFrame(region_rows)
    gene_df = pd.DataFrame(gene_rows)

    if region_df.empty:
        raise ValueError(f"No region features found in GFF3: {path}")
    if gene_df.empty:
        raise ValueError(f"No gene features found in GFF3: {path}")

    return region_df, gene_df


def parse_attributes(attr_text: str) -> Dict[str, str]:
    result = {}
    for item in attr_text.split(";"):
        if "=" in item:
            key, value = item.split("=", 1)
            result[key] = value
    return result


def build_label_dictionary(region_df: pd.DataFrame, subgenome_suffix: str) -> pd.DataFrame:
    df = region_df.copy()
    df["subgenome_suffix"] = subgenome_suffix
    df["label_normalized"] = df["chromosome_label"].astype(str).str.lower()
    df = df.sort_values(["label_normalized", "seqid"]).reset_index(drop=True)
    return df[["seqid", "chromosome_label", "label_normalized", "length", "subgenome_suffix"]]


def infer_label_from_source_contig(contig: str) -> Optional[str]:
    if not isinstance(contig, str):
        return None
    match = re.search(r"chr_([A-Z])_sg_([CE])", contig)
    if not match:
        return None
    letter = match.group(1).upper()
    subgenome = match.group(2).lower()
    if letter not in LETTER_TO_NUMBER:
        return None
    return f"{LETTER_TO_NUMBER[letter]}{subgenome}"


def make_crosswalk(
    existing_crosswalk: pd.DataFrame,
    label_dictionary: pd.DataFrame,
    subgenome_suffix: str,
) -> pd.DataFrame:
    lookup = {
        row["label_normalized"]: row
        for _, row in label_dictionary.iterrows()
    }

    rows = []
    for _, row in existing_crosswalk.iterrows():
        source_contig = row["source_contig"]
        inferred_label = infer_label_from_source_contig(source_contig)
        rescue_possible = inferred_label is not None and inferred_label in lookup

        if rescue_possible:
            target = lookup[inferred_label]
            mapping_method = "label_inference_from_source_contig"
            mapping_confidence = "medium"
            note = (
                "Source contig was matched to the annotation by inferred chromosome label "
                "derived from the contig name (A->1, B->2, ... with subgenome suffix)."
            )
            target_seqid = target["seqid"]
            target_length = target["length"]
            target_label = target["chromosome_label"]
            mapping_ready = True
        else:
            mapping_method = "unresolved"
            mapping_confidence = "unresolved"
            note = "No chromosome-label rescue was possible for this source contig."
            target_seqid = np.nan
            target_length = np.nan
            target_label = np.nan
            mapping_ready = False

        rows.append({
            "source_contig": source_contig,
            "source_subgenome": subgenome_suffix,
            "source_length": row.get("source_length"),
            "inferred_source_label": inferred_label,
            "target_seqid": target_seqid,
            "target_length": target_length,
            "target_chromosome_label": target_label,
            "mapping_method": mapping_method,
            "coordinate_method": "length_scaled" if mapping_ready else "unresolved",
            "mapping_confidence": mapping_confidence,
            "note": note,
            "mapping_ready": mapping_ready,
            "exact_length_match": False,
        })

    return pd.DataFrame(rows)


def remap_hotspots(
    hotspot_df: pd.DataFrame,
    crosswalk_df: pd.DataFrame,
) -> pd.DataFrame:
    crosswalk = crosswalk_df.set_index("source_contig").to_dict("index")
    rows = []

    for _, row in hotspot_df.iterrows():
        source_contig = row["source_contig"]
        source_length = row.get("source_length")
        mapping = crosswalk.get(source_contig)

        out = dict(row)

        if not mapping or not mapping["mapping_ready"]:
            out.update({
                "target_seqid": np.nan,
                "target_chromosome_label": np.nan,
                "mapping_method": mapping["mapping_method"] if mapping else "unresolved",
                "coordinate_method": "unresolved",
                "mapping_confidence": mapping["mapping_confidence"] if mapping else "unresolved",
                "target_start": np.nan,
                "target_end": np.nan,
                "target_span_bp": np.nan,
                "length_scale_factor": np.nan,
                "annotation_ready": False,
                "mapping_note": mapping["note"] if mapping else "No crosswalk could be derived for this hotspot.",
            })
            rows.append(out)
            continue

        if pd.isna(source_length) or source_length in (0, None):
            out.update({
                "target_seqid": mapping["target_seqid"],
                "target_chromosome_label": mapping["target_chromosome_label"],
                "mapping_method": mapping["mapping_method"],
                "coordinate_method": "unresolved",
                "mapping_confidence": "unresolved",
                "target_start": np.nan,
                "target_end": np.nan,
                "target_span_bp": np.nan,
                "length_scale_factor": np.nan,
                "annotation_ready": False,
                "mapping_note": "Source contig length was missing, preventing coordinate scaling.",
            })
            rows.append(out)
            continue

        target_length = mapping["target_length"]
        scale = float(target_length) / float(source_length)

        source_start = int(row["source_start"])
        source_end = int(row["source_end"])

        target_start = max(1, int(math.floor(source_start * scale)))
        target_end = min(int(target_length), int(math.ceil(source_end * scale)))

        if target_start > target_end:
            target_start, target_end = target_end, target_start

        out.update({
            "target_seqid": mapping["target_seqid"],
            "target_chromosome_label": mapping["target_chromosome_label"],
            "mapping_method": mapping["mapping_method"],
            "coordinate_method": "length_scaled",
            "mapping_confidence": mapping["mapping_confidence"],
            "target_start": target_start,
            "target_end": target_end,
            "target_span_bp": target_end - target_start + 1,
            "length_scale_factor": scale,
            "annotation_ready": True,
            "mapping_note": mapping["note"],
        })
        rows.append(out)

    return pd.DataFrame(rows)


def annotate_hotspots_with_genes(
    mapped_hotspots: pd.DataFrame,
    gene_df: pd.DataFrame,
) -> pd.DataFrame:
    rows = []

    genes_by_seqid = {
        seqid: sub_df.sort_values(["start", "end"]).reset_index(drop=True)
        for seqid, sub_df in gene_df.groupby("seqid", sort=False)
    }

    for _, hotspot in mapped_hotspots.iterrows():
        if not bool(hotspot.get("annotation_ready", False)):
            continue

        seqid = hotspot["target_seqid"]
        start = int(hotspot["target_start"])
        end = int(hotspot["target_end"])

        sub = genes_by_seqid.get(seqid)
        if sub is None or sub.empty:
            continue

        overlap = sub[(sub["start"] <= end) & (sub["end"] >= start)].copy()
        if overlap.empty:
            continue

        for _, gene in overlap.iterrows():
            description = gene.get("description")
            keyword_categories = detect_keyword_categories(description)
            rows.append({
                "hotspot_identifier": hotspot.get("hotspot_identifier"),
                "contrast": hotspot.get("contrast"),
                "source_contig": hotspot.get("source_contig"),
                "target_seqid": seqid,
                "target_chromosome_label": hotspot.get("target_chromosome_label"),
                "target_start": start,
                "target_end": end,
                "gene_seqid": gene["seqid"],
                "gene_start": gene["start"],
                "gene_end": gene["end"],
                "gene_id": gene.get("gene_id"),
                "gene_symbol": gene.get("gene_symbol"),
                "gene_name": gene.get("name"),
                "gene_biotype": gene.get("gene_biotype"),
                "strand": gene.get("strand"),
                "description": description,
                "keyword_categories": ", ".join(sorted(keyword_categories)) if keyword_categories else None,
                "keyword_supported": bool(keyword_categories),
            })

    return pd.DataFrame(rows)


def detect_keyword_categories(description: Optional[str]) -> List[str]:
    if not isinstance(description, str) or not description.strip():
        return []
    text = description.lower()
    categories = []
    for category, patterns in KEYWORD_PATTERNS.items():
        if any(re.search(pattern, text, re.IGNORECASE) for pattern in patterns):
            categories.append(category)
    return categories


def summarize_hotspot_gene_support(mapped_df: pd.DataFrame, gene_hits: pd.DataFrame) -> pd.DataFrame:
    out = mapped_df.copy()
    if out.empty:
        return out

    default_columns = {
        "n_overlapping_genes": 0,
        "n_protein_coding_genes": 0,
        "n_keyword_supported_genes": 0,
        "keyword_categories_detected": None,
    }

    for col, default_value in default_columns.items():
        if col not in out.columns:
            out[col] = default_value

    required_gene_columns = {"hotspot_identifier", "gene_id", "gene_biotype", "keyword_supported", "keyword_categories"}
    if gene_hits is None or gene_hits.empty or not required_gene_columns.issubset(set(gene_hits.columns)):
        out["n_overlapping_genes"] = pd.Series(out["n_overlapping_genes"]).fillna(0).astype(int)
        out["n_protein_coding_genes"] = pd.Series(out["n_protein_coding_genes"]).fillna(0).astype(int)
        out["n_keyword_supported_genes"] = pd.Series(out["n_keyword_supported_genes"]).fillna(0).astype(int)
        out["keyword_categories_detected"] = out["keyword_categories_detected"].replace("", np.nan)
        return out

    gene_summary = (
        gene_hits.groupby("hotspot_identifier", dropna=False)
        .agg(
            n_overlapping_genes=("gene_id", "count"),
            n_protein_coding_genes=("gene_biotype", lambda s: int((s == "protein_coding").sum())),
            n_keyword_supported_genes=("keyword_supported", lambda s: int(pd.Series(s).fillna(False).astype(bool).sum())),
            keyword_categories_detected=("keyword_categories", lambda s: ", ".join(sorted(set(x for x in s.dropna() if x)))),
        )
        .reset_index()
    )

    out = out.drop(columns=[c for c in default_columns if c in out.columns], errors="ignore").merge(
        gene_summary, how="left", on="hotspot_identifier"
    )

    for col, default_value in default_columns.items():
        if col not in out.columns:
            out[col] = default_value

    out["n_overlapping_genes"] = pd.Series(out["n_overlapping_genes"]).fillna(0).astype(int)
    out["n_protein_coding_genes"] = pd.Series(out["n_protein_coding_genes"]).fillna(0).astype(int)
    out["n_keyword_supported_genes"] = pd.Series(out["n_keyword_supported_genes"]).fillna(0).astype(int)
    out["keyword_categories_detected"] = out["keyword_categories_detected"].replace("", np.nan)
    return out


def build_candidate_shortlist(shared_sgc_genes: pd.DataFrame, shared_sge_genes: pd.DataFrame, subgenome_genes: pd.DataFrame) -> pd.DataFrame:
    all_hits = []
    for df, section in [
        (shared_sgc_genes, "primary_shared_sgC"),
        (shared_sge_genes, "primary_shared_sgE"),
        (subgenome_genes, "primary_subgenome_specific"),
    ]:
        if df is None or df.empty:
            continue
        tmp = df.copy()
        tmp["section"] = section
        all_hits.append(tmp)

    if not all_hits:
        return pd.DataFrame(columns=[
            "gene_id", "gene_symbol", "gene_name", "description", "gene_biotype",
            "n_supporting_hits", "n_unique_hotspots", "n_keyword_supported_hits",
            "sections_detected", "keyword_categories_detected"
        ])

    combined = pd.concat(all_hits, ignore_index=True)
    grouped = (
        combined.groupby(["gene_id", "gene_symbol", "gene_name", "description", "gene_biotype"], dropna=False)
        .agg(
            n_supporting_hits=("hotspot_identifier", "count"),
            n_unique_hotspots=("hotspot_identifier", "nunique"),
            n_keyword_supported_hits=("keyword_supported", lambda s: int(pd.Series(s).fillna(False).astype(bool).sum())),
            sections_detected=("section", lambda s: ", ".join(sorted(set(s)))),
            keyword_categories_detected=("keyword_categories", lambda s: ", ".join(sorted(set(x for x in s.dropna() if x)))),
        )
        .reset_index()
        .sort_values(
            ["n_keyword_supported_hits", "n_unique_hotspots", "n_supporting_hits", "gene_symbol", "gene_id"],
            ascending=[False, False, False, True, True],
            na_position="last",
        )
        .reset_index(drop=True)
    )
    return grouped


def build_readiness_table(sgc_crosswalk: pd.DataFrame, sge_crosswalk: pd.DataFrame, shared_sgc: pd.DataFrame, shared_sge: pd.DataFrame, subgenome_mapped: pd.DataFrame) -> pd.DataFrame:
    def summarize(section: str, df: pd.DataFrame) -> Dict[str, object]:
        if df is None or df.empty:
            return {
                "section": section,
                "n_source_contigs": 0,
                "n_mapping_ready": 0,
                "n_high_confidence": 0,
                "n_medium_confidence": 0,
                "n_low_confidence": 0,
                "n_unresolved": 0,
            }
        conf = df.get("mapping_confidence", pd.Series(dtype=object)).fillna("unresolved")
        return {
            "section": section,
            "n_source_contigs": int(len(df)),
            "n_mapping_ready": int(pd.Series(df.get("annotation_ready", df.get("mapping_ready", False))).fillna(False).astype(bool).sum()),
            "n_high_confidence": int((conf == "high").sum()),
            "n_medium_confidence": int((conf == "medium").sum()),
            "n_low_confidence": int((conf == "low").sum()),
            "n_unresolved": int((conf == "unresolved").sum()),
        }

    rows = [
        summarize("sgC_label_rescue", sgc_crosswalk),
        summarize("sgE_label_rescue", sge_crosswalk),
        summarize("primary_shared_hotspot_sgC", shared_sgc),
        summarize("primary_shared_hotspot_sgE", shared_sge),
        summarize("primary_subgenome_specific", subgenome_mapped),
    ]
    return pd.DataFrame(rows)


def build_input_inventory(args) -> pd.DataFrame:
    rows = []
    for label, path in [
        ("part5_workbook", args.part5_workbook),
        ("supplementary_workbook", args.supplementary_workbook),
        ("sgc_gff3", args.sgc_gff3),
        ("sge_gff3", args.sge_gff3),
        ("assembly_data_report_jsonl", args.assembly_data_report_jsonl),
    ]:
        p = Path(path)
        rows.append({
            "input_label": label,
            "path": str(p),
            "exists": p.exists(),
            "size_bytes": p.stat().st_size if p.exists() else np.nan,
        })
    return pd.DataFrame(rows)


def build_run_parameters(args) -> pd.DataFrame:
    return pd.DataFrame([{
        "run_utc": utc_now_string(),
        "script": Path(__file__).name,
        "part5_workbook": str(args.part5_workbook),
        "supplementary_workbook": str(args.supplementary_workbook),
        "sgc_gff3": str(args.sgc_gff3),
        "sge_gff3": str(args.sge_gff3),
        "assembly_data_report_jsonl": str(args.assembly_data_report_jsonl),
        "output_workbook": str(args.output_workbook),
        "supplementary_output_workbook": str(args.supplementary_output_workbook),
        "rescue_strategy": "Label inference from source contig names combined with RefSeq GFF3 chromosome labels",
    }])


def write_workbook(path: Path, sheets: Dict[str, pd.DataFrame]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with pd.ExcelWriter(path, engine="openpyxl") as writer:
        for sheet_name, df in sheets.items():
            clean_name = clean_sheet_name(sheet_name)
            out = df.copy()
            out.to_excel(writer, index=False, sheet_name=clean_name)


def prepare_gene_table(df: pd.DataFrame) -> pd.DataFrame:
    if df is None or df.empty:
        return pd.DataFrame(columns=[
            "hotspot_identifier", "contrast", "source_contig", "target_chromosome_label",
            "gene_id", "gene_symbol", "gene_name", "gene_biotype", "description", "keyword_categories"
        ])
    cols = [
        "hotspot_identifier", "contrast", "source_contig", "target_chromosome_label",
        "gene_id", "gene_symbol", "gene_name", "gene_biotype", "description", "keyword_categories"
    ]
    cols = [c for c in cols if c in df.columns]
    return df[cols].copy()


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Part 5c label-based hotspot annotation rescue for the coffee introgression analysis."
    )
    parser.add_argument("--part5-workbook", required=True, type=Path)
    parser.add_argument("--supplementary-workbook", required=True, type=Path)
    parser.add_argument("--sgc-gff3", required=True, type=Path)
    parser.add_argument("--sge-gff3", required=True, type=Path)
    parser.add_argument("--assembly-data-report-jsonl", required=True, type=Path)
    parser.add_argument("--output-workbook", required=True, type=Path)
    parser.add_argument("--supplementary-output-workbook", required=True, type=Path)
    return parser.parse_args()


def main() -> None:
    args = parse_args()

    print("Loading Part 5 workbooks")
    primary_shared_sgc = required_sheet(args.part5_workbook, PRIMARY_SHARED_SGC_SHEET)
    primary_shared_sge = required_sheet(args.part5_workbook, PRIMARY_SHARED_SGE_SHEET)
    primary_subgenome = required_sheet(args.part5_workbook, PRIMARY_SUBGENOME_SHEET)
    sgc_crosswalk_existing = required_sheet(args.part5_workbook, "04_sgC_contig_crosswalk")
    sge_crosswalk_existing = required_sheet(args.part5_workbook, "05_sgE_contig_crosswalk")

    print("Loading JSONL metadata context")
    jsonl_inventory = parse_jsonl_inventory(args.assembly_data_report_jsonl)

    print("Parsing GFF3 annotation files")
    sgc_regions, sgc_genes = parse_gff3_regions_and_genes(args.sgc_gff3)
    sge_regions, sge_genes = parse_gff3_regions_and_genes(args.sge_gff3)

    print("Building label dictionaries")
    sgc_label_dict = build_label_dictionary(sgc_regions, "C")
    sge_label_dict = build_label_dictionary(sge_regions, "E")

    print("Rescuing contig crosswalks by chromosome labels")
    sgc_crosswalk = make_crosswalk(sgc_crosswalk_existing, sgc_label_dict, "C")
    sge_crosswalk = make_crosswalk(sge_crosswalk_existing, sge_label_dict, "E")

    # Add source lengths for shared hotspot rows from the crosswalk tables if missing.
    sgc_lengths = sgc_crosswalk.set_index("source_contig")["source_length"].to_dict()
    sge_lengths = sge_crosswalk.set_index("source_contig")["source_length"].to_dict()

    if "source_length" not in primary_shared_sgc.columns:
        primary_shared_sgc = primary_shared_sgc.copy()
        primary_shared_sgc["source_length"] = primary_shared_sgc["source_contig"].map(sgc_lengths)

    if "source_length" not in primary_shared_sge.columns:
        primary_shared_sge = primary_shared_sge.copy()
        primary_shared_sge["source_length"] = primary_shared_sge["source_contig"].map(sge_lengths)

    if "source_length" not in primary_subgenome.columns:
        primary_subgenome = primary_subgenome.copy()
        combined_lengths = {}
        combined_lengths.update(sgc_lengths)
        combined_lengths.update(sge_lengths)
        primary_subgenome["source_length"] = primary_subgenome["source_contig"].map(combined_lengths)

    print("Transferring hotspot intervals into the annotation coordinate system")
    rescued_shared_sgc = remap_hotspots(primary_shared_sgc, sgc_crosswalk)
    rescued_shared_sge = remap_hotspots(primary_shared_sge, sge_crosswalk)

    # Route subgenome-specific rows to the matching crosswalk by hotspot subgenome.
    sub_rows = []
    for _, row in primary_subgenome.iterrows():
        sub = str(row.get("hotspot_subgenome", "")).lower()
        if sub == "sgc" or sub == "c":
            mapped = remap_hotspots(pd.DataFrame([row]), sgc_crosswalk)
        elif sub == "sge" or sub == "e":
            mapped = remap_hotspots(pd.DataFrame([row]), sge_crosswalk)
        else:
            mapped = remap_hotspots(pd.DataFrame([row]), pd.DataFrame(columns=sgc_crosswalk.columns))
        sub_rows.append(mapped)
    rescued_subgenome = pd.concat(sub_rows, ignore_index=True) if sub_rows else pd.DataFrame()

    print("Annotating rescued hotspots with overlapping genes")
    shared_sgc_genes = annotate_hotspots_with_genes(rescued_shared_sgc, sgc_genes)
    shared_sge_genes = annotate_hotspots_with_genes(rescued_shared_sge, sge_genes)

    sgc_sub = rescued_subgenome[rescued_subgenome.get("hotspot_subgenome", pd.Series(dtype=object)).astype(str).str.lower().eq("sgc")].copy()
    sge_sub = rescued_subgenome[rescued_subgenome.get("hotspot_subgenome", pd.Series(dtype=object)).astype(str).str.lower().eq("sge")].copy()

    subgenome_genes_sgc = annotate_hotspots_with_genes(sgc_sub, sgc_genes)
    subgenome_genes_sge = annotate_hotspots_with_genes(sge_sub, sge_genes)
    rescued_subgenome_genes = pd.concat([subgenome_genes_sgc, subgenome_genes_sge], ignore_index=True)

    print("Summarising support and preparing candidate tables")
    rescued_shared_sgc = summarize_hotspot_gene_support(rescued_shared_sgc, shared_sgc_genes)
    rescued_shared_sge = summarize_hotspot_gene_support(rescued_shared_sge, shared_sge_genes)
    rescued_subgenome = summarize_hotspot_gene_support(rescued_subgenome, rescued_subgenome_genes)
    candidate_shortlist = build_candidate_shortlist(shared_sgc_genes, shared_sge_genes, rescued_subgenome_genes)
    readiness = build_readiness_table(sgc_crosswalk, sge_crosswalk, rescued_shared_sgc, rescued_shared_sge, rescued_subgenome)

    diagnostics = pd.DataFrame([
        {
            "metric": "jsonl_has_sequence_alias_fields",
            "value": bool(jsonl_inventory["sequence_alias_fields_present"].fillna(False).iloc[0]),
            "note": "False is acceptable. This script can still rescue crosswalks by label inference.",
        },
        {
            "metric": "shared_sgC_annotation_ready",
            "value": int(pd.Series(rescued_shared_sgc.get("annotation_ready", False)).fillna(False).astype(bool).sum()),
            "note": "Number of primary shared sgC hotspots successfully transferred into the annotation coordinate system.",
        },
        {
            "metric": "shared_sgE_annotation_ready",
            "value": int(pd.Series(rescued_shared_sge.get("annotation_ready", False)).fillna(False).astype(bool).sum()),
            "note": "Number of primary shared sgE hotspots successfully transferred into the annotation coordinate system.",
        },
        {
            "metric": "shared_sgC_gene_hits",
            "value": int(len(shared_sgc_genes)),
            "note": "Overlapping gene rows recovered for rescued primary shared sgC hotspots.",
        },
        {
            "metric": "shared_sgE_gene_hits",
            "value": int(len(shared_sge_genes)),
            "note": "Overlapping gene rows recovered for rescued primary shared sgE hotspots.",
        },
        {
            "metric": "candidate_shortlist_rows",
            "value": int(len(candidate_shortlist)),
            "note": "Candidate genes supported by one or more rescued hotspot overlaps.",
        },
    ])

    output_sheets = {
        OUTPUT_SHEETS["run_parameters"]: build_run_parameters(args),
        OUTPUT_SHEETS["input_inventory"]: build_input_inventory(args),
        OUTPUT_SHEETS["jsonl_inventory"]: jsonl_inventory,
        OUTPUT_SHEETS["annotation_readiness"]: readiness,
        OUTPUT_SHEETS["sgc_labels"]: sgc_label_dict,
        OUTPUT_SHEETS["sge_labels"]: sge_label_dict,
        OUTPUT_SHEETS["shared_sgc_map"]: rescued_shared_sgc,
        OUTPUT_SHEETS["shared_sge_map"]: rescued_shared_sge,
        OUTPUT_SHEETS["shared_sgc_genes"]: shared_sgc_genes,
        OUTPUT_SHEETS["shared_sge_genes"]: shared_sge_genes,
        OUTPUT_SHEETS["subgenome_map"]: rescued_subgenome,
        OUTPUT_SHEETS["subgenome_genes"]: rescued_subgenome_genes,
        OUTPUT_SHEETS["candidate_shortlist"]: candidate_shortlist,
        OUTPUT_SHEETS["shared_table"]: pd.concat(
            [prepare_gene_table(shared_sgc_genes), prepare_gene_table(shared_sge_genes)],
            ignore_index=True,
        ),
        OUTPUT_SHEETS["candidate_table"]: candidate_shortlist.copy(),
        OUTPUT_SHEETS["diagnostics"]: diagnostics,
    }

    supplementary_sheets = {
        "S5c_jsonl_inventory": jsonl_inventory,
        "S5c_sgC_crosswalk": sgc_crosswalk,
        "S5c_sgE_crosswalk": sge_crosswalk,
        "S5c_rescued_shared_sgC": rescued_shared_sgc,
        "S5c_rescued_shared_sgE": rescued_shared_sge,
        "S5c_shared_sgC_genes": shared_sgc_genes,
        "S5c_shared_sgE_genes": shared_sge_genes,
        "S5c_rescued_subgenome": rescued_subgenome,
        "S5c_subgenome_genes": rescued_subgenome_genes,
        "S5c_candidate_shortlist": candidate_shortlist,
        "S5c_diagnostics": diagnostics,
    }

    print(f"Writing workbook: {args.output_workbook}")
    write_workbook(args.output_workbook, output_sheets)

    print(f"Writing workbook: {args.supplementary_output_workbook}")
    write_workbook(args.supplementary_output_workbook, supplementary_sheets)

    print("Part 5c label rescue finished successfully.")


if __name__ == "__main__":
    main()
