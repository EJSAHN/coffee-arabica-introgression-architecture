# coffee-arabica-introgression-architecture

Subgenome-resolved Python workflow for reanalyzing public *Coffea arabica* genomics resources to prioritize introgressed accessions, shared chromosome 4 hotspots, and annotation-supported candidate gene space.

## Overview

This repository contains the scripted analysis workflow used to derive the accession-level, hotspot-level, and annotation-supported evidence tables for a public reanalysis of modern *Coffea arabica* introgression architecture.

The workflow operates on public resources and does not include large input datasets. Users should download the required input files from the original repositories listed in `docs/DATA_SOURCES.md`.

## Required input resources

Place the Dryad files in one project directory:

- `Accession_info.xlsx`
- `Arabica_sgC.TIP.BB.vcf.gz`
- `Arabica_sgE.TIP.BB.vcf.gz`
- `Coffea_syntenic_alignments.tar.gz`
- `README.md`

For annotation projection, also place or reference:

- `C_arabica_ET39_sgC.gff3`
- `C_arabica_ET39_sgE.gff3`
- `assembly_data_report.jsonl` from NCBI Datasets for `GCA_036785775.1`

## Installation

```bash
conda create -n coffee_introgression python=3.10 pandas numpy openpyxl
conda activate coffee_introgression
```

or install from the requirements file:

```bash
pip install -r requirements.txt
```

## Workflow

Use `<PROJECT_DIR>` as the directory containing public input files and `<OUTPUT_DIR>` as the directory for derived outputs.

### 1. Prepare analysis-ready variant panels

```bash
python scripts/01_prepare_variant_panel.py ^
  --project-dir "<PROJECT_DIR>" ^
  --output-dir "<OUTPUT_DIR>"
```

### 2. Summarize population structure and accession affinities

```bash
python scripts/02_population_structure.py ^
  --analysis-workbook "<OUTPUT_DIR>/coffee_introgression_manuscript_analysis.xlsx" ^
  --supplementary-workbook "<OUTPUT_DIR>/coffee_introgression_supplementary_data_s1.xlsx" ^
  --output-workbook "<OUTPUT_DIR>/coffee_introgression_part2_manuscript.xlsx" ^
  --supplementary-output-workbook "<OUTPUT_DIR>/coffee_introgression_part2_supplementary_data_s2.xlsx"
```

### 3. Prioritize introgression-associated evidence

```bash
python scripts/03_prioritize_introgression_signals.py ^
  --analysis-workbook "<OUTPUT_DIR>/coffee_introgression_part2_manuscript.xlsx" ^
  --supplementary-workbook "<OUTPUT_DIR>/coffee_introgression_part2_supplementary_data_s2.xlsx" ^
  --output-workbook "<OUTPUT_DIR>/coffee_introgression_part3_manuscript.xlsx" ^
  --supplementary-output-workbook "<OUTPUT_DIR>/coffee_introgression_part3_supplementary_data_s3.xlsx"
```

### 4. Rank accessions and hotspot intervals

```bash
python scripts/04_rank_accessions_and_hotspots.py ^
  --analysis-workbook "<OUTPUT_DIR>/coffee_introgression_part3_manuscript.xlsx" ^
  --supplementary-workbook "<OUTPUT_DIR>/coffee_introgression_part3_supplementary_data_s3.xlsx" ^
  --output-workbook "<OUTPUT_DIR>/coffee_introgression_part4_manuscript.xlsx" ^
  --supplementary-output-workbook "<OUTPUT_DIR>/coffee_introgression_part4_supplementary_data_s4.xlsx"
```

### 5. Project prioritized hotspots into annotation space

```bash
python scripts/05_project_hotspots_to_annotation.py ^
  --part1-workbook "<OUTPUT_DIR>/coffee_introgression_manuscript_analysis.xlsx" ^
  --part4-workbook "<OUTPUT_DIR>/coffee_introgression_part4_manuscript.xlsx" ^
  --sgc-gff3 "<PROJECT_DIR>/C_arabica_ET39_sgC.gff3" ^
  --sge-gff3 "<PROJECT_DIR>/C_arabica_ET39_sgE.gff3" ^
  --output-workbook "<OUTPUT_DIR>/coffee_introgression_part5_manuscript.xlsx" ^
  --supplementary-output-workbook "<OUTPUT_DIR>/coffee_introgression_part5_supplementary_data_s5.xlsx"
```

### 6. Refine annotation mapping by chromosome-label reconciliation

```bash
python scripts/06_refine_annotation_mapping.py ^
  --part5-workbook "<OUTPUT_DIR>/coffee_introgression_part5_manuscript.xlsx" ^
  --supplementary-workbook "<OUTPUT_DIR>/coffee_introgression_part5_supplementary_data_s5.xlsx" ^
  --sgc-gff3 "<PROJECT_DIR>/C_arabica_ET39_sgC.gff3" ^
  --sge-gff3 "<PROJECT_DIR>/C_arabica_ET39_sgE.gff3" ^
  --assembly-data-report-jsonl "<PROJECT_DIR>/ncbi_dataset/data/assembly_data_report.jsonl" ^
  --output-workbook "<OUTPUT_DIR>/coffee_introgression_part5c_manuscript.xlsx" ^
  --supplementary-output-workbook "<OUTPUT_DIR>/coffee_introgression_part5c_supplementary_data_s5c.xlsx"
```

## Notes

- Large public input files are not included in this repository.
- Derived outputs should be written outside the repository or to `outputs/`, which is ignored by Git.
- The scripts are ordered by execution sequence and retain the internal workbook dependencies required for reproducibility.
