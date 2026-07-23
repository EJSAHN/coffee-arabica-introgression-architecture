# coffee-arabica-introgression-architecture

Reproducible analysis workflow for strict subgenome filtering, SNP-panel sensitivity analysis, genome-wide cultivated–introgressed differentiation scans, direct chromosome 4 sequence alignment, and alignment-supported annotation extraction from public *Coffea arabica* genomic resources.

The repository contains the analytical workflow and machine-readable source data underlying the manuscript figures and supplementary tables. Manuscript-specific graphic layout and figure-rendering code are intentionally maintained separately from the analysis repository.

## Analysis overview

The workflow addresses four reproducibility questions:

1. Are population structure and accession rankings stable across SNP-panel sizes and random seeds?
2. Does inclusion of diploid progenitor references alter within-Arabica ordination geometry?
3. Does chromosome 4 emerge from a genome-wide cultivated–introgressed scan after strict subgenome filtering?
4. Do the source intervals align to the expected ET-39 chromosome 4 pseudomolecules, and which annotations occur within directly supported blocks?

Both public VCFs use a combined Arabica reference. Every analysis therefore filters records by explicit sgC/sgE pseudomolecule label before PCA, ranking, interval, or permutation analysis.

## Repository structure

```text
scripts/
  01_subgenome_filtered_sensitivity.py
  02_direct_interval_alignment.py
  03_finalize_alignment_annotations.py
  04_export_supplementary_outputs.py
  analysis_common.py
  native_stack_preflight.py
  run_and_tee.py

config/
  analysis_intervals.csv
  prioritized_accessions.csv
  accession_context.csv

outputs/figure_source_data/
  machine-readable CSV files underlying all manuscript figures

legacy/submitted_workflow/
  scripts preserved from the initial submission for provenance
```

## Public inputs

Large public datasets are not redistributed. See `docs/DATA_SOURCES.md` and `data/raw_public_inputs/README.md` for the source repositories and expected local inputs. The direct-alignment script downloads source assembly `GCA_036785775.1` when it is not already cached.

## Installation

Python 3.10 or later is required.

```bash
python -m pip install -r requirements.txt
```

NCBI BLAST+ is an external dependency for direct interval alignment.

- Linux/macOS with Conda: `conda install -c bioconda blast`
- Debian/Ubuntu: `apt-get install ncbi-blast+`
- macOS with Homebrew: `brew install blast`
- Windows: the script first uses BLAST+ on `PATH`; when absent, it can download and cache the official Win64 release.

Explicit executable paths may be supplied with `--blastn-path` and `--makeblastdb-path`.

## Platform-neutral workflow

### 1. Strict subgenome-filtered sensitivity analysis

```bash
python scripts/01_subgenome_filtered_sensitivity.py \
  --project-dir "<PROJECT_DIR>" \
  --output-dir "<SENSITIVITY_OUTPUT>" \
  --panel-sizes "6000,12000,24000,48000" \
  --seeds "20250416,20250417,20250418,20250419,20250420" \
  --submitted-panel-size 12000 \
  --submitted-seed 20250416 \
  --n-permutations 500 \
  --contig-filter-mode matching_pseudomolecules_only
```

### 2. Direct chromosome 4 alignment

```bash
python scripts/02_direct_interval_alignment.py \
  --project-dir "<PROJECT_DIR>" \
  --output-dir "<ALIGNMENT_OUTPUT>" \
  --intervals-csv config/analysis_intervals.csv \
  --threads 4
```

### 3. Alignment-supported annotations

```bash
python scripts/03_finalize_alignment_annotations.py \
  --project-dir "<PROJECT_DIR>" \
  --alignment-dir "<ALIGNMENT_OUTPUT>" \
  --output-dir "<ALIGNMENT_OUTPUT>/final"
```

### 4. Supplementary workbooks and source-data export

```bash
python scripts/04_export_supplementary_outputs.py \
  --sensitivity-workbook "<SENSITIVITY_OUTPUT>/sampling_population_structure_validation.xlsx" \
  --alignment-workbook "<ALIGNMENT_OUTPUT>/final/alignment_validation_corrected.xlsx" \
  --priority-accessions config/prioritized_accessions.csv \
  --accession-context config/accession_context.csv \
  --output-dir "<EXPORT_OUTPUT>"
```

The export step creates Supplementary Tables, Supplementary Data S1, an analysis summary, and CSV source data underlying the manuscript figures. It does not render PDF or PNG figures.

## Windows convenience launchers

The files under `launch/` are optional wrappers. Supply the project directory as the first argument; if omitted, they use `data\raw_public_inputs` within the repository.

```bat
launch\RUN_FULL_SUBGENOME_VALIDATION.cmd "<PROJECT_DIR>"
launch\RUN_DIRECT_INTERVAL_ALIGNMENT.cmd "<PROJECT_DIR>"
launch\RUN_EXPORT_SUPPLEMENTARY_OUTPUTS.cmd "<PROJECT_DIR>"
```

## Interpretation boundaries

- sgC contains the stronger chromosome 4 differentiation signal.
- sgE contains a weaker chromosome 4-associated signal in the same accession contrast; it is not interpreted as equivalent alien introgression.
- Chromosome-level correspondence is alignment-supported, but exact base-pair interval boundaries are not treated as resolved.
- Candidate annotations are hypotheses for downstream validation, not causal assignments.
- Arabica-only PCA is used for within-Arabica inference; full-panel PCA is retained as progenitor context.

## Provenance

The scripts used for the initial manuscript submission are preserved unchanged under `legacy/submitted_workflow/`. They are retained for auditability and are not recommended for the corrected analysis. See `CHANGELOG.md` and `docs/PROVENANCE.md`.
