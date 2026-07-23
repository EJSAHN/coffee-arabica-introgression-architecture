# Output structure

## Subgenome-filtered sensitivity analysis

The full sensitivity workflow creates a timestamped directory containing:

- `sampling_population_structure_validation.xlsx`
- `sampling_validation_status.json`
- strict contig-filter audit files
- PCA scores and explained-variance tables
- genome-wide chromosome and 1-Mb-window scans
- permutation null distributions
- accession-ranking stability summaries

The public script writes tabular and JSON outputs only; it does not render manuscript figures.

## Direct interval alignment

The direct-alignment workflow creates:

- `submitted_interval_alignment_validation.xlsx`
- `submitted_interval_alignment_status.json`
- BLAST tabular outputs
- source interval chunk FASTA files

The final annotation step creates a `final/` subdirectory containing:

- `alignment_validation_corrected.xlsx`
- `alignment_validation_status.json`
- `ALIGNMENT_INTERPRETATION.md`

## Supplementary and source-data export

The export workflow consumes the completed sensitivity and alignment workbooks and creates:

- `tables/Supplementary_Tables.xlsx`
- `data/Supplementary_Data_S1.xlsx`
- `figure_source_data/*.csv`
- `ANALYSIS_SUMMARY.md`
- `EXPORT_MANIFEST.json`

The CSV files contain the machine-readable values underlying the manuscript and supplementary figures. Publication-specific PDF and PNG layout is outside the repository workflow.
