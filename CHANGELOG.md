# Changelog

## Validated peer-review workflow

The current workflow adds the following validation and correction steps:

- strict filtering to the 11 explicitly labelled pseudomolecules of each intended subgenome;
- 6k/12k/24k/48k SNP-panel and five-seed sensitivity analyses;
- Arabica-only PCA for within-Arabica inference, with the full panel retained as progenitor context;
- genome-wide 1-Mb cultivated–introgressed scans with 500 label permutations;
- direct sequence alignment of chromosome 4 source intervals to ET-39 pseudomolecules;
- annotation summaries restricted to direct alignment-supported blocks;
- relative-path Windows launchers and platform-neutral BLAST discovery; and
- export of supplementary workbooks and machine-readable figure source data without manuscript figure rendering.

The original six-script workflow is preserved unchanged under `legacy/submitted_workflow/` for provenance. It should not be used for the corrected analyses.
