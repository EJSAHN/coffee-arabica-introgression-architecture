# Analysis provenance

## Initial submission

The original staged workflow is archived under `legacy/submitted_workflow/`. It generated the analyses available at the time of initial submission.

## Corrected analysis

Peer review prompted two material corrections. First, audit of the public VCF headers and records showed that both VCFs use a combined Arabica reference; current analyses retain only pseudomolecules explicitly labelled for the intended subgenome. Second, length-scaled annotation projection was replaced by direct sequence alignment and annotation extraction from alignment-supported target blocks.

The root workflow supersedes the archived workflow. The archive is retained to make the analytical history transparent, not as an alternative recommended analysis.

The repository distributes analytical code, supplementary-data export code, and machine-readable source data underlying all figures. Manuscript-specific graphic layout and rendering are maintained separately from the analysis workflow.
