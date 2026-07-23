# Analysis asset summary

## Inputs

- Sensitivity workbook: `sampling_population_structure_validation.xlsx`
- Alignment workbook: `alignment_validation_v44_corrected.xlsx`

## Principal validation results

- Overall sensitivity status: **PASS**
- Median PCA distance correlation across sensitivity conditions: **0.9981**
- Median accession-rank Spearman correlation: **0.9729**
- Median top-six overlap: **6/6**
- Median full-panel vs. Arabica-only distance correlation: **0.860**

## Genome-wide chromosome 4 evidence

- sgC top chromosome mean |ΔAF|: **0.3878**
- sgE top chromosome mean |ΔAF|: **0.0843**
- sgC permutation: observed **0.5738**, null 99th **0.3482**, empirical P <= **0.0020**
- sgE permutation: observed **0.2232**, null 99th **0.1656**, empirical P <= **0.0040**

## Direct alignment validation

- sgC: **PASS** chromosome assignment; **70.3%** aligned source coverage; **99.643%** weighted identity; boundary resolution **WARN**.
- sgE: **PASS** chromosome assignment; **70.4%** aligned source coverage; **99.653%** weighted identity; boundary resolution **WARN**.

## Interpretation constraints

1. Arabica-only PCA should be used for within-Arabica inference; full-panel PCA should be retained only as progenitor context.
2. sgC carries the stronger chromosome 4 differentiation signal. sgE retains a weaker chromosome 4-associated signal in the same accession contrast.
3. Exact interval endpoints should remain operational and supplementary. Chromosome-level correspondence is validated, but exact base-pair liftover boundaries are not resolved.
4. Candidate annotations are restricted to direct alignment-supported blocks and remain provisional.
5. Earlier accession-level subgenome-asymmetry candidates derived from mixed-contig analyses should not be retained without a new corrected definition.
