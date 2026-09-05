# Mucosal inflammation and biologic response in ulcerative colitis

This repository contains the reproducible R workflow supporting the manuscript, **"Baseline Mucosal Inflammation Marks a Shared Adverse Prognostic State Across Anti-TNF, Ustekinumab and Vedolizumab in Ulcerative Colitis."**

The analysis asks whether a high baseline mucosal inflammatory state is associated with lower subsequent response across anti-TNF, ustekinumab and vedolizumab cohorts. Randomized-trial treatment-by-score interactions are used to distinguish a shared prognostic association from a drug-selection biomarker. Longitudinal analyses assess whether successful treatment is accompanied by a decline in the state. Single-cell analyses are restricted to state localization and do not independently validate clinical prediction.

## What is included

- `scripts/`: checkpointed R workflow.
- `config/`: immutable input manifest and frozen gene sets.
- `docs/`: analysis contract, data dictionary and reproducibility checklist.
- `environment/`: package-version snapshot and environment notes.
- `results/reproducibility/`: validation report and record of the last pipeline run.
- `results/source_data/`: machine-readable source data supporting manuscript figures.
- `results/supplementary_tables/TableS9_platform_mapping_and_score_standardization.csv`: revision Table S9 with accession-level platform mapping, signature coverage, and score-standardization parameters.
- `results/tables/common_state/table_covariate_availability_matrix.csv`: revision Supplementary Audit 1 for covariate availability and adjustment eligibility.
- `results/tables/common_state/table_singlecell_dataset_level_audit.csv` and `table_singlecell_longitudinal_pairing_audit.csv`: revision Supplementary Audit 2 inputs for sample eligibility, donor overlap, and longitudinal pairing.
- `results/tables/common_state/table_sample_exclusion_audit.csv`: revision exclusion audit with sample-level eligibility and exclusion reasons.

Raw datasets and large intermediate files are deliberately excluded. They are publicly available from the sources recorded in `config/data_manifest.csv` and are downloaded and checksum-verified by the workflow.

## Reproduce the analysis

```sh
Rscript scripts/install_dependencies.R
Rscript scripts/run_all.R
```

The first complete run downloads approximately 8 GB of public inputs. Allow at least 20 GB of disk space and 16 GB of RAM. To inspect available stages, run `Rscript scripts/run_all.R --list-stages`. For an existing download-only offline run, use `Rscript scripts/run_all.R --offline`.

## Interpretation boundary

The results support a shared adverse prognostic state and molecular reversal accompanying successful treatment. They do not establish a drug-selection biomarker, a clinical decision threshold, causal reversal, or independent single-cell response validation.

## Data sources

Bulk transcriptomic inputs are from NCBI GEO: GSE16879, GSE73661, GSE23597, GSE92415 and GSE206285. The single-cell input is GSE282122. Exact URLs, source roles and fixed MD5 checksums are listed in `config/data_manifest.csv`.

## Citation

If you use this workflow, please cite the associated manuscript and this repository. Citation metadata are provided in `CITATION.cff`. The revision-specific audit package is released as `v-1.1.0`; its corresponding versioned archive DOI will be added after Zenodo minting.
