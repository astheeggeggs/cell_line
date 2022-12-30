# Cell line PRS evaluation pipeline

The file structure of the resultant data (and input data) on the cluster is shown below.

```
PRS_cell_data
   |-data
   |---1000G
   |-----B37
   |-----B38
   |---PGS_score_files
   |-----pgs-catalog-20221123-hg19
   |-------scores
   |-----pgs-catalog-20221123-hg38
   |-------scores
   |---PGS_score_results
   |-----imputed
   |---TOPMed_variants
   |---UKB_imputed_subset
   |---celldataB36
   |---celldataB37
   |-----imputed
   |-------HRC
   |-------TOPMED
   |---celldataB38
   |-----imputed
   |-------TOPMED
   |---------B38
   |---------lifted_from_B37
   |---for_imputation
```

The scripts in this repository are located in main `PRS_cell_data` directory. Here's a brief description of what each of the scripts does, the required inputs, and where they were grabbed from.

## celldata_qc.r
The first step is to clean-up the cell-line data
## ukb_geno_qc.r
## ukb_imputation_preparation.sh
## ukb_geno_qc_combined.r
## ukb_imputation_preparation_combined.sh
## cell_line_imputation_preparation.sh
## grab_PGS_scorings.r
## pgs_1000G_submission.sh
## filter_cell_lines.sh
## pgs_cell_line_submission.sh
## aligning_1000G_and_cell_line.sh
## UKB_subset_vcf_creation.sh
## plot_cell_line_prs.r
