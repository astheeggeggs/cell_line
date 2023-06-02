library(data.table)
library(dplyr)

# export PATH="/well/lindgren/dpalmer/:$PATH"

# 1000 genomes population
dt_cell_pop <- fread(
  "/well/lindgren/UKBIOBANK/dpalmer/PRS_cell_data/superpopulation_assignments/superpopulation_labels.tsv") %>% 
  select(`sample.ID`, `classification_strict`, `super.population`)

# Cell-line ancestry assignment
directory <- "/well/lindgren/UKBIOBANK/dpalmer/PRS_cell_data/data/celldataB37/"
plink_cell <- paste0(directory, "celldataB37_dp_X")
system(paste("plink --bfile", paste0(directory, "celldataB37_dp"),
  "--out", plink_cell, "--geno 0.05 --make-bed"))
# Sample QC: call rate in cases or controls ≥ 0.97
system(paste("plink --bfile", plink_cell, "--out", plink_cell, "--mind 0.03 --make-bed"))
# SNP QC: call rate ≥ 0.98
system(paste("plink --bfile", plink_cell, "--chr 23 --out", plink_cell, "--geno 0.02 --make-bed"))

tmp_prune <- paste0(directory, "tmp")
plink_cell_pruned <- paste0(directory, "celldataB37_dp_X")
system(paste("plink", "--bfile", plink_cell, "--maf 0.01 --indep-pairwise 50 5 0.05 --out", tmp_prune))
system(paste("plink", "--bfile", plink_cell, "--extract", paste0(tmp_prune, ".prune.in"), "--make-bed --out", plink_cell_pruned))

# For ascertaining relatedness
plink_cell_auto <- paste0(directory, "celldataB37_dp_test-updated-autosomes-combined")
plink_cell_pruned_auto <- paste0(directory, "celldataB37_dp_test-updated-autosomes-combined_pruned")
system(paste("plink", "--bfile", plink_cell_auto, "--maf 0.01 --indep-pairwise 50 5 0.05 --out", tmp_prune))
system(paste("plink", "--bfile", plink_cell_auto, "--extract", paste0(tmp_prune, ".prune.in"), "--make-bed --out", plink_cell_pruned_auto))
