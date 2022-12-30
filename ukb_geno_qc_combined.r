library(data.table)
library(dplyr)

# Do it again, combining across the pops
pops <- c("EUR", "AFR", "EAS", "AMR", "SAS")
merge_list <- "/well/lindgren-ukbb/projects/ukbb-11867/dpalmer/PRS_cell_data/data/bfiles.txt"
fwrite(data.table(
    file = paste0("/well/lindgren-ukbb/projects/ukbb-11867/dpalmer/PRS_cell_data/data/", pops, "/UKB_subset_", pops, "_test")),
    file=merge_list, sep=" ", quote=FALSE, col.names=FALSE)

# Merge everything
system("mkdir data/combined")
system(paste("plink --merge-list", merge_list,
    "--out data/combined/UKB_subset_combined --make-bed"))