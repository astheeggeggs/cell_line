library(data.table)
library(dplyr)

# module purge
# module load BCFtools
# module load R

# Sanity checking the variant counts that go into the PRS evalution
imputed_path <- "/well/lindgren/UKBIOBANK/dpalmer/PRS_cell_data/data/celldataB37/HRC/imputed"
ukb_imputed_subset_dir <- "/well/lindgren/UKBIOBANK/dpalmer/PRS_cell_data/data/UKB_imputed_subset"

# Check that all UKB imputed variants restricted to are in the cell line data
for (chr in c(seq(1, 22), "X")) {
	system(paste0("bcftools query -f '%CHROM:%POS:%REF:%ALT{0}\n' ", imputed_path, "/chr", chr, ".dose.vcf.gz > ", imputed_path, "/tmp_", chr))
	dt_cell_line <- fread(paste0(imputed_path, "/tmp_", chr), header=FALSE, key="V1")
	if (chr != "X") {
		dt_ukb_imputed <- fread(paste0(ukb_imputed_subset_dir, "/UKB_imputed_subset_chr", chr, "_thresholded_recode.bim")) %>% mutate(varid=paste0(V1,":",V4, ":", V6, ":", V5))
	} else {
		dt_ukb_imputed <- fread(paste0(ukb_imputed_subset_dir, "/UKB_imputed_subset_chr", chr, "_thresholded_recode.bim")) %>% mutate(varid=paste0("X:", V4, ":", V6, ":", V5))
	}
	print(all(dt_ukb_imputed$varid %in% dt_cell_line$V1))
}

# Check that the number of variants present in the 'final' UKB vcf is the sum of the variants in the imputed data and the array data, minus the variants that were not phased.
for (chr in c(seq(1,22), "X")) {
	imputed <- fread(paste0(ukb_imputed_subset_dir, "/UKB_imputed_subset_chr", chr, "_thresholded_recode.bim"), header=FALSE)
	array <- fread(paste0(ukb_imputed_subset_dir, "/UKB_subset_combined-updated-chr", chr, "-no-missing-no-unphased.bim"), header=FALSE) 
	array_and_imputed <- fread(paste0(ukb_imputed_subset_dir, "/UKB_imputed_subset_chr", chr, "_thresholded_recode_merge_typed.bim"), header=FALSE)
	print(length(intersect((imputed %>% transmute(varid=paste(V1, V4, V6, V5, sep=":")))$varid, array$V2)) == 0)
}
