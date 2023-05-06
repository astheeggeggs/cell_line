library(data.table)
library(dplyr)

# Edge case! Some of the unphased variants would have had their alleles flipped
# by the HRC perl script - find those and add them to the list.

ukb_imputed_subset_dir <- "/well/lindgren/UKBIOBANK/dpalmer/PRS_cell_data/data/UKB_imputed_subset"
fread(paste0(ukb_imputed_subset_dir, "/array_variants.tsv"))
dt_array <- fread(paste0(ukb_imputed_subset_dir, "/array_variants.tsv"), header=FALSE, sep=":")
names(dt_array) <- c("chr", "pos", "ref_old", "alt_old")

# odd variants to fix
to_add <- c()
for (chr in seq(1,23)) {
	setkeyv(dt_array, c("chr", "pos"))
	dt_ukb_subset_array <- fread(paste0(ukb_imputed_subset_dir, "/UKB_subset_combined-updated-chr", chr, "-no-missing.bim"))
	names(dt_ukb_subset_array)[c(1,2,4,5,6)] <- c("chr", "varid_new", "pos", "alt_new", "ref_new")
	dt_ukb_subset_array <- dt_ukb_subset_array %>% mutate(varid_new = paste(chr, pos, ref_new, alt_new, sep=":"))
	setkeyv(dt_ukb_subset_array, c("chr", "pos"))
	dt_before <- merge(dt_ukb_subset_array, dt_array) %>% mutate(varid_old=paste(chr, pos, ref_old, alt_old, sep=":"))
	to_add <- c(to_add, (dt_before %>% filter(varid_old!=varid_new))$varid_new)
	print(to_add)
}

dt_unphased <- fread(paste0(ukb_imputed_subset_dir, "/array_variants.tsv"), header=FALSE)
dt_unphased <- rbind(dt_unphased, data.table(V1=to_add))
fwrite(dt_unphased, file=paste0(ukb_imputed_subset_dir, "/unphased_and_fixed_variants.tsv"))
