library(data.table)
library(dplyr)

# export PATH="/well/lindgren/dpalmer/:$PATH"

# UK Biobank genotype QC for phasing and imputation

# Filter down to the at collection of samples that are in EUR, SAS, EAS, AFR
# Use the superpopulation labelling file for that

dt_pop <- fread("/well/lindgren/UKBIOBANK/dpalmer/superpopulation_assignments/superpopulation_labels.tsv") %>% filter(sample.ID > 0)
sample_qc_file <- "/well/lindgren-ukbb/projects/ukbb-11867/DATA/QC/ukb_sqc_v2.txt"
UKB_fam <- "/well/lindgren/UKBIOBANK/DATA/SAMPLE_FAM/ukb11867_cal_chr1_v2_s488363.fam"
dt_sample_qc <- cbind(fread(UKB_fam) %>% rename(sample.ID=V1, sex=V5) %>% mutate(sample.ID=as.character(sample.ID)), fread(sample_qc_file))
dt_snp_qc <- fread("/well/lindgren-ukbb/projects/ukbb-11867/DATA/QC/ukb_snp_qc.txt")
setkey(dt_sample_qc, "sample.ID")
setkey(dt_pop, "sample.ID")
# Filtering to not in kinship table excludes any relateds up to third degree.
dt_pop <- merge(dt_pop, dt_sample_qc) 
dt_pop <- dt_pop %>% filter(in.kinship.table==0) %>% filter(Submitted.Gender==Inferred.Gender)

EUR_samples <- sample((dt_pop %>% filter(classification_strict == "EUR"))$sample.ID, 15000)
UKB_subset <- dt_pop %>% filter(classification_strict %in% c("AFR", "EAS", "SAS", "AMR"))
UKB_subset <- rbind(dt_pop %>% filter(dt_pop$sample.ID %in% EUR_samples), UKB_subset)

UKB_subset_filename <- "/well/lindgren-ukbb/projects/ukbb-11867/dpalmer/PRS_cell_data/data/UKB_subset.txt"
UKB_bed <- paste0("/well/lindgren/UKBIOBANK/DATA/CALLS/ukb_cal_chr",
	c(seq(1,22), "X", "Y", "MT"), "_v2.bed")
UKB_bim <- paste0("/well/lindgren/UKBIOBANK/DATA/CALLS/ukb_snp_chr",
	c(seq(1,22), "X", "Y", "MT"), "_v2.bim")
UKB_fam <- "/well/lindgren/UKBIOBANK/DATA/SAMPLE_FAM/ukb11867_cal_chr1_v2_s488363.fam"

fwrite(UKB_subset %>% mutate(fam.ID=sample.ID) %>% select(sample.ID, fam.ID), sep=" ",
	file=UKB_subset_filename, col.names=FALSE, quote=FALSE)

UKB_out <- paste0("/well/lindgren-ukbb/projects/ukbb-11867/dpalmer/PRS_cell_data/data/UKB_subset_chr",
	c(seq(1,22), "X", "Y", "MT"))
# Write the file, and the list of samples
for(chr in seq(1, 25)) {
	# SNP QC: call rate ≥ 0.95, filter to subset of samples
	system(paste("plink --bim", UKB_bim[chr], "--bed", UKB_bed[chr],
		"--fam", UKB_fam, "--keep", UKB_subset_filename,
		"--geno 0.05", "--out",  UKB_out[chr], "--make-bed"))
}

# Apply population specific HWE filters
dt_subset <- fread(UKB_subset_filename, header=FALSE) %>% rename(sample.ID=V1, fam.ID=V2) %>% mutate(sample.ID=as.character(sample.ID))
dt_pop_subset <- dt_pop %>% filter(sample.ID %in% dt_subset$sample.ID)
for (pop in c("EUR", "EAS", "SAS", "AFR", "AMR")) {
	fwrite(dt_pop_subset %>% filter(classification_strict == pop) %>% mutate(fam.ID=sample.ID) %>% select(sample.ID, fam.ID),
	sep=" ", file=paste0("/well/lindgren-ukbb/projects/ukbb-11867/dpalmer/PRS_cell_data/data/UKB_subset_", pop, ".txt"),
	col.names=FALSE, quote=FALSE)
}

# Autosomes
merge_list <- "/well/lindgren-ukbb/projects/ukbb-11867/dpalmer/PRS_cell_data/data/bfiles.txt"
fwrite(data.table(file = paste0("/well/lindgren-ukbb/projects/ukbb-11867/dpalmer/PRS_cell_data/data/UKB_subset_chr",
	seq(1,22))), file=merge_list, sep=" ", quote=FALSE, col.names=FALSE)

# Merge everything
# Sample QC: call rate ≥ 0.97, SNP QC: call rate ≥ 0.98
for (pop in c("EUR", "SAS", "EAS", "AFR", "AMR")) {
	system(paste("plink --merge-list", merge_list, "--keep",
	paste0("/well/lindgren-ukbb/projects/ukbb-11867/dpalmer/PRS_cell_data/data/UKB_subset_", pop, ".txt"),
	"--mind 0.03 --geno 0.02 --hwe 1e-10 --out", paste0("data/UKB_subset_", pop, "_aut"), "--make-bed"))
}

# That defines a set of samples to remove - use the output use keep-fam file.
# HWE filter for the X filter to the women - take the fam file and filter to the 2s
for (pop in c("EUR", "SAS", "EAS", "AFR", "AMR")) {
	fwrite(fread(paste0("data/UKB_subset_", pop, "_aut.fam")) %>% filter(V5 == 2),
		file=paste0("data/UKB_subset_", pop, "_X.fam"), sep=" ", quote=FALSE, col.names=FALSE)
	# Filter the X to the snplist
	system(paste("plink --bfile data/UKB_subset_chrX --keep-fam", paste0("data/UKB_subset_", pop, "_X.fam"),
	"--hwe 1e-10 --out", paste0("data/UKB_subset_", pop, "_X"), "--write-snplist"))
	# Now, filter the full set to those SNPs
	system(paste("plink --bfile data/UKB_subset_chrX --keep-fam", paste0("data/UKB_subset_", pop, "_aut.fam"),
	"--extract", paste0("data/UKB_subset_", pop, "_X.snplist"), "--out", paste0("data/UKB_subset_", pop, "_X"), "--make-bed"))
}

for (pop in c("EUR", "SAS", "AFR", "EAS", "AMR")) {
	bfile <- paste0(
		"/well/lindgren-ukbb/projects/ukbb-11867/dpalmer/PRS_cell_data/data/UKB_subset_",
		paste0(pop, c("_aut", "_X"))
		)
	bfiles <- data.table(files=bfile)
	merge_list <- "/well/lindgren-ukbb/projects/ukbb-11867/dpalmer/PRS_cell_data/data/bfiles.txt"
	fwrite(data.table(file = bfiles), file=merge_list, sep=" ", quote=FALSE, col.names=FALSE)
	system(paste("plink --merge-list", merge_list, "--out", paste0("data/UKB_subset_", pop, "_test"), "--make-bed"))
}

# Merge everything, and cleanup the plink-files split by chromosome
system("rm data/UKB_subset_chr*bed")
system("rm data/UKB_subset_chr*bim")
system("rm data/UKB_subset_chr*fam")
system("rm data/UKB_subset_chr*log")
system("rm data/UKB_subset_chr*nosex")
system("rm data/UKB_subset_*aut*")
system("rm data/UKB_subset_*X*")
system("rm data/bfile.txt")

