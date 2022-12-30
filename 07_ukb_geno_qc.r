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

fwrite(UKB_subset %>% mutate(fam.ID=sample.ID) %>% select(sample.ID, fam.ID), sep=" ",
	file=UKB_subset_filename, col.names=FALSE, quote=FALSE)
