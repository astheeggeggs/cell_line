library(ggplot2)
library(data.table)
library(dplyr)
library(jsonlite)

Sys.setlocale(category = "LC_ALL", locale = "en_US.UTF-8")

# Grab all of the score files and merge them
reference_scores <- grep("scores.txt",
	dir("/well/lindgren/UKBIOBANK/dpalmer/PRS_cell_data/data/PGS_score_results",
		recursive=TRUE, full.names=TRUE), value=TRUE)

# Get all the data together
# First: the UK Biobank data
ukb_reference_scores <- grep("ukb_job_", reference_scores, value=TRUE)
# Intersected with the cell line and 1000G
ukb_reference_scores_fully_intersected <- grep("intersect", ukb_reference_scores, value=TRUE)
# Just intersected with the cell line
ukb_reference_scores_intersected <- setdiff(ukb_reference_scores, ukb_reference_scores_fully_intersected)

# Sanity checking - check that the number of samples, and the number of PGS is the same in all versions
root <- "/well/lindgren/UKBIOBANK/dpalmer/PRS_cell_data/data/PGS_score_results/ukb_job_"
for (i in 1:33) {
	fully_intersected <- paste0(root, i, "/ukb_intersect_1000G_cell_line_job_", i, "_scores.txt")
	intersected <- paste0(root, i, "/ukb_job_", i, "_scores.txt")
	if (file.exists(fully_intersected) & file.exists(intersected)) {
		dt1 <- fread(fully_intersected)
		dt2 <- fread(intersected)
		if (all(dim(dt1) == dim(dt2))) {
			cat("Dimensions match...")
			if (all(names(dt1) == names(dt2))) {
				cat("PGS names match...")
				for (PGS in names(dt1)[2:length(names(dt1))]) {
					print(cor(dt1[[PGS]], dt2[[PGS]]))
				}
			} else {
				cat("\nPGS names don't match!")
			}
		} else {
			cat("\nDimensions don't match!")
		}
	} else {
		cat(paste(fully_intersected, ifelse(file.exists(fully_intersected), "exists\n", "does NOT exist\n")))
		cat(paste(intersected, ifelse(file.exists(intersected), "exists\n", "does NOT exist\n")))
	}
}

# Second: the 1000G data
reference_scores_1000G <- grep("/job_[0-9]+/", reference_scores, value=TRUE)
# Intersected with the cell line and 1000G
reference_scores_fully_intersected_1000G <- grep("intersect", reference_scores_1000G, value=TRUE)
# Just intersected with the cell line
reference_scores_intersected_1000G <- setdiff(reference_scores_1000G, reference_scores_fully_intersected_1000G)

# Sanity checking - check that the number of samples, and the number of PGS is the same in all versions
root <- "/well/lindgren/UKBIOBANK/dpalmer/PRS_cell_data/data/PGS_score_results/job_"
for (i in 1:33) {
	fully_intersected <- paste0(root, i, "/1000G_intersect_ukb_cell_line_job_", i, "_scores.txt")
	intersected <- paste0(root, i, "/job_", i, "_scores.txt")
	if (file.exists(fully_intersected) & file.exists(intersected)) {
		dt1 <- fread(fully_intersected)
		dt2 <- fread(intersected)
		if (all(dim(dt1) == dim(dt2))) {
			cat("Dimensions match...")
			if (all(names(dt1) == names(dt2))) {
				cat("PGS names match...")
				for (PGS in names(dt1)[2:length(names(dt1))]) {
					print(cor(dt1[[PGS]], dt2[[PGS]]))
				}
			} else {
				cat("\nPGS names don't match!")
			}
		} else {
			cat("\nDimensions don't match!")
		}
	} else {
		cat(paste(fully_intersected, ifelse(file.exists(fully_intersected), "exists\n", "does NOT exist\n")))
		cat(paste(intersected, ifelse(file.exists(intersected), "exists\n", "does NOT exist\n")))
	}
}

# Finally: the cell line data
reference_scores_cell_line <- grep("/cell_line_job", reference_scores, value=TRUE)
# Intersected with 1000G and UKB
reference_scores_cell_line_fully_intersected <- grep("cell_line_intersect_UKB_1000G_", reference_scores_cell_line, value=TRUE)
# Intersected with UKB
reference_scores_cell_line_ukb_intersected <- setdiff(grep("cell_line_intersect_UKB_", reference_scores_cell_line, value=TRUE), cell_line_fully_intersected)
# Intersected with 1000G
reference_scores_cell_line_1000G_intersected <- grep("cell_line_job_[0-9]+/cell_line_job_", reference_scores_cell_line, value=TRUE)

root <- "/well/lindgren/UKBIOBANK/dpalmer/PRS_cell_data/data/PGS_score_results/cell_line_job_"
for (i in 1:33) {
	dt1 <- fread(reference_scores_cell_line_fully_intersected[i])
	dt2 <- fread(reference_scores_cell_line_1000G_intersected[i])
	dt3 <- fread(reference_scores_cell_line_ukb_intersected[i])
	if (all(dim(dt1) == dim(dt2)) & all(dim(dt2) == dim(dt3))) {
		cat("Dimensions all match...")
		if (all(names(dt1) == names(dt2)) & all(names(dt1) == names(dt3))) {
			cat("PGS names all match...")
			for (PGS in names(dt1)[2:length(names(dt1))]) {
				print(cor(data.table(PGS1 = dt1[[PGS]], PGS2 = dt2[[PGS]], PGS3 = dt3[[PGS]])))
			}
		} else {
			cat("\nPGS names don't match!")
		}
	} else {
		cat("\nDimensions don't match!")
	}
}

create_interleaved_hist <- function(dt)
{
	p <- ggplot(dt, aes(x=coverage, color=type)) +
	geom_histogram(fill="white", position="dodge") +
	theme(legend.position="top")
	return(p)
}

check_report_agreement <- function(string1, string2, plot=FALSE)
{
	for (i in 1:33) {
		cat(paste("job", i, "\n"))
		path1 <- gsub("@", i, string1)
		path2 <- gsub("@", i, string2)
		if ((file.exists(path1)) & (file.exists(path2))) {
			# Check if md5sum is the same
			md51 <- strsplit(system(paste("md5sum", path1), intern=TRUE), split=" ")[[1]][1]
			md52 <- strsplit(system(paste("md5sum", path2), intern=TRUE), split=" ")[[1]][1]
			dt1 <- fread(path1)
			dt2 <- fread(path2)
			dt1 <- dt1 %>% mutate(type='1')
			dt2 <- dt2 %>% mutate(type='2')
			setkeyv(dt1, c("score", "trait", "trait_efo", "type", "coverage"))
			setkeyv(dt2, c("score", "trait", "trait_efo", "type", "coverage"))
			dt_hist <- merge(dt1, dt2, all=TRUE) %>% select(coverage, type)
			if (plot) {
				p_hist <- create_interleaved_hist(dt_hist)
				print(p_hist)
			}
			if (md51 == md52) {
				cat("Reports are identical, great!\n")
			} else {
				cat("There's a discrepancy in the reports\n")
				dt1 <- fread(path1)
				dt2 <- fread(path2)
				dt1 <- setkeyv(dt1, c("score", "trait", "trait_efo"))
				dt2 <- setkeyv(dt2, c("score", "trait", "trait_efo"))
				differences_1 <- dt1[which(dt1 != dt2, arr.ind=TRUE)[,1],]
				differences_2 <- dt2[which(dt1 != dt2, arr.ind=TRUE)[,1],]
				print(differences_1)
				print(differences_2)
				cat("smallest PGS (by number of SNPs) with a difference")
				print(differences_1 %>% filter(variants_used == min(differences_1$variants_used)))
				print(differences_2 %>% filter(variants_used == min(differences_2$variants_used)))
			}
		}
	}
}

determine_variants_used_in_pgs <- function(pgs_variant_file, PGS_id) {
	dt <- fread(pgs_variant_file) %>% filter(score==PGS_id)
	dt <- dt %>% mutate(variant=paste0(chr_name, ":", chr_position))
	print(dt)
}

# Further checks - ensure that the same number of SNPs go into the scores (where they should). Check the reports align
# where there should be intersections.

# First, check agreement between:
# UKB and cell-line (single intersection)
cell_line_string <- "/well/lindgren/UKBIOBANK/dpalmer/PRS_cell_data/data/PGS_score_results/cell_line_job_@/cell_line_interesect_UKB_job_@_scores_report.csv"
ukb_string <- "/well/lindgren/UKBIOBANK/dpalmer/PRS_cell_data/data/PGS_score_results/ukb_job_@/ukb_job_@_scores_report.csv"

check_report_agreement(cell_line_string, ukb_string)

# Next, check agreement between 
# 1000G and cell-line (single intersection)
cell_line_string <- "/well/lindgren/UKBIOBANK/dpalmer/PRS_cell_data/data/PGS_score_results/cell_line_job_@/cell_line_job_@_scores_report.csv"
string_1000 <- "/well/lindgren/UKBIOBANK/dpalmer/PRS_cell_data/data/PGS_score_results/job_@/job_@_scores_report.csv"

check_report_agreement(cell_line_string, string_1000)
# In 1000G
pgs_variant_file <- "/well/lindgren/UKBIOBANK/dpalmer/PRS_cell_data/data/PGS_score_results/job_31/job_31_scores_variants_used.csv"
determine_variants_used_in_pgs(pgs_variant_file, "PGS003175")
# In cell line
pgs_variant_file <- "/well/lindgren/UKBIOBANK/dpalmer/PRS_cell_data/data/PGS_score_results/cell_line_job_31/cell_line_job_31_scores_variants_used.csv"
determine_variants_used_in_pgs(pgs_variant_file, "PGS003175")

# How to check for variants in indexed vcf
# system("bcftools view -rX:82139752-82139952 /well/lindgren/UKBIOBANK/dpalmer/PRS_cell_data/data/1000G/B37/all_phase3_r2_0.3_chrX.vcf.gz")

# Fully intersected
string_1000 <- "/well/lindgren/UKBIOBANK/dpalmer/PRS_cell_data/data/PGS_score_results/job_@/1000G_intersect_ukb_cell_line_job_@_scores_report.csv"
ukb_string <- "/well/lindgren/UKBIOBANK/dpalmer/PRS_cell_data/data/PGS_score_results/ukb_job_@/ukb_intersect_1000G_cell_line_job_@_scores_report.csv"
cell_line_string <- "/well/lindgren/UKBIOBANK/dpalmer/PRS_cell_data/data/PGS_score_results/cell_line_job_@/cell_line_intersect_UKB_1000G_job_@_scores_report.csv" 

check_report_agreement(cell_line_string, string_1000)
check_report_agreement(cell_line_string, ukb_string)
check_report_agreement(ukb_string, string_1000)

# Now, create a series of plots
# We can do this for 1000G first, because we know everything overlaps
refs <- reference_scores_intersected_1000G
dt_ref <- fread(refs[1], key="sample")
for (ref_score in refs[-1]) {
	dt_ref <- merge(dt_ref, fread(ref_score, key="sample"))
}

refs <- reference_scores_cell_line_1000G_intersected
dt_cell_line <- fread(refs[1], key="sample")
for (cell_score in refs[-1]) {
	dt_cell_line <- merge(dt_cell_line, fread(cell_score, key="sample"))
}

# Merge in the superpopulation information
dt_cell_line_pop <- fread("/well/lindgren/UKBIOBANK/dpalmer/PRS_cell_data/superpopulation_assignments/superpopulation_labels.tsv") %>% rename(sample=sample.ID)
dt_cell_line <- dt_cell_line %>% mutate(sample = gsub("Genotype_.*", "Genotype", sample))
setkey(dt_cell_line_pop, "sample")
setkey(dt_cell_line, "sample")
dt_cell_line <- merge(dt_cell_line_pop, dt_cell_line)

# Merge this information with the 1000G population information
# system("wget -P /well/lindgren/UKBIOBANK/dpalmer/PRS_cell_data/data/1000G/ ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/working/20130606_sample_info/20130606_g1k.ped")
# system("wget -P /well/lindgren/UKBIOBANK/dpalmer/PRS_cell_data/data/1000G/ ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/20131219.populations.tsv")

dt_pop <- fread("/well/lindgren/UKBIOBANK/dpalmer/PRS_cell_data/data/1000G/20131219.populations.tsv")
dt_sample_pop <- fread("/well/lindgren/UKBIOBANK/dpalmer/PRS_cell_data/data/1000G/20130606_g1k.ped")

dt_sample_pop <- left_join(dt_sample_pop, dt_pop[, 1:3], by = c("Population" = "Population Code"))
dt_ref <- dt_ref %>% mutate(sample=gsub("0_", "", sample))
dt_ref <- merge(dt_sample_pop %>% rename(sample = `Individual ID`, super_population=`Super Population`), dt_ref, by="sample")

dt <- rbind(
	dt_ref %>% mutate(source="1000 Genomes"),
	dt_cell_line %>% mutate(super_population=classification_strict, source="Cell line"),
	fill=TRUE)

system("mkdir /well/lindgren/UKBIOBANK/dpalmer/PRS_cell_data/final_output/")
fwrite(dt, file="/well/lindgren/UKBIOBANK/dpalmer/PRS_cell_data/final_output/PRS_1000G_cell_line.tsv.gz", sep='\t', quote=FALSE)
