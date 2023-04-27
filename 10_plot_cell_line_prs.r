library(ggplot2)
library(data.table)
library(dplyr)
library(jsonlite)

Sys.setlocale(category = "LC_ALL", locale = "en_US.UTF-8")

# Grab all of the score files and merge them
reference_scores <- grep("scores.txt", dir("/well/lindgren/UKBIOBANK/dpalmer/PRS_cell_data/data/PGS_score_results", recursive=TRUE, full.names=TRUE), value=TRUE)
cell_line_scores <- grep("cell_line", reference_scores, value=TRUE)
ukb_scores <- grep("ukb", reference_scores, value=TRUE)
reference_scores <- setdiff(reference_scores, c(ukb_scores, cell_line_scores))

dt_ref <- fread(reference_scores[1], key="sample")
for (ref_score in reference_scores[-1]) {
	dt_ref <- merge(dt_ref, fread(ref_score, key="sample"))
}

dt_cell_line <- fread(cell_line_scores[1], key="sample")
for (cell_score in cell_line_scores[-1]) {
	dt_cell_line <- merge(dt_cell_line, fread(cell_score, key="sample"))
}

dt_ukb <- fread(ukb_scores[1], key="sample")
for (ukb in ukb_scores[-1]) {
	dt_ukb <- merge(dt_ukb, fread(ukb, key="sample"))
}
dt_ukb <- dt_ukb %>% mutate(sample = gsub("_.*", "", sample))

# Merge in the superpopulation information
dt_cell_line_pop <- fread("/well/lindgren/UKBIOBANK/dpalmer/PRS_cell_data/superpopulation_assignments/superpopulation_labels.tsv") %>% rename(sample=sample.ID)
dt_cell_line <- dt_cell_line %>% mutate(sample = gsub("Genotype_.*", "Genotype", sample))
setkey(dt_cell_line_pop, "sample")
setkey(dt_cell_line, "sample")
dt_cell_line <- merge(dt_cell_line_pop, dt_cell_line)
dt_cell_line <- dt_cell_line %>% filter(classification_strict != "unsure")

# Merge this information with the 1000G population information
# system("wget -P /well/lindgren/UKBIOBANK/dpalmer/PRS_cell_data/data/1000G/ ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/working/20130606_sample_info/20130606_g1k.ped")
# system("wget -P /well/lindgren/UKBIOBANK/dpalmer/PRS_cell_data/data/1000G/ ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/20131219.populations.tsv")

dt_pop <- fread("/well/lindgren/UKBIOBANK/dpalmer/PRS_cell_data/data/1000G/20131219.populations.tsv")
dt_sample_pop <- fread("/well/lindgren/UKBIOBANK/dpalmer/PRS_cell_data/data/1000G/20130606_g1k.ped")

dt_sample_pop <- left_join(dt_sample_pop, dt_pop[, 1:3], by = c("Population" = "Population Code"))
dt_ref <- dt_ref %>% mutate(sample=gsub("0_", "", sample))
dt_ref <- merge(dt_sample_pop %>% rename(sample = `Individual ID`, super_population=`Super Population`), dt_ref, by="sample")

meta <- fromJSON("/well/lindgren/UKBIOBANK/dpalmer/PRS_cell_data/data/PGS_score_files/pgs-catalog-20221123-hg19/scores.meta.json")
pdf("/well/lindgren-ukbb/projects/ukbb-11867/dpalmer/PRS_cell_data/data/PGS_score_results/PRS_plots.pdf", width=6, height=4)
i <- 1
for (PRS in grep("PGS", sort(intersect(names(dt_ref), names(dt_cell_line))), value=TRUE)) {
	if ((i %% 100) == 0) print(i)
	p <- ggplot() + 
	geom_density(data = dt_ref, aes_string(x = PRS, col="super_population")) +
	geom_rug(data = dt_cell_line, aes_string(x = PRS, col="classification_strict")) + theme_classic() + theme(legend.title=element_blank())# + theme(legend.position="none")
	p <- p + ggtitle(paste(meta[[PRS]]$trait)) + #theme(plot.title = element_textbox_simple()) +
	labs(caption = paste(
		meta[[PRS]]$link,
		paste0(meta[[PRS]]$publication$firstauthor, " et al., ", gsub("-.*", "", meta[[PRS]]$publication$date), ", ", meta[[PRS]]$publication$journal),
		meta[[PRS]]$publication$doi,
		paste(meta[[PRS]]$variants, "variants,", meta[[PRS]]$samples, "samples"), sep="\n"))
	print(p)
	i <- i+1
}
dev.off()

dt <- rbind(
	dt_ref %>% mutate(source="1000 Genomes"),
	dt_cell_line %>% mutate(super_population=classification_strict, source="Cell line"),
	fill=TRUE)

pdf("/well/lindgren-ukbb/projects/ukbb-11867/dpalmer/PRS_cell_data/data/PGS_score_results/PRS_plots_boxplots.pdf", width=7, height=4)
i <- 1
for (PRS in grep("PGS", sort(intersect(names(dt_ref), names(dt_cell_line))), value=TRUE)) {
	if ((i %% 100) == 0) print(i)

	p <- ggplot(data=dt, aes_string(x="super_population", y=PRS, fill="source")) + 
	geom_boxplot() +
	theme_classic()
	p <- p + ggtitle(paste(meta[[PRS]]$trait)) + theme(legend.title=element_blank())#theme(plot.title = element_textbox_simple()) +
	p <- p + xlab("Super Population") +
	labs(caption = paste(
		meta[[PRS]]$link,
		paste0(meta[[PRS]]$publication$firstauthor, " et al., ", gsub("-.*", "", meta[[PRS]]$publication$date), ", ", meta[[PRS]]$publication$journal),
		meta[[PRS]]$publication$doi,
		paste(meta[[PRS]]$variants, "variants,", meta[[PRS]]$samples, "samples"), sep="\n"))
	print(p)
	i <- i+1
}
dev.off()


# # Extract the population labels for the UK Biobank data
# pops <- c("EUR", "EAS", "SAS", "AMR", "AFR")
# dt_loc <- paste0("/well/lindgren/UKBIOBANK/dpalmer/PRS_cell_data/data/", pops, "/imputed/scores.txt")
# dt_UKB <- rbind(
# 	fread(dt_loc[1]) %>% mutate(super_population=pops[1]),
# 	fread(dt_loc[2]) %>% mutate(super_population=pops[2]),
# 	fread(dt_loc[3]) %>% mutate(super_population=pops[3]),
# 	fread(dt_loc[4]) %>% mutate(super_population=pops[4]),
# 	fread(dt_loc[5]) %>% mutate(super_population=pops[5])
# 	)
# setkey(dt_UKB, "sample")

ukb_subset <- fread("/well/lindgren-ukbb/projects/ukbb-11867/dpalmer/PRS_cell_data/data/UKB_subset.txt")
dt_ukb_pop <- fread("/well/lindgren/UKBIOBANK/dpalmer/superpopulation_assignments/superpopulation_labels.tsv") %>% filter(sample.ID > 0)
ukb_subset <- ukb_subset %>% rename(sample.ID=V1)
ukb_subset <- ukb_subset %>% mutate(sample.ID=as.character(sample.ID))
setkey(ukb_subset, "sample.ID")
setkey(dt_ukb_pop, "sample.ID")
dt_ukb_pop <- merge(dt_ukb_pop, ukb_subset)
dt_ukb <- merge(dt_ukb_pop %>% rename(sample=sample.ID, super_population=classification_strict) %>% setkey("sample"), dt_ukb %>% setkey("sample"))
# Now, read in and merge all of the PGS for UK Biobank samples

# Now, read in all the UK Biobank information and create violin plots
dt <- rbind(
	dt_ref %>% mutate(source="1000 Genomes"),
	dt_ukb %>% mutate(source="UK Biobank"),
	dt_cell_line %>% mutate(super_population=classification_strict, source="Cell line"),
	fill=TRUE)
fwrite(dt, file="/well/lindgren-ukbb/projects/ukbb-11867/dpalmer/PRS_cell_data/data/PGS_score_results/combined_results.tsv.gz", sep="\t")

pdf("/well/lindgren-ukbb/projects/ukbb-11867/dpalmer/PRS_cell_data/data/PGS_score_results/PRS_plots_boxplots_and_UKB.pdf", width=7, height=4)
i <- 1
for (PRS in grep("PGS", sort(intersect(names(dt_ukb), names(dt_cell_line))), value=TRUE)) {
	if ((i %% 100) == 0) print(i)
	p <- ggplot(data=dt, aes_string(x="super_population", y=PRS, fill="source")) + 
	# geom_violin(trim=FALSE) + 
	geom_boxplot() +
	theme_classic() #+ theme(legend.position="none")
	p <- p + ggtitle(paste(meta[[PRS]]$trait)) + theme(legend.title=element_blank())#theme(plot.title = element_textbox_simple()) +
	p <- p + xlab("Super Population") +
	labs(caption = paste(
		meta[[PRS]]$link,
		paste0(meta[[PRS]]$publication$firstauthor, " et al., ", gsub("-.*", "", meta[[PRS]]$publication$date), ", ", meta[[PRS]]$publication$journal),
		meta[[PRS]]$publication$doi,
		paste(meta[[PRS]]$variants, "variants,", meta[[PRS]]$samples, "samples"), sep="\n"))
	print(p)
	i <- i+1
}
dev.off()
