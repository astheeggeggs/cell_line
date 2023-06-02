library(data.table)
library(dplyr)
library(ggplot2)

# Determine the passed sex
# Determine the estimated sex, using F-stat
# Determine the estimated sex, using the number of Y calls

# IMPUTESEX_FILE <- '/well/lindgren/UKBIOBANK/dpalmer/PRS_cell_data/data/celldataB37/celldataB37_dp_impute_sex.tsv'
# Y_NCALLED <- '/well/lindgren/UKBIOBANK/dpalmer/PRS_cell_data/data/celldataB37/celldataB37_dp_y_called.tsv'
IMPUTESEX_FILE <- 'celldataB37_dp_impute_sex.tsv'
Y_NCALLED <- 'celldataB37_dp_y_called.tsv'

dt_y <- fread(Y_NCALLED)
dt_x <- fread(IMPUTESEX_FILE)

setkey(dt_x, "s")
setkey(dt_y, "s")

dt <- merge(dt_x, dt_y)
library(stringr)
dt_meta <- data.table(s=dt$s, str_split_fixed(dt$s, "_", 5))
names(dt_meta) <- c("s", "array_abbreviation", "cell", "treatment", "replicate", "obtained_by")
setkey(dt_meta, "s")
dt <- merge(dt_meta, dt)
dt <- dt %>% mutate(Replicate = as.integer(gsub("rep", "", replicate)))

# https://genome.ucsc.edu/ENCODE/cellTypes.html
dt_cell_line_info <- fread("cell_line_information.tsv")

# Rename cell lines
dt$cell[which(dt$cell == "SK-N-SH-RA")] <- "SK-N-SH_RA"
dt$cell[which(dt$cell == "BC-Jejunum-H12817N")] <- "BC_Jejunum_H12817N"
dt$cell[which(dt$cell == "BC-SmallIntestine-01-11002")] <- "BC_Small_Intestine_01-11002"

setkey(dt_cell_line_info, "cell")
setkey(dt, "cell")

dt <- merge(dt, dt_cell_line_info)

# Plot things
ggplot(dt, aes(x=impute_sex.f_stat, y=n_called, col=as.factor(Sex))) + 
	geom_point() + theme_bw() + 
	ylab("Number of calls on the Y") + xlab("F statistic") + scale_colour_discrete("Reported sex")

dt %>% filter((Sex == "M" & impute_sex.f_stat < 0.5) | (Sex == "F" & dt$impute_sex.f_stat > 0.5))
dt <- dt %>% mutate(reported.super.population = ifelse(grepl("[C,c]auc", Description) | grepl("[E,e]uropean", Description), "EUR", ifelse(grepl("Yoruba", Description), "AFR", NA)))

# Append this additional information to the cell line information used for the plotting
setkey(dt, "s")
dt_prs <- fread(cmd="gzcat PRS_1000G_cell_line.tsv.gz")
dt_prs[, s:=sample]
dt_prs[, sample:=NULL]
setkey(dt_prs, "s")

dt <- merge(dt, dt_prs, all.y=TRUE)

# Sanity check
dt <- dt %>% mutate(sample = gsub("rep[0-9]_", "", s))
dt_summary_check <- dt %>% filter(source == "Cell line") %>% group_by(sample) %>% summarise(check=length(unique((classification_strict))))

# Repeated samples
replicates <- dt %>% filter(sample %in% names(which(table(dt$sample)>1)))
ggplot(data=replicates, aes(x=PC1, y=PC2, col=factor(sample))) + geom_point()

# Remove samples where the Sex doesn't match the imputed sex
dt <- dt %>% filter(
	((Sex == "M") & (!impute_sex.is_female)) | 
	((Sex == "F") & (impute_sex.is_female)) |
	(Sex == "U") | (is.na(Sex)))
dt <- dt %>% filter(
	(`reported.super.population` == classification_strict) |
	(is.na(reported.super.population)),
	classification_strict != "unsure"
	)

dt <- dt %>% filter((Karyotype != "cancer") | (is.na(Karyotype)))
# Remove replicate from the sample ID so we can use it as a facet.
dt <- dt %>% mutate(sample = gsub("rep[0-9]_", "", s))
# Ensure that samples with the same sample ID have the same estimated population
dt %>% filter(source == "Cell line") %>% group_by(sample) %>% summarise(check=length(unique((classification_strict))))
# Use the above data.table in combination with the .json below to explore the data.
meta <- fromJSON("/well/lindgren/UKBIOBANK/dpalmer/PRS_cell_data/data/PGS_score_files/pgs-catalog-20221123-hg19/scores.meta.json")




dt_ukb <- fread(ukb_scores[1], key="sample")
for (ukb in ukb_scores[-1]) {
	dt_ukb <- merge(dt_ukb, fread(ukb, key="sample"))
}
dt_ukb <- dt_ukb %>% mutate(sample = gsub("_.*", "", sample))




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

