library(ggplot2)
library(data.table)
library(dplyr)
library(jsonlite)

Sys.setlocale(category = "LC_ALL", locale = "en_US.UTF-8")

# Grab all of the score files and merge them
reference_scores <- grep("scores.txt", dir("/well/lindgren/UKBIOBANK/dpalmer/PRS_cell_data/data/PGS_score_results", recursive=TRUE, full.names=TRUE), value=TRUE)
cell_line_scores <- grep("cell_line", reference_scores, value=TRUE)
reference_scores <- setdiff(reference_scores, cell_line_scores)

dt_ref <- fread(reference_scores[1], key="sample")
for (ref_score in reference_scores[-1]) {
	dt_ref <- merge(dt_ref, fread(ref_score, key="sample"))
}

dt_cell_line <- fread(cell_line_scores[1], key="sample")
for (cell_score in cell_line_scores[-1]) {
	dt_cell_line <- merge(dt_cell_line, fread(cell_score, key="sample"))
}

pops <- c("EUR", "EAS", "SAS", "AMR", "AFR")
dt_loc <- paste0("/well/lindgren/UKBIOBANK/dpalmer/PRS_cell_data/data/", pops, "/imputed/scores.txt")
dt_UKB <- rbind(
	fread(dt_loc[1]) %>% mutate(super_population=pops[1]),
	fread(dt_loc[2]) %>% mutate(super_population=pops[2]),
	fread(dt_loc[3]) %>% mutate(super_population=pops[3]),
	fread(dt_loc[4]) %>% mutate(super_population=pops[4]),
	fread(dt_loc[5]) %>% mutate(super_population=pops[5])
	)
setkey(dt_UKB, "sample")

# Merge this information with the 1000G population information
system("wget -P /well/lindgren/UKBIOBANK/dpalmer/PRS_cell_data/data/1000G/ ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/working/20130606_sample_info/20130606_g1k.ped")
system("wget -P /well/lindgren/UKBIOBANK/dpalmer/PRS_cell_data/data/1000G/ ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/20131219.populations.tsv")

dt_sample_pop <- fread("/well/lindgren/UKBIOBANK/dpalmer/PRS_cell_data/data/1000G/20130606_g1k.ped")
dt_pop <- fread("/well/lindgren/UKBIOBANK/dpalmer/PRS_cell_data/data/1000G/20131219.populations.tsv")

dt_pop <- fread("/well/lindgren/UKBIOBANK/dpalmer/PRS_cell_data/data/1000G/20131219.populations.tsv")
dt_sample_pop <- fread("/well/lindgren/UKBIOBANK/dpalmer/PRS_cell_data/data/1000G/20130606_g1k.ped")

dt_sample_pop <- left_join(dt_sample_pop, dt_pop[, 1:3], by = c("Population" = "Population Code"))
dt_ref <- dt_ref %>% mutate(sample=gsub("0_", "", sample))
dt_ref <- merge(dt_sample_pop %>% rename(sample = `Individual ID`, super_population=`Super Population`), dt_ref, by="sample")

meta <- fromJSON("/well/lindgren/UKBIOBANK/dpalmer/PRS_cell_data/data/PGS_score_files/scores.meta.json")
pdf("PRS_plots.pdf", width=3, height=3)
i <- 1
for (PRS in grep("PGS", intersect(names(dt_ref), names(dt_cell_line)), value=TRUE)) {
	if ((i %% 100) == 0) print(i)
	p <- ggplot() + 
	geom_density(data = dt_ref, aes_string(x = PRS, col="super_population")) +
	geom_rug(data = dt_cell_line, aes_string(x = PRS)) + theme_classic() + theme(legend.position="none")
	p <- p + ggtitle(paste(meta[[PRS]]$trait)) + #theme(plot.title = element_textbox_simple()) +
	labs(caption = paste(
		meta[[PRS]]$link,
		paste0(meta[[PRS]]$publication$doi),
		paste(meta[[PRS]]$variants, "variants"), sep="\n"))
	print(p)
	i <- i+1
}
dev.off()

dt <- rbind(
	dt_UKB %>% mutate(source="UK Biobank"),
	dt_cell_line %>% mutate(source="Cell line"),
	fill=TRUE)
pdf("PRS_plots_violin.pdf", width=6, height=3)
i <- 1
for (PRS in grep("PGS", intersect(names(dt_UKB), names(dt_cell_line)), value=TRUE)[1:10]) {
	if ((i %% 100) == 0) print(i)
	p <- ggplot() + 
	geom_violin(data=dt, aes_string(x="super_population", y=PRS, fill="source"), trim=FALSE) + 
	theme_classic() + theme(legend.position="none")
	p <- p + ggtitle(meta[[PRS]]$trait) +
	labs(caption = paste(
		meta[[PRS]]$link,
		paste0(meta[[PRS]]$publication$doi),
		paste(meta[[PRS]]$variants, "variants"), sep="\n"))
	print(p)
	i <- i+1
}
dev.off()



