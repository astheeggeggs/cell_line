library(bigsnpr)
library(data.table)
library(dplyr)
library(ggplot2)
library(randomForest)
library(ggrastr)

# Plotting functions:
T_RF <- 0.95

# Function to sanity check the random matrix theory projections
project_onto_ref_PCs <- function(bed.ref, obj.bed, n_PCs = 20, n_PC_plot = 8,
  outdir = "/well/lindgren/UKBIOBANK/dpalmer/PRS_cell_data/superpopulation_assignments", filename = "output", strand_flip=FALSE) {

  test <- bed_projectPCA(bed.ref, obj.bed, k = n_PCs, strand_flip = FALSE)#, ncores = nb_cores())

  plot(test$obj.svd.ref)
  ggsave(paste0(outdir, filename, "_scree.pdf"), width = 13, height = 7)
  plot(test$obj.svd, type = "loadings", loadings = 1:n_PCs, coeff = 0.4)
  ggsave(paste0(outdir, filename, "_loadings.pdf"), width = 13, height = 7)

  PC.ref <- predict(test$obj.svd.ref)
  proj1 <- test$simple_proj
  proj2 <- test$OADP_proj
  save("test", "PC.ref", "proj1", "proj2", file=paste0(outdir, filename, "_data.Rdata"))

  # shrinkage coefficients
  shrinkage <- unname(sapply(1:n_PCs, function(k) {
    MASS::rlm(proj2[, k] ~ proj1[, k] + 0, maxit=100)$coef
  }))
  print(round(shrinkage, 2))

  ind <- seq(1, min(20e3, nrow(proj1)))

  plot_grid(plotlist = lapply(1:floor(n_PC_plot/2), function(k) {
    k1 <- 2 * k - 1
    k2 <- 2 * k
    plot_grid(
      qplot(PC.ref[, k1], PC.ref[, k2], size = I(2)) +
        geom_point(aes(proj1[ind, k1], proj1[ind, k2]), color = "red") +
        theme_bigstatsr(0.5) +
        labs(x = paste0("PC", k1), y = paste0("PC", k2)),
      qplot(PC.ref[, k1], PC.ref[, k2], size = I(2)) +
        geom_point(aes(proj2[ind, k1], proj2[ind, k2]), color = "blue") +
        theme_bigstatsr(0.5) +
        labs(x = paste0("PC", k1), y = paste0("PC", k2)),
      scale = 0.95
    )
  }), nrow = floor(sqrt(n_PCs)))
  ggsave(file=paste0(outdir, filename, "_shrinkage.pdf"), width=15, height=15)
  return(round(shrinkage, 2))

}

# First, generate cleaned up genotype data
write(sapply(1:22, function(chr) {
  c(paste0("/well/lindgren/UKBIOBANK/dpalmer/PRS_cell_data/data/celldataB37/celldataB37_dp_test-updated-chr", chr, ".bed"),
    paste0("/well/lindgren/UKBIOBANK/dpalmer/PRS_cell_data/data/celldataB37/celldataB37_dp_test-updated-chr", chr, ".bim"),
    paste0("/well/lindgren/UKBIOBANK/dpalmer/PRS_cell_data/data/celldataB37/celldataB37_dp_test-updated-chr", chr, ".fam"))
}), '/well/lindgren/UKBIOBANK/dpalmer/PRS_cell_data/data/celldataB37/plink_get_pop_filepaths.txt', ncolumns = 3)

# If this dies, it's likely due to memory. Requires at least 10 slots on qe.
snp_plinkQC(
  plink.path = '/well/lindgren/dpalmer/plink',
  prefix.in = '/well/lindgren/UKBIOBANK/dpalmer/PRS_cell_data/data/celldataB37/plink_get_pop_filepaths.txt',
  file.type = "--merge-list",
  prefix.out = "/well/lindgren/UKBIOBANK/dpalmer/PRS_cell_data/data/celldataB37/celldataB37_dp_test-updated-autosomes-combined",
  geno = 0.01,
  autosome.only = TRUE
)

obj.bed <- bed('/well/lindgren/UKBIOBANK/dpalmer/PRS_cell_data/data/celldataB37/celldataB37_dp_test-updated-autosomes-combined.bed')
bed.ref <- bed(download_1000G("/well/lindgren/dpalmer/ukb_get_EUR/data"))

system("mkdir /well/lindgren/UKBIOBANK/dpalmer/PRS_cell_data/superpopulation_assignments")
outdir <- '/well/lindgren/UKBIOBANK/dpalmer/PRS_cell_data/superpopulation_assignments/'
filename <- 'cell_line_projected_to_1kg_PCs'
n_PCs <- 10
save_figures <- TRUE
perform_plotting <- TRUE
creating_new_EUR_def <- TRUE

project_onto_ref_PCs(bed.ref, obj.bed, filename=filename, outdir=outdir, n_PCs=n_PCs)

# Grab 1000G information containing populations
# download.file("ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/working/20130606_sample_info/20130606_g1k.ped",
#               destfile = "/well/lindgren/UKBIOBANK/dpalmer/1kg_for_EUR_assign/20130606_g1k.ped")
dt_pop <- fread("/well/lindgren/UKBIOBANK/dpalmer/1kg_for_EUR_assign/20130606_g1k.ped") %>% 
  mutate(sample.ID = `Individual ID`, population=Population) %>% 
  select(sample.ID, population) %>%
  # Determine the superpopulations
  mutate(super.population = case_when(
    population %in% c("CHB", "JPT", "CHS", "CDX", "KHV") ~ "EAS",
    population %in% c("CEU", "TSI", "FIN", "GBR", "IBS") ~ "EUR",
    population %in% c("YRI", "LWK", "MAG", "MSL", "ESN", "ASW", "ACB", "GWD") ~ "AFR",
    population %in% c("AMR", "PUR", "CLM", "PEL", "MXL") ~ "AMR",
    population %in% c("GIH", "PJL","BEB", "STU", "ITU") ~ "SAS",
    TRUE ~ "other"
  )
)

load(paste0(outdir, filename, "_data.Rdata"))
PC.ref <- data.table(PC.ref)
names(PC.ref) <- paste0("PC", seq(1,n_PCs))

dt <- cbind(data.table(bed.ref$fam), PC.ref)
setkey(dt, "sample.ID")
setkey(dt_pop, "sample.ID")
dt_1kg <- merge(dt, dt_pop)

proj2 <- data.table(proj2)
names(proj2) <- paste0("PC", seq(1,n_PCs))
dt_ukb <- cbind(data.table(obj.bed$fam), proj2)

# Now, need to grab the population labels to train the random forest.
dt_train = dt_1kg %>%
  select(c(super.population, population, PC1, PC2, PC3, PC4, PC5, PC6, PC7, PC8, PC9, PC10))

dt_predict = dt_ukb %>%
  select(c(sample.ID, PC1, PC2, PC3, PC4, PC5, PC6, PC7, PC8, PC9, PC10))

PCs_to_use <- paste0('PC', seq(1,4))

# Determine a classifier.
set.seed(160487)
rf <- randomForest(x=dt_train %>% select(PCs_to_use), y=as.factor(as.character(dt_train$super.population)), ntree=10000)
rf_probs <- predict(rf, dt_predict %>% select(PCs_to_use), type='prob')

check_thres <- function(row, threshold) {
  return(!any(row > threshold))
}

unsure <- apply(rf_probs, 1, check_thres, T_RF)
classification <- as.character(predict(rf, dt_predict %>% select(PCs_to_use)))
dt_predict$classification_loose <- as.factor(classification)
classification[unsure] <- 'unsure'
dt_predict$classification_strict <- as.factor(classification)

dt_predict <- dt_predict %>% mutate(sample.ID = as.character(sample.ID)) %>% select(sample.ID, classification_strict, classification_loose, starts_with("PC"))
dt_predict <- data.table(dt_predict)
setkeyv(dt_predict, c("sample.ID", paste0("PC", seq(1,10))))
setkeyv(dt_1kg, c("sample.ID", paste0("PC", seq(1,10))))
dt_classify <- merge(dt_predict, dt_1kg, all=TRUE) %>% select(sample.ID, classification_strict, classification_loose, starts_with("PC"), super.population, population)

# Write the result to file.
fwrite(dt_classify, file = "/well/lindgren/UKBIOBANK/dpalmer/PRS_cell_data/superpopulation_assignments/superpopulation_labels.tsv", sep='\t')
