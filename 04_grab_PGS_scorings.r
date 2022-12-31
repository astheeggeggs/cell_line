library(data.table)
library(dplyr)

# Install the pgs-calc tool
# https://github.com/lukfor/pgs-calc

# Download .json containing the locations of the harmonised scoring files
# wget https://raw.githubusercontent.com/genepi/imputationserver/master/files/pgs-catalog.json

# Load the package required to read JSON files.
library("rjson")

# Give the input file name to the function.
PGS_score_file_location <- "/well/lindgren/UKBIOBANK/dpalmer/PRS_cell_data/data/PGS_score_files/"
PGS <- fromJSON(file = paste0(PGS_score_file_location, "pgs-catalog.json"))$results
# These are the harmonized scoring files provided by Cambridge
for (i in 1:length(PGS)) {
	system(paste("wget -P", PGS_score_file_location, PGS[[i]]$ftp_harmonized_scoring_files$GRCh37$positions))
	system(paste("wget -P", PGS_score_file_location, PGS[[i]]$ftp_harmonized_scoring_files$GRCh38$positions))
}

# These are the harmonized scoring files provided by the Michigan Imputation server people
# https://imputationserver.sph.umich.edu/resources/pgs-catalog/
# See https://github.com/lukfor/pgs-repository-builder - this is what we should use, as these 
# files can be directly input into pgs-calc.
system("wget https://imputationserver.sph.umich.edu/resources/pgs-catalog/pgs-catalog-20221123-hg19.zip")
system("wget https://imputationserver.sph.umich.edu/resources/pgs-catalog/pgs-catalog-20221123-hg38.zip")

# Create a set of files to run the PGS on, split into 20 separate jobs
# List the set of scoring files on B37
score_files <- dir("/well/lindgren/UKBIOBANK/dpalmer/PRS_cell_data/data/PGS_score_files/pgs-catalog-20221123-hg19/scores", full.names=FALSE)
# List the set of scoring files of B37 and B38
dt_B37 <- data.table(scores=score_files) %>% filter(grepl("txt.gz", scores))
dt_B37 <- dt_B37 %>% mutate(job = ceiling(seq(1,nrow(dt_B37))/100))

for (j in unique(dt_B37$job)) {
	fwrite(dt_B37 %>% filter(job == j) %>% select(scores),
		file=paste0("/well/lindgren/UKBIOBANK/dpalmer/PRS_cell_data/data/PGS_score_files/pgs-catalog-20221123-hg19/scores/job_", j, ".txt"),
		quote=FALSE)
}
