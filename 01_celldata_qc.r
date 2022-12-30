library(data.table)
library(dplyr)

reverse_complement <- function(ACGT)
{
	if (ACGT == "A") {
		ACGT <- "T"
	} else if (ACGT == "C") {
		ACGT <- "G"
	} else if (ACGT == "G") {
		ACGT <- "C"
	} else if (ACGT == "T") {
		ACGT <- "A"
	}
	return(ACGT)
}

# export PATH="/well/lindgren/dpalmer/:$PATH"
setwd("/well/lindgren/UKBIOBANK/dpalmer/PRS_cell_data/data")

dt <- fread("/well/lindgren/ferreira/cell_lines/ENCFF735RWY.tsv")
dt <- dt %>% filter(rowSums(dt[,4:ncol(dt)] == "NC") < 80)
dt <- dt %>% rename(B36_Pos = `Position (hg18)`, B36_Chr = Chr)
setkeyv(dt, c("Name", "B36_Chr", "B36_Pos"))

# Restrict to the non-missing data first
mapping <- fread("/well/lindgren/ferreira/cell_lines/ENCFF509XPS.tsv", skip=5, key="Name")
mapping <- mapping[, -c(14)]
mapping <- mapping %>% rename(B36_Chr = Chr, B36_Pos = MapInfo)
setkeyv(mapping, c("Name", "B36_Chr", "B36_Pos"))

dt <- merge(dt, mapping)

# Now, we need to assign to the ref and alt to the correct strand
# Given that we're using the illumina strand information, the encoding should be preserved across genome builds. Let's see if that's true.
ilmn_b36 <- fread("../Human1M-Duov3_B-b36.Ilmn.strand")
names(ilmn_b36) <- c("Name", paste0("B36_", c("Chr", "Pos", "Match", "Orientation", "Strand")))
ilmn_b37 <- fread("../Human1M-Duov3_B-b37.Ilmn.strand")
names(ilmn_b37) <- c("Name", paste0("B37_", c("Chr", "Pos", "Match", "Orientation", "Strand")))
ilmn <- merge(ilmn_b36, ilmn_b37, by="Name")
ilmn_b38 <- fread("../Human1M-Duov3_B-b38.Ilmn.strand")
names(ilmn_b38) <- c("Name", paste0("B38_", c("Chr", "Pos", "Match", "Orientation", "Strand")))
ilmn <- merge(ilmn_b38, ilmn, by="Name")
setkeyv(ilmn, c("Name", "B36_Chr", "B36_Pos"))

# Now merge this information with dt, to remove inconsistensies
dt <- merge(dt, ilmn)

# Now, determine the A replacement for everything
# dt <- dt %>% mutate(
# 	SNP_1 = gsub("\\[([A-Z])/.*", "\\1", dt$SNP),
# 	SNP_2 = gsub("\\[[A-Z]/([A-Z]).*", "\\1", dt$SNP)
# 	)

# dt <- dt %>% mutate(
# 	A = ifelse(IlmnStrand == "TOP", SNP_1, SNP_2),
# 	B = ifelse(IlmnStrand == "BOT", SNP_1, SNP_2)
# 	)

dt <- dt %>% mutate(
	A = gsub("\\[([A-Z])/.*", "\\1", dt$SNP),
	B = gsub("\\[[A-Z]/([A-Z]).*", "\\1", dt$SNP)
)

# Now, replace all instances of A in the data with the column A, and all instances of B in the data with the column B
# Now, define 'chrom', 'pos', 'A', 'B'

cols <- grep("HAIB", names(dt), value=TRUE)
data <- matrix(nrow=length(cols), ncol=nrow(dt))

# Not efficient, but should do the right thing.
for (pos in seq(1,nrow(dt))) {
	if (pos %% 1000 == 0) print(pos)
	row <- gsub("A", dt$A[pos], dt[pos, ..cols])
	row <- gsub("B", dt$B[pos], row)
	# Remove indels
	row <- gsub("I", "0", row)
	row <- gsub("D", "0", row)
	# Set no calls to missing
	row <- gsub("NC", "0 0", row)
	row <- gsub("([A-Z])([A-Z])", "\\1 \\2", row)
	data[,pos] <- row
}

ped <- cbind(cols, cols, 0, 0, 0, 0, data)

fwrite(ped, sep=" ", file="celldataB36.ped", quote=FALSE, col.names=FALSE)

map <- data.table()
map$Chr <- paste0("chr", dt$B36_Chr)
map$Name <- dt$Name
map$Genomic_distance <- 0
map$Pos <- dt$B36_Pos

# Now, we have the non-binary plink files
fwrite(map, sep=" ", file="celldataB36.map", quote=FALSE, col.names=FALSE)
# Convert to plink binary
system("plink --file celldataB36 --out celldataB36")

# We immediately have a build 37 version
fwrite(ped, sep=" ", file="celldataB37.ped", quote=FALSE, col.names=FALSE)

map <- data.table()
map$Chr <- paste0("chr", dt$B37_Chr)
map$Name <- dt$Name
map$Genomic_distance <- 0
map$Pos <- dt$B37_Pos

fwrite(map, sep=" ", file="celldataB37.map", quote=FALSE, col.names=FALSE)
# Convert to plink binary
system("plink --file celldataB37 --out celldataB37_dp")

# We immediately have a build 38 version
fwrite(ped, sep=" ", file="celldataB38.ped", quote=FALSE, col.names=FALSE)

map <- data.table()
map$Chr <- paste0("chr", dt$B38_Chr)
map$Name <- dt$Name
map$Genomic_distance <- 0
map$Pos <- dt$B38_Pos

fwrite(map, sep=" ", file="celldataB38.map", quote=FALSE, col.names=FALSE)
# Convert to plink binary
system("plink --file celldataB38 --out celldataB38")

# Now, all fixed (100% AF) variants in the dataset will not have the ref or alt allele - fill that information in in the .bim file.
bim <- fread("celldataB37_dp.bim")
# Convert chromosomes from strings to ints where required
bim <- bim %>% mutate(V1 = ifelse(V1 == 23, "X", ifelse(V1 == 24, "Y", ifelse(V1 == 26, "MT", V1))))
# Merge the bim with the dt

dt_bim <- dt %>% select(B37_Chr, B37_Pos, Name, A, B, B37_Orientation)
bim <- bim %>% mutate(order = seq(1, nrow(bim)), Name=V2, B37_Chr=V1, B37_Pos=V4)

setkeyv(dt_bim, c("Name", "B37_Chr", "B37_Pos"))
setkeyv(bim, c("Name", "B37_Chr", "B37_Pos"))

dt_bim <- merge(dt_bim, bim) %>% select(V1, V2, V3, V4, V5, V6, A, B, order, B37_Orientation)
dt_bim <- dt_bim %>% mutate(V5 = ifelse(V5 == "0", ifelse(V6 == A, B, A), V5)) %>% select(-c("A", "B"))

# Finally, use the ilmn strand information to perform reverse complements where required
dt_bim <- dt_bim %>% mutate(
	V5 = ifelse(B37_Orientation == "-", sapply(V5, reverse_complement), V5),
	V6 = ifelse(B37_Orientation == "-", sapply(V6, reverse_complement), V6),
)
dt_bim <- dt_bim[order(order),]

fwrite(dt_bim %>% select(V1, V2, V3, V4, V5, V6), file="celldataB37_dp.bim", sep=" ", col.names=FALSE, quote=FALSE)

# Use the reference information to ensure that the ref/alt is correct
refalt <- fread("../Human1M-Duov3_B-b37.strand.RefAlt", header=FALSE)
fwrite(refalt %>% filter(nchar(V2) == 1), file="../Human1M-Duov3_B-b37.strand.RefAlt.snv", sep=" ", col.names=FALSE, quote=FALSE)

# Ensure that all indels are removed
fwrite(dt_bim %>% filter(!(dt_bim$V5 %in% c("A", "C", "G", "T"))) %>% select(V2), file="variants_to_exclude.txt", sep=" ", col.names=FALSE)
system("plink --bfile celldataB37_dp --exclude variants_to_exclude.txt --out celldataB37_dp --make-bed")

# Following (loosely) the Ricopili suggestions for pre-imputation QC
# SNP QC: call rate ≥ 0.95
system("plink --bfile celldataB37_dp --out celldataB37_dp_test --geno 0.05 --make-bed")
# Sample QC: call rate in cases or controls ≥ 0.97
system("plink --bfile celldataB37_dp_test --out celldataB37_dp_test --mind 0.03 --make-bed")
# SNP QC: call rate ≥ 0.98
system("plink --bfile celldataB37_dp_test --out celldataB37_dp_test --geno 0.02 --make-bed")
system("plink --bfile celldataB37_dp_test --out celldataB37_dp_test --hwe 1e-6 --make-bed")
# Set the correct reference allele
# Don't do this - the perl script will break things if you do!
# system("plink --bfile celldataB37_dp_test --a1-allele ../Human1M-Duov3_B-b37.strand.RefAlt.snv --out celldataB37_dp_test --make-bed")

# DO THE SAME THING FOR BUILD 38
# Now, all fixed (100% AF) variants in the dataset will not have the ref or alt allele - fill that information in in the .bim file.
bim <- fread("celldataB38.bim")
# Convert chromosomes from strings to ints where required
bim <- bim %>% mutate(V1 = ifelse(V1 == 23, "X", ifelse(V1 == 24, "Y", ifelse(V1 == 26, "MT", V1))))
# Merge the bim with the dt

dt_bim <- dt %>% select(B38_Chr, B38_Pos, Name, A, B, B38_Orientation)
bim <- bim %>% mutate(order = seq(1, nrow(bim)), Name=V2, B38_Chr=V1, B38_Pos=V4)

setkeyv(dt_bim, c("Name", "B38_Chr", "B38_Pos"))
setkeyv(bim, c("Name", "B38_Chr", "B38_Pos"))

dt_bim <- merge(dt_bim, bim) %>% select(V1, V2, V3, V4, V5, V6, A, B, order, B38_Orientation)
dt_bim <- dt_bim %>% mutate(V5 = ifelse(V5 == "0", ifelse(V6 == A, B, A), V5)) %>% select(-c("A", "B"))

# Finally, use the ilmn strand information to perform reverse complements where required
dt_bim <- dt_bim %>% mutate(
	V5 = ifelse(B38_Orientation == "-", sapply(V5, reverse_complement), V5),
	V6 = ifelse(B38_Orientation == "-", sapply(V6, reverse_complement), V6),
)
dt_bim <- dt_bim[order(order),]

fwrite(dt_bim %>% select(V1, V2, V3, V4, V5, V6), file="celldataB38.bim", sep=" ", col.names=FALSE, quote=FALSE)

# DO THE SAME THING FOR BUILD 38
# Use the reference information to ensure that the ref/alt is correct
refalt <- fread("../Human1M-Duov3_B-b38.strand.RefAlt", header=FALSE)
fwrite(refalt %>% filter(nchar(V2) == 1), file="../Human1M-Duov3_B-b38.strand.RefAlt.snv", sep=" ", col.names=FALSE, quote=FALSE)

# Ensure that all indels are removed
fwrite(dt_bim %>% filter(!(dt_bim$V5 %in% c("A", "C", "G", "T"))) %>% select(V2), file="variants_to_exclude.txt", sep=" ", col.names=FALSE)
system("plink --bfile celldataB38 --exclude variants_to_exclude.txt --out celldataB38 --make-bed")

# Following (loosely) the Ricopili suggestions for pre-imputation QC
# SNP QC: call rate ≥ 0.95
system("plink --bfile celldataB38 --out celldataB38_test --geno 0.05 --make-bed")
# Sample QC: call rate in cases or controls ≥ 0.97
system("plink --bfile celldataB38_test --out celldataB38_test --mind 0.03 --make-bed")
# SNP QC: call rate ≥ 0.98
system("plink --bfile celldataB38_test --out celldataB38_test --geno 0.02 --make-bed")
system("plink --bfile celldataB38_test --out celldataB38_test --hwe 1e-6 --make-bed")
# Set the correct reference allele
# Don't do this - the perl script will break things if you do!
# system("plink --bfile celldataB38_test --a1-allele ../Human1M-Duov3_B-b38.strand.RefAlt.snv --out celldataB38_test --make-bed")
