# Cell line PRS evaluation pipeline

The file structure of the resultant data (and input data) on the cluster is shown below.

```
PRS_cell_data
   |-data
   |---1000G
   |-----B37
   |-----B38
   |---PGS_score_files
   |-----pgs-catalog-20221123-hg19
   |-------scores
   |-----pgs-catalog-20221123-hg38
   |-------scores
   |---PGS_score_results
   |-----imputed
   |---TOPMed_variants
   |---UKB_imputed_subset
   |---celldataB36
   |---celldataB37
   |-----imputed
   |-------HRC
   |-------TOPMED
   |---celldataB38
   |-----imputed
   |-------TOPMED
   |---------B38
   |---------lifted_from_B37
   |---for_imputation
```

The scripts in this repository are located in main `PRS_cell_data` directory. Here's a brief description of what each of the scripts does, the required inputs, and where they were grabbed from.

## 01_celldata_qc.r
The first step is to clean-up the cell-line data. This script grabs the raw data that was downloaded (`ENCFF509XPS.tsv` and `ENCFF735RWY.tsv`) and cleans it up. We use the strand files for each of the genome builds, taken from https://www.well.ox.ac.uk/~wrayner/strand/ for `Human1M-Duov3_B` on the Illumina strand. A technical note about the Illumina strand and how A/B are mapped is [here](https://www.illumina.com/documents/products/technotes/technote_topbot.pdf). Note that the first table is misleading - the text is correct (re: ordering of A/B). Using the information in the strand file, we create `.ped` and `.map` files (the human readable original plink format), and then convert to binary plink files using plink. As the Illumina strand allocation of A/B is invariant to genome build, it's straightforward to generate plink files for each build. We need to then apply reverse complements to align to the forward strand under each build to match to the reference.

We then follow (loosely) the Ricopili suggestions for pre-imputation QC.
* SNP QC: call rate ≥ 0.95
* Sample QC: call rate in cases or controls ≥ 0.97
* SNP QC: call rate ≥ 0.98
* HWE p > 1e-6

The resultant plink files are then munged a little more to get them ready to pass to the Michigan imputation server.

## 02_cell_line_imputation_preparation.sh

For this portion, we make use of some great tools from [Will Rayner](https://www.well.ox.ac.uk/~wrayner/tools/). For these, we first need dbSNP submitted sites for TOPMED (e.g. https://bravo.sph.umich.edu/freeze5/hg38/). On the cluster, these have been downloaded to `data/TOPMed_variants/bravo-dbsnp-all.GRCh37.vcf.gz` and `data/TOPMed_variants/bravo-dbsnp-all.GRCh38.vcf.gz`. Note, this first step isn't required for HRC. To check allele frequency differences (and the possibility of allele flips), we evaluate the frequency of variant sites using plink. This (`.freq` file), together with the `.bim` is passed to `HRC-1000G-check-bim.pl` to perform the relevant checks. The perl script (`HRC-1000G-check-bim.pl`) provides a series of plink commands in a single shell script (`Run-plink.sh`). Finally, we block gzip (`bgzip`) the output chromosome specific `.vcf` files. The data is now ready for imputation!

I ran imputation on both the Michigan imputation server (against the HRC reference panel), and the TOPMED imputation server. In each case, I used the web-based GUI, but it's worth noting that there is software that lets you do this via the terminal which is super easy to use and very helpful if you're a power user! Everything you need to know about it is [here](https://imputationbot.readthedocs.io/en/latest/).

I downloaded the resultant imputed files to:

```
PRS_cell_data
   |-data
   |---celldataB37
   |-----imputed
   |-------HRC
   |-------TOPMED
   |---celldataB38
   |-----imputed
   |-------TOPMED
   |---------B38
   |---------lifted_from_B37
```
This was overkill - we imputed to HRC and TOPMED. For the remainder of the pipeline, I use the cell line data imputed to HRC using build 37 as this was the reference panel used for imputing UK Biobank (which we'll be contrasting to) and so we'd like similar error modes.

Importantly, we use the R2 > 0.3 option (the highest option) to ensure that the sites included in the output were imputed to a reasonably high accuracy (this will be important later on).

## 03_aligning_1000G_and_cell_line.sh

Grab the cleaned up 1000 genomes phase 3 data. See this [overview](https://alanaw1.github.io/post/2021/03/03/visualizing-1000-genomes-data/), and these instructions for download of against different [reference builds](https://www.cog-genomics.org/plink/2.0/resources).

I convert to plink (`.bed/.bim/.fam`) format, and restrict to the set of sites that were imputed with an R2 > 0.3 in the cell-line data.

I then send the resultant files to `.vcf` format, as this is the format required to use `pgs-calc` (see later).

I then also restrict the imputed data to the set of sites that the 1000 genomes data is now defined on. This ensure that the resultant set of sites for which both datasets are defined on is identical.

This last step is slow, so we created a submission script to run each chromosomes concurrently (`03_filter_cell_lines.sh`).

## 04_grab_PGS_scorings.r

Others have done tons of work to organise and curate PGS/PRS. It turns out there are two teams that have worked to harmonise PGS on the [PGS catalog](http://www.pgscatalog.org/) across genome builds. Code to download both sets is in `04_grab_PGS_scorings.r`. However, I found [`pgs_calc`](https://github.com/lukfor/pgs-calc) very straightforward to use, and their curated PGS scoring files work directly with data straight off the Michigan imputation server very easily together with `pgs-calc`, so I went with those scoring files on build 37 for the remainder of the pipeline.

The all important location of freezes of PGS-catalog are:
```
wget https://imputationserver.sph.umich.edu/resources/pgs-catalog/pgs-catalog-20221123-hg19.zip
wget https://imputationserver.sph.umich.edu/resources/pgs-catalog/pgs-catalog-20221123-hg38.zip
```
I then create subsets of scoring files to parallelise over when submitting jobs to the cluster.

Importantly, the authors also provide software [pgs-repository-builder](https://github.com/lukfor/pgs-repository-builder) to create your own harmonised scoring files if there are specific PGS scoring files that are not present in one of the freezes.

## 05_pgs_1000G_submission.sh

I then download and install `pgs-calc` and apply it to all of the scoring files in the PGS catalog for the 1000 genomes dataset. Note that this requires a dbSNP file for those scoring files that just have rsids.

The software creates a set of scores across samples, as well as an `.html` report.

## 06_pgs_cell_line_submission.sh

I do the same thing for the cell lines.

## ukb_geno_qc.r

A simple R script to read in merge and extract a set of samples in the UK Biobank estimated to be unrelated members of superpopulation ancestry labels. We extract 15,000 estimated to be in EUR, and as many as possible from AFR, AMR, EAS, SAS. Population labels were taken from `/well/lindgren/UKBIOBANK/dpalmer/superpopulation_assignments/superpopulation_labels.tsv`.

## ukb_subset_vcf_creation.sh

## plot_cell_line_prs.r
