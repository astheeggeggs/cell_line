#!/bin/bash 
#SBATCH -J extract-vcf-ukb-subset
#SBATCH -A lindgren.prj 
#SBATCH -o output-extract-vcf.out 
#SBATCH -e error-extract-vcf.err 
#SBATCH -c 4 
#SBATCH -p short
#SBATCH --array 1-22:1 
#SBATCH --requeue

export PATH="/well/lindgren/dpalmer/:$PATH"

chr=${SLURM_ARRAY_TASK_ID}

awk '{print $1":"$4"_"$6"_"$5}' /well/lindgren/UKBIOBANK/dpalmer/PRS_cell_data/data/1000G/B37/all_phase3_r2_0.3_chr${chr}.bim > /well/lindgren/UKBIOBANK/dpalmer/PRS_cell_data/data/1000G/B37/all_phase3_r2_0.3_chr${chr}_variants.txt

/well/lindgren/dpalmer/qctool2/qctool/build/release/apps/qctool_v2.2.0 -g /well/lindgren-ukbb/projects/ukbb-11867/DATA/IMPUTATION/ukb_imp_chr${chr}_v3.bgen -s /well/lindgren-ukbb/projects/ukbb-11867/DATA/SAMPLE_FAM/ukb11867_imp_chr1_v3_s487395.sample -incl-snpids /well/lindgren/UKBIOBANK/dpalmer/PRS_cell_data/data/1000G/B37/all_phase3_r2_0.3_chr${chr}_variants.txt -incl-samples /well/lindgren/UKBIOBANK/dpalmer/PRS_cell_data/data/UKB_subset_one_col.txt -og /well/lindgren/UKBIOBANK/dpalmer/PRS_cell_data/data/UKB_imputed_subset/UKB_imputed_subset_chr${chr}.vcf.gz
/well/lindgren/dpalmer/qctool2/qctool/build/release/apps/qctool_v2.2.0 -g /well/lindgren/UKBIOBANK/dpalmer/PRS_cell_data/data/UKB_imputed_subset/UKB_imputed_subset_chr${chr}.vcf.gz -threshold 0.5001 -og /well/lindgren/UKBIOBANK/dpalmer/PRS_cell_data/data/UKB_imputed_subset/UKB_imputed_subset_chr${chr}_thresholded_05.vcf.gz

# Then, convert them all to plink to perform hard calls, then convert back from plink to vcf, then bgzip, then run then through the pgs-calculator.

plink --vcf /well/lindgren/UKBIOBANK/dpalmer/PRS_cell_data/data/UKB_imputed_subset/UKB_imputed_subset_chr${chr}_thresholded_05.vcf.gz --out /well/lindgren/UKBIOBANK/dpalmer/PRS_cell_data/data/UKB_imputed_subset/UKB_imputed_subset_chr${chr}_thresholded_recode --real-ref-alleles --make-bed

# Recode the bim files to be in the right format
awk -v OFS='\t' '{print $1,$1":"$4":"$6":"$5,$3,$4,$5,$6}' /well/lindgren/UKBIOBANK/dpalmer/PRS_cell_data/data/UKB_imputed_subset/UKB_imputed_subset_chr${chr}_thresholded_recode.bim > /well/lindgren/UKBIOBANK/dpalmer/PRS_cell_data/data/UKB_imputed_subset/UKB_imputed_subset_chr${chr}_thresholded_recode.bim.tmp

mv /well/lindgren/UKBIOBANK/dpalmer/PRS_cell_data/data/UKB_imputed_subset/UKB_imputed_subset_chr${chr}_thresholded_recode.bim.tmp /well/lindgren/UKBIOBANK/dpalmer/PRS_cell_data/data/UKB_imputed_subset/UKB_imputed_subset_chr${chr}_thresholded_recode.bim

# Then recode to vcf
plink --double-id --bfile /well/lindgren/UKBIOBANK/dpalmer/PRS_cell_data/data/UKB_imputed_subset/UKB_imputed_subset_chr${chr}_thresholded_recode --real-ref-alleles --recode vcf --out /well/lindgren/UKBIOBANK/dpalmer/PRS_cell_data/data/UKB_imputed_subset/UKB_imputed_subset_chr${chr}_thresholded_recode_final

module purge
module load BCFtools
bgzip /well/lindgren/UKBIOBANK/dpalmer/PRS_cell_data/data/UKB_imputed_subset/UKB_imputed_subset_chr${chr}_thresholded_recode_final.vcf
