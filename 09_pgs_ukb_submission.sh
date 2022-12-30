#!/bin/bash 
#SBATCH -J 1000G-prs
#SBATCH -A lindgren.prj 
#SBATCH -o output.out 
#SBATCH -e error.err 
#SBATCH -c 4 
#SBATCH -p short
#SBATCH --array 1-33:1 
#SBATCH --requeue

export PATH="/well/lindgren/dpalmer/:$PATH"

mkdir /well/lindgren/UKBIOBANK/dpalmer/PRS_cell_data/data/PGS_score_results/job_ukb_${SLURM_ARRAY_TASK_ID}

pgs-calc apply --ref /well/lindgren/UKBIOBANK/dpalmer/PRS_cell_data/data/PGS_score_files/scores/job_${SLURM_ARRAY_TASK_ID}.txt /well/lindgren/UKBIOBANK/dpalmer/PRS_cell_data/data/UKB_imputed_subset/UKB_imputed_subset_chr*_thresholded_recode_final.vcf.gz --out /well/lindgren/UKBIOBANK/dpalmer/PRS_cell_data/data/PGS_score_results/job_ukb_${SLURM_ARRAY_TASK_ID}/job_ukb_${SLURM_ARRAY_TASK_ID}_scores.txt --dbsnp /well/lindgren/UKBIOBANK/dpalmer/PRS_cell_data/data/dbsnp154_hg19.txt.gz --report-html /well/lindgren/UKBIOBANK/dpalmer/PRS_cell_data/data/PGS_score_results/job_ukb_${SLURM_ARRAY_TASK_ID}/job_ukb_${SLURM_ARRAY_TASK_ID}_scores_report.html --genotypes GT
