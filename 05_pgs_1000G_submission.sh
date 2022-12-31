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

mkdir /well/lindgren/UKBIOBANK/dpalmer/PRS_cell_data/data/PGS_score_results/job_${SLURM_ARRAY_TASK_ID}

pgs-calc apply --ref /well/lindgren/UKBIOBANK/dpalmer/PRS_cell_data/data/PGS_score_files/pgs-catalog-20221123-hg19/scores/job_${SLURM_ARRAY_TASK_ID}.txt /well/lindgren/UKBIOBANK/dpalmer/PRS_cell_data/data/1000G/B37/all_phase3_r2_0.3_chr*.vcf.gz --out /well/lindgren/UKBIOBANK/dpalmer/PRS_cell_data/data/PGS_score_results/job_${SLURM_ARRAY_TASK_ID}/job_${SLURM_ARRAY_TASK_ID}_scores.txt --dbsnp /well/lindgren/UKBIOBANK/dpalmer/PRS_cell_data/data/dbsnp154_hg19.txt.gz --report-html /well/lindgren/UKBIOBANK/dpalmer/PRS_cell_data/data/PGS_score_results/job_${SLURM_ARRAY_TASK_ID}/job_${SLURM_ARRAY_TASK_ID}_scores_report.html --genotypes GT --meta /well/lindgren/UKBIOBANK/dpalmer/PRS_cell_data/data/PGS_score_files/pgs-catalog-20221123-hg19/scores.meta.json --report-csv /well/lindgren/UKBIOBANK/dpalmer/PRS_cell_data/data/PGS_score_results/job_${SLURM_ARRAY_TASK_ID}/job_${SLURM_ARRAY_TASK_ID}_scores_report.csv
