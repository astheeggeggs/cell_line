#!/bin/bash 
#SBATCH -J cell-line-prs
#SBATCH -A lindgren.prj 
#SBATCH -o output.out 
#SBATCH -e error.err 
#SBATCH -c 4 
#SBATCH -p short
#SBATCH --array 1-33:1 
#SBATCH --requeue

export PATH="/well/lindgren/dpalmer/:$PATH"

mkdir /well/lindgren/UKBIOBANK/dpalmer/PRS_cell_data/data/PGS_score_results/cell_line_job_${SLURM_ARRAY_TASK_ID}

pgs-calc apply --ref /well/lindgren/UKBIOBANK/dpalmer/PRS_cell_data/data/PGS_score_files/scores/job_${SLURM_ARRAY_TASK_ID}.txt /well/lindgren/UKBIOBANK/dpalmer/PRS_cell_data/data/celldataB37/imputed/HRC/chr*.dose.subset.vcf.gz --out /well/lindgren/UKBIOBANK/dpalmer/PRS_cell_data/data/PGS_score_results/cell_line_job_${SLURM_ARRAY_TASK_ID}/cell_line_job_${SLURM_ARRAY_TASK_ID}_scores.txt --dbsnp /well/lindgren/UKBIOBANK/dpalmer/PRS_cell_data/data/dbsnp154_hg19.txt.gz --report-html /well/lindgren/UKBIOBANK/dpalmer/PRS_cell_data/data/PGS_score_results/cell_line_job_${SLURM_ARRAY_TASK_ID}/cell_line_job_${SLURM_ARRAY_TASK_ID}_scores_report.html --genotypes GT --meta /well/lindgren/UKBIOBANK/dpalmer/PRS_cell_data/data/PGS_score_files/scores.meta.json --report-csv /well/lindgren/UKBIOBANK/dpalmer/PRS_cell_data/data/PGS_score_results/cell_line_job_${SLURM_ARRAY_TASK_ID}/cell_line_job_${SLURM_ARRAY_TASK_ID}_scores_report.csv
