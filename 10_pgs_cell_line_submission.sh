#!/bin/bash 
#SBATCH -J cell-line-prs
#SBATCH -A lindgren.prj 
#SBATCH -o cell_line.ukb.int.prs.output.%A_%a.out 
#SBATCH -e cell_line.ukb.int.prs.error.%A_%a.err 
#SBATCH -c 4 
#SBATCH -p short
#SBATCH --array 1-33:1 
#SBATCH --requeue

export PATH="/well/lindgren/dpalmer/:$PATH"
pgs_results_dir="/well/lindgren/UKBIOBANK/dpalmer/PRS_cell_data/data/PGS_score_results"
pgs_score_dir="/well/lindgren/UKBIOBANK/dpalmer/PRS_cell_data/data/PGS_score_files"
imputed_path="/well/lindgren/UKBIOBANK/dpalmer/PRS_cell_data/data/celldataB37/HRC/imputed"
mkdir -p ${pgs_results_dir}/cell_line_job_${SLURM_ARRAY_TASK_ID}

pgs-calc apply \
	--ref ${pgs_score_dir}/pgs-catalog-20221123-hg19/scores/job_${SLURM_ARRAY_TASK_ID}.txt \
	${imputed_path}/chr*_UKB_intersected.dose.vcf.gz \
	--out ${pgs_results_dir}/cell_line_job_${SLURM_ARRAY_TASK_ID}/cell_line_intersect_UKB_job_${SLURM_ARRAY_TASK_ID}_scores.txt \
	--dbsnp /well/lindgren/UKBIOBANK/dpalmer/PRS_cell_data/data/dbsnp154_hg19.txt.gz \
	--report-html ${pgs_results_dir}/cell_line_job_${SLURM_ARRAY_TASK_ID}/cell_line_intersect_UKB_job_${SLURM_ARRAY_TASK_ID}_scores_report.html \
	--genotypes GT \
	--meta ${pgs_score_dir}/pgs-catalog-20221123-hg19/scores.meta.json \
	--report-csv ${pgs_results_dir}/cell_line_job_${SLURM_ARRAY_TASK_ID}/cell_line_intersect_UKB_job_${SLURM_ARRAY_TASK_ID}_scores_report.csv \
	--write-variants ${pgs_results_dir}/cell_line_job_${SLURM_ARRAY_TASK_ID}/cell_line_intersect_UKB_job_${SLURM_ARRAY_TASK_ID}_scores_variants_used.csv
