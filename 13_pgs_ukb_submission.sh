#!/bin/bash 
#SBATCH -J ukb-prs
#SBATCH -A lindgren.prj 
#SBATCH -o ukb.full.int.prs.output.%A_%a.out
#SBATCH -e ukb.full.int.prs.error.%A_%a.err 
#SBATCH -c 4 
#SBATCH -p long
#SBATCH --array 1-33:1 
#SBATCH --requeue

export PATH="/well/lindgren/dpalmer/:$PATH"
pgs_results_dir="/well/lindgren/UKBIOBANK/dpalmer/PRS_cell_data/data/PGS_score_results"
pgs_score_dir="/well/lindgren/UKBIOBANK/dpalmer/PRS_cell_data/data/PGS_score_files"
ukb_imputed_intersected_1000G_subset_dir="/well/lindgren/UKBIOBANK/dpalmer/PRS_cell_data/data/UKB_imputed_subset/UKB_cell_line_1000G_intersection"
mkdir -p ${pgs_results_dir}/ukb_job_${SLURM_ARRAY_TASK_ID}/fully_intersected

pgs-calc apply \
	--ref ${pgs_score_dir}/pgs-catalog-20221123-hg19/scores/job_${SLURM_ARRAY_TASK_ID}.txt \
	${ukb_imputed_intersected_1000G_subset_dir}/UKB_imputed_subset_chr*_thresholded_recode_merge_typed_final_UKB_cell_line_1000G_intersected.vcf.gz \
	--out ${pgs_results_dir}/ukb_job_${SLURM_ARRAY_TASK_ID}/fully_intersected/ukb_intersect_1000G_cell_line_job_${SLURM_ARRAY_TASK_ID}_scores.txt \
	--dbsnp /well/lindgren/UKBIOBANK/dpalmer/PRS_cell_data/data/dbsnp154_hg19.txt.gz \
	--report-html ${pgs_results_dir}/ukb_job_${SLURM_ARRAY_TASK_ID}/fully_intersected/ukb_intersect_1000G_cell_line_job_${SLURM_ARRAY_TASK_ID}_scores_report.html \
	--genotypes GT \
	--meta ${pgs_score_dir}/pgs-catalog-20221123-hg19/scores.meta.json \
	--report-csv ${pgs_results_dir}/ukb_job_${SLURM_ARRAY_TASK_ID}/fully_intersected/ukb_intersect_1000G_cell_line_job_${SLURM_ARRAY_TASK_ID}_scores_report.csv \
	--write-variants ${pgs_results_dir}/ukb_job_${SLURM_ARRAY_TASK_ID}/fully_intersected/ukb_intersect_1000G_cell_line_job_${SLURM_ARRAY_TASK_ID}_scores_variants_used.csv

gzip -f ${pgs_results_dir}/ukb_job_${SLURM_ARRAY_TASK_ID}/fully_intersected/ukb_intersect_1000G_cell_line_job_${SLURM_ARRAY_TASK_ID}_scores_variants_used.csv
