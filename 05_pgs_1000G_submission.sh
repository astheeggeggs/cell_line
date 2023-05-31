#!/bin/bash 
#SBATCH -J 1000G-prs
#SBATCH -A lindgren.prj 
#SBATCH -o 1000G.prs.output.%A_%a.out 
#SBATCH -e 1000G.prs.error.%A_%a.err 
#SBATCH -c 4 
#SBATCH -p short
#SBATCH --array 1-33:1 
#SBATCH --requeue

export PATH="/well/lindgren/dpalmer/:$PATH"
pgs_results_dir="/well/lindgren/UKBIOBANK/dpalmer/PRS_cell_data/data/PGS_score_results"
pgs_score_dir="/well/lindgren/UKBIOBANK/dpalmer/PRS_cell_data/data/PGS_score_files"
phase3_1kg_path="/well/lindgren/UKBIOBANK/dpalmer/PRS_cell_data/data/1000G/B37"
mkdir ${pgs_results_dir}/job_${SLURM_ARRAY_TASK_ID}

pgs-calc apply \
	--ref ${pgs_score_dir}/pgs-catalog-20221123-hg19/scores/job_${SLURM_ARRAY_TASK_ID}.txt \
	${phase3_1kg_path}/all_phase3_r2_0.3_chr*.vcf.gz \
	--out ${pgs_results_dir}/job_${SLURM_ARRAY_TASK_ID}/job_${SLURM_ARRAY_TASK_ID}_scores.txt \
	--dbsnp /well/lindgren/UKBIOBANK/dpalmer/PRS_cell_data/data/dbsnp154_hg19.txt.gz \
	--report-html ${pgs_results_dir}/job_${SLURM_ARRAY_TASK_ID}/job_${SLURM_ARRAY_TASK_ID}_scores_report.html \
	--genotypes GT \
	--meta ${pgs_score_dir}/pgs-catalog-20221123-hg19/scores.meta.json \
	--report-csv ${pgs_results_dir}/job_${SLURM_ARRAY_TASK_ID}/job_${SLURM_ARRAY_TASK_ID}_scores_report.csv \
	--write-variants ${pgs_results_dir}/job_${SLURM_ARRAY_TASK_ID}/job_${SLURM_ARRAY_TASK_ID}_scores_variants_used.csv
