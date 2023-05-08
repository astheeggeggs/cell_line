#!/bin/bash 
#SBATCH -J intersect-vcf
#SBATCH -A lindgren.prj 
#SBATCH -o intersect-output.out 
#SBATCH -e intersect-error.err 
#SBATCH -c 4 
#SBATCH -p short
#SBATCH --array 1-23:1 
#SBATCH --requeue

module purge
module load BCFtools

imputed_path="/well/lindgren/UKBIOBANK/dpalmer/PRS_cell_data/data/celldataB37/HRC/imputed"
phase3_1kg_path="/well/lindgren/UKBIOBANK/dpalmer/PRS_cell_data/data/1000G/B37"
ukb_imputed_subset_dir="/well/lindgren/UKBIOBANK/dpalmer/PRS_cell_data/data/UKB_imputed_subset"

chr=${SLURM_ARRAY_TASK_ID}

# if [ $chr -eq 23 ]; then
# 	chr="X"
# fi

# bcftools isec -c none -n=3 \
# -p ${ukb_imputed_subset_dir}/UKB_1000G_cell_line_intersection_chr${chr} -O z \
# ${ukb_imputed_subset_dir}/UKB_imputed_subset_chr${chr}_thresholded_recode_merge_typed_final.vcf.gz \
# ${phase3_1kg_path}/all_phase3_r2_0.3_chr${chr}.vcf.gz \
# ${imputed_path}/chr${chr}.dose.vcf.gz

mkdir -p ${ukb_imputed_subset_dir}/UKB_cell_line_1000G_intersection
# Move the vcfs
mv ${ukb_imputed_subset_dir}/UKB_1000G_cell_line_intersection_chr${chr}/0000.vcf.gz \
${ukb_imputed_subset_dir}/UKB_cell_line_1000G_intersection/UKB_imputed_subset_chr${chr}_thresholded_recode_merge_typed_final_UKB_cell_line_1000G_intersected.vcf.gz
mv ${ukb_imputed_subset_dir}/UKB_1000G_cell_line_intersection_chr${chr}/0001.vcf.gz \
${phase3_1kg_path}/all_phase3_r2_0.3_chr${chr}_UKB_cell_line_1000G_intersected.vcf.gz
mv ${ukb_imputed_subset_dir}/UKB_1000G_cell_line_intersection_chr${chr}/0002.vcf.gz \
${imputed_path}/chr${chr}_UKB_cell_line_1000G_intersected.dose.vcf.gz

# Move the index files
mv ${ukb_imputed_subset_dir}/UKB_1000G_cell_line_intersection_chr${chr}/0000.vcf.gz.tbi \
${ukb_imputed_subset_dir}/UKB_cell_line_1000G_intersection/UKB_imputed_subset_chr${chr}_thresholded_recode_merge_typed_final_UKB_cell_line_1000G_intersected.vcf.gz.tbi
mv ${ukb_imputed_subset_dir}/UKB_1000G_cell_line_intersection_chr${chr}/0001.vcf.gz.tbi \
${phase3_1kg_path}/all_phase3_r2_0.3_chr${chr}_UKB_cell_line_1000G_intersected.vcf.gz.tbi
mv ${ukb_imputed_subset_dir}/UKB_1000G_cell_line_intersection_chr${chr}/0002.vcf.gz.tbi \
${imputed_path}/chr${chr}_UKB_cell_line_1000G_intersected.dose.vcf.gz.tbi
