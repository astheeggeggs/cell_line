#!/bin/bash 
#SBATCH -J extract-vcf-ukb-subset
#SBATCH -A lindgren.prj 
#SBATCH -o output-extract-vcf.out 
#SBATCH -e error-extract-vcf.err 
#SBATCH -c 4 
#SBATCH -p short
#SBATCH --array 1-24:1 
#SBATCH --requeue

export PATH="/well/lindgren/dpalmer/:$PATH"
phase3_1kg_path="/well/lindgren/UKBIOBANK/dpalmer/PRS_cell_data/data/1000G/B37"
ukb_imputed_subset_dir="/well/lindgren/UKBIOBANK/dpalmer/PRS_cell_data/data/UKB_imputed_subset"
chr=${SLURM_ARRAY_TASK_ID}
sample_file="/well/lindgren-ukbb/projects/ukbb-11867/DATA/SAMPLE_FAM/ukb11867_imp_chr1_v3_s487395.sample"

if [ $chr -eq 23 ]; then
	chr="X"
	sample_file="/well/lindgren-ukbb/projects/ukbb-11867/DATA/SAMPLE_FAM/ukb11867_imp_chrX_v3_s486743.sample"
fi

if [ $chr -eq 24 ]; then
	chr="XY"
	sample_file="/well/lindgren-ukbb/projects/ukbb-11867/DATA/SAMPLE_FAM/ukb11867_imp_chrXY_v3_s486429.sample"
fi

module purge
module load BCFtools

awk '{print $1":"$4"_"$6"_"$5}' \
	${phase3_1kg_path}/all_phase3_r2_0.3_chr${chr}.bim \
	> ${phase3_1kg_path}/all_phase3_r2_0.3_chr${chr}_variants.txt

# Filter to subset of UK Biobank using qctool
/well/lindgren/dpalmer/qctool2/qctool/build/release/apps/qctool_v2.2.0 \
	-g /well/lindgren-ukbb/projects/ukbb-11867/DATA/IMPUTATION/ukb_imp_chr${chr}_v3.bgen \
	-s ${sample_file} \
	-incl-snpids ${phase3_1kg_path}/all_phase3_r2_0.3_chr${chr}_variants.txt \
	-incl-samples /well/lindgren/UKBIOBANK/dpalmer/PRS_cell_data/data/UKB_subset_one_col.txt \
	-og ${ukb_imputed_subset_dir}/UKB_imputed_subset_chr${chr}.vcf.gz

# Force hard-calls
/well/lindgren/dpalmer/qctool2/qctool/build/release/apps/qctool_v2.2.0 \
	-g ${ukb_imputed_subset_dir}/UKB_imputed_subset_chr${chr}.vcf.gz \
	-threshold 0.5001 \
	-og ${ukb_imputed_subset_dir}/UKB_imputed_subset_chr${chr}_thresholded_05.vcf.gz

# Convert to plink, and fill in missing calls
plink \
	--vcf ${ukb_imputed_subset_dir}/UKB_imputed_subset_chr${chr}_thresholded_05.vcf.gz \
	--real-ref-alleles --fill-missing-a2 --make-bed \
	--out ${ukb_imputed_subset_dir}/UKB_imputed_subset_chr${chr}_thresholded_recode

# Fill in missing calls in the 'UKB_subset_combined' 
# (these are the genotyped files that are not present in the UKB imputed data)
plink \
	--bfile /well/lindgren/UKBIOBANK/dpalmer/PRS_cell_data/data/combined/UKB_subset_combined-updated-chr${chr} \
	--make-bed --fill-missing-a2 --real-ref-alleles \
	--out ${ukb_imputed_subset_dir}/UKB_subset_combined-updated-chr${chr}-no-missing

# Send to compressed vcf
plink \
	--bfile ${ukb_imputed_subset_dir}/UKB_subset_combined-updated-chr${chr}-no-missing \
	--recode vcf bgz --real-ref-alleles \
	--out ${ukb_imputed_subset_dir}/UKB_subset_combined-updated-chr${chr}-no-missing

# Merge the imputed data with the genotype (array) data
plink \
	--bfile ${ukb_imputed_subset_dir}/UKB_imputed_subset_chr${chr}_thresholded_recode \
	--bmerge ${ukb_imputed_subset_dir}/UKB_subset_combined-updated-chr${chr}-no-missing \
	--real-ref-alleles --make-bed \
	--out ${ukb_imputed_subset_dir}/UKB_imputed_subset_chr${chr}_thresholded_recode_merge_typed

# Remove missing data
plink \
	--bfile ${ukb_imputed_subset_dir}/UKB_imputed_subset_chr${chr}_thresholded_recode_merge_typed \
	--fill-missing-a2 --real-ref-alleles --make-bed \
	--out ${ukb_imputed_subset_dir}/UKB_imputed_subset_chr${chr}_thresholded_recode_merge_typed

# Recode the bim files to be in the right format
awk -v OFS='\t' '{print $1,$1":"$4":"$6":"$5,$3,$4,$5,$6}' \
${ukb_imputed_subset_dir}/UKB_imputed_subset_chr${chr}_thresholded_recode_merge_typed.bim \
> ${ukb_imputed_subset_dir}/UKB_imputed_subset_chr${chr}_thresholded_recode_merge_typed.bim.tmp

mv ${ukb_imputed_subset_dir}/UKB_imputed_subset_chr${chr}_thresholded_recode_merge_typed.bim.tmp \
	${ukb_imputed_subset_dir}/UKB_imputed_subset_chr${chr}_thresholded_recode_merge_typed.bim

# Then recode to vcf
plink \
	--double-id \
	--bfile ${ukb_imputed_subset_dir}/UKB_imputed_subset_chr${chr}_thresholded_recode_merge_typed \
	--real-ref-alleles --recode vcf bgz \
	--out ${ukb_imputed_subset_dir}/UKB_imputed_subset_chr${chr}_thresholded_recode_merge_typed_final
