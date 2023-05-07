#!/bin/bash 
#SBATCH -J extract-vcf-ukb-subset
#SBATCH -A lindgren.prj 
#SBATCH -o output-extract-vcf.out 
#SBATCH -e error-extract-vcf.err 
#SBATCH -c 4 
#SBATCH -p short
#SBATCH --array 1-23:1 
#SBATCH --requeue

export PATH="/well/lindgren/dpalmer/:$PATH"
imputed_path="/well/lindgren/UKBIOBANK/dpalmer/PRS_cell_data/data/celldataB37/HRC/imputed"
ukb_imputed_subset_dir="/well/lindgren/UKBIOBANK/dpalmer/PRS_cell_data/data/UKB_imputed_subset"
chr=${SLURM_ARRAY_TASK_ID}
chr_dose=${SLURM_ARRAY_TASK_ID}
chr_plink=${SLURM_ARRAY_TASK_ID}
sample_file="/well/lindgren-ukbb/projects/ukbb-11867/DATA/SAMPLE_FAM/ukb11867_imp_chr1_v3_s487395.sample"

# if [ $chr_plink -eq 23 ]; then
# 	chr="X"
# 	chr_dose="X"
# 	sample_file="/well/lindgren-ukbb/projects/ukbb-11867/DATA/SAMPLE_FAM/ukb11867_imp_chrX_v3_s486743.sample"
# fi

# module purge
# module load BCFtools

# # Extract the variant ID from vcf files to restrict the plink file
# bcftools query -f '%CHROM:%POS\_%REF\_%ALT{0}\n' ${imputed_path}/chr${chr_dose}.dose.vcf.gz > ${imputed_path}/tmp_${chr}
# module purge
# module load OpenBLAS/0.3.12-GCC-10.2.0 

# # Filter to subset of UK Biobank using qctool
# /well/lindgren/dpalmer/qctool2/qctool/build/release/apps/qctool_v2.2.0 \
# -g /well/lindgren-ukbb/projects/ukbb-11867/DATA/IMPUTATION/ukb_imp_chr${chr}_v3.bgen \
# -s ${sample_file} \
# -incl-snpids ${imputed_path}/tmp_${chr} \
# -incl-samples /well/lindgren/UKBIOBANK/dpalmer/PRS_cell_data/data/UKB_subset_one_col.txt \
# -og ${ukb_imputed_subset_dir}/UKB_imputed_subset_chr${chr}.vcf.gz

# rm ${imputed_path}/tmp_${chr}

# # Force hard-calls
# /well/lindgren/dpalmer/qctool2/qctool/build/release/apps/qctool_v2.2.0 \
# -g ${ukb_imputed_subset_dir}/UKB_imputed_subset_chr${chr}.vcf.gz \
# -threshold 0.5001 \
# -og ${ukb_imputed_subset_dir}/UKB_imputed_subset_chr${chr}_thresholded_05.vcf.gz

# # Convert to plink, and fill in missing calls
# plink \
# --vcf ${ukb_imputed_subset_dir}/UKB_imputed_subset_chr${chr}_thresholded_05.vcf.gz \
# --real-ref-alleles --fill-missing-a2 --make-bed \
# --out ${ukb_imputed_subset_dir}/UKB_imputed_subset_chr${chr}_thresholded_recode

# # Fill in missing calls in the 'UKB_subset_combined' - need to create this for the X
# # (these are the genotyped files that are not present in the UKB imputed data)
# plink \
# --bfile /well/lindgren/UKBIOBANK/dpalmer/PRS_cell_data/data/combined/UKB_subset_combined-updated-chr${chr_plink} \
# --make-bed --fill-missing-a2 --real-ref-alleles \
# --out ${ukb_imputed_subset_dir}/UKB_subset_combined-updated-chr${chr}-no-missing

# if [ ! -f $unphased ]; then
# 	# Remove the variants that are not phased, and the three strange chromosome 9 variants.
# 	awk '{ if ($NF == 0) { print $4":"$5":"$6":"$7} }' /well/lindgren-ukbb/projects/ukbb-11867/DATA/QC/ukb_snp_qc.txt > ${ukb_imputed_subset_dir}/unphased_variants.tsv
# 	awk '{ print $4":"$5":"$6":"$7 }' /well/lindgren-ukbb/projects/ukbb-11867/DATA/QC/ukb_snp_qc.txt | tail -n +2 > ${ukb_imputed_subset_dir}/array_variants.tsv
# 	# Add the three chromosome 9 variants
# 	echo -e "9:139888121:C:T\n9:140086963:G:A\n9:140225101:C:T\n" >> ${ukb_imputed_subset_dir}/unphased_variants.tsv
# 	echo -e "9:139888121:C:T\n9:140086963:G:A\n9:140225101:C:T\n" >> ${ukb_imputed_subset_dir}/array_variants.tsv
# fi

# Rscript 08_ukb_fix_edge_case.r

# # Recode the bim files to be in the right format
# awk -v OFS='\t' '{print $1,$1":"$4":"$6":"$5,$3,$4,$5,$6}' \
# ${ukb_imputed_subset_dir}/UKB_subset_combined-updated-chr${chr}-no-missing.bim \
# > ${ukb_imputed_subset_dir}/UKB_subset_combined-updated-chr${chr}-no-missing.tmp

# mv ${ukb_imputed_subset_dir}/UKB_subset_combined-updated-chr${chr}-no-missing.tmp \
# ${ukb_imputed_subset_dir}/UKB_subset_combined-updated-chr${chr}-no-missing.bim

# plink \
# --bfile ${ukb_imputed_subset_dir}/UKB_subset_combined-updated-chr${chr}-no-missing \
# --exclude ${ukb_imputed_subset_dir}/unphased_and_fixed_variants.tsv \
# --real-ref-alleles --make-bed \
# --out ${ukb_imputed_subset_dir}/UKB_subset_combined-updated-chr${chr}-no-missing-no-unphased

# # Merge the imputed data with the genotype (array) data
# plink \
# --bfile ${ukb_imputed_subset_dir}/UKB_imputed_subset_chr${chr}_thresholded_recode \
# --bmerge ${ukb_imputed_subset_dir}/UKB_subset_combined-updated-chr${chr}-no-missing-no-unphased \
# --real-ref-alleles --make-bed \
# --out ${ukb_imputed_subset_dir}/UKB_imputed_subset_chr${chr}_thresholded_recode_merge_typed

# # Remove missing data
# plink \
# --bfile ${ukb_imputed_subset_dir}/UKB_imputed_subset_chr${chr}_thresholded_recode_merge_typed \
# --fill-missing-a2 --real-ref-alleles --make-bed \
# --out ${ukb_imputed_subset_dir}/UKB_imputed_subset_chr${chr}_thresholded_recode_merge_typed

# if [ $chr_plink -eq 23 ]; then
# 	# Recode the bim files to be in the right format
# 	awk -v OFS='\t' '{print "X","X:"$4":"$6":"$5,$3,$4,$5,$6}' \
# 	${ukb_imputed_subset_dir}/UKB_imputed_subset_chr${chr}_thresholded_recode_merge_typed.bim \
# 	> ${ukb_imputed_subset_dir}/UKB_imputed_subset_chr${chr}_thresholded_recode_merge_typed.bim.tmp
# else
# 	# Recode the bim files to be in the right format
# 	awk -v OFS='\t' '{print $1,$1":"$4":"$6":"$5,$3,$4,$5,$6}' \
# 	${ukb_imputed_subset_dir}/UKB_imputed_subset_chr${chr}_thresholded_recode_merge_typed.bim \
# 	> ${ukb_imputed_subset_dir}/UKB_imputed_subset_chr${chr}_thresholded_recode_merge_typed.bim.tmp
# fi

# mv ${ukb_imputed_subset_dir}/UKB_imputed_subset_chr${chr}_thresholded_recode_merge_typed.bim.tmp \
# ${ukb_imputed_subset_dir}/UKB_imputed_subset_chr${chr}_thresholded_recode_merge_typed.bim

# # Then recode to vcf
# plink \
# --double-id \
# --bfile ${ukb_imputed_subset_dir}/UKB_imputed_subset_chr${chr}_thresholded_recode_merge_typed \
# --real-ref-alleles --recode vcf bgz \
# --out ${ukb_imputed_subset_dir}/UKB_imputed_subset_chr${chr}_thresholded_recode_merge_typed_final

# # These are a lot smaller but there will still be a lot of differences between the
# # variants in each of the VCFs. There will be a subset of variants present in the 
# # UKB data because they are present in the array data for UK-Biobank, but not present in the cell-line data

# # Finally, we need to intersect the UK Biobank VCFs and the cell-lines VCFs
# # These are then passed for PGS scoring.

# module purge
# module load BCFtools

imputed_path="/well/lindgren/UKBIOBANK/dpalmer/PRS_cell_data/data/celldataB37/HRC/imputed"

# if [ $chr_plink -eq 23 ]; then
# 	chr="X"
# 	echo "23 ${chr}" >> chr_name_conv.txt
# 	bcftools annotate --rename-chrs chr_name_conv.txt ${ukb_imputed_subset_dir}/UKB_imputed_subset_chr${chr}_thresholded_recode_merge_typed_final.vcf.gz | bgzip > test.vcf.gz
# 	mv test.vcf.gz ${ukb_imputed_subset_dir}/UKB_imputed_subset_chr${chr}_thresholded_recode_merge_typed_final.vcf.gz
# 	bcftools index -f ${ukb_imputed_subset_dir}/UKB_imputed_subset_chr${chr}_thresholded_recode_merge_typed_final.vcf.gz
# 	bcftools index -f ${imputed_path}/chr${chr}.dose.vcf.gz
# 	bcftools isec -c none -n=2 \
# 	-p ${ukb_imputed_subset_dir}/UKB_cell_line_intersection_chr${chr} -O z  \
# 	${ukb_imputed_subset_dir}/UKB_imputed_subset_chr${chr}_thresholded_recode_merge_typed_final.vcf.gz \
# 	${imputed_path}/chr${chr}.dose.vcf.gz
# else
# 	bcftools index -f ${ukb_imputed_subset_dir}/UKB_imputed_subset_chr${chr}_thresholded_recode_merge_typed_final.vcf.gz
# 	bcftools index -f ${imputed_path}/chr${chr}.dose.vcf.gz
# 	bcftools isec -c none -n=2 \
# 	-p ${ukb_imputed_subset_dir}/UKB_cell_line_intersection_chr${chr} -O z \
# 	${ukb_imputed_subset_dir}/UKB_imputed_subset_chr${chr}_thresholded_recode_merge_typed_final.vcf.gz \
# 	${imputed_path}/chr${chr}.dose.vcf.gz
# fi

# mkdir -p ${ukb_imputed_subset_dir}/UKB_cell_line_intersection
# # Move the vcfs
# mv ${ukb_imputed_subset_dir}/UKB_cell_line_intersection_chr${chr}/0000.vcf.gz  \
# ${ukb_imputed_subset_dir}/UKB_cell_line_intersection/UKB_imputed_subset_chr${chr}_thresholded_recode_merge_typed_final_intersected.vcf.gz
mv ${ukb_imputed_subset_dir}/UKB_cell_line_intersection_chr${chr}/0001.vcf.gz \
${imputed_path}/chr${chr}_UKB_intersected.dose.vcf.gz

# Move the index files
# mv ${ukb_imputed_subset_dir}/UKB_cell_line_intersection_chr${chr}/0000.vcf.gz.tbi  \
# ${ukb_imputed_subset_dir}/UKB_cell_line_intersection/UKB_imputed_subset_chr${chr}_thresholded_recode_merge_typed_final_intersected.vcf.gz.tbi
mv ${ukb_imputed_subset_dir}/UKB_cell_line_intersection_chr${chr}/0001.vcf.gz.tbi \
${imputed_path}/chr${chr}_UKB_intersected.dose.vcf.gz.tbi


