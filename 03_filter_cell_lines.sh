#!/bin/bash 
#SBATCH -J filter-cell-line
#SBATCH -A lindgren.prj 
#SBATCH -o output.out 
#SBATCH -e error.err 
#SBATCH -c 4 
#SBATCH -p short
#SBATCH --array 1-23:1 
#SBATCH --requeue

module purge
module load BCFtools

phase3_1kg_path="/well/lindgren/UKBIOBANK/dpalmer/PRS_cell_data/data/1000G/B37"
input="${phase3_1kg_path}/all_phase3_r2_0.3"
imputed_path="/well/lindgren/UKBIOBANK/dpalmer/PRS_cell_data/data/celldataB37/HRC/imputed"

chr=${SLURM_ARRAY_TASK_ID}
echo "chromosome ${chr}"

if [ $chr -eq 23 ]; then
	chr="X"
	echo "23 ${chr}" >> chr_name_conv.txt
	# bgzip -f ${input}_chr${chr}.vcf
	bcftools annotate --rename-chrs chr_name_conv.txt ${input}_chr${chr}.vcf.gz | bgzip > test.vcf.gz
	mv test.vcf.gz ${input}_chr${chr}.vcf.gz
	bcftools index -f ${input}_chr${chr}.vcf.gz
	bcftools index -f ${imputed_path}/chr${chr}.dose.vcf.gz
	bcftools isec -c none -n=2 -w1 \
	-O z -o ${imputed_path}/chr${chr}.dose.subset.vcf.gz \
	${imputed_path}/chr${chr}.dose.vcf.gz ${input}_chr${chr}.vcf.gz
	bcftools convert -h ${input}_test --haploid2diploid ${input}_chr${chr}.vcf.gz
	mkdir ${phase3_1kg_path}/original_X
	cp ${input}_chr${chr}.vcf.gz ${phase3_1kg_path}/original_X/all_phase3_r2_0.3_chr${chr}_hap_and_dip.vcf.gz
	cp ${input}_chr${chr}.vcf.gz.csi ${phase3_1kg_path}/original_X/all_phase3_r2_0.3_chr${chr}_hap_and_dip.vcf.gz.csi
	bcftools convert -H ${input}_test -o test.vcf.gz -O z
	bcftools annotate --set-id +'%CHROM:%POS:%REF:%FIRST_ALT' test.vcf.gz -o ${input}_chr${chr}.vcf.gz -O z
	bcftools index -f ${input}_chr${chr}.vcf.gz
	bcftools index -f ${imputed_path}/chr${chr}.dose.subset.vcf.gz
else
	# bgzip -f ${input}_chr${chr}.vcf
	bcftools index -f ${input}_chr${chr}.vcf.gz
	bcftools index -f ${imputed_path}/chr${chr}.dose.vcf.gz
	bcftools isec -c none -n=2 -w1 \
	-O z -o ${imputed_path}/chr${chr}.dose.subset.vcf.gz \
	${imputed_path}/chr${chr}.dose.vcf.gz ${input}_chr${chr}.vcf.gz
	bcftools index -f ${imputed_path}/chr${chr}.dose.subset.vcf.gz
fi
