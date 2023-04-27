#!/bin/bash 
#SBATCH -J filter-cell-line
#SBATCH -A lindgren.prj 
#SBATCH -o output.out 
#SBATCH -e error.err 
#SBATCH -c 4 
#SBATCH -p short
#SBATCH --array 13-22:1 
#SBATCH --requeue

module purge
module load BCFtools

phase3_1kg_path="/well/lindgren/UKBIOBANK/dpalmer/PRS_cell_data/data/1000G/B37"
input="${phase3_1kg_path}/all_phase3_r2_0.3"
imputed_path="/well/lindgren/UKBIOBANK/dpalmer/PRS_cell_data/data/celldataB37/imputed/HRC"

chr=${SLURM_ARRAY_TASK_ID}
echo "chromosome ${chr}"

bgzip -f ${input}_chr${chr}.vcf
bcftools index -f ${input}_chr${chr}.vcf.gz
bcftools index -f ${imputed_path}/chr${chr}.dose.vcf.gz
bcftools isec -c none -n=2 -w1 \
-O z -o ${imputed_path}/chr${chr}.dose.subset.vcf.gz \
${input}_chr${chr}.vcf.gz ${imputed_path}/chr${chr}.dose.vcf.gz
