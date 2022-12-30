#!/bin/bash 
#SBATCH -J filter-cell-line
#SBATCH -A lindgren.prj 
#SBATCH -o output.out 
#SBATCH -e error.err 
#SBATCH -c 4 
#SBATCH -p short
#SBATCH --array 1-22:1 
#SBATCH --requeue

module purge
module load BCFtools

input="/well/lindgren/UKBIOBANK/dpalmer/PRS_cell_data/data/1000G/B37/all_phase3_r2_0.3"
chr=${SLURM_ARRAY_TASK_ID}
echo "chromosome ${chr}"
bcftools filter -R ${input}_variants_bcftools.txt /well/lindgren/UKBIOBANK/dpalmer/PRS_cell_data/data/celldataB37/imputed/HRC/chr${chr}.dose.vcf.gz -o /well/lindgren/UKBIOBANK/dpalmer/PRS_cell_data/data/celldataB37/imputed/HRC/chr${chr}.dose.subset.vcf.gz
