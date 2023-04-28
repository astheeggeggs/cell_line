# Instructions

export PATH="/well/lindgren/dpalmer/:$PATH"
imputed_path="/well/lindgren/UKBIOBANK/dpalmer/PRS_cell_data/data/celldataB37/HRC/imputed"
phase3_1kg_path="/well/lindgren/UKBIOBANK/dpalmer/PRS_cell_data/data/1000G/B37"
module purge
module load BCFtools

# Get the 1000G data and clean it up
# https://alanaw1.github.io/post/2021/03/03/visualizing-1000-genomes-data/

# Build 38 and Build 37
# https://www.cog-genomics.org/plink/2.0/resources

# Follow the instructions
# Build 37

cd /well/lindgren/UKBIOBANK/dpalmer/PRS_cell_data/data/1000G/B37
mv phase3_corrected.psam all_phase3.psam
plink2 --zst-decompress all_phase3.pgen.zst > all_phase3.pgen
plink2 --zst-decompress all_phase3.pvar.zst > all_phase3.pvar
plink2 --make-bed --pfile all_phase3 --max-alleles 2 --out all_phase3
rm all_phase3.pgen
rm all_phase3.pvar

# Determine the intersection of the 1000G files with the cell line files for variants with R2 > 0.3
# To merge with plink, get chr:pos:ref:alt for the bim file
awk '$2=$1":"$4":"$6":"$5' all_phase3.bim > all_phase3_chr_pos_ref_alt.bim

# Extract the variant ID from vcf files to restrict the plink file
bcftools query -f '%CHROM:%POS:%REF:%ALT{0}\n' ${imputed_path}/chr1.dose.vcf.gz > ${imputed_path}/tmp
# Need to unzip the X chr
for chr in {2..22}; do
	bcftools query -f '%CHROM:%POS:%REF:%ALT{0}\n' ${imputed_path}/chr${chr}.dose.vcf.gz >> ${imputed_path}/tmp
done

bcftools query -f '%CHROM:%POS:%REF:%ALT{0}\n' ${imputed_path}/chrX.dose.vcf.gz >> ${imputed_path}/tmp

plink --extract ${imputed_path}/tmp --make-bed \
	--bed ${phase3_1kg_path}/all_phase3.bed \
	--bim ${phase3_1kg_path}/all_phase3_chr_pos_ref_alt.bim \
	--fam ${phase3_1kg_path}/all_phase3.fam \
	--out ${phase3_1kg_path}/all_phase3_r2_0.3 \
	--allow-extra-chr --keep-allele-order

# Now, send these files to vcfs, I can pilfer the commands output from HRC-1000G-check-bim.pl to ensure this will be correct
input="${phase3_1kg_path}/all_phase3_r2_0.3"

for chr in {1..22}; do
	plink --bfile ${input} --real-ref-alleles --chr ${chr} --out ${input}_chr${chr} --make-bed
	plink --bfile ${input} --real-ref-alleles --recode vcf --chr ${chr} --out ${input}_chr${chr}
done

chr="X"
plink --bfile ${input} --real-ref-alleles --recode vcf --chr ${chr} --out ${input}_chr${chr}
plink --bfile ${input} --real-ref-alleles --chr ${chr} --out ${input}_chr${chr} --make-bed

# for chr in {1..22}; do
# 	echo "chromosome ${chr}"
# 	bgzip -f ${input}_chr${chr}.vcf
# 	bcftools index -f ${input}_chr${chr}.vcf.gz
# 	bcftools index -f ${imputed_path}/chr${chr}.dose.vcf.gz
# 	bcftools isec -c none -n=2 -w1 \
# 	-O z -o ${imputed_path}/chr${chr}.dose.subset.vcf.gz \
# 	${input}_chr${chr}.vcf.gz ${imputed_path}/chr${chr}.dose.vcf.gz
# done

# chr="X"
# echo "23 ${chr}" >> chr_name_conv.txt
# bgzip -f ${input}_chr${chr}.vcf
# bcftools annotate --rename-chrs chr_name_conv.txt ${input}_chr${chr}.vcf.gz | bgzip > test.vcf.gz
# mv test.vcf.gz ${input}_chr${chr}.vcf.gz
# bcftools index -f ${input}_chr${chr}.vcf.gz
# bcftools index -f ${imputed_path}/chr${chr}.dose.vcf.gz
# bcftools isec -c none -n=2 -w1 \
# -O z -o ${imputed_path}/chr${chr}.dose.subset.vcf.gz \
# ${input}_chr${chr}.vcf.gz ${imputed_path}/chr${chr}.dose.vcf.gz

