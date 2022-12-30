# Instructions

export PATH="/well/lindgren/dpalmer/:$PATH"

# Get the 1000G data and clean it up
# https://alanaw1.github.io/post/2021/03/03/visualizing-1000-genomes-data/

# Build 38 and Build 37
# https://www.cog-genomics.org/plink/2.0/resources

# Follow the instructions
# Build 37
cd B37
mv phase3_corrected.psam all_phase3.psam
plink2 --zst-decompress all_phase3.pgen.zst > all_phase3.pgen
plink2 --zst-decompress all_phase3.pvar.zst > all_phase3.pvar
plink2 --make-bed --pfile all_phase3 --max-alleles 2 --out all_phase3
rm all_phase3.pgen
rm all_phase3.pvar
cd ..

# Determine the intersection of the 1000G files with the cell line files for variants with R2 > 0.3
# To merge with plink, get chr:pos:ref:alt for the bim file
awk '$2=$1":"$4":"$6":"$5' all_phase3.bim > all_phase3_chr_pos_ref_alt.bim

# Extract the variant ID from vcf files to restrict the plink file
bcftools query -f '%CHROM:%POS:%REF:%ALT{0}\n' /well/lindgren/UKBIOBANK/dpalmer/PRS_cell_data/data/celldataB37/imputed/HRC/chr1.dose.vcf.gz > /well/lindgren/UKBIOBANK/dpalmer/PRS_cell_data/data/celldataB37/imputed/HRC/tmp
# Need to unzip the X chr
for chr in {2..22}
do
	bcftools query -f '%CHROM:%POS:%REF:%ALT{0}\n' /well/lindgren/UKBIOBANK/dpalmer/PRS_cell_data/data/celldataB37/imputed/HRC/chr${chr}.dose.vcf.gz >> /well/lindgren/UKBIOBANK/dpalmer/PRS_cell_data/data/celldataB37/imputed/HRC/tmp
done

bcftools query -f '%CHROM:%POS:%REF:%ALT{0}\n' /well/lindgren/UKBIOBANK/dpalmer/PRS_cell_data/data/celldataB37/imputed/HRC/chrX.dose.vcf.gz >> /well/lindgren/UKBIOBANK/dpalmer/PRS_cell_data/data/celldataB37/imputed/HRC/tmp

plink --extract /well/lindgren/UKBIOBANK/dpalmer/PRS_cell_data/data/celldataB37/imputed/HRC/tmp --make-bed --bed /well/lindgren/UKBIOBANK/dpalmer/PRS_cell_data/data/1000G/B37/all_phase3.bed --bim /well/lindgren/UKBIOBANK/dpalmer/PRS_cell_data/data/1000G/B37/all_phase3_chr_pos_ref_alt.bim --fam /well/lindgren/UKBIOBANK/dpalmer/PRS_cell_data/data/1000G/B37/all_phase3.fam --out /well/lindgren/UKBIOBANK/dpalmer/PRS_cell_data/data/1000G/B37/all_phase3_r2_0.3 --allow-extra-chr --keep-allele-order

# Now, send these files to vcfs, I can pilfer the commands output from HRC-1000G-check-bim.pl to ensure this will be correct
input="/well/lindgren/UKBIOBANK/dpalmer/PRS_cell_data/data/1000G/B37/all_phase3_r2_0.3"
for chr in {1..22}
do
	plink --bfile ${input} --real-ref-alleles --chr ${chr} --out ${input}_chr${chr} --make-bed
	plink --bfile ${input} --real-ref-alleles --recode vcf --chr ${chr} --out ${input}_chr${chr}
done

chr="X"
plink --bfile ${input} --real-ref-alleles --recode vcf --chr ${chr} --out ${input}_chr${chr}
plink --bfile ${input} --real-ref-alleles --chr ${chr} --out ${input}_chr${chr} --make-bed

for chr in {1..22}
do
    bgzip ${input}_chr${chr}.vcf
done

chr="X"
bgzip ${input}_chr${chr}.vcf

# Combine the bim files
chr="1"
cat ${input}_chr${chr}.bim > ${input}_variants.txt
for chr in {2..22}
do
	cat ${input}_chr${chr}.bim >> ${input}_variants.txt
done

chr="X"
awk -v OFS='\t' '{print "X\t"$2,$3,$4,$5, $6}' ${input}_chr${chr}.bim >> ${input}_variants.txt

awk -v OFS='\t' '{print $1, $4}' ${input}_variants.txt > ${input}_variants_bcftools.txt

for chr in {1..22}
do
	echo "chromosome ${chr}"
	bcftools index /well/lindgren/UKBIOBANK/dpalmer/PRS_cell_data/data/celldataB37/imputed/HRC/chr${chr}.dose.vcf.gz
done

chr="X"
echo "chromosome ${chr}"
bcftools index /well/lindgren/UKBIOBANK/dpalmer/PRS_cell_data/data/celldataB37/imputed/HRC/chr${chr}.dose.vcf.gz

for chr in {1..22}
do
	echo "chromosome ${chr}"
	bcftools filter -R ${input}_variants_bcftools.txt /well/lindgren/UKBIOBANK/dpalmer/PRS_cell_data/data/celldataB37/imputed/HRC/chr${chr}.dose.vcf.gz -o /well/lindgren/UKBIOBANK/dpalmer/PRS_cell_data/data/celldataB37/imputed/HRC/chr${chr}.dose.subset.vcf.gz
done

chr="X"
echo "chromosome ${chr}"
bcftools filter -R ${input}_variants_bcftools.txt /well/lindgren/UKBIOBANK/dpalmer/PRS_cell_data/data/celldataB37/imputed/HRC/chr${chr}.dose.vcf.gz -o /well/lindgren/UKBIOBANK/dpalmer/PRS_cell_data/data/celldataB37/imputed/HRC/chr${chr}.dose.subset.vcf.gz
