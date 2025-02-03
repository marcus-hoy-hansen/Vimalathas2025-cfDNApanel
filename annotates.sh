#!/bin/bash

fn=${1%.vcf}  # Remove .vcf extension from input

# Index input VCF file
gatk-4.3.0.0/gatk IndexFeatureFile -I "$1"

# Filter out common variants
java -jar snpEff/SnpSift.jar annotate -info COMMON 00-common_all.vcf.gz "$1" | \
grep -v ";COMMON=1" | grep -v "=chrUn_" | grep -v "=HLA-" | grep -v "VQSRTrancheSNP" > "${fn}.uncommon.vcf"

# Hard filter variants below absolute minimum depth (< 30)
bcftools view -i 'INFO/DP > 29' "${fn}.uncommon.vcf" > "${fn}.dbsnp.vcf"

# Index filtered VCF
gatk-4.3.0.0/gatk IndexFeatureFile -I "${fn}.dbsnp.vcf"

# Split by chromosome (1-22)
for counter in {1..22}
do
    echo "Processing chr$counter"
    gatk-4.3.0.0/gatk SelectVariants \
        -R Homo_sapiens_assembly38.fasta \
        -V "${fn}.dbsnp.vcf" \
        -L "chr$counter" \
        -O "${fn}.chr$counter.vcf" &
done

# Process chrX and chrY in parallel
gatk-4.3.0.0/gatk SelectVariants -R Homo_sapiens_assembly38.fasta -V "${fn}.dbsnp.vcf" -L "chrX" -O "${fn}.chrX.vcf" &
gatk-4.3.0.0/gatk SelectVariants -R Homo_sapiens_assembly38.fasta -V "${fn}.dbsnp.vcf" -L "chrY" -O "${fn}.chrY.vcf" &
wait  # Ensure all chromosome files are created before proceeding

# Annotate with gnomAD
for counter in {1..22}
do
    echo "Annotating chr$counter"
    java -jar snpEff/SnpSift.jar annotate -name gnomAD. -info AF "v4/gnomad.exomes.v4.0.sites.chr$counter.vcf.bgz" \
        "${fn}.chr$counter.vcf" > "${fn}.chr$counter.gnomad.vcf" &
done

# Annotate chrX and chrY
java -jar snpEff/SnpSift.jar annotate -name gnomAD. -info AF "v4/gnomad.exomes.v4.0.sites.chrX.vcf.bgz" "${fn}.chrX.vcf" > "${fn}.chrX.gnomad.vcf" &
java -jar snpEff/SnpSift.jar annotate -name gnomAD. -info AF "v4/gnomad.exomes.v4.0.sites.chrY.vcf.bgz" "${fn}.chrY.vcf" > "${fn}.chrY.gnomad.vcf" &
wait

# Merge all annotated VCFs
gatk-4.3.0.0/gatk GatherVcfs -O "${fn}.gnomad.vcf" $(for i in {1..22} X Y; do echo "-I ${fn}.chr$i.gnomad.vcf"; done)

# Filter out common variants (gnomAD AF >= 0.01)
bcftools view -e 'INFO/gnomAD.AF >= 0.01' "${fn}.gnomad.vcf" > "${fn}.gnomad.nonSNP.vcf"

# Add ClinVar annotations
java -jar snpEff/SnpSift.jar annotate -info CLNSIG hg38/clinvar_20221211.vcf.gz "${fn}.gnomad.nonSNP.vcf" > "${fn}.gnomad.nonSNP.clinvar.vcf"

# Run snpEff
java -jar snpEff/snpEff.jar GRCh38.p13 "${fn}.gnomad.nonSNP.clinvar.vcf" > "${fn}.gnomad.nonSNP.clinvar.snpeff.vcf"

# Remove Panel of Normal (PON) variants
java -jar snpEff/SnpSift.jar annotate -info PON hg38/somatic-hg38-1000g_pon.hg38.mod.vcf.gz "${fn}.gnomad.nonSNP.clinvar.snpeff.vcf" | \
grep -v "PON=1" > "${fn}.gnomad.nonSNP.clinvar.snpeff.nonpon.vcf"

# Annotate with COSMIC
java -jar snpEff/SnpSift.jar annotate -info CNT hg38/CosmicCodingMuts.vcf.gz "${fn}.gnomad.nonSNP.clinvar.snpeff.nonpon.vcf" > "${fn}.gnomad.nonSNP.clinvar.snpeff.nonpon.cosmic.vcf"

# Annotate with dbNSFP
java -jar snpEff/SnpSift.jar dbnsfp -v -db hg38/dbNSFP4.1a.txt.gz "${fn}.gnomad.nonSNP.clinvar.snpeff.nonpon.cosmic.vcf" > "${fn}.gnomad.nonSNP.clinvar.snpeff.nonpon.cosmic.dbnsfp.vcf"

# Index final VCF
gatk-4.3.0.0/gatk IndexFeatureFile -I "${fn}.gnomad.nonSNP.clinvar.snpeff.nonpon.cosmic.dbnsfp.vcf"

# Convert to table
gatk-4.3.0.0/gatk VariantsToTable --show-filtered -V "${fn}.gnomad.nonSNP.clinvar.snpeff.nonpon.cosmic.dbnsfp.vcf" \
    -F CHROM -F POS -F REF -F ALT -F ID -F TYPE -F FILTER -F CNT -GF AD -GF AF -F ANN -O "${fn}.out.tsv"

# Filter B-cell activation genes
grep -f hg38/QuickGO-annotations_B-cell_symbols.txt "${fn}.out.tsv" > "${fn}.annotated.table.B-cellactivation.tsv"

# Organize output files
mkdir -p tmp
mv "${fn}"*.vcf tmp/
mv "${fn}"*.vcf.idx tmp/
mv "${fn}"*.tsv tmp/
mv "${fn}"*.table tmp/
mv "tmp/${fn}.gnomad.nonSNP.clinvar.snpeff.nonpon.cosmic.dbnsfp.vcf" .
