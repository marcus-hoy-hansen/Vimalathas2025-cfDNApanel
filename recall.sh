#!/bin/bash
set -x

fn=$1

#sudo apt update
#sudo apt -y install samtools

samtools index $fn'.bam'

panel=Probes_merged_ok_OHU_Bcell_clonality_v1_1X_TE-94124764_hg38.bed

gatk-4.3.0.0/gatk BaseRecalibrator \
-I $fn'.bam' \
-R hg38/resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta \
-L 'hg38/'$panel \
--known-sites hg38/resources_broad_hg38_v0_1000G_phase1.snps.high_confidence.hg38.vcf.gz \
--known-sites hg38/resources_broad_hg38_v0_Homo_sapiens_assembly38.dbsnp138.vcf \
--known-sites hg38/resources_broad_hg38_v0_Homo_sapiens_assembly38.known_indels.vcf.gz \
--known-sites hg38/resources_broad_hg38_v0_Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
-O $fn.recal.table

gatk-4.3.0.0/gatk ApplyBQSR \
-I $fn'.bam' \
-R hg38/resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta \
-L 'hg38'/$panel \
--bqsr-recal-file $fn.recal.table \
-O $fn'.recal.bam'

echo "Running Mutect2..."
gatk-4.3.0.0/gatk Mutect2 -R hg38/resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta -L 'hg38/'$panel -I $fn'.recal.bam' --germline-resource hg38/somatic-hg38_af-only-gnomad.hg38.vcf.gz --panel-of-normals hg38/somatic-hg38_1000g_pon.hg38.vcf.gz -O $fn'_somatic.vcf.gz'

gatk-4.3.0.0/gatk FilterMutectCalls --min-allele-fraction 0.02 --min-reads-per-strand 2 --unique-alt-read-count 5 -V $fn'_somatic.vcf.gz' -O $fn'_somatic_filtered.vcf' -R hg38/resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta -L 'hg38/'$panel

