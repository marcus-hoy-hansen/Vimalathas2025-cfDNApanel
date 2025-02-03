#!/bin/bash
set -x  # Debugging mode (prints each command before execution)

fn=$1

# Define paths
PANEL="hg38/Probes_merged_ok_OUH_Bcell_clonality_v1_1X_TE-94124764_hg38.bed"
REF_GENOME="hg38/resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta"
GATK="gatk-4.3.0.0/gatk"

# Ensure samtools is installed
if ! command -v samtools &> /dev/null; then
    echo "Installing samtools"
    sudo apt update && sudo apt -y install samtools
fi

# Index BAM file
samtools index "${fn}.bam"

# Base Recalibration
$GATK BaseRecalibrator \
    -I "${fn}.bam" \
    -R "$REF_GENOME" \
    -L "$PANEL" \
    --known-sites hg38/resources_broad_hg38_v0_1000G_phase1.snps.high_confidence.hg38.vcf.gz \
    --known-sites hg38/resources_broad_hg38_v0_Homo_sapiens_assembly38.dbsnp138.vcf \
    --known-sites hg38/resources_broad_hg38_v0_Homo_sapiens_assembly38.known_indels.vcf.gz \
    --known-sites hg38/resources_broad_hg38_v0_Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
    -O "${fn}.recal.table"

$GATK ApplyBQSR \
    -I "${fn}.bam" \
    -R "$REF_GENOME" \
    -L "$PANEL" \
    --bqsr-recal-file "${fn}.recal.table" \
    -O "${fn}.recal.bam"

# Run Mutect2
echo "Running Mutect2..."
$GATK Mutect2 \
    -R "$REF_GENOME" \
    -L "$PANEL" \
    -I "${fn}.recal.bam" \
    --germline-resource hg38/somatic-hg38_af-only-gnomad.hg38.vcf.gz \
    --panel-of-normals hg38/somatic-hg38_1000g_pon.hg38.vcf.gz \
    -O "${fn}_somatic.vcf.gz" | tee "${fn}_mutect2.log"

# Filter Mutect Calls
$GATK FilterMutectCalls \
    --min-allele-fraction 0.01 \
    --min-reads-per-strand 2 \
    --unique-alt-read-count 5 \
    -V "${fn}_somatic.vcf.gz" \
    -O "${fn}_somatic_filtered.vcf.gz" \
    -R "$REF_GENOME" \
    -L "$PANEL"
