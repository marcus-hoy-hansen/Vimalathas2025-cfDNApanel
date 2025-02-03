#!/bin/bash

# Define number of threads
threads=8
SAMPLE=$1
REF_GENOME="hg38/resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta"

# Convert FASTQ to BAM
gatk-4.3.0.0/gatk FastqToSam \
-F1 "${SAMPLE}_R1_001.fastq.gz" \
-F2 "${SAMPLE}_R2_001.fastq.gz" \
-O "${SAMPLE}_unmapped.bam" \
-SM "$SAMPLE"

# Extract UMIs
fgbio ExtractUmisFromBam \
-i "${SAMPLE}_unmapped.bam" \
-o "${SAMPLE}_unmapped_umi_extracted.bam" \
-r 5M2S+T 5M2S+T \
-t RX \
-a true

# Convert BAM to FASTQ
gatk-4.3.0.0/gatk SamToFastq \
-I "${SAMPLE}_unmapped_umi_extracted.bam" \
-F "${SAMPLE}_umi_extracted_R1.fastq" \
-F2 "${SAMPLE}_umi_extracted_R2.fastq" \
--CLIPPING_ATTRIBUTE XT \
--CLIPPING_ACTION 2

# Trim reads using fastp
fastp \
-i "${SAMPLE}_umi_extracted_R1.fastq" \
-o "${SAMPLE}_umi_extracted_trimmed_R1.fastq" \
-I "${SAMPLE}_umi_extracted_R2.fastq" \
-O "${SAMPLE}_umi_extracted_trimmed_R2.fastq" \
-g -W 5 -q 30 -u 40 -x -3 -l 75 -c -e 30 \
-j "${SAMPLE}.fastp.json" \
-h "${SAMPLE}_fastp.html" \
-w "$threads"

# Align reads with BWA MEM
bwa mem \
-R "@RG\tID:A\tDS:GAYA_UMI\tPL:ILLUMINA\tLB:$SAMPLE\tSM:$SAMPLE" \
-t "$threads" -M \
"$REF_GENOME" \
"${SAMPLE}_umi_extracted_trimmed_R1.fastq" \
"${SAMPLE}_umi_extracted_trimmed_R2.fastq" | \
samtools view -@ "$threads" -Sb > "${SAMPLE}_umi_aligned.bam"

# Merge BAM alignment
gatk-4.3.0.0/gatk MergeBamAlignment \
--ATTRIBUTES_TO_RETAIN X0 \
--ATTRIBUTES_TO_REMOVE NM \
--ATTRIBUTES_TO_REMOVE MD \
--ALIGNED_BAM "${SAMPLE}_umi_aligned.bam" \
--UNMAPPED_BAM "${SAMPLE}_unmapped_umi_extracted.bam" \
--OUTPUT "${SAMPLE}_umi_extracted_aligned_merged.bam" \
--REFERENCE_SEQUENCE "$REF_GENOME" \
--SORT_ORDER queryname \
--ALIGNED_READS_ONLY true \
--PRIMARY_ALIGNMENT_STRATEGY MostDistant \
--ALIGNER_PROPER_PAIR_FLAGS true \
--CLIP_OVERLAPPING_READS false

# Group reads by UMI
fgbio GroupReadsByUmi \
--input="${SAMPLE}_umi_extracted_aligned_merged.bam" \
--output="${SAMPLE}_umi_grouped.bam" \
--strategy=paired \
--edits=1 \
-t RX \
-f "${SAMPLE}_umi_group_data.tsv"

# Call Duplex Consensus Reads
fgbio CallDuplexConsensusReads \
--input="${SAMPLE}_umi_grouped.bam" \
--output="${SAMPLE}_umi_consensus_unmapped.bam" \
--error-rate-post-umi 40 \
--error-rate-pre-umi 45 \
--min-reads 3 1 1 \
--max-reads 50 \
--min-input-base-quality 30 \
--threads "$threads" \
--read-name-prefix='consensus'

# Convert consensus BAM to FASTQ
gatk-4.3.0.0/gatk SamToFastq \
-I "${SAMPLE}_umi_consensus_unmapped.bam" \
-F "${SAMPLE}_umi_consensus_unmapped_R1.fastq" \
-F2 "${SAMPLE}_umi_consensus_unmapped_R2.fastq" \
--CLIPPING_ATTRIBUTE XT \
--CLIPPING_ACTION 2

# Align consensus reads
bwa mem \
-R "@RG\tID:A\tDS:GAYA_UMI\tPL:ILLUMINA\tLB:$SAMPLE\tSM:$SAMPLE" \
-v 3 -Y -M \
-t "$threads" \
"$REF_GENOME" \
"${SAMPLE}_umi_consensus_unmapped_R1.fastq" \
"${SAMPLE}_umi_consensus_unmapped_R2.fastq" | \
samtools view -@ "$threads" -bh - > "${SAMPLE}_umi_consensus_mapped.bam"

# Sort BAMs
gatk-4.3.0.0/gatk SortSam \
-R "$REF_GENOME" \
-I "${SAMPLE}_umi_consensus_mapped.bam" \
-O "${SAMPLE}_umi_consensus_mappedB.bam" \
-SO queryname

gatk-4.3.0.0/gatk SortSam \
-R "$REF_GENOME" \
-I "${SAMPLE}_umi_consensus_unmapped.bam" \
-O "${SAMPLE}_umi_consensus_unmappedB.bam" \
-SO queryname

# Merge BAM files
picard MergeBamAlignment \
UNMAPPED="${SAMPLE}_umi_consensus_unmappedB.bam" \
ALIGNED="${SAMPLE}_umi_consensus_mappedB.bam" \
O="${SAMPLE}_final.bam" \
R="$REF_GENOME" \
CLIP_ADAPTERS=false \
VALIDATION_STRINGENCY=SILENT \
CREATE_INDEX=true \
EXPECTED_ORIENTATIONS=FR \
MAX_GAPS=-1 \
SO=coordinate \
ALIGNER_PROPER_PAIR_FLAGS=false

# Add Read Groups
picard AddOrReplaceReadGroups \
I="${SAMPLE}_final.bam" \
O="${SAMPLE}.bam" \
RGID=1 \
RGLB=lib1 \
RGPL=ILLUMINA \
RGPU=unit1 \
RGSM="$SAMPLE"

# Run recalibration script
./recal.sh "$SAMPLE"
