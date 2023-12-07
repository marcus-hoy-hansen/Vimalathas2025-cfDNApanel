#!/bin/bash

threads=8
SAMPLE=$1

gatk-4.3.0.0/gatk FastqToSam \
-F1 $SAMPLE'_R1_001.fastq.gz' \
-F2 $SAMPLE'_R2_001.fastq.gz' \
-O $SAMPLE'_unmapped.bam' \
-SM $SAMPLE

fgbio ExtractUmisFromBam \
-i $SAMPLE'_unmapped.bam' \
-o $SAMPLE'_unmapped_umi_extracted.bam' \
-r 5M2S+T 5M2S+T \
-t RX \
-a true

gatk-4.3.0.0/gatk SamToFastq \
-I $SAMPLE'_unmapped_umi_extracted.bam' \
-F $SAMPLE'_umi_extracted_R1.fastq' \
-F2 $SAMPLE'_umi_extracted_R2.fastq' \
--CLIPPING_ATTRIBUTE XT \
--CLIPPING_ACTION 2

fastp \
-i $SAMPLE'_umi_extracted_R1.fastq' \
-o $SAMPLE'_umi_extracted_trimmed_R1.fastq' \
-I $SAMPLE'_umi_extracted_R2.fastq' \
-O $SAMPLE'_umi_extracted_trimmed_R2.fastq' \
-g -W 5 -q 30 -u 40 -x -3 -l 75 -c -e 30 \
-j $SAMPLE'.fastp.json' \
-h $SAMPLE'fastp.html' \
-w $threads


bwa mem \
-R "@RG\tID:A\tDS:GAYA_UMI\tPL:ILLUMINA\tLB:$SAMPLE\tSM:$SAMPLE" \
-t $threads -M \
hg38/resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta \
$SAMPLE'_umi_extracted_trimmed_R1.fastq' \
$SAMPLE'_umi_extracted_trimmed_R2.fastq' | \
samtools view -Sb > $SAMPLE'_umi_aligned.bam'

gatk-4.3.0.0/gatk MergeBamAlignment \
--ATTRIBUTES_TO_RETAIN X0 \
--ATTRIBUTES_TO_REMOVE NM \
--ATTRIBUTES_TO_REMOVE MD \
--ALIGNED_BAM $SAMPLE'_umi_aligned.bam' \
--UNMAPPED_BAM $SAMPLE'_unmapped_umi_extracted.bam' \
--OUTPUT $SAMPLE'_umi_extracted_aligned_merged.bam' \
--REFERENCE_SEQUENCE hg38/resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta \
--SORT_ORDER queryname \
--ALIGNED_READS_ONLY true \
--MAX_INSERTIONS_OR_DELETIONS -1 \
--PRIMARY_ALIGNMENT_STRATEGY MostDistant \
--ALIGNER_PROPER_PAIR_FLAGS true \
--CLIP_OVERLAPPING_READS false

#samtools view -@ $threads -f 2 -bh $SAMPLE'_umi_extracted_aligned_merged.bam' > $SAMPLE'_umi_extracted_aligned_merged_filtered.bam'

fgbio GroupReadsByUmi \
--input=$SAMPLE'_umi_extracted_aligned_merged.bam' \
--output=$SAMPLE'_umi_grouped.bam' \
--strategy=paired \
--edits=1 \
-t RX \
-f $SAMPLE'_umi_group_data.tsv'

fgbio CallDuplexConsensusReads \
--input=$SAMPLE'_umi_grouped.bam' \
--output=$SAMPLE'_umi_consensus_unmapped.bam' \
--error-rate-post-umi 40 \
--error-rate-pre-umi 45 \
--min-reads 3 1 1 \
--max-reads 50 \
--min-input-base-quality 30 \
--threads $threads \
--read-name-prefix='consensus'


gatk-4.3.0.0/gatk SamToFastq \
-I $SAMPLE'_umi_consensus_unmapped.bam' \
-F $SAMPLE'_umi_consensus_unmapped_R1.fastq' \
-F2 $SAMPLE'_umi_consensus_unmapped_R2.fastq' \
--CLIPPING_ATTRIBUTE XT \
--CLIPPING_ACTION 2

bwa mem \
-R "@RG\tID:A\tDS:GAYA_UMI\tPL:ILLUMINA\tLB:$SAMPLE\tSM:$SAMPLE" \
-v 3 -Y -M \
-t $threads \
hg38/resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta \
$SAMPLE'_umi_consensus_unmapped_R1.fastq' \
$SAMPLE'_umi_consensus_unmapped_R2.fastq' \
| \
samtools view -@ $threads -bh - > $SAMPLE'_umi_consensus_mapped.bam'

gatk-4.3.0.0/gatk SortSam \
-R hg38/resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta \
-I $SAMPLE'_umi_consensus_mapped.bam' \
-O $SAMPLE'_umi_consensus_mappedB.bam' \
-SO queryname

gatk-4.3.0.0/gatk SortSam \
-R hg38/resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta \
-I $SAMPLE'_umi_consensus_unmapped.bam' \
-O $SAMPLE'_umi_consensus_unmappedB.bam' \
-SO queryname

picard MergeBamAlignment UNMAPPED=$SAMPLE'_umi_consensus_unmappedB.bam' ALIGNED=$SAMPLE'_umi_consensus_mappedB.bam' O=$SAMPLE'_final.bam' R=hg38/resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta CLIP_ADAPTERS=false VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true EXPECTED_ORIENTATIONS=FR MAX_GAPS=-1 SO=coordinate ALIGNER_PROPER_PAIR_FLAGS=false

picard AddOrReplaceReadGroups \
I=$SAMPLE'_final.bam' \
O=$SAMPLE'.bam' \
RGID=1 \
RGLB=lib1 \
RGPL=ILLUMINA \
RGPU=unit1 \
RGSM=20

./recal.sh $SAMPLE
