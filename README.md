# Supporting information on targeted cfDNA gene panel for investigating circulating tumor DNA (ctDNA) in DLBCL, Vimalathas et al. 2025
The provided [panel](Probes_merged_ok_OUH_Bcell_clonality_v1_1X_TE-94124764_hg38.bed) is implemented in a UMI-based targeted sequencing pipeline based on cfDNA from DLBCL. Workflow consists of [UMI processing](umi_analyze.sh) and variant calling followed by variant [recalling](recall.sh), and [annotation](annotates.sh).

## Manuscript title
Liquid Biopsy for Enhanced Specificity in Identifying Somatic Mutations in Aggressive Non-Hodgkin Large B-cell Lymphoma: A Comparative Study of Cell-free DNA and Formalin-Fixed Paraffin-Embedded Tissue 

## Abstract
Introduction: Formalin-fixed paraffin-embedded (FFPE) tumor biopsy is the current mainstay of genotyping, but is limited by its invasiveness and tumor heterogeneity. Plasma cell-free DNA (cfDNA) constitutes a minimally invasive alternative that may better capture tumor-derived profiles from circulating tumor DNA (ctDNA). This study compares the performance and genomic concordance of cfDNA and FFPE tumor DNA in aggressive non-Hodgkin large B-cell lymphoma. 

Methods: Paired diagnostic FFPE tissue and plasma samples from 15 patients were sequenced with a custom 53-gene panel. 

Results: Detection thresholds were empirically guided at 1% variant allele frequency (VAF) for cfDNA and 10% for unpaired FFPE DNA. The median number of cfDNA variants was 6 (interquartile range (IQR): 2–11) versus 63 (IQR: 15–250) in FFPE DNA at 1% VAF. Collectively, 102 somatic variants were shared between cfDNA and FFPE DNA with a median of 5 (range: 0–23). cfDNA showed a five-fold lower median VAF for shared variants than FFPE DNA (7% vs. 36%, p < 0.0001).  Eighty percent of patients harbored at least one cfDNA variant. A maximum cfDNA recall rate of 83% was observed at FFPE DNA VAF >50%. COSMIC database overlap was twice as high for cfDNA compared to FFPE DNA (22% vs. 11%) at 10% VAF.

Conclusion: cfDNA has superior specificity for somatic mutation detection but lower sensitivity than FFPE DNA. Modest concordance was demonstrated between the two compartments. Our results support a complementary role of ctDNA in mutational profiling at a 1% VAF threshold in a clinically applicable set-up.



## Non-immunoglobulin genes included (50 genes)
*ARID1A*, *ACTB*, *B2M*, *BCL2*, *BCL6*, *BCL10*, *BTG1*, *CARD11*, *CIITA*, *CCND3*, *CD70*, *CD79A*, *CD79B*, *CDKN2A*, *CREBBP*, *EP300*, *ETV6*, *EZH2*, *FOXO1*, *GNA13*, *H1-2*, *H1-4*, *IGLL5*, *IRF4*, *ID3*, *ITPKB*, *KMT2D*, *MEF2B*, *MYC*, *MYD88*, *NFKBIE*, *NOTCH1*, *NOTCH2*, *PIM1*, *PRDM1*, *P2RY8*, *SOCS1*, *STAT3*, *STAT6*, *SPEN*, *SGK1*, *TNFAIP3*, *TNFRSF14*, *TP53*, *TET2*, *CD58*, *UBEA2*, *PIM2*, *DTX1*, *IRF8*

## Immunoglobulin genes included (3 loci: IGH, IGL, IGK, 137 individual Ig genes)
*IGHV6-1*, *IGHV1-2*, *IGHV1-3*, *IGHV4-4*, *IGHV7-4-1*, *IGHV2-5*, *IGHV3-7*, *IGHV1-8*, *IGHV3-9*, *IGHV3-11*, *IGHV3-13*, *IGHV3-15*, *IGHV1-18*, *IGHV3-20*, *IGHV3-21*, *IGHV3-23*, *IGHV1-24*, *IGHV2-26*, *IGHV4-28*, *IGHV3-30*, *IGHV4-31*, *IGHV3-33*, *IGHV4-34*, *IGHV4-39*, *IGHV3-43*, *IGHV1-45*, *IGHV1-46*, *IGHV3-48*, *IGHV3-49*, *IGHV5-51*, *IGHV3-53*, *IGHV1-58*, *IGHV4-59*, *IGHV4-61*, *IGHV3-64*, *IGHV3-66*, *IGHV1-69*, *IGHV2-70*, *IGHV3-72*, *IGHV3-73*, *IGHV3-74*, *IGHD1-26*, *IGHD6-25*, *IGHD3-22*, *IGHD2-21*, *IGHD1-20*, *IGHD6-19*, *IGHD5-18*, *IGHD4-17*, *IGHD3-16*, *IGHD2-15*, *IGHD6-13*, *IGHD5-12*, *IGHD3-10*, *IGHD3-9*, *IGHD2-8*, *IGHD1-7*, *IGHD6-6*, *IGHD5-5*, *IGHD4-4*, *IGHD3-3*, *IGHD2-2*, *IGHD1-1*, *IGHD7-27*, *IGHJ1*, *IGHJ2*, *IGHJ3*, *IGHJ4*, *IGHJ5*, *IGHJ6*, *IGKV4-1*, *IGKV5-2*, *IGKV1-5*, *IGKV1-6*, *IGKV1-8*, *IGKV1-9*, *IGKV3-11*, *IGKV1-12*, *IGKV1-13*, *IGKV3-15*, *IGKV1-16*, *IGKV1-17*, *IGKV3-20*, *IGKV6-21*, *IGKV2-24*, *IGKV1-27*, *IGKV2-28*, *IGKV2-29*, *IGKV2-30*, *IGKV1-33*, *IGKV1-39*, *IGKV2-40*, *IGKJ1*, *IGKJ2*, *IGKJ3*, *IGKJ4*, *IGKJ5*, *IGLV4-69*, *IGLV8-61*, *IGLV4-60*, *IGLV6-57*, *IGLV10-54*, *IGLVpreB*, *IGLV5-52*, *IGLV1-51*, *IGLV9-49*, *IGLV1-47*, *IGLV7-46*, *IGLV5-45*, *IGLV1-44*, *IGLV7-43*, *IGLV1-40*, *IGLV5-37*, *IGLV1-36*, *IGLV3-27*, *IGLV3-25*, *IGLV2-23*, *IGLV3-22*, *IGLV3-21*, *IGLV3-19*, *IGLV2-18*, *IGLV3-16*, *IGLV2-14*, *IGLV3-12*, *IGLV2-11*, *IGLV3-10*, *IGLV3-9*, *IGLV2-8*, *IGLV4-3*, *IGLV3-1*, *IGLV5-39*, *IGLV5-37*, *IGLJ1*, *IGLJ2*, *IGLJ3*, *IGLJ6*, *IGLJ7*


| **Design Information**              | **Value**                   |
|-------------------------------------|-----------------------------|
| **Customer**                        | OHU                         |
| **Design Name**                     | Bcell_clonality_v1_1X       |
| **Design ID**                       | TE-94124764                 |
| **Genome**                          | hg38                        |
| **Probe Length**                    | 120                         |
| **Mean GC%**                        | 55.93%                      |
| **Stringency**                      | Medium                      |
| **Target Size (bp)**                 | 169,605                     |
| **Total Target Regions**            | 816                         |
| **Target Regions with Probes**      | 792                         |
| **Number of Probes**                | 1,710                       |
| **Probes Removed (Repeats)**        | 61                          |
| **Probes Removed (Ambiguous Bases)** | 0                           |
| **Non-Covered Whole Target Regions (#)** | 24                     |
| **Non-Covered Whole Target Regions (%)** | 2.94%                  |
| **Non-Covered Whole Target Regions (bp)** | 2,665                 |
| **Covered Size (bp)**               | 164,551                     |
| **Overall Coverage (%)**            | 97.02%                      |
| **Design Size (bp)**                | 183,255                     |
