# pub-PABPN1_APA-shift_Raz

This is the code repository of scripts and code used in the data analysis of article: "PABPN1 loss-of-function causes APA-shift in oculopharyngeal muscular dystrophy", by Milad Shademan, Hailiang Mei, Yavuz Ariyurek, Susan Kloet, Vered Raz

The directories are organized based on the actual steps taken during the data analysis.

## step1_proximal_distal_annotations
This the script used to create GTF file annotating proximal and distal regions of all 3' UTR. Because motif based annotation file (see create_motif_based_gtf.py) does not give reliable annotations, we opt for a simpler and a bit naive approach by dividing 3' UTR into just one proximal and one distal region. Since this division is a bit abtrary, we allow a parameter to specify a shift of the division point between these two regions to investigate the robustness of this approach, SHIFT_RATIO_3UTR

## step2_run_biowdl_RNAseq_pipeline
All mouse and human reads in FASTQ format are processed using [BioWDL RNAseq pipeline version 5.0.0](https://github.com/biowdl/RNA-seq). It first filtered using Cutadapt (v2.10) to remove all remaining adapter sequences. Then the reads were aligned to the Ensembl transcriptome version 104 using STAR (v2.7.5a) including UMI based deduplication using UMI-Tools (v1.1.1) to generate transcriptome-based alignment in BAM format.

## step3_quantification_APA
With this APA proximal and distal annotation GTF file generated in step 1 and the human and mouse transcriptome-based alignment file generated in step2 as input, we quantified the APA enrichment signal at both proximal and distal regions using featureCounts (v2.0.1) with option “-M” to include multiple mapped reads. 

## step4_differential_analysis_APA

## step5_visualization
