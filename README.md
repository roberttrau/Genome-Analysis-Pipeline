# ChIP-seq Analysis Pipeline
Applied Genomics master’s coursework at UConn. This automates sequencing read quality control, quality and adapter trimming, alignment to the hg38 reference genome, peak calling, and DNA motif analysis.

The script integrates FastQC, Cutadapt, Bowtie2, SAMtools, BEDTools, MACS2, and MEME Suite tools. It processes original and trimmed reads, records processing logs, and generates coverage and peak tracks for genome browser visualization.

This script uses file paths and reference resources specific to the original course environment.
