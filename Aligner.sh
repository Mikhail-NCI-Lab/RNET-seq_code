#!/bin/sh

# Remove adaptor from the original FastQ file.
cutadapt -a GATCGGAAGAGCACACGTCTGAACTCCAGTCAC -o original_data_R1_001_cutadapt.fastq original_data_R1_001.fastq

# Remove PCR duplicates from the processed FastQ file.
clumpify.sh in=original_data_R1_001_cutadapt.fastq out=original_data_R1_001_cutadapt_remD.fastq subs=0 dedupe=t

# Remove barcode from the processed FastQ file.
cutadapt --cut 6 original_data_R1_001_cutadapt_remD.fastq -o original_data_R1_001_cutadapt_remD_6Ntrim.fastq

# Align reads in the processed FastQ file to E. coli MG1655 genome NC_000913.2.
bowtie ./Bowtie_index/NC_000913.2 -q original_data_R1_001_cutadapt_remD_6Ntrim.fastq --best --strata -v 1 -m 1 -p 10 -S original_data_processed.sam

# Sort and transform the Sam file to Bam file.
samtools sort original_data_processed.sam -o original_data_processed_sort.sam
samtools view -F 4 original_data_processed_sort.sam -b -o original_data_processed_sort.bam

# Generate the Wig files containing the count of reads where 5' end locate on the positive and negative strands.
bedtools genomecov -d -5 -strand - -ibam original_data_processed_sort.bam > original_data_processed_sort_Positive.wig
bedtools genomecov -d -5 -strand + -ibam original_data_processed_sort.bam > original_data_processed_sort_Negative.wig
