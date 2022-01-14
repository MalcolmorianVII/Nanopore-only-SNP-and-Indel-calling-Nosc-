#! /usr/bin/bash

ref=/home/ubuntu/data/belson/guppy5_guppy3_comparison/reference_genomes/Acinetobacter_baumannii.fna
res=/home/ubuntu/data/belson/guppy5_guppy3_comparison/nosc/nosc_gatk/results/2021.11.25
reads=/home/ubuntu/data/belson/illumina_reads/Acinetobacter_baumannii_J9

bwa mem ${ref}  ${reads}/J9_S142_L001_R1_001.fastq.gz ${reads}/J9_S142_L001_R2_001.fastq.gz > ${res}/Acinetobacter_bwa-mem_original-mapped.sam

samtools sort -@ 8 -o ${res}/Acinetobacter_bwa-mem_original-mapped.bam ${res}/Acinetobacter_bwa-mem_original-mapped.sam


samtools view -h ${res}/Acinetobacter_bwa-mem_original-mapped.bam | grep -v -e 'XA:Z:' -e 'SA:Z:' | samtools view -b > ${res}/Acinetobacter_bwa-mem_original-mapped_uniquely_mapped.bam

samtools depth -aa ${res}/Acinetobacter_bwa-mem_original-mapped_uniquely_mapped.bam | awk '{if ($3 > 5) print $0}' > ${res}/more5XcoverageAcinetobacter_mapped_original.txt
