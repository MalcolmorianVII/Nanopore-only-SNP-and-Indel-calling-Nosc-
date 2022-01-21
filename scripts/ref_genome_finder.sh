#! /usr/bin/bash
eval "$(conda shell.bash hook)"
conda activate "snakemake"

db_path='/home/ubuntu/data/belson/bacteria-refseq/'
samples='/home/ubuntu/data/belson/bioinformatics/projects_2021/napa/results/2021.08.02/guppy3'
ref_genomes=/home/ubuntu/data/belson/bioinformatics/projects_2021/nosc/results/Refseeker/${name}Refseeker.txt
for sample in $samples/*
do
	name=`basename $sample`
	referenceseeker $db_path $sample/*Illumina/contigs.fasta | tee $ref_genomes
done

echo "Reference sequence finding complete! Result is in $ref_genomes "


