#! /usr/bin/bash

set -e 
eval "$(conda shell.bash hook)"
conda activate "snakemake"

species=(barcode01_Acinetobacter_baumannii_J9 barcode02_Citrobacter_koseri_MINF_9D barcode03_Enterobacter_kobei_MSB1_1B barcode04_Haemophilus_unknown_M1C132_1 barcode05_Klebsiella_oxytoca_MSB1_2C barcode06_CHF10J1 barcode07_Klebsiella_variicola_INF345 barcode08_Serratia_marcescens_17-147-1671)
BASE=/home/ubuntu/data/belson/bioinformatics/projects_2021/nosc/results/nosc_clair/2022.01.02
model=/home/belson/Clair3/models/ont_guppy5
ref=/home/ubuntu/data/belson/bioinformatics/projects_2021/nosc/reference_genomes
claired=/home/ubuntu/data/belson/bioinformatics/projects_2021/nosc/results/nosc_clair/2022.01.02
THREADS="8" 
MODEL_NAME="ont_guppy5"
sort_index() {
	# $1=BASE,$2=SAM,$3=BAM
	samtools sort -@ 8 -o $1/$3 $1/$2 && samtools index $1/$3
}

minimap() {
	g5=/home/ubuntu/data/belson/bioinformatics/projects_2021/napa/input_data/guppy_5
	for i in ${species[@]}
	do	
	#Map the nanopore reads to the ref genome
	out=${claired}/${i}_clair_results
	if [ -f ${out}/${i}.bam ]
	then
		echo "Minimap not run!!!!"
		continue
	else
		if [ -f ${out}/${i}.sam ]
		then
			sort_index ${out} ${i}.sam ${i}.bam
		else
			mkdir -p ${out}
			minimap2 -ax map-ont ${ref}/${i}_ref.fna ${g5}/${i}.fastq.gz > ${out}/${i}.sam
			sort_index ${out} ${i}.sam ${i}.bam
		fi
	fi
	done
}


faidx() {
	# $1=ref_path $2=ref_genome
	samtools faidx $1/$2	#Somehow giving an error of not finding fai file for consensus.fasta
}

echo "<<<< STARTING WITH ALL ABOUT MAPPINGS I.E. MINIMAPING,SORTING & INDEXING >>>>"
minimap
echo "<<<< STARTED Clair3 >>>>"
conda deactivate
for dir in ${claired}/barcode*
do
	# mkdir -p ${dir}/pepper_results
	conda activate "clair3"
	slashed=$(basename $dir | sed 's/_clair_results//')
	BAM=${slashed}.bam
	ref_genome=${slashed}_ref.fna
	if [ ! -f ${ref}/${ref_genome}.fai ]
	then
		faidx ${ref} ${ref_genome}
	fi

	if [ -f ${claired}/${slashed}_clair_results/${slashed}_clair.norm.vcf.gz ]
	then
		continue
	else
		mkdir -p ${claired}/${slashed}_clair_results
		/home/belson/Clair3/run_clair3.sh \
  --bam_fn=${dir}/${BAM}\
  --ref_fn=${ref}/${ref_genome} \
  --threads=${THREADS} \
  --platform="ont" \
  --model_path="${CONDA_PREFIX}/bin/models/${MODEL_NAME}" \
  --output=${claired}/${slashed}_clair_results \
  --include_all_ctgs \
  --haploid_precise  
		# bash ~/Clair3/run_clair3.sh  --model_path=${model} --ref_fn=${ref}/${ref_genome} \
	 # --bam_fn=${dir}/${BAM} --platform="ont" --haploid_precise \
	 # --output=${claired}/${slashed}_clair_results --threads=8 --include_all_ctgs
	 # proceeding with normalization
	 normVCF=${claired}/${slashed}_clair_results/${slashed}_clair.norm.vcf.gz
	 conda activate "bcf"
	 bcftools norm  -f ${ref}/${ref_genome} ${claired}/${slashed}_clair_results/merge_output.vcf.gz  -Oz -o $normVCF
	 tabix -p vcf $normVCF
	fi

done

echo "<<<<DONE WITH CLAIR>>>>"
