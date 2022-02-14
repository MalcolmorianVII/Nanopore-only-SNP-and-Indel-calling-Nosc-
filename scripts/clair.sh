#! /usr/bin/bash

set -e 
eval "$(conda shell.bash hook)"
conda activate "snakemake"

#model=/home/belson/Clair3/models/ont_guppy5
ref_path=/home/ubuntu/data/belson/bioinformatics/projects_2021/nosc/reference_genomes
claired=/home/ubuntu/data/belson/bioinformatics/projects_2021/nosc/results/nosc_clair/2022.01.02
THREADS="8" 
MODEL_NAME="ont_guppy5"

sort_index() {
	# This function takes 2 arguments:
	# $1=SAM  i.e the sam file to sort and convert to bam file
	# $2=BAM  i.e the name of the bam file to be produced
	# And it produces a bam file which is sorted and indexed
	samtools sort -@ 8 -o $2 $1 && samtools index $2
}

species=(barcode01_Acinetobacter_baumannii_J9 barcode02_Citrobacter_koseri_MINF_9D barcode03_Enterobacter_kobei_MSB1_1B barcode04_Haemophilus_unknown_M1C132_1 barcode05_Klebsiella_oxytoca_MSB1_2C barcode06_CHF10J1 barcode07_Klebsiella_variicola_INF345 barcode08_Serratia_marcescens_17-147-1671)
minimap() {
	# This function does mapping using minimap2
	g5=/home/ubuntu/data/belson/bioinformatics/projects_2021/napa/input_data/guppy_5
	for specie in ${species[@]}
	do
		#Map the nanopore reads to the ref genome
		out=${claired}/${specie}_clair_results
		sam_file=${out}/${specie}.sam
		bam_file=${out}/${specie}.bam
		if [ -f ${bam_file} ]
		then
			echo "Sorted and indexed BAM file already exists!!!!"
			echo "Therefore skipping minimap step"
			continue
		else
			if [ -f ${sam_file} && ! -f ${bam_file} ]
			then
				sort_index ${sam_file} ${bam_file}
			else
				mkdir -p ${out}
				minimap2 -ax map-ont ${ref_path}/${specie}_ref.fna ${g5}/${specie}.fastq.gz > ${out}/${sam_file}
				sort_index ${sam_file} ${bam_file}
			fi
		fi
	done
}


faidx() {
	# This function takes two argument:
	# $1=ref_path
	# $2=ref_genome
	samtools faidx $1/$2
}
	 	
run_clair () {
	conda activate "clair3"

	base=$(basename ${1} | sed 's/_clair_results//')
	BAM=${base}.bam
	ref_genome=${base}_ref.fna

	clair_vcf=${claired}/${base}_clair_results/merge_output.vcf.gz
	# norm_vcf=${claired}/${base}_clair_results/${base}_clair.norm.vcf.gz

	# Check if the reference genome is indexed
	if [ ! -f ${ref_path}/${ref_genome}.fai ]
	then
		faidx ${ref_path} ${ref_genome}
	fi

	if [ -f ${clair_vcf} ]
	then
		echo "Skipp running clair"
		return
	else
		# First check if the vcf from clair is already found
		mkdir -p ${claired}/${base}_clair_results
		# Run clair3 first
		/home/belson/Clair3/run_clair3.sh \
  		--bam_fn=${dir}/${BAM}\
  		--ref_fn=${ref_path}/${ref_genome} \
  		--threads=${THREADS} \
  		--platform="ont" \
  		--model_path="${CONDA_PREFIX}/bin/models/${MODEL_NAME}" \
  		--output=${claired}/${base}_clair_results \
  		--include_all_ctgs \
  		--haploid_precise 

		
	fi
	
}

unique=/home/ubuntu/data/belson/bioinformatics/projects_2021/nosc/results/uniquely_mapped/2022.02.07
filter_vcfs () {
	# This function takes a clair vcf & filters those regions that <5x
	# Has 3 arguments
	# $1 clair_vcf
	# $2 bed file
	# $3 filtered vcf
	conda activate "snakemake"
	bedtools intersect -a ${1} -b ${2} -header > ${3}
}

skip_indels () {
	# This function aims to return a vcf file with indels skipped
	# Has 2 arguments
	# $1 vcf 
	# $2 output
	echo "Skipping indels"
	echo ${2},"using",${1}
	/home/belson/bin/vcftools/src/cpp/vcftools --gzvcf ${1} --remove-indels  --recode --recode-INFO-all --stdout | gzip -c > ${2}

}

normalize () {
	# This function performs vcf normalization
	# Takes in 4 arguments
	# $1 = ref_path
	# $2 = ref_genome
	# $3 = indel_skipped_vcf
	# $4 = normalized vcf
	conda deactivate
	conda activate "bcf"
	echo "Normalization started"
	# echo ${1},${2},${3},${4}
	bcftools norm  -f ${1}/${2} ${3}  -Oz -o ${4}
	tabix -p vcf ${4}
	echo "Normalization ended"

}

echo "----- Started mapping reads to reference genome -----"
minimap
echo "----- Started Clair3 -----"
for dir in ${claired}/barcode*
do
	run_clair ${dir}
	# Filter clair-vcf for <5x regions
	bed_file=${unique}/${base}_uniquely_mapped/${base}_more5x.bed
	filtered_vcf=${claired}/${base}_clair_results/${base}_clair_filtered.vcf.gz
	echo "echo vars!!!!!!!!!!!!!!"
	echo ${clair_vcf},${bed_file},${filtered_vcf}
	echo "vars echoed!!!!!!!!!!!!!!"
	echo 

	echo "Filtering vcfs"
	filter_vcfs ${clair_vcf} ${bed_file} ${filtered_vcf}
	echo "Filtering vcfs done"
	# Skip indels
	indel_skiped_vcf=${claired}/${base}_clair_results/${base}_clair_indel_skiped.vcf.gz
	 
	skip_indels ${filtered_vcf}  ${indel_skiped_vcf}

	# Then proceed with normalization
	#normVCF=${claired}/${base}_clair_results/${base}_clair.norm.vcf.gz
	
	norm_vcf=${claired}/${base}_clair_results/${base}_clair.norm.vcf.gz
	normalize  ${ref_path} ${ref_genome} ${indel_skiped_vcf} ${norm_vcf}
	

done

echo " Clair pipeline run successfully "
