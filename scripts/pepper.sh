#! /usr/bin/bash
eval "$(conda shell.bash hook)"
conda activate "snakemake"

species=(barcode01_Acinetobacter_baumannii_J9 barcode02_Citrobacter_koseri_MINF_9D barcode03_Enterobacter_kobei_MSB1_1B barcode04_Haemophilus_unknown_M1C132_1 barcode05_Klebsiella_oxytoca_MSB1_2C barcode06_CHF10J1 barcode07_Klebsiella_variicola_INF345 barcode08_Serratia_marcescens_17-147-1671)
ref_path=/home/ubuntu/data/belson/bioinformatics/projects_2021/nosc/reference_genomes
mapped=/home/ubuntu/data/belson/bioinformatics/projects_2021/nosc/results/nosc_pepper/2022.02.07
# Set the number of CPUs to use
THREADS="8"

sort_index() {
	# This function takes 2 arguments:
	# $1=SAM  i.e the sam file to sort and convert to bam file
	# $2=BAM  i.e the name of the bam file to be produced
	# And it produces a bam file which is sorted and indexed
	samtools sort -@ 8 -o "$2" "$1" && samtools index "$2"
}

species=(barcode01_Acinetobacter_baumannii_J9 barcode02_Citrobacter_koseri_MINF_9D barcode03_Enterobacter_kobei_MSB1_1B barcode04_Haemophilus_unknown_M1C132_1 barcode05_Klebsiella_oxytoca_MSB1_2C barcode06_CHF10J1 barcode07_Klebsiella_variicola_INF345 barcode08_Serratia_marcescens_17-147-1671)
minimap() {
	# This function does mapping using minimap2
	g5=/home/ubuntu/data/belson/bioinformatics/projects_2021/napa/input_data/guppy_5
	for specie in ${species[@]}
	do
		#Map the nanopore reads to the ref genome
		out=${mapped}/${specie}_pepper_results
		sam_file=${out}/${specie}.sam
		bam_file=${out}/${specie}.bam
		if [ -f "${bam_file}" ]
		then
			echo "Sorted and indexed BAM file already exists!!!!"
			echo "Therefore skipping minimap step"
			continue
		else
			if [ -f "${sam_file}" ]  && [ ! -f "${bam_file}" ]
			then
				sort_index "${sam_file}" "${bam_file}"
			else
				mkdir -p "${out}"
				minimap2 -ax map-ont ${ref_path}/"${specie}"_ref.fna ${g5}/"${specie}".fastq.gz > "${out}"/"${sam_file}"
				sort_index "${sam_file}" "${bam_file}"
			fi
		fi
	done
}


pepper(){
	# Run PEPPER-Margin-DeepVariant
	# $1=bam_path,$2=ref_path,$3=outpath,$3=BAM,$4=REF,$5=OUTPUT_PREFIX,$6=THREADS
	# Check if the reference genome is indexed
	pepper_vcf=${dir}/${base}.vcf.gz
	if [ ! -f ${ref_path}/"${ref_genome}".fai ]
	then
		faidx ${ref_path} "${ref_genome}"
	fi
	if [[ -f ${pepper_vcf} && $(stat -c %s "${pepper_vcf}") -ne 0 ]]
	then
		echo "PREVIOUS PEPPER RUN SUCCESSFUL"
		echo "EXITING....."
		return
	else
		docker run --ipc=host \
		-v "$1":"$1"\
		-v "$2":"$2"\
		kishwars/pepper_deepvariant:r0.7 \
		run_pepper_margin_deepvariant call_variant \
		-b "${1}"/"${3}" \
		-f "${2}"/"${4}" \
		-o "${1}" \
		-p "$5" \
		-t "$6" \
		--ont_r9_guppy5_sup
	fi
}

unique=/home/ubuntu/data/belson/bioinformatics/projects_2021/nosc/results/uniquely_mapped/2022.02.07
filter_vcfs () {
	# This function takes a pepper vcf & filters those regions that <5x
	# Has 3 arguments
	# $1 pepper_vcf
	# $2 bed file
	# $3 filtered vcf
	conda activate "snakemake"
	bedtools intersect -a "${1}" -b "${2}" -header > "${3}"
}

skip_indels () {
	# This function aims to return a vcf file with indels skipped
	# Has 2 arguments
	# $1 vcf 
	# $2 output
	echo "Skipping indels"
	echo "${2}","using","${1}"
	/home/belson/bin/vcftools/src/cpp/vcftools --gzvcf "${1}" --remove-indels  --recode --recode-INFO-all --stdout | gzip -c > "${2}"

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
	echo "----- Normalization started -----"
	echo "----- Checking if normalized file already exists -----"
	if [[ -f ${dir}/${base}_pepper.norm.vcf.gz && $(stat -c %s "${dir}"/"${base}"_pepper.norm.vcf.gz) -ne 0 ]]
	then
		echo "----- No need for normalization -----"
		return 
	else
		echo "----- Normalied vcf not found ----- "
		echo "----- Proceed with normalization -----"
		# echo ${1},${2},${3},${4}
		bcftools norm  -f "${1}"/"${2}" "${3}"  -Oz -o "${4}"
		tabix -p vcf "${4}"
		echo "----- Normalization ended -----"
	fi

}


minimap 

for dir in ${mapped}/barcode*
do
	# mkdir -p ${dir}/pepper_results
	base=$(basename "$dir" | sed 's/_pepper_results//')
	BAM=${base}.bam
	ref_genome=${base}_ref.fna

	pepper "${dir}" ${ref_path} "${BAM}" "${ref_genome}" "${base}" ${THREADS}
	#Filter vcfs for <5x regions

	bed_file=${unique}/${base}_uniquely_mapped/${base}_more5x.bed
	filtered_vcf=${mapped}/${base}_pepper_results/${base}_pepper_filtered.vcf.gz
	
	echo "echo vars!!!!!!!!!!!!!!"
	echo "${pepper_vcf}","${bed_file}","${filtered_vcf}"
	echo "vars echoed!!!!!!!!!!!!!!"
	echo 

	echo "Filtering vcfs"
	filter_vcfs "${pepper_vcf}" "${bed_file}" "${filtered_vcf}"
	echo "Filtering vcfs done"

	# Skip indels
	indel_skiped_vcf=${mapped}/${base}_pepper_results/${base}_pepper_indel_skiped.vcf.gz
	skip_indels "${filtered_vcf}"  "${indel_skiped_vcf}"

	# Then proceed with normalization
	norm_vcf=${mapped}/${base}_pepper_results/${base}_pepper.norm.vcf.gz
	normalize  ${ref_path} "${ref_genome}" "${indel_skiped_vcf}" "${norm_vcf}"
done


