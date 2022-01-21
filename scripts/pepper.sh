#! /usr/bin/bash
species=(barcode01_Acinetobacter_baumannii_J9 barcode02_Citrobacter_koseri_MINF_9D barcode03_Enterobacter_kobei_MSB1_1B barcode04_Haemophilus_unknown_M1C132_1 barcode05_Klebsiella_oxytoca_MSB1_2C barcode06_CHF10J1 barcode07_Klebsiella_variicola_INF345 barcode08_Serratia_marcescens_17-147-1671)
ref=/home/ubuntu/data/belson/bioinformatics/projects_2021/nosc/reference_genomes
BASE="/home/ubuntu/data/belson/bioinformatics/projects_2021/nosc/results/nosc_pepper/2022.01.02"
mapped=/home/ubuntu/data/belson/bioinformatics/projects_2021/nosc/results/nosc_pepper/2022.01.02
# Set the number of CPUs to use
THREADS="8"

sort_index() {
	# $1=BASE,$2=SAM,$3=BAM
	samtools sort -@ 8 -o $1/$3 $1/$2 && samtools index $1/$3
}

minimap() {
	g5=/home/ubuntu/data/belson/bioinformatics/projects_2021/napa/input_data/guppy_5
	for i in ${species[@]}
	do	
	#Map the nanopore reads to the ref genome
	out=${BASE}/${i}_pepper_results
	if [ -f ${out}/${i}.bam ]
	then
		echo "No need for minimapping"
		continue
	fi

	if [ -f ${out}/${i}.sam ]
	then
		sort_index ${out} ${i}.sam ${i}.bam
	else
		mkdir -p ${out}
		minimap2 -ax map-ont ${ref}/${i}_ref.fna ${g5}/${i}.fastq.gz > ${out}/${i}.sam
		sort_index ${out} ${i}.sam ${i}.bam
	fi
	done
}
pepper(){
	# Run PEPPER-Margin-DeepVariant
	# $1=bam_path,$2=ref_path,$3=outpath,$3=BAM,$4=REF,$5=OUTPUT_PREFIX,$6=THREADS
	docker run --ipc=host \
	-v $1:$1\
	-v $2:$2\
	kishwars/pepper_deepvariant:r0.5 \
	run_pepper_margin_deepvariant call_variant \
	-b ${1}/${3} \
	-f ${2}/${4} \
	-o ${1} \
	-p "$5" \
	-t $6 \
	--ont
}

minimap 

for dir in ${mapped}/barcode*
do
	# mkdir -p ${dir}/pepper_results
	slashed=$(basename $dir | sed 's/_pepper_results//')
	BAM=${slashed}.bam
	REF=${slashed}_ref.fna
	if [ -d ${dir}/{slashed}.vcf.gz ]
	then
		echo "PREVIOUS PEPPER RUN WAS SUCCESSFUL"
	else
		pepper ${dir} ${ref} ${BAM} ${REF} ${slashed} ${THREADS}
	fi
done

