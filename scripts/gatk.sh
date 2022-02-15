#! /usr/bin/bash
set -e
eval "$(conda shell.bash hook)"

reads=/home/ubuntu/data/belson/bioinformatics/projects_2021/napa/input_data/illumina_reads
gatk=/home/ubuntu/data/belson/bioinformatics/projects_2021/nosc/results/nosc_gatk/2022.01.05
ref_path=/home/ubuntu/data/belson/bioinformatics/projects_2021/nosc/reference_genomes
species=(barcode01_Acinetobacter_baumannii_J9
 barcode02_Citrobacter_koseri_MINF_9D 
 barcode03_Enterobacter_kobei_MSB1_1B 
 barcode04_Haemophilus_unknown_M1C132_1 
 barcode05_Klebsiella_oxytoca_MSB1_2C 
 barcode06_CHF10J1 
 barcode07_Klebsiella_variicola_INF345 
 barcode08_Serratia_marcescens_17-147-1671)


bwa_map () {
     #  Create index file
     # Arguments
     # $1 = ${ref_genome}
     # $2 = ${read1}
     # $3 = ${read2}
     # $4 = ${SAM}
     echo "Started  BWA indexing"
     if [ ! -f ${1} ] 
     then
          bwa index "${1}"
     fi
     if [ -f ${4} ]
     then
          echo "Sam file found...quiting mapping"
          return
     else
          echo "Starting bwa mapping"
          # Ensure SAM/BAM files have read groups
          bwa mem -R '@RG\tID:dip1\tSM:dip1' "${1}" "${2}" "${3}" > "${4}"
          echo "BWA mapping done"
     fi
}


sort_index () {
     # Arguments
     # $1 = ${SAM}
     # $2 = ${BAM}
     echo "Started sorting and indexing"
     samtools sort -@ 8 -o "${2}" "${1}" && samtools index "${2}"
}

CreateSequenceDictionary() {
     # Arguments
     # $1 = ${ref_genome}
     # $2 = ${ref_path}/${specie}_ref.dict
     # First check if sequence dictionary already exists
     if [ -f "$2" ]
     then
          echo "Sequence dictionary already exists"
          return
     else
          java -jar ~/picard.jar CreateSequenceDictionary -R "${1}"
          echo "Sequence dictionary created"
     fi

}

run_gatk () {
     # First check the existence of gatk_vcf file
     # Arguments :
     # $1 = ${BAM}
     # $2 = ${ref_genome}
     # $3 = ${gatk_dir}/${specie}_illumina_original.vcf.gz
     conda activate "gatk"
     echo "----- Started running gatk -----"
     echo "------ Checking if gatk vcf file exists"
     if [ -f "${3}" ]
     then 
          echo "Gatk vcf found ..... Skipping the gatk step"
          return
     else
          gatk3 -T HaplotypeCaller -ploidy 1 -I "${1}" -R "${2}" -o "${3}"
     fi
}

skip_indels () {
     # This function aims to return a vcf file with indels skipped
     # Has 2 arguments
     # $1 gatk_vcf 
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
     echo "Normalization started"
     if [ -f ${3} ]
     then
          echo "Skip normalization"
          return
     else
          bcftools norm  -f "${1}" "${2}"  -Oz -o "${3}"
          tabix -p vcf "${3}"
          echo "Normalization ended"
     fi

}

          
for specie in ${species[@]}
do
     base=$(echo "$specie" | sed 's/^barcode0[1-8]_//')
     specie=$(echo "$specie" | sed 's/_gatk_results//')
     gatk_dir=${gatk}/${specie}_gatk_results
     BAM=${gatk_dir}/${specie}.bam
     SAM=${gatk_dir}/${specie}.sam

     ref_genome=${ref_path}/${specie}_ref.fna
     
     gatk_vcf=${gatk_dir}/${specie}_illumina_original.vcf.gz
     
     if [ ${specie} = 'barcode06_CHF10J1' ]
     then
        #base="Salmonella_isangi_17-762-33892"
        read1="${reads}/Salmonella_isangi_17-762-33892/17762-33892_1_71_bbduk_1.fastq.gz"
        read2="${reads}/Salmonella_isangi_17-762-33892/17762-33892_1_71_bbduk_2.fastq.gz"
     else
        read1=${reads}/${base}/$(ls ${reads}/${base} | grep "R1_001.fastq.gz")
        read2=${reads}/${base}/$(ls ${reads}/${base} | grep "R2_001.fastq.gz")
     fi
     if [ ! -d "${gatk_dir}" ]
     then
          mkdir -p "${gatk_dir}"
     fi
     # First BWA mapping
     bwa_map "${ref_genome}" ${read1} ${read2} "${SAM}"

     # Sort & index 
     sort_index "${SAM}" "${BAM}"

     # Check sequence dictionary
     CreateSequenceDictionary "${ref_genome}" ${ref_path}/"${specie}"_ref.dict
     # Run gatk
     run_gatk "${BAM}" "${ref_genome}" "${gatk_vcf}"

     # Skip indels
     indel_skiped_vcf=${gatk}/${specie}_gatk_results/${specie}_gatk_indel_skiped.vcf.gz
      
     skip_indels "${gatk_vcf}"  ${indel_skiped_vcf}

     # Then proceed with normalization
     
     norm_vcf=${gatk}/${specie}_gatk_results/${specie}_gatk.norm.vcf.gz
     normalize  "${ref_genome}" ${indel_skiped_vcf} ${norm_vcf}

     #echo "----- Removing sam file ${SAM} -----"
     #rm "${SAM}"

done


