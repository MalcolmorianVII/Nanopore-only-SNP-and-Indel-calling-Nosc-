#! /usr/bin/bash
set -e
eval "$(conda shell.bash hook)"
conda activate "gatk"
reads=/home/ubuntu/data/belson/bioinformatics/projects_2021/napa/input_data/clean_illumina_reads
out=/home/ubuntu/data/belson/bioinformatics/projects_2021/nosc/results/nosc_gatk/2022.01.05
refs=/home/ubuntu/data/belson/bioinformatics/projects_2021/nosc/reference_genomes
species=(barcode01_Acinetobacter_baumannii_J9
 barcode02_Citrobacter_koseri_MINF_9D 
 barcode03_Enterobacter_kobei_MSB1_1B 
 barcode04_Haemophilus_unknown_M1C132_1 
 barcode05_Klebsiella_oxytoca_MSB1_2C 
 barcode06_CHF10J1 
 barcode07_Klebsiella_variicola_INF345 
 barcode08_Serratia_marcescens_17-147-1671)

for sp in ${species[@]}
do
     base=$(echo $sp | sed 's/^barcode0[1-8]_//')
     out_dir=${out}/${sp}_gatk_results
     BAM=${out_dir}/${sp}.bam
     SAM=${out_dir}/${sp}.sam
     ref_genome=${refs}/${sp}_ref.fna
     read1=${reads}/${base}/${base}_R1.fastq.gz
     read2=${reads}/${base}/${base}_R2.fastq.gz
     if [ $base = 'CHF10J1' ]
     then
          base="Salmonella_isangi_17-762-33892"
          read1="${reads}/${base}/17762-33892_1_71_bbduk_1.fastq.gz"
          read2="${reads}/${base}/17762-33892_1_71_bbduk_2.fastq.gz"
     fi
     
     if [ -f ${out_dir}/${sp}_illumina_original.vcf.gz ]
     then 
          echo "PREVIOUS GATK RUN SUCCESSFUL"
          echo "SKIPPING GATK RUN FOR ${base}"
     else
          mkdir -p ${out_dir}
          #  Create index file
          echo "STARTED WITH BWA INDEXING"
          bwa index ${ref_genome}
          # Ensure SAM/BAM files have read groups
          echo "BWA INDEXING DONE"
          echo "STARTING BWA MAPPING"
          bwa mem -R '@RG\tID:dip1\tSM:dip1' ${ref_genome} ${read1} ${read2} > ${SAM}
          echo "BWA MAPPING DONE"
          echo "STARTING SORTING AND INDEXING"
          samtools sort -@ 8 -o ${BAM} ${SAM} && samtools index ${BAM}
          if [ -f ${refs}/${sp}_ref.dict ]
          then
               echo "SEQUENCE DICTIONARY ALREADY EXISTS"
          else
               echo "CREATING SEQUENCE DICTIONARY"
               java -jar /home/belson/picard.jar CreateSequenceDictionary -R ${ref_genome}
          fi
          echo "STARTED RUNNING GATK3"
          gatk3 -T HaplotypeCaller -ploidy 1 -I ${BAM} -R ${ref_genome} -o ${out_dir}/${sp}_illumina_original.vcf.gz

          echo "<<< REMOVING SAM FILE: ${SAM} >>>"
          rm ${SAM}
     fi

done

