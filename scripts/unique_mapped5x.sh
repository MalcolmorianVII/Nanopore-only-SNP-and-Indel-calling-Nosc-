#! /usr/bin/bash
set -e
eval "$(conda shell.bash hook)"
conda activate "snakemake"

reads=/home/ubuntu/data/belson/bioinformatics/projects_2021/napa/input_data/illumina_reads
out=/home/ubuntu/data/belson/bioinformatics/projects_2021/nosc/results/uniquely_mapped/2022.02.07
refs=/home/ubuntu/data/belson/bioinformatics/projects_2021/nosc/reference_genomes
species=(barcode01_Acinetobacter_baumannii_J9
 barcode02_Citrobacter_koseri_MINF_9D 
 barcode03_Enterobacter_kobei_MSB1_1B 
 barcode04_Haemophilus_unknown_M1C132_1 
 barcode05_Klebsiella_oxytoca_MSB1_2C 
 barcode06_CHF10J1 
 barcode07_Klebsiella_variicola_INF345 
 barcode08_Serratia_marcescens_17-147-1671)

for sp in "${species[@]}"
do
    base=$(echo "$sp" | sed 's/^barcode0[1-8]_//')
    out_dir=${out}/${sp}_uniquely_mapped
    BAM=${out_dir}/${sp}.bam
    unique_bam=${out_dir}/${sp}_unique.bam
    less5X_mapped=${out_dir}/${sp}_less5xcoverage.txt
    SAM=${out_dir}/${sp}.sam
    less5_bed=${out_dir}/${sp}_less5x.bed
    ref_genome=${refs}/${sp}_ref.fna
    if [ "$base" = 'CHF10J1' ]
    then
        #base="Salmonella_isangi_17-762-33892"
        read1="${reads}/Salmonella_isangi_17-762-33892/17762-33892_1_71_bbduk_1.fastq.gz"
        read2="${reads}/Salmonella_isangi_17-762-33892/17762-33892_1_71_bbduk_2.fastq.gz"
    else
        read1=${reads}/${base}/$(ls ${reads}/${base} | grep "R1_001.fastq.gz")
        read2=${reads}/${base}/$(ls ${reads}/${base} | grep "R2_001.fastq.gz")
    fi
     
    if [ -f "${out_dir}"/success.txt ]
    then 
        echo "NO NEED TO DO THIS"
     
    else
        mkdir -p "${out_dir}"
        #  Create index file
        echo "STARTED WITH BWA INDEXING"
        bwa index "${ref_genome}"
        # Ensure SAM/BAM files have read groups
        echo "BWA INDEXING DONE"
        echo "STARTING BWA MAPPING"
        bwa mem -R '@RG\tID:dip1\tSM:dip1' "${ref_genome}" ${read1} ${read2} > "${SAM}"
        echo "BWA MAPPING DONE"
        echo "STARTING SORTING AND INDEXING"
        samtools sort -@ 8 -o "${BAM}" "${SAM}" && samtools index "${BAM}"

        echo "STARTED UNIQUE MAPPING"
        samtools view -h "${BAM}" | grep -v -e 'XA:Z:' -e 'SA:Z:' | samtools view -b > "${unique_bam}"
          
        echo "Generating less < 5x coverage depth"
        samtools depth -aa "${unique_bam}" | awk '{if ($3 < 5) print $0}' > "${less5X_mapped}"

        echo "Filtering the less < 5x coverage depth"
        cut -f 1 "${less5X_mapped}" | sort| uniq | while read X; do awk -v X="$X" '($1==X) { printf("%s\t%d\t%d\n",$1,$2,int($2)+1);}' "${less5X_mapped}" | sort -k1,1 -k2,2n | ~/bedtools merge -i - | sed "s/\$/\t${X}/" ; done | cut -f 1,2,3 > "${less5_bed}"

        touch "${out_dir}"/success.txt
    fi
    if  [ -f ${SAM} ]
    then
 	rm ${SAM}
    fi
done

