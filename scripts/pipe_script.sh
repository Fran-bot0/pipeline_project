#!/bin/bash

# Exporting bin path for the symbolic links to work
export PATH="$HOME/pipeline_project/bin:$PATH"

# Asking the user which organism the pipeline will run with.
read -p  'Please insert of which organism to download the primary_assembly (ex: Homo sapiens): ' input_organism

# Making the input_organism lowercase and separated by an underscore.
organism=$(tr -s ' ' '_' <<< ${input_organism,,})

# Creating the necesary directories if none exist.
mkdir -p /home/francisco/pipeline_project/pipeline_data/$organism/input_files
mkdir -p /home/francisco/pipeline_project/pipeline_data/$organism/output_files
cd /home/francisco/pipeline_project/pipeline_data/$organism/input_files

# Downloads the primary assembly of the specified organism from Ensembl. If the primary assembly doesn't exist, it 
# downloads the toplevel file. Also downloads the gtf and the vcf files.
python /home/francisco/pipeline_project/scripts/downloader.py $organism || exit

# Extracting gtf file
yes n | gunzip -k gtf_$organism.gtf.gz

# Sort the unsorted gtf file
[ -f gtf_$organism.sorted.gtf.gz ] || (grep ^"#" gtf_$organism.gtf; grep -v ^"#" gtf_$organism.gtf | sort -k1,1 -k4,4n) | bgzip  > gtf_$organism.sorted.gtf.gz

# Converting primary_assembly from gzip to bgzip to be used with samtools
[ -f primary_assembly_toindex_$organism.fa.gz ] || gunzip -c primary_assembly_$organism.fa.gz | bgzip > primary_assembly_toindex_$organism.fa.gz

# Unziping primary_assembly to be used by STAR
yes n | gunzip -k primary_assembly_$organism.fa.gz

# Indexing the genome using samtools
samtools faidx primary_assembly_toindex_$organism.fa.gz

# Indexing the genome/transcriptome using STAR
[[ -d STAR ]] || STAR --runThreadN 2 \
		      --runMode genomeGenerate \
    		      --genomeDir STAR \
       		      --sjdbGTFfile gtf_$organism.gtf \
    		      --genomeFastaFiles primary_assembly_$organism.fa \
     		      --sjdbOverhang 100 \
     	 	      --genomeSAindexNbases 10 \
		      --limitBAMsortRAM 1362848454

# Indexing the vcf and sorted gtf files using tabix
tabix vcf_$organism.vcf.gz
tabix gtf_$organism.sorted.gtf.gz

# Creating a sequence dictionary using picard
java -Xmx4g -jar /home/francisco/pipeline_project/software/picard/picard.jar CreateSequenceDictionary -R primary_assembly_toindex_$organism.fa.gz

cd /home/francisco/pipeline_project/pipeline_data/$organism/output_files

# Takes the accessions the user wants to work with and returns the bam files and bai files necessary for the next step
python /home/francisco/pipeline_project/scripts/seq_treat_script.py $organism || exit
