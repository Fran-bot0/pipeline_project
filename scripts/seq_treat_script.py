#!/usr/bin/env python

import os
import subprocess
import sys

TRIM_LEAD = '3'
TRIM_TRAIL = '3'
TRIM_MIN_LEN = '36'
TRIM_SLIDE = '4:15'
JAVA_CMD = ['java', '-Xmx4g','-jar']

organism = sys.argv[1]
adapter = '/home/francisco/pipeline_project/software/Trimmomatic-0.39/adapters/TruSeq3-PE-2.fa'
path_gatk = '/home/francisco/pipeline_project/software/gatk-4.2.6.1/gatk-package-4.2.6.1-local.jar'
path_genome_index = f'/home/francisco/pipeline_project/pipeline_data/{organism}/input_files/STAR'
path_picard = '/home/francisco/pipeline_project/software/picard/picard.jar'
path_primary_assebly = f'/home/francisco/pipeline_project/pipeline_data/{organism}/input_files/primary_assembly_toindex_{organism}.fa.gz'
path_snps = f'/home/francisco/pipeline_project/pipeline_data/{organism}/input_files/vcf_{organism}.vcf.gz'
path_trimmomatic = '/home/francisco/pipeline_project/software/Trimmomatic-0.39/trimmomatic-0.39.jar'


def run(operation, output_file=None):
    '''Takes a command and an output file and runs the command in the shell if the output file does not yet exist.'''
    if output_file and os.path.exists(output_file):
        print(f'The {output_file} file is already in this directory. Proceeding to next step...')
    else:
        try:
            subprocess.run(operation, shell=False, check=True)
        except subprocess.CalledProcessError as e:
            print(e)
            exit()


nr_accessions = int(input('How many accession do you want to work with? '))
accessions = []
for num in range(1, nr_accessions+1):
    accessions.append(input(f'Please insert acession number {num} '))


for accession in accessions:

    output_path = f'/home/francisco/pipeline_project/pipeline_data/{organism}/output_files/{accession}/'
    if not os.path.exists(output_path):
        os.makedirs(output_path)
    os.chdir(output_path)

    # Download the reads from SRA
    run(['fasterq-dump', '-p', '--split-files', accession], f'{accession}_1.fastq.gz')


    # bgzip the downloaded right and left reads
    run(['bgzip', f'{accession}_1.fastq'], f'{accession}_1.fastq.gz')
    run(['bgzip', f'{accession}_2.fastq'], f'{accession}_2.fastq.gz')


    # Trim the reads using Trimmomatic

    trim_input = [f'{accession}_1.fastq.gz',
                  f'{accession}_2.fastq.gz']

    trim_output =  [f'{accession}_1.P.fastq.gz',
                    f'{accession}_1.U.fastq.gz',
                    f'{accession}_2.P.fastq.gz',
                    f'{accession}_2.U.fastq.gz']

    trim_opt = [f'ILLUMINACLIP:{adapter}:2:30:10',
                f'LEADING:{TRIM_LEAD}',
                f'TRAILING:{TRIM_TRAIL}',
                f'SLIDINGWINDOW:{TRIM_SLIDE}',
                f'MINLEN:{TRIM_MIN_LEN}']

    trim_command = JAVA_CMD + [path_trimmomatic, 'PE'] + trim_input + trim_output + trim_opt

    run(trim_command, trim_output[0])


    # Unzip the trimmed left and right reads to be used in STAR
    run(['gunzip', '-k', trim_output[0]], trim_output[0][:-3])
    run(['gunzip', '-k', trim_output[2]], trim_output[2][:-3])


    # Align the reads using STAR
    aligned_file = f'{accession}_Aligned.sortedByCoord.out.bam'

    star_command = ['STAR',
                    '--runMode', 'alignReads',
                    '--limitBAMsortRAM', '1362848454',
                    '--readFilesIn', f'{trim_output[0][:-3]}', f'{trim_output[2][:-3]}', 
                    '--genomeDir', path_genome_index,
                    '--outSAMtype', 'BAM', 'SortedByCoordinate',
                    '--outFileNamePrefix', f'{accession}_']

    run(star_command, aligned_file)


    # Indexing the aligned reads with samtools
    run(['samtools', 'index', aligned_file, f'{aligned_file}.bai'])


    # Remove PCR and optical duplicates with picard
    removed_duplicates_file = f'{accession}_aligned_rd.bam'

    pic_command = JAVA_CMD + [path_picard, 'MarkDuplicates', '--REMOVE_DUPLICATES', 'true', 
                            '-I', aligned_file, 
                            '-O', removed_duplicates_file, 
                            '-M', f'{accession}_metrics.txt']

    run(pic_command, removed_duplicates_file)


    # Indexing the aligned reads without duplicates with samtools
    run(['samtools', 'index', removed_duplicates_file, f'{removed_duplicates_file}.bai'])


    # Add more read groups to the bam file (name, sample run, and library)
    read_groups_file = f'{accession}_aligned_rd_rg.bam'
    read_groups_command = JAVA_CMD + [path_picard, 'AddOrReplaceReadGroups', 
                                    '-I', removed_duplicates_file, 
                                    '-O', read_groups_file,
                                    '-PL', 'ILLUMINA', '-PU', 'run', '-LB', accession[3:], '-SM', accession]
    run(read_groups_command, read_groups_file)


    # Indexing the aligned reads with more info
    run(['samtools', 'index', read_groups_file, f'{read_groups_file}.bai'])


    # Making and checking recalibration table
    cov_file1 = f'{accession}_cov1.txt'

    recalibration_command1 = JAVA_CMD + [path_gatk, 'BaseRecalibrator', 
                                            '-R', path_primary_assebly, '--known-sites', path_snps, 
                                            '-I', read_groups_file, 
                                            '-O', cov_file1]

    run(recalibration_command1, cov_file1)


    # Base recalibration using GATK
    recal_file = f'{accession}_aligned_rd_rg_rc.bam'

    bqsr_command = JAVA_CMD + [path_gatk, 'ApplyBQSR', 
                            '-R', path_primary_assebly, '-bqsr', cov_file1, 
                            '-I', read_groups_file, 
                            '-O', recal_file] 

    run(bqsr_command, recal_file)


    # Making a new recalibration table to confirm results after base recalibration
    cov_file2 = f'{accession}_cov2.txt'

    recalibration_command2 = JAVA_CMD + [path_gatk, 'BaseRecalibrator', 
                                        '-R', path_primary_assebly, '--known-sites', path_snps, 
                                        '-I', recal_file, 
                                        '-O', cov_file2]

    run(recalibration_command2, cov_file2)


    # Analyzing covariates using GATK
    pdf_output = f'{accession}_analyze_covariants.pdf'
    cov_command = JAVA_CMD + [path_gatk, 'AnalyzeCovariates', 
                            '-before', cov_file1, 
                            '-after', cov_file2, 
                            '-plots', pdf_output]

    run(cov_command, pdf_output)


    # Indexing recalibrated bam files
    run(['samtools', 'index', recal_file], f'{recal_file}.bai')


    # Mapping statistics
    run(['samtools', 'flagstat', recal_file])
