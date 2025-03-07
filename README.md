# RNA-seq-pipeline
### How it works
The pipe_script.sh is a bash script which uses 2 python scripts (downloader.py and seq_treat_script.py) to take user input and return the final bam files with the raw counts ready to be analysed. The graphs.py script takes the bam files and returns a clustermap and a volcano plot.
1. In the first step of the pipeline, the downloader.py script asks the user which organism they want to work with. Then, it downloads the primary assembly, the vcf, and the gtf files from ensembl.
2. The pipe_script.sh does some file preparations to the downloaded files, such as unziping sorting and indexing.
3. The seq_treat_script.py starts by asking the user with how many and which accession from NCBI they want to work with and downloads them using fasterq-dump.
4. The reads are trimmed using Trimmomatic.
5. The alignment is made using STAR
6. The duplicates are removed using pacard.
7. Read groups are added to the bam file using picard.
8. Base recalibration using GATK.
9. Indexing the recalibrated bam files using samtools.
10. Finally, the graph.py takes the bam files, calculates the RPKM values, performs quantile normalization, executes a t-test between the samples, calculates the fold change, and returns a clustermap and a volcano plot.

Currently, the pipeline up until returning the bam files for analysis (Step 9) should work with any user input as long as the organism exists in ensembl (with available gtf and vcf files) and the accession(s) exist in NCBI. However, the graphs.py scrip is still hard-coded to work with the bam files originated by running the pipeline with the organism and accessions used in this [study](https://doi.org/10.3389/fmicb.2024.1394880). You can check the resulting clustermap and volcano plot in the img directory.

(Organism = Saccharomyces cerevisiae; Accessions = SRR28797006, SRR28797007, SRR28797008, SRR28797012, SRR28797013, SRR28797014)

# Requirements
* Ubuntu 24.04
* Python 3.12.3
* STAR aligner 2.7.10b
* Trimmomatic 0.39
* bcftools 1.9
* htslib 1.9
* samtools 1.9
* sratoolkit 3.0.0
* gatk 4.2.6.1
* picard 2.26.2

The script software_installation.sh has the commands to install all of the required software. This script is inside the software folder.

# Python Modules, Packages and Libraries
* BeautifulSoup
* HTSeq
* matplotlib
* numpy
* os
* pandas
* pysam
* requests
* scipy
* subprocess
* sys
* urllib

# To make the code work...
You will need to at least edit the user in the paths used in the python scripts and pipe_script.sh.

You'll need to give execute permissions to pipe_script.sh using chmod before running the script.

Inside the software folder there is a script which has the commands to install the software needed to run the pipeline. However, for some softwares, intermediate steps are necessary:
* For the installation of Trimmomatic to be successful, after cloning the repo, you'll need to edit the build.xml file and change the source and target values from 1.5 to 1.8 or above. After this, the ant command will work.
* After installing the sratoolkit, to be able to use it, you'll need to create or select a folder to be used as cache. The command needed to perform this action is in software_installation.sh.

Finally, for the pipe_script.sh to be able to run the different software, you'll have to create a bin directory on the same level as the software directory. This bin directory should contian symbolic links leading to the desired softwares. The code needed to achieve this is also in software_installation.sh.

To better understand how the directories should be organized, you can follow this scheme:


![file_order](https://github.com/Fran-bot0/pipeline_project/blob/main/zhow_directories_should_be.png)
