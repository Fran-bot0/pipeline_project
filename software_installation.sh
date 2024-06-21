#!/bin/bash

mkdir software
cd software

# Install updates and required packages
sudo apt-get update
sudo apt-get install gcc
sudo apt-get install make
sudo apt-get install libbz2-dev
sudo apt-get install zlib1g-dev
sudo apt-get install libncurses5-dev
sudo apt-get install libncursesw5-dev
sudo apt-get install liblzma-dev


# Install HTSLIB
wget https://github.com/samtools/htslib/releases/download/1.9/htslib-1.9.tar.bz2
tar -vxjf htslib-1.9.tar.bz2
cd htslib-1.9
make

# Install Samtools
cd ..
wget https://github.com/samtools/samtools/releases/download/1.9/samtools-1.9.tar.bz2
tar -vxjf samtools-1.9.tar.bz2
cd samtools-1.9
make

# Install BCFTools
cd ..
wget https://github.com/samtools/bcftools/releases/download/1.9/bcftools-1.9.tar.bz2
tar -vxjf bcftools-1.9.tar.bz2
cd bcftools-1.9
make

# Install Picard
cd ..
mkdir picard
cd picard
wget -O picard.jar https://github.com/broadinstitute/picard/releases/download/2.26.2/picard.jar

# Install STAR
cd ..
git clone https://github.com/alexdobin/STAR.git
cd STAR/source
git checkout ee50484
make
cd ..

# Install GATK
cd ..
wget https://github.com/broadinstitute/gatk/releases/download/4.2.6.1/gatk-4.2.6.1.zip
unzip gatk-4.2.6.1.zip

# Install SRAToolkit
wget https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/3.0.0/sratoolkit.3.0.0-centos_linux64-cloud.tar.gz
mkdir sratoolkit
tar -xzf sratoolkit.3.0.0-centos_linux64-cloud.tar.gz -C sratoolkit
# For SRAToolkit to be configured run the folloing command:
# ./sratoolkit/usr/local/ncbi/sra-tools/bin/vdb-config --interactive
# and select/create a direcotory to be used as the cache

# Install Trimmomatic
# git clone https://github.com/usadellab/Trimmomatic.git
# cd Trimmomatic
# For the build of Trimmomatic to work, go inside the build.xml file and replace source=1.5 and target=1.5 by source=1.8 and target=1.8 and run ant command
# Use the command 'nano build.xml' to edit the build.xml file
# ant


# Create symbolic links (Make sure to change the paths of the symlinks)
# cd ..
# cd ..
# mkdir bin
# cd bin
# ln -s /home/user/pipeline_project/software/STAR/source/STAR STAR
# ln -s /home/user/pipeline_project/software/bcftools-1.9/bcftools bcftools
# ln -s /home/user/pipeline_project/software/htslib-1.9/bgzip bgzip
# ln -s /home/user/pipeline_project/software/sratoolkit/sratoolkit.3.1.0-ubuntu64/bin/sratools.3.1.0 fasterq-dump
# ln -s /home/user/pipeline_project/software/samtools-1.9/samtools samtools
# ln -s /home/user/pipeline_project/software/htslib-1.9/tabix tabix
