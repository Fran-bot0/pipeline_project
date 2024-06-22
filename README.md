# RNA-seq-pipeline
Description.

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

The script software_installation.sh, has the commands to install all of the required software.This script is inside the software folder.

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
You will need to atleast edit the user in the paths used in the python scripts and the pipe_script.sh, which can be found in the scripts folder.

Inside the software folder there is a script whcih has the commands to install the software needed to run the pipeline. However, for some softwares, some intermediate steps are necessary:
* For the installation of Trimmomatic to be successful, after cloning the repo, you'll need to edit the build.xml file and change the source and target values from 1.5 to 1.8 or above. After this, the ant command will work.
* After installing the sratoolkit, to be able to use it, you'll need to create or select a folder to be used as cache. The command needed to perform this action is in software_installation.sh.

Finally, for the pipe_script.sh to be able to run the different software, you'll have to create a bin directory on the same level as the scripts and software directories. This bin directory should contian symbolic links leading to the desired softwares. The code needed to achive this is also in software_installation.sh




