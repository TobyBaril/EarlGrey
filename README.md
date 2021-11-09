![earlGreyIcon](https://user-images.githubusercontent.com/46785187/136248346-21e980ee-1154-48c2-9398-70938bbe2404.png)

![GitHub all releases](https://img.shields.io/github/downloads/TobyBaril/EarlGrey/total) [![DOI](https://zenodo.org/badge/412126708.svg)](https://zenodo.org/badge/latestdoi/412126708) ![Twitter URL](https://img.shields.io/twitter/url?style=social&url=https%3A%2F%2Fgithub.com%2FTobyBaril%2FEarlGrey%2F) 

# Earl Grey

Earl Grey is a full-automated transposable element (TE) annotation pipeline, leveraging the most widely-used tools and combining these with a consensus elongation process to better define _de novo_ consensus sequences when annotating new genome assemblies.

# Contents

[Example Run](#example)

[References and Acknowledgements](#references-and-acknowledgements)

[Installation](#installation)


<!-- toc -->

# Example

Given an input genome, Earl Grey will run through numerous steps to identify and mask transposable elements (TEs). We recommend running earlGrey within a tmux or screen session, so that you can log off and leave Earl Grey running.

```
# remember to activate the conda environment before running

conda activate earlGrey

# run earl grey

earlGrey -g [genome.fasta] -s [speciesName] -r [repeatmasker-search-term] -o [outputDirectory] -t [numberThreads]

# e.g

earlGrey -g myzusPersicae.fasta -s myzusPersicae -r arthropoda -o ./earlGreyOutputs -t 8

```

The runtime of Earl Grey will depend on the repeat content of your input genome. Once finished, you will notice that a number of directories have been created by Earl Grey. The most important results are found within the "summaryFiles" folder, however intermediate results are kept in case you wish to use alignments for further manual curation or investigation, for example. 

Directories created by earl grey:

- [speciesName]EarlGrey/
  * [speciesName]_RepeatMasker/
    + Results of the initial RepeatMasker run used to mask previously characterised TEs
  * [speciesName]_Database/
    + Database created from the masked genome output of the initial RepeatMasker run. Required for RepeatModeler.
  * [speciesName]_RepeatModeler/
    + Results of the RepeatModeler _de novo_ TE identification step
  * [speciesName]_BLASTN/
    + Results of the first BLASTn search of _de novo_ TEs against the input genome (Step 1 of the "BLAST, Extract, Extend" process)
  * [speciesName]_ExtractAlign/
    + Results of all iterations of the "BLAST, Extract, Extend" process
  * [speciesName]_Masked_de_novo_Repeats/
    + Intermediate results for the reclassification of _de novo_ TEs following the iterative "BLAST, Extract, Extend" process
  * [speciesName]_Curated_Library/
    + Contains the _de novo_ repeat libraries generated with each iteration of the "BLAST, Extract, Extend" process, the library of known repeats used by RepeatMasker, and a combined library containing both sets of repeats (used in the final RepeatMasker run)
  * [speciesName]_RepeatMasker_Against_Custom_Library/
    + Results of the ReepatMasker run using known repeats AND _de novo_ identified repeats
  * [speciesName]_RepeatLandscape/
    + Intermediate files for the generation of RepeatLandscapes (RepeatMasker .divsum files)
  * [speciesName]_mergedRepeats/
    + Intermediate files and results of TE defragmentation step using RepeatCraft.
  * [speciesName]_summaryFiles/
    + Results and plots from Earl Grey:
    + TE annotations in GFF3 and BED format
    + Repeat Landscape showing TE activity (PDF)
    + Pie chart of genome repeat content (PDF)
    + _de novo_ repeat library in FASTA format
    + Combined repeat library in FASTA format (used by final repeat mask)
  * [speciesName]_clusTErs/
    + BED file containing coordinates of TE clusters across genome scaffolds

### Example Outputs (NOTE: example data has been used here):

- Pie chart summarising TE content in input genome

![image](https://user-images.githubusercontent.com/46785187/140897482-fc80c30a-3e5f-4bf6-99b0-0f18b11b33d1.png)

- RepeatLandscape summarising relative TE activity (recent activity towards the RHS)

![image](https://user-images.githubusercontent.com/46785187/140897874-5a2cbc9b-808e-4f51-b781-3b3a22b1b223.png)

- TE annotations

```
# BED format
NC_045808.1     4964941 4965925 LINE/Penelope   5073    +
NC_045808.1     7291353 7291525 LINE/L2 1279    +
NC_045808.1     8922477 8923791 DNA/TcMar-Tc1   11957   +
NC_045808.1     9848807 9849591 LTR/Copia       734     -


# GFF3 format
NC_045808.1     RepeatMasker    LINE/Penelope   4964942 4965925 5073    +       .       ID="RND-1_FAMILY-48";Tend="2677";Tstart="1700";shortTE="F";uid="cb016329-11ac-492a-b11c-db2dc50e9d89"
NC_045808.1     RepeatMasker    LINE/L2 7291354 7291525 1279    +       .       ID="RND-5_FAMILY-151";Tend="1398";Tstart="1226";shortTE="F";uid="b1ffaec9-e362-4ffe-89c2-5830cea00ab7"
NC_045808.1     RepeatMasker    DNA/TcMar-Tc1   8922478 8923791 11957   +       .       ID="RND-1_FAMILY-124";Tend="3292";Tstart="1979";shortTE="F";uid="1e507e68-39f3-45f0-b62b-65bd9b7919c8"
NC_045808.1     RepeatMasker    LTR/Copia       9848808 9849591 734     -       .       ID="COPIA-3_BTA-I";Tend="2057";Tstart="1246";shortTE="F";uid="2a786d4d-4d0d-4013-824e-92f97c60ecfe"
```

# References and Acknowledgements

This pipeline has been designed to be used and shared openly by the community.

### When using Earl Grey, please cite:

Baril, Tobias, Imrie, Ryan, and Hayward, Alexander. (2021) Earl Grey. Zenodo [doi:10.5281/zenodo.5654616](https://doi.org/10.5281/zenodo.5654616)

### This pipeline makes use of scripts from:

[Bioinfo Tools (David Ray)](https://github.com/davidaray/bioinfo_tools) - which provides the extraction script for the automated version of the "BLAST, Extract, Extend" desrcibed in: Platt RN, Blanco-Berdugo L, Ray DA. Accurate transposable element annotation is vital when analyzing new genome assemblies. Genome Biol Evol 2016;8:403–10. https://doi.org/10.1093/gbe/evw009.

[RepeatCraft](https://github.com/niccw/repeatcraftp) - Wong WY, Simakov O. RepeatCraft: a meta-pipeline for repetitive element de-fragmentation and annotation. Bioinformatics 2018;35:1051–2. https://doi.org/10.1093/bioinformatics/bty745.

### Please also cite the following open source software utilised by this pipeline:

Smit AFA, Hubley RR, Green PR. RepeatMasker Open-4.0. Http://RepeatmaskerOrg 2013.

Smit A, Hubley R. RepeatModeler Open-1.0. Http://RepeatmaskerOrg 2015.

Bao Z, Eddy SR. Automated De Novo Identification of Repeat Sequence Families in Sequenced Genomes. Genome Res 2002;12:1269–76. https://doi.org/10.1101/gr.88502.

Price AL, Jones NC, Pevzner PA. De novo identification of repeat families in large genomes. Bioinformatics 2005;21:i351–8.

Camacho C, Coulouris G, Avagyan V, Ma N, Papadopoulos J, Bealer K, et al. BLAST+: Architecture and applications. BMC Bioinformatics 2009;10: 1–9. https://doi.org/10.1186/1471-2105-10-421.

Katoh K, Standley DM. MAFFT multiple sequence alignment software version 7: Improvements in performance and usability. Mol Biol Evol 2013;30:772–80. https://doi.org/10.1093/molbev/mst010.

Capella-Gutiérrez S, Silla-Martínez JM, Gabaldón T. trimAl: A tool for automated alignment trimming in large-scale phylogenetic analyses. Bioinformatics 2009;25:1972–3. https://doi.org/10.1093/bioinformatics/btp348.

Rice P, Longden L, Bleasby A. EMBOSS: The European Molecular Biology Open Software Suite. Trends Genet 2000;16:276–7. https://doi.org/10.1016/S0168-9525(00)02024-2.

Xu Z, Wang H. LTR_FINDER: an efficient tool for the prediction of full-length LTR retrotransposons. Nucleic Acids Res 2007;35:W265–8. https://doi.org/10.1093/nar/gkm286.

Ou S, Jiang N. LTR_FINDER_parallel: parallelization of LTR_FINDER enabling rapid identification of long terminal repeat retrotransposons. BioRxiv 2019:2–6.

Rubino F, Creevey CJ. MGkit: Metagenomic framework for the study of microbial communities 2014.


# Installation

Before using Earl Grey, please ensure RepeatMasker (version 4.1.2) and RepeatModeler (version 2.0.2) are installed and configured. If these are not, please follow the instructions below to install these before continuing with Earl Grey Installation. 
NOTE: These instructions are provided to install RepeatMasker, RepeatModeler and related programs with sudo priveleges. If you are working on a shared cluster, please request installation of RepeatMasker and RepeatModeler by your sysadmin before working with Earl Grey. Earl Grey will function with RepeatMasker and RepeatModeler installed in the local path environment, or when the modules are loaded on a HPC cluster.


#==============================================================================================================================================================================#

## Earl Grey Installation and Configuration (If you already have RepeatMasker and RepeatModeler)

If you do not currently have RepeatMasker and RepeatModeler installed, the instructions are provided further down this page. If you do have them installed, **ensure the executables are in your PATH environment, including the RepeatMasker/util/ directory!**

All of the scripts and associated modules are contained within this github repository. Earl Grey runs inside an anaconda environment to ensure all required packages are present and are the correct version. Therefore to run Earl Grey, you will require anaconda to be installed on your system. 

**If anaconda is NOT installed on your system, please install it following these instructions:**

```
# Change to /tmp directory as we won't need the script after running it
cd /tmp

# Download the anaconda installation script
curl https://repo.anaconda.com/archive/Anaconda3-2020.02-Linux-x86_64.sh --output anaconda.sh

# Run the script to install anaconda
bash anaconda.sh
# answer yes when asked, and install anaconda3 in the specified location (recommended) unless you want it to be installed elsewhere.
# When asked "Do you wish the installer to initialize Anaconda3 by running conda init?", answer "yes" for ease of use.
# Activate conda by refreshing terminal session

source ~/.bashrc

# If successful, you should now see (base) on the left of your username on the command line
```

**Now that Anaconda is installed, we can install Earl Grey**

Clone the Earl Grey github repo

```
# Clone into a home directory, or somewhere you want to install Earl Grey
git clone https://github.com/TobyBaril/EarlGrey
```

Enter the Earl Grey directory and configure the program

```
cd ./EarlGrey
chmod +x ./configure
./configure
```

Once this is complete, remember to activate the earlGrey conda environment before attempting to run the Earl Grey pipeline

```
conda activate earlGrey

earlGrey -g genome.fasta -s speciesName -o outputDirectory -r repeatMaskerTerm -t threads
```

For suggestions or questions, please use the discussion and issues functions in github.

Thank you for trying Earl Grey!

#==============================================================================================================================================================================#

## Earl Grey Installation and Configuration (If you DO NOT have RepeatMasker and RepeatModeler) - WITH SUDO PRIVELAGES

These instructions will guide you through configuring all required programs and scripts to run Earl Grey.

All of the scripts and associated modules are contained within this github repository. Earl Grey runs inside an anaconda environment to ensure all required packages are present and are the correct version. Therefore to run Earl Grey, you will require anaconda to be installed on your system. 

**If anaconda is NOT installed on your system, please install it following these instructions:**

```
# Change to /tmp directory as we won't need the script after running it
cd /tmp

# Download the anaconda installation script
curl https://repo.anaconda.com/archive/Anaconda3-2020.02-Linux-x86_64.sh --output anaconda.sh

# Run the script to install anaconda
bash anaconda.sh
# answer yes when asked, and install anaconda3 in the specified location (recommended) unless you want it to be installed elsewhere.
# When asked "Do you wish the installer to initialize Anaconda3 by running conda init?", answer "yes" for ease of use.
# Activate conda by refreshing terminal session

source ~/.bashrc

# If successful, you should now see (base) on the left of your username on the command line
```

**Now that Anaconda is installed, we can install Earl Grey**

### To install RepeatMasker

RepeatMasker can be downloaded from: http://www.repeatmasker.org/RepeatMasker/. Installation instructions can be found on the website. Alternatively, please use the code below:

You will need to download and install a couple of programs for RepeatMasker to work. 

Download rmblast

```wget http://www.repeatmasker.org/rmblast-2.11.0+-x64-linux.tar.gz```

Extract the rmblast package

```tar -zxvf rmblast-2.11.0+-x64-linux.tar.gz```

Make a note of the full path to rmblast-2.11.0/bin/ as you will need this in the RepeatMasker configuration
If you are not certain of the full path, please run the following command

```realpath ./rmblast-2.11.0/bin/```

Download RepeatMasker (this will download it to the current directory)

```wget http://www.repeatmasker.org/RepeatMasker/RepeatMasker-4.1.2-p1.tar.gz```

Copy the RepeatMasker package to /usr/local/, or somewhere that all users will be able to access the installation. 
Copying to /usr/local/ might require sudo priveleges

```sudo cp RepeatMasker-4.1.2-p1.tar.gz /usr/local/```

Change directory to /usr/local, and extract the RepeatMasker package. This might require sudo priveleges.

```cd /usr/local/```

```sudo tar -zxvf RepeatMasker-4.1.2-p1.tar.gz```

Install the required RepeatMasker libraries - Earl Grey has been tested with Dfam 3.3 and RepBase. 
Unfortunately, RepBase is now behind a paywall, but to ensure Earl Grey remains open it does not rely on RepBase, although inclusion of RepBase can improve classification of repeats by RepeatModeler. If you have access to this database, please include it in your configuration of RepeatMasker. 
We recommend that you download Dfam 3.3 as a minimum before using Earl Grey. The Dfam library is large - this could take a while!

We recommend downloading Dfam into your home directory (~/) or a subdirectory of home

Change directory to home

```cd ~/```

Download lastest Dfam release - This may take a while

```wget https://www.dfam.org/releases/current/families/Dfam.h5.gz```

Unzip the Dfam release - This may take a while with no indication that anything is happening, please be patient!

```gunzip Dfam.h5.gz```

Move the Dfam library to the RepeatMasker library folder
NOTE, a warning might come up that this will overwrite the existing file, allow this by pressing "y" then Enter

```mv Dfam.h5 /usr/local/RepeatMasker/Libraries/```

DO NOT configure RepeatMasker just yet...

Add RepeatMasker and util directory to your path environment

```echo 'export PATH=$PATH:/usr/local/RepeatMasker:/usr/local/RepeatMasker/util/' >> ~/.bashrc```

### To install RepeatModeler

RepeatModeler can be downloaded from: http://www.repeatmasker.org/RepeatModeler/. Installation instructions can be found on the website. Alternatively, please use the code below:

You will need to download and install a couple of programs for RepeatModeler to work. 

Download RECON and RepeatScout

```
wget http://www.repeatmasker.org/RepeatModeler/RECON-1.08.tar.gz
wget http://www.repeatmasker.org/RepeatScout-1.0.6.tar.gz
```

Extract the RECON and RepeatScout packages

```
tar -zxvf RECON-1.08.tar.gz
tar -zxvf RepeatScout-1.0.6.tar.gz
```

Enter the RECON directory and make from source

```
cd ./RECON-1.08/src/
make
make install
```

Enter the RepeatScout directory and make from source

```
cd ../../RepeatScout-1.0.6/
make
```

Make a note of the full path to RECON-1.08/bin/ and RepeatScout-1.0.6/ as you will need these in the RepeatModeler configuration
If you are not certain of the full paths, please run the following commands

```
realpath ./RECON-1.08/bin/
realpath ./RepeatScout-1.0.6/
```

Install CD-Hit

```
sudo apt install cd-hit
```

Make a note of the full path to CD-Hit as you will need these in the RepeatModeler configuration
If you are not certain of the full path, please run the following command

```
which cd-hit
```

Install UCSC TwoBit Tools. Make a note of the full path to the directory where these have been installed as you will need these in the RepeatModeler configuration

```
mkdir ~/ucscTwoBitTools/
cd ~/ucscTwoBitTools/
wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/faToTwoBit
wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/twoBitToFa
wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/twoBitInfo
wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/twoBitMask
wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/twoBitDup

# make scripts executable
chmod +x ~/ucscTwoBitTools/*
```

For Earl Grey, the LTR structural pipeline of RepeatModeler is not required, as this is run as part of Earl Grey's defragmentation process.

Download RepeatModeler

```wget http://www.repeatmasker.org/RepeatModeler/RepeatModeler-2.0.2a.tar.gz```

Copy the RepeatModeler package to /usr/local/, or somewhere that all users will be able to access the installation. 
Copying to /usr/local/ might require sudo priveleges

```sudo cp RepeatModeler-2.0.2a.tar.gz /usr/local/```

Change directory to /usr/local, and extract the RepeatModeler package. This might require sudo priveleges.

```cd /usr/local/```

```sudo tar -zxvf RepeatModeler-2.0.2a.tar.gz```

DO NOT configure RepeatModeler just yet...

Add RepeatModeler directory to your path environment

```echo 'export PATH=$PATH:/usr/local/RepeatModeler-2.0.2a/' >> ~/.bashrc```

### To Install Earl Grey and Configure all Programs

Clone the Earl Grey github repo

```
# Clone into a home directory, or somewhere you want to install Earl Grey
git clone https://github.com/TobyBaril/EarlGrey
```

Enter the Earl Grey directory and configure the program

```
cd ./EarlGrey
chmod +x ./configure
./configure
```

Activate the conda environment

```
conda activate earlGrey
```

For configuring RepeatMasker with sudo, h5py is required in the root path. The easiest way to do this is to run the following two lines:

```
sudo apt install python3-pip
sudo pip install h5py
```

The linux64 version of TRF is included in the modules directory of Earl Grey. Make a note of the full path to this file as this will be needed in the RepeatMasker configuration.

```realpath ./modules/trf409*```

Now to finally configure RepeatMasker! - Run these lines and then follow the on-screen prompts from RepeatMasker. Again, this might require sudo priveleges.

```
cd /usr/local/RepeatMasker
# check which perl interpreter you should use 
which perl
# replace [perl] in the below command with the path printed from the command above and use this as your perl interpreter (this should be the anaconda perl interpreter)
sudo [perl] ./configure
```

Time to configure RepeatModeler! You need to enter the paths to lots of the directories of programs we have installed, please note them down before running the configuration script!! 


```
cd /usr/local/RepeatModeler-2.0.2a/
# check which perl interpreter you should use - COPY THIS
which perl
# replace perl in the below command with the path printed from the command above and use this as your perl interpreter
sudo perl ./configure
```

Once this is complete, remember to activate the earlGrey conda environment before attempting to run the Earl Grey pipeline

```
conda activate earlGrey

earlGrey -g genome.fasta -s speciesName -o outputDirectory -r repeatMaskerTerm -t threads
```

#==============================================================================================================================================================================#

## Earl Grey Installation and Configuration (If you DO NOT have RepeatMasker and RepeatModeler) - WITHOUT SUDO PRIVELAGES

These instructions will guide you through configuring all required programs and scripts to run Earl Grey.

All of the scripts and associated modules are contained within this github repository. Earl Grey runs inside an anaconda environment to ensure all required packages are present and are the correct version. Therefore to run Earl Grey, you will require anaconda to be installed on your system. 

**If anaconda is NOT installed on your system, please install it following these instructions:**

```
# Change to /tmp directory as we won't need the script after running it
cd /tmp

# Download the anaconda installation script
curl https://repo.anaconda.com/archive/Anaconda3-2020.02-Linux-x86_64.sh --output anaconda.sh

# Run the script to install anaconda
bash anaconda.sh
# answer yes when asked, and install anaconda3 in the specified location (recommended) unless you want it to be installed elsewhere.
# When asked "Do you wish the installer to initialize Anaconda3 by running conda init?", answer "yes" for ease of use.
# Activate conda by refreshing terminal session

source ~/.bashrc

# If successful, you should now see (base) on the left of your username on the command line
```

**Now that Anaconda is installed, we can install Earl Grey**

### To install RepeatMasker

RepeatMasker can be downloaded from: http://www.repeatmasker.org/RepeatMasker/. Installation instructions can be found on the website. Alternatively, please use the code below:

You will need to download and install a couple of programs for RepeatMasker to work. 

Download rmblast

```wget http://www.repeatmasker.org/rmblast-2.11.0+-x64-linux.tar.gz```

Extract the rmblast package

```tar -zxvf rmblast-2.11.0+-x64-linux.tar.gz```

Make a note of the full path to rmblast-2.11.0/bin/ as you will need this in the RepeatMasker configuration
If you are not certain of the full path, please run the following command

```realpath ./rmblast-2.11.0/bin/```

Download RepeatMasker (this will download it to the current directory)

```wget http://www.repeatmasker.org/RepeatMasker/RepeatMasker-4.1.2-p1.tar.gz```

Extract the RepeatMasker package.

```
tar -zxvf RepeatMasker-4.1.2-p1.tar.gz
```

Install the required RepeatMasker libraries - Earl Grey has been tested with Dfam 3.3 and RepBase. 
Unfortunately, RepBase is now behind a paywall, but to ensure Earl Grey remains open it does not rely on RepBase, although inclusion of RepBase can improve classification of repeats by RepeatModeler. If you have access to this database, please include it in your configuration of RepeatMasker. 
We recommend that you download Dfam 3.3 as a minimum before using Earl Grey. The Dfam library is large - this could take a while!


Change directory to RepeatMasker Libraries

```cd ./RepeatMasker/Libraries/```

Download lastest Dfam release - This may take a while

```wget https://www.dfam.org/releases/current/families/Dfam.h5.gz```

Unzip the Dfam release - This may take a while with no indication that anything is happening, please be patient!

```gunzip Dfam.h5.gz```

NOTE, a warning might come up that this will overwrite the existing file, allow this by pressing "y" then Enter

DO NOT configure RepeatMasker just yet...

Add RepeatMasker and util directory to your path environment (Replace /path/to/ with the full path to your installation directory)

```echo 'export PATH=$PATH:/path/to/RepeatMasker/:/path/to/RepeatMasker/util/' >> ~/.bashrc```

### To install RepeatModeler

RepeatModeler can be downloaded from: http://www.repeatmasker.org/RepeatModeler/. Installation instructions can be found on the website. Alternatively, please use the code below:

You will need to download and install a couple of programs for RepeatModeler to work. 

Download RECON and RepeatScout (Make sure you are no longer inside your RepeatMasker directory!)

```
wget http://www.repeatmasker.org/RepeatModeler/RECON-1.08.tar.gz
wget http://www.repeatmasker.org/RepeatScout-1.0.6.tar.gz
```

Extract the RECON and RepeatScout packages

```
tar -zxvf RECON-1.08.tar.gz
tar -zxvf RepeatScout-1.0.6.tar.gz
```

Enter the RECON directory and make from source

```
cd ./RECON-1.08/src/
make
make install
```

Enter the RepeatScout directory and make from source

```
cd ../../RepeatScout-1.0.6/
make
```

Make a note of the full path to RECON-1.08/bin/ and RepeatScout-1.0.6/ as you will need these in the RepeatModeler configuration
If you are not certain of the full paths, please run the following commands

```
realpath ./RECON-1.08/bin/
realpath ./RepeatScout-1.0.6/
```

Install CD-Hit - This step will require SUDO privelages, if you are working on a cluster, this might already be installed. Please check with your sysadmin.

```
sudo apt install cd-hit
```

Make a note of the full path to CD-Hit as you will need these in the RepeatModeler configuration
If you are not certain of the full path, please run the following command

```
which cd-hit
```

Install UCSC TwoBit Tools. Make a note of the full path to the directory where these have been installed as you will need these in the RepeatModeler configuration

```
mkdir ./ucscTwoBitTools/
cd ./ucscTwoBitTools/
wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/faToTwoBit
wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/twoBitToFa
wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/twoBitInfo
wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/twoBitMask
wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/twoBitDup

# make scripts executable
chmod +x ./ucscTwoBitTools/*
```

For Earl Grey, the LTR structural pipeline of RepeatModeler is not required, as this is run as part of Earl Grey's defragmentation process.

Download RepeatModeler

```wget http://www.repeatmasker.org/RepeatModeler/RepeatModeler-2.0.2a.tar.gz```

Unpack RepeatModeler

```sudo tar -zxvf RepeatModeler-2.0.2a.tar.gz```

DO NOT configure RepeatModeler just yet...

Add RepeatModeler directory to your path environment (Replace /path/to/ with the full path to your installation directory)


```echo 'export PATH=$PATH:/path/to/RepeatModeler-2.0.2a/' >> ~/.bashrc```

### To Install Earl Grey and Configure all Programs

Clone the Earl Grey github repo

```
# Clone into a home directory, or somewhere you want to install Earl Grey
git clone https://github.com/TobyBaril/EarlGrey
```

Enter the Earl Grey directory and configure the program

```
cd ./EarlGrey
chmod +x ./configure
./configure
```

Activate the conda environment

```
conda activate earlGrey
```

The linux64 version of TRF is included in the modules directory of Earl Grey. Make a note of the full path to this file as this will be needed in the RepeatMasker configuration.

```realpath ./modules/trf409*```

Now to finally configure RepeatMasker! - Run these lines and then follow the on-screen prompts from RepeatMasker. Again, this might require sudo priveleges.

```
cd /path/to/RepeatMasker
perl ./configure
```

Time to configure RepeatModeler! You need to enter the paths to lots of the directories of programs we have installed, please note them down before running the configuration script!! 


```
cd /path/to/RepeatModeler-2.0.2a/
perl ./configure
```

Once this is complete, remember to activate the earlGrey conda environment before attempting to run the Earl Grey pipeline

```
conda activate earlGrey

earlGrey -g genome.fasta -s speciesName -o outputDirectory -r repeatMaskerTerm -t threads
```

#==============================================================================================================================================================================#



