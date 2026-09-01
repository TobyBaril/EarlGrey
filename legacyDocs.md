# This file contains installation instructions for legacy versions of Earl Grey

## We recommend installing via conda or docker for the latest stable releases of Earl Grey, the below are retained as a reference

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

earlGrey -g genome.fasta -s speciesName -o outputDirectory -t threads
```

For suggestions or questions, please use the discussion and issues functions in github.

Thank you for trying Earl Grey!

#==============================================================================================================================================================================#

## Earl Grey Installation and Configuration (If you DO NOT have RepeatMasker and RepeatModeler) - WITH SUDO PRIVILEGES

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

```wget https://www.repeatmasker.org/rmblast/rmblast-2.13.0+-x64-linux.tar.gz```

Extract the rmblast package

```tar -zxvf rmblast-2.13.0+-x64-linux.tar.gz```

Make a note of the full path to rmblast-2.13.0/bin/ as you will need this in the RepeatMasker configuration
If you are not certain of the full path, please run the following command

```realpath ./rmblast-2.13.0/bin/```

Download RepeatMasker (this will download it to the current directory). NOTE some extra steps are required for Dfam 3.8, please refer to repeatmasker.org/RepeatMasker.

```wget https://www.repeatmasker.org/RepeatMasker/RepeatMasker-4.1.6.tar.gz```

Copy the RepeatMasker package to /usr/local/, or somewhere that all users will be able to access the installation. 
Copying to /usr/local/ might require sudo privileges

```sudo cp RepeatMasker-4.1.6.tar.gz /usr/local/```

Change directory to /usr/local, and extract the RepeatMasker package. This might require sudo privileges.

```cd /usr/local/```

```sudo tar -zxvf RepeatMasker-4.1.6.tar.gz```

Install the required RepeatMasker libraries - Earl Grey has been tested with Dfam 3.4 onwards and RepBase. 
Unfortunately, RepBase is now behind a paywall, but to ensure Earl Grey remains open it does not rely on RepBase, although inclusion of RepBase can improve classification of repeats by RepeatModeler. If you have access to this database, please include it in your configuration of RepeatMasker. 
We recommend that you download Dfam 3.8 before using Earl Grey, and this is required for RepeatMasker 4.1.6. The Dfam library is large - this could take a while!

We recommend downloading Dfam into your home directory (~/) or a subdirectory of home

Change directory to home

```cd ~/```

Download lastest Dfam release - This may take a while
*** NOTE: There are many partitions for Dfam, each containing curated and uncurated TE consensi split by taxonomy. As a minimum, you need to configure RepeatMasker with the root partition, but you can also include others. Note: There may be erroneous annotations in the uncurated Dfam database, use at your own risk***

For whole Dfam library run all of the below, or run each part as required (minimum the root (0) partition):
```
wget https://dfam.org/releases/Dfam_3.8/families/FamDB/dfam38_full.0.h5.gz
wget https://dfam.org/releases/Dfam_3.8/families/FamDB/dfam38_full.1.h5.gz
wget https://dfam.org/releases/Dfam_3.8/families/FamDB/dfam38_full.2.h5.gz
wget https://dfam.org/releases/Dfam_3.8/families/FamDB/dfam38_full.3.h5.gz
wget https://dfam.org/releases/Dfam_3.8/families/FamDB/dfam38_full.4.h5.gz
wget https://dfam.org/releases/Dfam_3.8/families/FamDB/dfam38_full.5.h5.gz
wget https://dfam.org/releases/Dfam_3.8/families/FamDB/dfam38_full.6.h5.gz
wget https://dfam.org/releases/Dfam_3.8/families/FamDB/dfam38_full.7.h5.gz
wget https://dfam.org/releases/Dfam_3.8/families/FamDB/dfam38_full.8.h5.gz
```

Unzip the Dfam release - This may take a while with no indication that anything is happening, please be patient!

For whole Dfam library (including uncurated elements):
```gunzip *.h5.gz```

Move the Dfam library to the RepeatMasker library folder
NOTE, a warning might come up that this will overwrite the existing file, allow this by pressing "y" then Enter

```mv *.h5 /usr/local/RepeatMasker/Libraries/famdb/```

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

```wget http://www.repeatmasker.org/RepeatModeler/RepeatModeler-2.0.4.tar.gz```

Copy the RepeatModeler package to /usr/local/, or somewhere that all users will be able to access the installation. 
Copying to /usr/local/ might require sudo privileges

```sudo cp RepeatModeler-2.0.4.tar.gz /usr/local/```

Change directory to /usr/local, and extract the RepeatModeler package. This might require sudo privileges.

```cd /usr/local/```

```sudo tar -zxvf RepeatModeler-2.0.4.tar.gz```

DO NOT configure RepeatModeler just yet...

Add RepeatModeler directory to your path environment

```echo 'export PATH=$PATH:/usr/local/RepeatModeler-2.0.4/' >> ~/.bashrc```

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

NOTE: For the latest version of RepeatModeler, this must be called trf. We can make a symlinked version: ```ln -s ./modules/trf409* ./modules/trf```

Now to finally configure RepeatMasker! - Run these lines and then follow the on-screen prompts from RepeatMasker. Again, this might require sudo privileges.

```
cd /usr/local/RepeatMasker
# check which perl interpreter you should use 
which perl
# replace [perl] in the below command with the path printed from the command above and use this as your perl interpreter (this should be the anaconda perl interpreter)
sudo [perl] ./configure
```

Time to configure RepeatModeler! You need to enter the paths to lots of the directories of programs we have installed, please note them down before running the configuration script!! 


```
cd /usr/local/RepeatModeler-2.0.4/
# check which perl interpreter you should use - COPY THIS
which perl
# replace perl in the below command with the path printed from the command above and use this as your perl interpreter
sudo perl ./configure
```

Once this is complete, remember to activate the earlGrey conda environment before attempting to run the Earl Grey pipeline

```
conda activate earlGrey

earlGrey -g genome.fasta -s speciesName -o outputDirectory -t threads
```

#==============================================================================================================================================================================#

## Earl Grey Installation and Configuration (If you DO NOT have RepeatMasker and RepeatModeler) - WITHOUT SUDO PRIVILEGES

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

```wget https://www.repeatmasker.org/rmblast/rmblast-2.13.0+-x64-linux.tar.gz```

Extract the rmblast package

```tar -zxvf rmblast-2.13.0+-x64-linux.tar.gz```

Make a note of the full path to rmblast-2.13.0/bin/ as you will need this in the RepeatMasker configuration
If you are not certain of the full path, please run the following command

```realpath ./rmblast-2.13.0/bin/```

Download RepeatMasker (this will download it to the current directory). NOTE some extra steps are required if you wish to use Dfam 3.7, please refer to repeatmasker.org/RepeatMasker for these extra steps if required.

```wget https://www.repeatmasker.org/RepeatMasker/RepeatMasker-4.1.6.tar.gz```

Extract the RepeatMasker package.

```
tar -zxvf RepeatMasker-4.1.6.tar.gz
```

Install the required RepeatMasker libraries - Earl Grey has been tested with Dfam 3.3 and RepBase. 
Unfortunately, RepBase is now behind a paywall, but to ensure Earl Grey remains open it does not rely on RepBase, although inclusion of RepBase can improve classification of repeats by RepeatModeler. If you have access to this database, please include it in your configuration of RepeatMasker. 
We recommend that you download Dfam 3.8 before using Earl Grey, and this is required for RepeatMasker 4.1.6. The Dfam library is large - this could take a while!

Change directory to RepeatMasker Libraries

```cd ./RepeatMasker/Libraries/famdb/```

Download lastest Dfam release - This may take a while
*** NOTE: There are many partitions for Dfam, each containing curated and uncurated TE consensi split by taxonomy. As a minimum, you need to configure RepeatMasker with the root partition, but you can also include others. Note: There may be erroneous annotations in the uncurated Dfam database, use at your own risk***

For whole Dfam library run all of the below, or run each part as required (minimum the root (0) partition):
```
wget https://dfam.org/releases/Dfam_3.8/families/FamDB/dfam38_full.0.h5.gz
wget https://dfam.org/releases/Dfam_3.8/families/FamDB/dfam38_full.1.h5.gz
wget https://dfam.org/releases/Dfam_3.8/families/FamDB/dfam38_full.2.h5.gz
wget https://dfam.org/releases/Dfam_3.8/families/FamDB/dfam38_full.3.h5.gz
wget https://dfam.org/releases/Dfam_3.8/families/FamDB/dfam38_full.4.h5.gz
wget https://dfam.org/releases/Dfam_3.8/families/FamDB/dfam38_full.5.h5.gz
wget https://dfam.org/releases/Dfam_3.8/families/FamDB/dfam38_full.6.h5.gz
wget https://dfam.org/releases/Dfam_3.8/families/FamDB/dfam38_full.7.h5.gz
wget https://dfam.org/releases/Dfam_3.8/families/FamDB/dfam38_full.8.h5.gz
```

Unzip the Dfam release - This may take a while with no indication that anything is happening, please be patient!

For whole Dfam library (including uncurated elements):
```gunzip *.h5.gz```

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

```wget http://www.repeatmasker.org/RepeatModeler/RepeatModeler-2.0.4.tar.gz```

Unpack RepeatModeler

```sudo tar -zxvf RepeatModeler-2.0.4.tar.gz```

DO NOT configure RepeatModeler just yet...

Add RepeatModeler directory to your path environment (Replace /path/to/ with the full path to your installation directory)


```echo 'export PATH=$PATH:/path/to/RepeatModeler-2.0.4/' >> ~/.bashrc```

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

NOTE: For the latest version of RepeatModeler, this must be called trf. We can make a symlinked version: ```ln -s ./modules/trf409* ./modules/trf```

Now to finally configure RepeatMasker! - Run these lines and then follow the on-screen prompts from RepeatMasker. Again, this might require sudo privileges.

```
cd /path/to/RepeatMasker
perl ./configure
```

Time to configure RepeatModeler! You need to enter the paths to lots of the directories of programs we have installed, please note them down before running the configuration script!! 


```
cd /path/to/RepeatModeler-2.0.4/
perl ./configure
```

Once this is complete, remember to activate the earlGrey conda environment before attempting to run the Earl Grey pipeline

```
conda activate earlGrey

earlGrey -g genome.fasta -s speciesName -o outputDirectory -t threads
```

# Alternative installation methods

Before using Earl Grey, please ensure RepeatMasker (>=version 4.1.4) and RepeatModeler (>=version 2.0.4) are installed and configured. If these are not, please follow the instructions below to install these before continuing with Earl Grey Installation. 
NOTE: These instructions are provided to install RepeatMasker, RepeatModeler and related programs with sudo privileges. If you are working on a shared cluster, please request installation of RepeatMasker and RepeatModeler by your sysadmin before working with Earl Grey. Earl Grey will function with RepeatMasker and RepeatModeler installed in the local path environment, or when the modules are loaded on a HPC cluster.


#==============================================================================================================================================================================#

## Earl Grey is available in a [Docker container](./Docker) preconfigured with Dfam version 3.9 or empty, but with the option to install required partitions

To use this container, please make sure Docker is installed and configured on your system. All instructions relating to the docker installation are found within the [Docker directory](./Docker) in this repository. Please consult the README in the Docker directory for installation instructions.

#==============================================================================================================================================================================#

## Earl Grey is available for [Singularity](./Singularity) preconfigured with Dfam version 3.9 or empty, but with the option to install required partitions

To use this container, please make sure Singularity is installed and configured on your system. All instructions relating to the singularity installation are found within the [Singularity directory](./Singularity) in this repository. Please consult the README in the Singularity directory for installation instructions.
