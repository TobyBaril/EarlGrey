# Earl Grey

Earl Grey is a full-automated transposable element (TE) annotation pipeline, leveraging the most widely-used tools and combining these with a consensus elongation process to better define _de novo_ consensus sequences when annotating new genome assemblies.

# Installation

Before using Earl Grey, please ensure RepeatMasker (version 4.1.2) and RepeatModeler (version 2.0.2) are installed and configured. If these are not, please follow the instructions below to install these before continuing with Earl Grey Installation. 
NOTE: These instructions are provided to install RepeatMasker, RepeatModeler and related programs with sudo priveleges. If you are working on a shared cluster, please request installation of RepeatMasker and RepeatModeler by your sysadmin before working with Earl Grey. Earl Grey will function with RepeatMasker and RepeatModeler installed in the local path environment, or when the modules are loaded on a HPC cluster.

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

Download TRF from https://tandem.bu.edu/trf/trf.download.html - In most cases the linux 64-bit will be approprite, but please check and download the version applicable to your system. Once TRF has been downloaded, make a note of the full path to the trf409 file. If you are not certain of the full path, please run the following command

```realpath ./trf409*```

Download python3-pip and h5py - Note for RepeatMasker configuration in /usr/local/ root priveleges are required for h5py, hence installation with sudo!

```
sudo apt update
sudo apt install python3-pip -y
sudo pip install h5py
```


Download RepeatMasker (this will download it to the current directory)

```wget http://www.repeatmasker.org/RepeatMasker/RepeatMasker-4.1.2-p1.tar.gz```

Copy the RepeatMasker package to /usr/local/, or somewhere that all users will be able to access the installation. 
Copying to /usr/local/ might require sudo priveleges

```sudo cp RepeatMasker-4.1.2-p1.tar.gz /usr/local/```

Change directory to /usr/local, and extract the RepeatMasker package. This might require sudo priveleges.

```cd /usr/local/```

```sudo tar -zxvf RepeatMasker-4.1.2-p1.tar.gz```

Install the required RepeatMasker libraries - Earl Grey has been tested with Dfam 3.3 and RepBase. 
Unfortunately, RepBase is now behind a paywall, but to ensure Earl Grey remains open it does not rely on RepBase. 
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

Now to finally configure RepeatMasker! - Run these lines and then follow the on-screen prompts from RepeatMasker. Again, this might require sudo priveleges.

```cd /usr/local/RepeatMasker```

```sudo perl ./configure```

Add RepeatMasker and util directory to your path environment

```export PATH=$PATH:/usr/local/RepeatMasker:/usr/local/RepeatMasker/util/```

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

Some perl modules are required by RepeatModeler. Install these:

```
sudo cpan install JSON File::Which URI LWP::UserAgent Devel::Size
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

Time to configure RepeatModeler! You need to enter the paths to lots of the directories of programs we have installed, please note them down before running the configuration script!!

```
cd /usr/local/RepeatModeler-2.0.2a/
sudo perl ./configure
```

Add RepeatModeler directory to your path environment

```export PATH=$PATH:/usr/local/RepeatModeler-2.0.2a/```

## Earl Grey Installation and Configuration

All of the scripts and associated modules are contained within this github repository. Earl Grey runs inside an anaconda environment to ensure all required packages are present and are the correct version. Therefore, to run Earl Grey, you will require anaconda to be installed on your system. 

If anaconda is NOT installed on your system, please install it following these instructions:

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
