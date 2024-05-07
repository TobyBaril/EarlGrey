# Get Earl Grey to run in Docker

## A [biocontainer](https://biocontainers.pro/tools/earlgrey) has been generated from the Earl Grey Bioconda package.

## Get the preconfigured docker container from biocontainers
## This image will be configured with Dfam curated elements only
In this case, we need to bind a system directory to the docker container. In the line below, we are binding a directory call `host_data` that is found on our current path to `/data/` in the docker container. Please replace the file path before `:` to the directory you wish to bind to `/data/` in the container. This container must be run in interactive mode the first time you use it.

```
docker run -it -v `pwd`/host_data/:/data/ quay.io/biocontainers/earlgrey:4.2.4--h4ac6f70_0
```

## If you are running the container for the first time, you need to enable Earl Grey to configure the Dfam libraries correctly in interactive mode.
This is done by running Earl Grey for the first time. You will be prompted if Dfam libraries have not been configured properly. Answer `y` to configure Dfam.
Ensure the genome you wish to annotate is found in the directory you want to bind to the container. It will then be found in `/data` in the container. Alternatively, if you source a genome assembly from a database, download it to `/data/`. Always output results to the bound directory. In our case `/data/`. The results will then be found in your system directory after exiting the container.

```
earlGrey -g /data/genome.fasta -s test_genome -t 8 -o /data/
```

## If you need the container to run offline and/or without interactive mode
I try to keep an up-to-date container in docker hub, but this might not always be the case depending on if I have had time to build and upload a new image. Currently, there is an image with Dfam 3.7 curated elements only, and this is version 4.2.4. You can use this image by pulling the container:

```
# Interactive mode
docker run -it -v `pwd`:/data/ tobybaril/earlgrey_dfam3.7:latest

# Non interactive mode example:
docker run -v `pwd`:/data/ tobybaril/earlgrey_dfam3.7:latest earlGrey -g /data/NC_045808_EarlWorkshop.fasta -s nonInteractiveTest -o /data/ -t 8
```

# OR Build the container (not recommended)...

## Note, the docker container has been tested with Dfam 3.7 without RepBase. If you would like to use with earlier versions of Dfam, please hash out lines 130-132 in the Dockerfile

Run the following from inside the EarlGrey/Docker/ directory

## Download dependencies

First, run the getFiles.sh script to download the required packages:

```
chmod +x getFiles.sh
./getFiles.sh
```

NOTE: This script will download the latest curated TEs from Dfam. If you would rather download the whole dataset, including uncurated TEs, please see below. 
We strongly recommend upgrading the RepeatMasker contained within Docker with the Dfam library v3.7 curated elements ONLY (which is sourced by `getFiles.sh`).

If you need to download these,
Dfam.h5.gz (NOTE ONLY TESTED ON DFAM 3.7) can be downloaded by running:

```
# all elements including uncurated
cd ./src/
wget https://www.dfam.org/releases/Dfam_3.7/families/Dfam.h5.gz 
# THIS IS A BIG FILE!
```

Repbase is now behind a paywall, The following lines are hashed out due to incompatabilities with the latest Dfam (lines 105-107). You are free to reconfigure the container if you are confident in configuring these programs.
```
# && cd /opt/RepeatMasker \
# && cp /opt/src/RepBaseRepeatMaskerEdition-20181026.tar.gz . \
# && tar -zxf RepBaseRepeatMaskerEdition-20181026.tar.gz \
```

## Build a docker container (run from inside the directory where the Dockerfile for EarlGrey is stored)
```
docker build -t earlgrey --build-arg USER_ID=$(id -u) --build-arg GROUP_ID=$(id -g) .
```

## start the docker container

```
docker run -it --init --mount type=bind,source="$(pwd)",target=/work --user "$(id -u):$(id -g)" --workdir "/work" --env "HOME=/work" earlgrey "$@"
```

## configure the earlGrey conda environment in the docker container (ie after starting the container)
```
eval "$(/anaconda3/bin/conda shell.bash  hook)"
conda env create -f /home/user/EarlGrey/earlGrey.yml
conda activate earlGrey
```
## Run these commands to activate the conda environment for EarlGrey when starting the container

```
eval "$(/anaconda3/bin/conda shell.bash  hook)"
conda activate earlGrey
```


