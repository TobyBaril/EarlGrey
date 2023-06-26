# Get Earl Grey to run in Docker

## An image for Earl Grey configured with Dfam 3.7 curated elements can be found [here](https://hub.docker.com/repository/docker/tobybaril/earlgrey/general)

# Get the preconfigured docker container from docker hub
## This image has been configured with Dfam v3.7 curated elements only

```
docker run -it --init --mount type=bind,source="$(pwd)",target=/work --user "$(id -u):$(id -g)" --workdir "/work" --env "HOME=/work" tobybaril/earlgrey "$@"
```

## configure the earlGrey conda environment in the docker container (ie after starting the container)
```
eval "$(/anaconda3/bin/conda shell.bash  hook)"
conda env create -f /home/user/EarlGrey/earlGrey.yml
conda activate earlGrey
Rscript /home/user/EarlGrey/scripts/install_r_packages.R
```
## Run these commands to activate the conda environment for EarlGrey when starting the container

```
eval "$(/anaconda3/bin/conda shell.bash  hook)"
conda activate earlGrey
```

# OR Build the container...

## Note, the docker container has been tested with Dfam 3.7 without RepBase. If you would like to use with earlier versions of Dfam, please hash out lines 130-132 in the Dockerfile

Run the following from inside the EarlGrey/Docker/ directory

## Download dependencies

First, run the getFiles.sh script to download the required packages:

```
chmod +x getFiles.sh
./getFiles.sh
```

NOTE: This script will download the latest curated TEs from Dfam. If you would rather download the whole dataset, including uncurated TEs, please see below. 
We strongly recommend upgrading the RepeatMasker contained within Docker with the Dfam library v3.7 curated elements ONLY, information on this is found: https://github.com/Dfam-consortium/TETools
note there are only ~6000 seqs in the basic RepeatMasker library in this container, real Dfam is much larger!

If you need to download these,
Dfam.h5.gz (NOTE ONLY TESTED ON DFAM 3.7) can be downloaded by running:

```
# only curated elements
wget https://www.dfam.org/releases/current/families/Dfam_curatedonly.h5.gz && mv Dfam_curatedonly.h5.gz Dfam.h5.gz

# all elements including uncurated
wget https://www.dfam.org/releases/current/families/Dfam.h5.gz 
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
Rscript /home/user/EarlGrey/scripts/install_r_packages.R
```
## Run these commands to activate the conda environment for EarlGrey when starting the container

```
eval "$(/anaconda3/bin/conda shell.bash  hook)"
conda activate earlGrey
```


