# Get Earl Grey to run in Docker

Run the following from inside the EarlGrey/Docker/ directory

## Download dependencies

First, run the getFiles.sh script to download the required packages:

```
chmod +x getFiles.sh
./getFiles.sh
```

NOTE: please make sure that Dfam.h5.gz and RepBase libraries tar.gz are inside the ./src directory after running ./getFiles.sh . 
we strongly recommend upgrading the RepeatMasker contained within Docker with the Dfam libraries and RepBase, information on this is found: https://github.com/Dfam-consortium/TETools
note there are only ~6000 seqs in the basic RepeatMasker library in this container, real Dfam is much larger!

If you need to download these,
Dfam.h5.gz (NOTE ONLY TESTED UP TO DFAM 3.6) can be downloaded by running:

```
wget https://www.dfam.org/releases/Dfam_3.6/families/Dfam.h5.gz 
# THIS IS A BIG FILE!
```

Repbase is now behind a paywall, if you do not have access, please comment out the following lines in the Dockerfile (lines 105-107)
```
&& cd /opt/RepeatMasker \
&& cp /opt/src/RepBaseRepeatMaskerEdition-20181026.tar.gz . \
&& tar -zxf RepBaseRepeatMaskerEdition-20181026.tar.gz \
```

## Build a docker container (run from inside the directory where the Dockerfile for EarlGrey is stored)
```
docker build . -t earlgrey
```

## start the docker container

```
docker run -it --init --mount type=bind,source="$(pwd)",target=/work --user "$(id -u):$(id -g)" --workdir "/work" --env "HOME=/work" earlgrey "$@"
```

## configure the earlGrey conda environment in the docker container (ie after starting the container)
```
eval "$(/anaconda3/bin/conda shell.bash  hook)"
conda env create -f /opt/EarlGrey/earlGrey.yml
conda activate earlGrey
Rscript /opt/EarlGrey/scripts/install_r_packages.R
```
## Run these commands to activate the conda environment for EarlGrey when starting the container

```
eval "$(/anaconda3/bin/conda shell.bash  hook)"
conda activate earlGrey
```
