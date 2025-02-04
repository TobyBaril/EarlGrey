# Get Earl Grey to run in Docker

## A [biocontainer](https://biocontainers.pro/tools/earlgrey) has been generated from the Earl Grey Bioconda package.

## Get the preconfigured docker container from biocontainers
## This image will be configured with Dfam curated elements only
In this case, we need to bind a system directory to the docker container. In the line below, we are binding a directory call `host_data` that is found on our current path to `/data/` in the docker container. Please replace the file path before `:` to the directory you wish to bind to `/data/` in the container. This container must be run in interactive mode the first time you use it.

```
docker run -it -v `pwd`/host_data/:/data/ quay.io/biocontainers/earlgrey:5.1.1--h9948957_0
```

## If you are running the container for the first time, you need to enable Earl Grey to configure the Dfam libraries correctly in interactive mode.
This is done by running Earl Grey for the first time. You will be prompted if Dfam libraries have not been configured properly. Answer `y` to configure Dfam.
Ensure the genome you wish to annotate is found in the directory you want to bind to the container. It will then be found in `/data` in the container. Alternatively, if you source a genome assembly from a database, download it to `/data/`. Always output results to the bound directory. In our case `/data/`. The results will then be found in your system directory after exiting the container.

```
earlGrey -g /data/genome.fasta -s test_genome -t 8 -o /data/
```

## If you need the container to run offline and/or without interactive mode
I try to keep an up-to-date container in docker hub, but this might not always be the case depending on if I have had time to build and upload a new image. Currently, there is an image with Dfam 3.7 curated elements only, and this is version 5.1.1. You can use this image by pulling the container:

```
# Interactive mode
docker run -it -v `pwd`:/data/ tobybaril/earlgrey_dfam3.7:latest

# Non interactive mode example:
docker run -v `pwd`:/data/ tobybaril/earlgrey_dfam3.7:latest earlGrey -g /data/NC_045808_EarlWorkshop.fasta -s nonInteractiveTest -o /data/ -t 8
```

