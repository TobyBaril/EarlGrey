# Singularity image for Earl Grey

A Singularity image can be built from the Docker image hosted on Docker Hub. The image contains no pre-configured Dfam partitions, but will generate a configuration script on first run of `earlGrey`. You will need to configure the Dfam libraries once, then you can use the configured sandbox for all subsequent analyses.

Wherever you start this is where your `/data` directory will point to inside the container. Make sure your genome assembly is found here, or change `$(pwd)` in the below commands to the directory where you want to perform your analysis.

**Write access is required inside the container to download and configure the Dfam partitions.** The recommended approach is to build a writable sandbox, configure Dfam, and then optionally convert the sandbox to a portable SIF file.

## First-time setup: building and configuring the sandbox

```bash
# Build a writable sandbox from the Docker image
singularity build --sandbox earlgrey_sandbox/ docker://tobybaril/earlgrey:latest-nodfam

# Run the sandbox interactively with write access, binding your working directory to /data
singularity shell -C --bind $(pwd):/data --writable earlgrey_sandbox/

# Inside the container: change to the data directory and run earlGrey to generate the configuration script
cd /data/
earlGrey

# Modify the configuration script to select the required Dfam 3.9 partitions
# (https://www.dfam.org/releases/Dfam_3.9/families/FamDB/)
# At minimum, partition 0 is required. Change 0-16 to whichever partitions you need.
##### e.g. for partitions 0-5:
sed -i '/^curl/ s/0-16/0-5/g' configure_dfam39.sh
##### e.g. for specific partitions 1, 3, and 5:
sed -i '/^curl/ s/0-16/1,3,5/g' configure_dfam39.sh

# Run the configuration script to download and install the selected Dfam partitions
bash configure_dfam39.sh

# Return to your data directory
cd /data/

# Exit the sandbox
exit
```

## Running Earl Grey after configuration

```bash
# Run the configured sandbox interactively
singularity shell -C --bind $(pwd):/data --writable earlgrey_sandbox/

# Or run non-interactively
singularity exec -C --bind $(pwd):/data --writable earlgrey_sandbox/ earlGrey -g /data/GENOME.fasta -s mySpecies -o /data/ -t 8
```

## Optional: convert the configured sandbox to a portable SIF

Once Dfam is configured, you can convert the sandbox to a read-only SIF file for easier distribution and deployment. Note that the resulting SIF will be large depending on which Dfam partitions were downloaded.

```bash
singularity build earlgrey_configured.sif earlgrey_sandbox/

# Run non-interactively with the configured SIF
singularity exec -C --bind $(pwd):/data earlgrey_configured.sif earlGrey -g /data/GENOME.fasta -s mySpecies -o /data/ -t 8
```
