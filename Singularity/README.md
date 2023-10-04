# Singularity image for Earl Grey

This directory contains the singularity image for Earl Grey preconfigured with Dfam 3.7 curated elements.

Wherever you start this is where your /work directory will point to. Make sure your genome assembly is found here, or change `$(pwd)` in the below command to the directory where you want to perform your analysis

To start:
```
singularity shell -C -H $(pwd):/work --writable-tmpfs -u earlgrey.sif
```

Then run the following before starting Earl Grey:

```
eval "$(/anaconda3/bin/conda shell.bash  hook)"
conda env create -f /home/user/EarlGrey/earlGrey.yml
conda activate earlGrey
Rscript /home/user/EarlGrey/scripts/install_r_packages.R
```

You are now ready to run Earl Grey!
