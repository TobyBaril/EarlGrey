# Singularity image for Earl Grey

These instructions will help you to generate the singularity image for Earl Grey preconfigured with Dfam 3.7 curated elements.

Wherever you start this is where your /work directory will point to. Make sure your genome assembly is found here, or change `$(pwd)` in the below command to the directory where you want to perform your analysis

To start:
```
# Build the image from the Docker image
singularity build earlgrey.sif docker://tobybaril/earlgrey:latest-nodfam

# Run the sandbox
singularity shell -C -H $(pwd):/work --writable-tmpfs -u earlgrey.sif

# Don't forget the first time you run earlGrey, it will guide you through getting the Dfam partitions! See main README Docker section for more info!
```

You are now ready to run Earl Grey!
