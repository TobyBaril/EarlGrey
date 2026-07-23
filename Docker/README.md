# Docker Container 

A Docker container is provided without pre-downloaded Dfam partitions.

I try to keep an up-to-date container in Docker Hub, but uploads may lag behind code releases. Currently, the recommended image ready for use is the `-nodfam` version. When running `earlGrey` interactively, instructions are printed to `stdout` and a helper script (`configure_dfam40.sh`) is generated in your working directory.

```
# Interactive mode
# Version 7.3.0 with no preconfigured partitions (RECOMMENDED!) - bind a directory, in my case the current directory using pwd
docker run -it -v 'pwd':/data/ tobybaril/earlgrey:latest-nodfam
# change to library directory
cd /data/
# run earlGrey to make the configuration script
earlGrey

# run the configuration script (uses interactive FamDB downloader)
bash configure_dfam40.sh

# return to your data directory
cd /data/
```

