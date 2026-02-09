# Docker Container 

A Docker container has been generated with none of Dfam 3.9, but with script generation to source required partitions

I try to keep an up-to-date container in docker hub, but this might not always be the case depending on if I have had time to build and upload a new image. Currently, the recommended image ready for use is `-nodfam` version. Upon running the container interactively and running the command `earlGrey`, instructions will print to `stdout` and a script that you can use will be placed in `/usr/local/share/RepeatMasker/Libraries/famdb/` when the container is running.

```
# Interactive mode
# Version 7.0.1 with no preconfigured partitions (RECOMMENDED!) - bind a directory, in my case the current directory using pwd
docker run -it -v 'pwd':/data/ tobybaril/earlgrey:latest-nodfam
# change to library directory
cd /data/
# run earlGrey to make the configuration script
earlGrey

# then alter script with required partitions and run the configuration script
# change 0-16 to whichever you require, but at least 0. This relates to the partitions of Dfam 3.9 (https://www.dfam.org/releases/Dfam_3.9/families/FamDB/)
##### e.g. for 0-5:
sed -i '/^curl/ s/0-16/0-5/g' configure_dfam39.sh
##### e.g for 1,3,5:
sed -i '/^curl/ s/0-16/1,3,5/g' configure_dfam39.sh

# run the configuration script
bash configure_dfam39.sh

# return to your data directory
cd /data/
```

