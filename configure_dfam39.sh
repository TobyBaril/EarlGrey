# first, change directory to the famdb library location
cd /data/toby/miniforge3/envs/earlgrey-pan-dev/share/RepeatMasker/Libraries/famdb/

# download the partitions you require from Dfam 3.9. In the below, change the numbers or range inside the square brackets to choose your subsets.
# e.g. to download partitions 0 to 10: [0-10]; or to download partitions 3,5, and 7: [3,5,7]; [0-16] is ALL PARTITIONS
curl -o 'dfam39_full.#1.h5.gz' 'https://dfam.org/releases/current/families/FamDB/dfam39_full.[0].h5.gz'

# decompress Dfam 3.9 paritions
gunzip *.gz

# move up to RepeatMasker main directory
cd /data/toby/miniforge3/envs/earlgrey-pan-dev/share/RepeatMasker/

# save the min_init partition as a backup, just in case!
mv /data/toby/miniforge3/envs/earlgrey-pan-dev/share/RepeatMasker/Libraries/famdb/min_init.0.h5 /data/toby/miniforge3/envs/earlgrey-pan-dev/share/RepeatMasker/Libraries/famdb/min_init.0.h5.bak

# Rerun RepeatMasker configuration
perl ./configure -libdir /data/toby/miniforge3/envs/earlgrey-pan-dev/share/RepeatMasker/Libraries/ -trf_prgm /data/toby/miniforge3/envs/earlgrey-pan-dev/bin/trf -rmblast_dir /data/toby/miniforge3/envs/earlgrey-pan-dev/bin -hmmer_dir /data/toby/miniforge3/envs/earlgrey-pan-dev/bin -abblast_dir /data/toby/miniforge3/envs/earlgrey-pan-dev/bin -crossmatch_dir /data/toby/miniforge3/envs/earlgrey-pan-dev/bin -default_search_engine rmblast

