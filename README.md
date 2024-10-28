![earlGreyIcon](https://user-images.githubusercontent.com/46785187/136248346-21e980ee-1154-48c2-9398-70938bbe2404.png)

[![DOI](https://zenodo.org/badge/412126708.svg)](https://zenodo.org/badge/latestdoi/412126708) [![Anaconda_version](https://anaconda.org/bioconda/earlgrey/badges/version.svg)](https://anaconda.org/bioconda/earlgrey) [![Anaconda_downloads](https://anaconda.org/bioconda/earlgrey/badges/downloads.svg)](
https://anaconda.org/bioconda/earlgrey) [![Anaconda_platforms](https://anaconda.org/bioconda/earlgrey/badges/platforms.svg)](https://anaconda.org/bioconda/earlgrey) ![Twitter URL](https://img.shields.io/twitter/url?style=social&url=https%3A%2F%2Fgithub.com%2FTobyBaril%2FEarlGrey%2F)

# Earl Grey

Earl Grey is a full-automated transposable element (TE) annotation pipeline, leveraging the most widely-used tools and combining these with a consensus elongation process to better define _de novo_ consensus sequences when annotating new genome assemblies.

# Contents

[Changes in Latest Release](#changes-in-latest-release)

[Example Run](#example)

[References and Acknowledgements](#references-and-acknowledgements)

[Usage Without Installation](#usage-without-installation)

[Recommended Installation](#recommended-installation-with-conda-or-mamba)

[Workaround for Dfam 3.8](#dfam-workaround-for-latest-release)

[Alternative Installation Methods](#alternative-installation-methods)

<!-- toc -->

# Changes in Latest Release
Big changes in the latest release!

*Earl Grey v5.0.0 is here!* 

This release incorporates the incremental improvements made throughout the life of version 4. 

It is now possible to run some subroutines in Earl Grey (run either of these new commands with `-h` to see a list of options):
- `earlGreyLibConstruct` can be used to run Earl Grey for _de novo_ TE detection, consensus generation, and improvement through the BEAT process. The output will be the strained TE consensus sequences, which can then be used for subsequent annotation. This is useful when you want to make a combined library from the libraries of several different genomes, where it is no longer required to waste time running annotations. Once the libraries are generated and you have curated them, you can then run the next step in isolation (next point!).
- `earlGreyAnnotationOnly` can be used to run the final annotation and defragmentation steps in Earl Grey. This is useful if you have already run the BEAT process and have a library of _de novo_ TE consensus sequences that you would like to use to annotate a given genome. This script is also compatible with the `-r` flag to take known repeats from the databases used to configure RepeatMasker in addition to the custom repeat library.
- *EXPERIMENTAL FEATURE:* I have also added an option to run [HELIANO](https://github.com/Zhenlisme/heliano) for improved detection of Helitrons, which are notoriously difficult to detect and classify using homology methods. This can be implemented by adding `-e yes` to the command line options after upgrading to v5.0.0. Currently, HELIANO annotations replace those which they overlap following the RepeatMasker run, which is performed during defragmentation (in a similar way to full-length LTRs being dealt with in `RepeatCraft`). Feedback is welcomed on this implementation, and I am continuing to test and improve the implementation of HELIANO within Earl Grey.
- The settings used for HELIANO are: `--nearest -dn 6000 -flank_sim 0.5 -w 10000`. These can be modified in the earlGrey script of your specific installation.

Thank you for your continued support and enthusiasm for Earl Grey!

# Example

Given an input genome, Earl Grey will run through numerous steps to identify, curate, and annotate transposable elements (TEs). We recommend running earlGrey within a tmux or screen session, so that you can log off and leave Earl Grey running.

There are several required and optional parameters for your Earl Grey run:
```
Required Parameters:
		-g == genome.fasta
		-s == species name
		-o == output directory

	Optional Parameters:
		-t == Number of Threads (DO NOT specify more than are available)
		-r == RepeatMasker search term (e.g arthropoda/eukarya)
		-l == Starting consensus library for an inital mask (in fasta format)
		-i == Number of Iterations to BLAST, Extract, Extend (Default: 10)
		-f == Number flanking basepairs to extract (Default: 1000)
		-c == Cluster TE library to reduce redundancy? (yes/no, Default: no)
		-m == Remove putative spurious TE annotations <100bp? (yes/no, Default: no)
		-d == Create soft-masked genome at the end? (yes/no, Default: no)
		-n == Max number of sequences used to generate consensus sequences (Default: 20)
		-a == minimum number of sequences required to build a consensus sequence (Default: 3)
                -e == Optional: Run HELIANO for detection of Helitrons (yes/no, Default: no)
		-h == Show help
```

```
# remember to activate the conda environment before running (NOTE: depending on your install, the environment name might vary)

conda activate earlgrey

# run earl grey with minimum command options

earlGrey -g [genome.fasta] -s [speciesName] -o [outputDirectory]

# e.g

earlGrey -g myzusPersicae.fasta -s myzusPersicae -o ./earlGreyOutputs

```

Following this, Earl Grey will run through several processes depending on the options selected. The general pipeline is illustrated below:

![Figure1_earlGreyFlow](https://github.com/TobyBaril/EarlGrey/assets/46785187/ba039063-9111-4264-84b6-e3580dc340ae)

For a more in-depth description of Earl Grey's steps, please refer to the implementation section in the [manuscript](https://academic.oup.com/mbe/article/41/4/msae068/7635926).

The runtime of Earl Grey will depend on the repeat content of your input genome. Once finished, you will notice that a number of directories have been created by Earl Grey. The most important results are found within the "summaryFiles" folder, however intermediate results are kept in case you wish to use alignments for further manual curation or investigation, for example. NOTE: RepeatModeler2 remains the rate-limiting step, with runtimes leading into several days, or even weeks, with large repeat-rich genomes. This is normal. 

Directories created by earl grey:
```
[speciesName]EarlGrey/
    |
    |---[speciesName]_RepeatMasker/
    |		+ Results of the optional initial RepeatMasker run used to mask previously characterised TEs
    |---[speciesName]_Database/
    |		+ Database created from the masked genome output of the initial RepeatMasker run. Required for RepeatModeler.
    |---[speciesName]_RepeatModeler/
    |		 + Results of the RepeatModeler2 _de novo_ TE identification step
    |---[speciesName]_strainer/
    |		+ Results of the "BLAST, Extract, Align, Trim" process
    |---[speciesName]_Curated_Library/
    |		+ Contains the _de novo_ repeat library generated with the "BLAST, Extract, Align, Trim" process, the library of known repeats used by RepeatMasker (OPTIONAL), and a combined library containing both sets of repeats (OPTIONAL)
    |---[speciesName]_RepeatMasker_Against_Custom_Library/
    |		+ Results of the RepeatMasker run using the final curated library
    |---[speciesName]_RepeatLandscape/
    |		+ Intermediate files for the generation of RepeatLandscapes (RepeatMasker .divsum files)
    |---[speciesName]_mergedRepeats/
    |		+ Intermediate files and results of TE defragmentation step using RepeatCraft.
    |---[speciesName]_summaryFiles/
    		+ Results and plots from Earl Grey:
    		+ TE annotations in GFF3 and BED format
    		+ High-level TE quantification table (tab delimited)
    		+ Family-level TE quantification table (tab delimited)
    		+ Repeat Landscapes showing TE activity (PDF)
    		+ Pie chart of genome repeat content (PDF)
    		+ _de novo_ repeat library in FASTA format
    		+ Combined repeat library in FASTA format (OPTIONAL)
```

As of `v4.4.5`, there is an option to generate _de novo_ TE libraries without running subsequent annotation. To run this option, use `earlGreyLibConstruct` with the same command-line options as `earlGrey`. This will run everything up to the end of TEstrainer, and leave you a `families-fa.strained` library file in the `summaryFiles` directory, which you can then use for manual curation, or for pangenome studies.

### Example Outputs (NOTE: example data has been used here):

- Pie chart summarising TE content in input genome

![image](https://user-images.githubusercontent.com/46785187/140897482-fc80c30a-3e5f-4bf6-99b0-0f18b11b33d1.png)

- RepeatLandscapes summarising relative TE activity using Kimura 2-Parameter Divergence (recent activity towards the RHS)

<img width="849" alt="Screenshot 2024-08-12 at 13 55 06" src="https://github.com/user-attachments/assets/5e1b18b6-a84a-47e9-bc6d-6f341a2297ec">

- TE annotations - These are in standard genomic information formats to be compatible with downstream analyses. NOTE: TE divergence calculated using Kimura 2-Parameter distance is now supplied for each insertion in column 9 of the GFF3 file:

```
# BED format
NC45808.1     4964941 4965925 LINE/Penelope   5073    +
NC45808.1     7291353 7291525 LINE/L2 1279    +
NC_045808.1     8922477 8923791 DNA/TcMar-Tc1   11957   +


# GFF3 format
scaffold_1	Earl_Grey	DNA/Mariner	71618	71814	892	+	.	TSTART=13;TEND=233;ID=RND-1_FAMILY-1;SHORTTE=F;KIMURA80=0.2141
scaffold_1	Earl_Grey	Unknown	81757	81927	785	-	.	TSTART=16;TEND=194;ID=RND-1_FAMILY-0;SHORTTE=F;KIMURA80=0.2037
```

# References and Acknowledgements

This pipeline has been designed to be used and shared openly by the community.

### When using Earl Grey, please cite:

Baril, T., Galbraith, J.G., and Hayward, A., Earl Grey: A Fully Automated User-Friendly Transposable Element Annotation and Analysis Pipeline, Molecular Biology and Evolution, Volume 41, Issue 4, April 2024, msae068 [doi:10.1093/molbev/msae068](https://doi.org/10.1093/molbev/msae068)

Baril, Tobias., Galbraith, James., and Hayward, Alexander. (2023) Earl Grey. Zenodo [doi:10.5281/zenodo.5654615](https://doi.org/10.5281/zenodo.5654615)

### This pipeline makes use of scripts from:

[RepeatCraft](https://github.com/niccw/repeatcraftp) - Wong WY, Simakov O. RepeatCraft: a meta-pipeline for repetitive element de-fragmentation and annotation. Bioinformatics 2018;35:1051–2. https://doi.org/10.1093/bioinformatics/bty745.

### The following open source software are utilised by this pipeline:

Smit AFA, Hubley RR, Green PR. RepeatMasker Open-4.0. Http://RepeatmaskerOrg 2013.

Flynn, J.M., Hubley, R., Goubert, C., Rosen, J., Clark, A.G., Feschotte, C. and Smit, A.F. RepeatModeler2 for automated genomic discovery of transposable element families. Proceedings of the National Academy of Sciences 2020;17;9451-9457. https://doi.org/10.1073/pnas.1921046117

Bao Z, Eddy SR. Automated De Novo Identification of Repeat Sequence Families in Sequenced Genomes. Genome Res 2002;12:1269–76. https://doi.org/10.1101/gr.88502.

Price AL, Jones NC, Pevzner PA. De novo identification of repeat families in large genomes. Bioinformatics 2005;21:i351–8.

Camacho C, Coulouris G, Avagyan V, Ma N, Papadopoulos J, Bealer K, et al. BLAST+: Architecture and applications. BMC Bioinformatics 2009;10: 1–9. https://doi.org/10.1186/1471-2105-10-421.

Katoh K, Standley DM. MAFFT multiple sequence alignment software version 7: Improvements in performance and usability. Mol Biol Evol 2013;30:772–80. https://doi.org/10.1093/molbev/mst010.

Capella-Gutiérrez S, Silla-Martínez JM, Gabaldón T. trimAl: A tool for automated alignment trimming in large-scale phylogenetic analyses. Bioinformatics 2009;25:1972–3. https://doi.org/10.1093/bioinformatics/btp348.

Rice P, Longden L, Bleasby A. EMBOSS: The European Molecular Biology Open Software Suite. Trends Genet 2000;16:276–7. https://doi.org/10.1016/S0168-9525(00)02024-2.

Xu Z, Wang H. LTR_FINDER: an efficient tool for the prediction of full-length LTR retrotransposons. Nucleic Acids Res 2007;35:W265–8. https://doi.org/10.1093/nar/gkm286.

Ou S, Jiang N. LTR_FINDER_parallel: parallelization of LTR_FINDER enabling rapid identification of long terminal repeat retrotransposons. BioRxiv 2019:2–6.

# Usage without installation
If you would like to try Earl Grey, or prefer to use it in a browser, you can do this through [gitpod](https://gitpod.io). You can get 50 hours free per month and use Earl Grey within a preconfigured environment. Simply select this repository by pasting the repository URL, and the system will be automatically configured for you. You can then upload your genome of interest, and run Earl Grey as you would on the command line. 

<img width="1919" alt="Screenshot 2023-09-29 at 13 38 43" src="https://github.com/TobyBaril/EarlGrey/assets/46785187/7dd2f2de-3c13-4553-b13a-007fdd8d94d6">

# Recommended Installation with Conda or Mamba

NOTE: This pipeline is currently running with Dfam 3.7 curated elements only. We are working on updating to Dfam 3.8 for a future release. If required, you can modify the conda installation of RepeatMasker within the conda environment (do at your own risk!)

Earl Grey version 5.0.0 (latest stable release) with all required and configured dependencies is found in the `toby_baril_bio` and `biooconda` conda channels. To install, simply run the following depending on your installation:
```
# With conda
conda create -n earlgrey -c conda-forge -c bioconda earlgrey=5.0.0

# With mamba
mamba create -n earlgrey -c conda-forge -c bioconda earlgrey=5.0.0
```

# Recommended Installation with Conda or Mamba on ARM-based Mac Systems (M chips)
## NOTE: This is currently experimental and requires altering several scripts and parameters. Not recommended unless you are confident with the terminal environment.

It is possible to install x86-based conda environments on M-chip Mac systems using Rosetta. The below form a guide that should work, but please reach out if you have any trouble!

First, we will run separate installations of conda for ARM and x86 architectures.

We want a simple way to activate the x86 architecture. We can do this by adding the following to `~/.zshrc`:
```
alias arm="env /usr/bin/arch -arm64 /bin/zsh --login"
alias intel="env /usr/bin/arch -x86_64 /bin/zsh --login"
```

Then close the terminal and open a new one.

To activate the intel environment, run the following in a new terminal:
```
intel
```

Next, we need to install the separate conda environments. This guide is based on [this article:](https://taylorreiter.github.io/2022-04-05-Managing-multiple-architecture-specific-installations-of-conda-on-apple-M1/)

First, check you are using the arm64 processor:
```
uname -m
```

If this returns `arm64`, all good. If not, run:
```
arm
```

Then, install miniforge for arm64:
```
curl -L https://github.com/conda-forge/miniforge/releases/download/23.3.1-1/Mambaforge-23.3.1-1-MacOSX-arm64.sh > miniforge_arm64.sh
sh miniforge_arm64.sh
```

Follow the prompts to accept the license, install, and initialise conda. The initialisation script, `conda init` will add something like the following block to your `~/.zshrc` file:
```
# >>> conda initialize >>>
# !! Contents within this block are managed by 'conda init' !!
__conda_setup="$('/Users/toby/mambaforge/bin/conda' 'shell.zsh' 'hook' 2> /dev/null)"
if [ $? -eq 0 ]; then
    eval "$__conda_setup"
else
    if [ -f "/Users/toby/mambaforge/etc/profile.d/conda.sh" ]; then
        . "/Users/toby/mambaforge/etc/profile.d/conda.sh"
    else
        export PATH="/Users/toby/mambaforge/bin:$PATH"
    fi
fi
unset __conda_setup

if [ -f "/Users/toby/mambaforge/etc/profile.d/mamba.sh" ]; then
    . "/Users/toby/mambaforge/etc/profile.d/mamba.sh"
fi
# <<< conda initialize <<<
```

We want to _cut_ this from the `~/.zshrc` file and place it in a new file:
```
# open a new file called ~/.start_miniforge3.sh and paste the text you cut from ~/.zshrc
nano ~/.start_miniforge3.sh
```

Next, we will install miniconda in the Rosetta, or _intel_, terminal.

First, check you are using the x86_64 processor:
```
uname -m
```

If this doesn't return `x86_64`, run:
```
intel
```

Then, install miniconda for x86_64:
```
curl -L https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-x86_64.sh > miniconda_x86_64.sh
sh miniconda_x86_64.sh
```

Again, you should follow the prompts to accept the license, install, and initialise conda. The initialisation script, `conda init` will add something like the following block to your `~/.zshrc` file:
```
# >>> conda initialize >>>
# !! Contents within this block are managed by 'conda init' !!
__conda_setup="$('/Users/toby/miniconda3/bin/conda' 'shell.zsh' 'hook' 2> /dev/null)"
if [ $? -eq 0 ]; then
    eval "$__conda_setup"
else
    if [ -f "/Users/toby/miniconda3/etc/profile.d/conda.sh" ]; then
        . "/Users/toby/miniconda3/etc/profile.d/conda.sh"
    else
        export PATH="/Users/toby/miniconda3/bin:$PATH"
    fi
fi
unset __conda_setup
# <<< conda initialize <<<
```

Again, we want to cut this from the `~/.zshrc` file and place it in a new file:
```
# open a new file called ~/.start_miniconda3.sh and paste the text you cut from ~/.zshrc
nano ~/.start_miniconda3.sh
```

Now we want to edit `~/.zshrc` to automatically initialise the correct conda for each architecture. To do this, add the following to your `~/.zshrc file`:
```
arch_name="$(uname -m)"

if [ "$arch_name" = "arm64" ]; then
    echo "Running on ARM64 using miniforge3"
    source ~/.start_miniforge3.sh
elif [ "$arch_name" = "x86_64" ]; then
    echo "Running on Rosetta using miniconda3"
    source ~/.start_miniconda3.sh
else
    echo "Unknown architecture: $arch_name"
fi
```

Exit all terminals. When you reopen them, it should all be ready.

Now, open a new terminal and get Earl Grey running.

First, activate the intel environment:
```
intel
```

Then, create an environment for Earl Grey
```
# I like mamba. This is optional but good
conda install -c conda-forge mamba

# create the environment with correct subdirectory
CONDA_SUBDIR=osx-64 conda create -n earlgrey

# activate conda environment
conda activate earlgrey

# install the packages - CHECK THIS IS USING THE MINICONDA3 MAMBA EXECUTABLE
mamba install -c bioconda earlgrey
```

Now we need to change some files to get them to behave with zsh!

Make sure you are in your intel environment!

Install a couple of utilities using homebrew.
```
# Make sure you are in the intel environemt and that the earlgrey conda environment is active

# install homebrew for the intel environment (if you don't already have it, it must be separate to the ARM installation!
/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"

# add the following to your ~/.zshrc
alias brewIntel="/usr/local/bin/brew"

# source ~/.zshrc
source ~/.zshrc
intel

brewIntel install gnu-sed
brewIntel install coreutils
```

Change TEstrainer_for_earlGrey.sh for the macOS version:
```
nano $(which earlGrey | gsed 's|bin.*|share/earlgrey-5.0.0-0/scripts/TEstrainer/TEstrainer_for_earlGrey.sh|g')

# delete everything in this file.
```

Then paste ALL of the following:
```
#!/bin/bash

usage() { echo "Usage: [-l Repeat library] [-g Genome ] [-t Threads (default 4) ] [-f Flank (default 1000) ] [-m Minimum number of sequences required in multiple sequence alignments during BEAT] [-r Number of iterations of BEAT to run (deafult 10)] [-d Out directory, if not specified wil be created ] [-h Print this help] [-M Ammount of memory TEstrainer needs to keep free]" 1>&2; exit 1; }


# default values

STRAIN_SCRIPTS=INSERT_FILENAME_HERE
FLANK=1000
THREADS=4
RUNS=10
NO_SEQ=20
# for potential folder name
TIME=$(date +"%s")
TIME=${TIME: -4}
MEM_FREE="200M"
MIN_SEQ=3

# parsing
while getopts l:g:t:f:m:r:d:h:n:M flag; do
  case "${flag}" in
    l) RM_LIBRARY_PATH=${OPTARG};;
    g) GENOME=${OPTARG};;
    t) THREADS=${OPTARG};;
    f) FLANK=${OPTARG};;
    m) MIN_SEQ=${OPTARG};;
    r) RUNS=${OPTARG};;
    d) DATA_DIR=${OPTARG};;
    n) NO_SEQ=${OPTARG};;
    M) MEM_FREE=${OPTARG};;
    h | *)
      print_usage
      exit_script
  esac
done

# determine further variables and check files exist
if [ ! -f ${RM_LIBRARY_PATH} ]; then echo "Library not found"; usage; fi
if [ -z ${RM_LIBRARY_PATH} ]; then echo "Library must be supplied"; usage; else RM_LIBRARY=$(echo $RM_LIBRARY_PATH | sed 's/.*\///'); fi
if [[ $RUNS -gt 0 ]]; then 
  if [ -z ${GENOME} ]; then echo "If refining genome must be supplied"; usage; fi
  if [ ! -f ${GENOME} ]; then echo "Refining genome not found"; usage; fi
fi
if [ -z "$DATA_DIR" ]; then DATA_DIR=$(echo "TS_"${RM_LIBRARY}"_"${TIME}); fi

# create data dir if missing
if [ ! -d ${DATA_DIR} ] 
then
    mkdir ${DATA_DIR}
fi


if [[ $THREADS -gt 4 ]]; then MAFFT_THREADS=$(($(($THREADS / 4)))); else MAFFT_THREADS=1; fi

# make directories
mkdir -p ${DATA_DIR}/run_0/

# initial copy
cp ${RM_LIBRARY_PATH} ${DATA_DIR}/${RM_LIBRARY}

# cp starting seq to starting directory
mkdir -p ${DATA_DIR}/run_0/og
cp ${RM_LIBRARY_PATH} ${DATA_DIR}/run_0/further_${RM_LIBRARY}
# create reference of original sequences
python ${STRAIN_SCRIPTS}/splitter.py -i ${DATA_DIR}/run_0/further_${RM_LIBRARY} -o ${DATA_DIR}/run_0/og

# runs
python ${STRAIN_SCRIPTS}/indexer.py -g ${GENOME}
if [ ! -f "${GENOME}".nsq ]; then
  makeblastdb -in ${GENOME} -dbtype nucl -out ${GENOME} # makeblastb if needed
fi

# curation
RUN_NO=1
while  [ $RUN_NO -le $RUNS ]
do

  # make directories
  mkdir -p ${DATA_DIR}/run_${RUN_NO}/raw \
           ${DATA_DIR}/run_${RUN_NO}/TEtrim_complete \
           ${DATA_DIR}/run_${RUN_NO}/initial_blast \
           ${DATA_DIR}/run_${RUN_NO}/self_search \
           ${DATA_DIR}/run_${RUN_NO}/to_align \
           ${DATA_DIR}/run_${RUN_NO}/mafft \
           ${DATA_DIR}/run_${RUN_NO}/TEtrim_con \
           ${DATA_DIR}/run_${RUN_NO}/TEtrim_unaln \
           ${DATA_DIR}/run_${RUN_NO}/TEtrim_blast \
           ${DATA_DIR}/run_${RUN_NO}/TEtrim_mafft \
           ${DATA_DIR}/run_${RUN_NO}/TEtrim_further \
           ${DATA_DIR}/run_${RUN_NO}/TEtrim_bp
  
  # split
  cp ${DATA_DIR}/run_$(expr $RUN_NO - 1)/further_${RM_LIBRARY} ${DATA_DIR}/run_${RUN_NO}/${RM_LIBRARY}
  echo "Splitting run "${RUN_NO}
  python ${STRAIN_SCRIPTS}/splitter.py -i ${DATA_DIR}/run_${RUN_NO}/${RM_LIBRARY} -o ${DATA_DIR}/run_${RUN_NO}/raw

  # run trf to determine if sequence is tandem repeat
  echo "Initial trf check for "${RUN_NO}
  parallel --bar --jobs ${THREADS} -a ${DATA_DIR}/run_${RUN_NO}/raw/${RM_LIBRARY}_split.txt trf ${DATA_DIR}/run_${RUN_NO}/raw/{} 2 7 7 80 10 50 500 -d -h -ngs ">" ${DATA_DIR}/run_${RUN_NO}/raw/{}.trf
  echo "Initial blast and preparation for MSA "${RUN_NO}
  # initial blast and extention
  parallel --bar --jobs ${THREADS} -a ${DATA_DIR}/run_${RUN_NO}/raw/${RM_LIBRARY}_split.txt python ${STRAIN_SCRIPTS}/initial_mafft_setup.py -d ${DATA_DIR} -r ${RUN_NO} -s {} -g ${GENOME} -f ${FLANK} -D -n ${NO_SEQ}
  
  ## first mafft alignment
  find ${DATA_DIR}/run_${RUN_NO}/to_align -type f | sed 's/.*\///' > ${DATA_DIR}/run_${RUN_NO}/to_align.txt
  echo "Primary alignment run "${RUN_NO}
  parallel --bar --jobs $MAFFT_THREADS -a ${DATA_DIR}/run_${RUN_NO}/to_align.txt mafft --thread 4 --quiet --localpair --adjustdirectionaccurately ${DATA_DIR}/run_${RUN_NO}/to_align/{} ">" ${DATA_DIR}/run_${RUN_NO}/mafft/{}

  # trim
  echo "Trimming run "${RUN_NO}
  parallel --bar --jobs $MAFFT_THREADS -a ${DATA_DIR}/run_${RUN_NO}/to_align.txt python ${STRAIN_SCRIPTS}/TEtrim.py -i ${DATA_DIR}/run_${RUN_NO}/mafft/{} -t 4 -f ${FLANK} -n ${RUN_NO} -d ${DATA_DIR} -m ${MIN_SEQ}
  
  # compile completed curations
  if [ -n "$(ls -A ${DATA_DIR}/run_${RUN_NO}/TEtrim_complete/ 2>/dev/null)" ]; then
     cat ${DATA_DIR}/run_${RUN_NO}/TEtrim_complete/*fasta > ${DATA_DIR}/run_${RUN_NO}/complete_${RM_LIBRARY}
  fi
  
  # Either compile consensuses for further curation
  if [ -n "$(ls -A ${DATA_DIR}/run_${RUN_NO}/TEtrim_further/ 2>/dev/null)" ]; then
    cat ${DATA_DIR}/run_${RUN_NO}/TEtrim_further/*fasta > ${DATA_DIR}/run_${RUN_NO}/further_${RM_LIBRARY}
    echo "Ready for " $(( RUN_NO++ ))
  else
  # or exit loop if no more curation needed
    echo "Finished extension"
    break
  fi

done

# Compile all completed
cat ${DATA_DIR}/run_*/complete_${RM_LIBRARY} > ${DATA_DIR}/${RM_LIBRARY}
# if more could have been extended append to compilation
if [ -s ${DATA_DIR}/run_${RUNS}/further_${RM_LIBRARY} ]; then
  cat ${DATA_DIR}/run_${RUNS}/further_${RM_LIBRARY} >> ${DATA_DIR}/${RM_LIBRARY}
fi

# add any families which went missing along the way due to alignment issues
find ${DATA_DIR}/run_*/mafft/ -empty -type f | sed 's/mafft/raw/' > ${DATA_DIR}/missing_consensi.txt
if [[ -s ${DATA_DIR}/missing_consensi.txt ]]; then
  while read file_name; do
    cat $file_name >> ${DATA_DIR}/${RM_LIBRARY}
  done < ${DATA_DIR}/missing_consensi.txt
fi
  
gsed -i 's/ .*//' ${DATA_DIR}/${RM_LIBRARY}

# Identify simple repeats and satellites, trim ends of LINEs/SINEs
echo "Splitting for simple/satellite packages"
mkdir -p ${DATA_DIR}/trf/split
python ${STRAIN_SCRIPTS}/splitter.py -i ${DATA_DIR}/${RM_LIBRARY} -o ${DATA_DIR}/trf/split
cp ${DATA_DIR}/${RM_LIBRARY} ${DATA_DIR}/trf/
# Run and parse TRF
echo "Running TRF"
parallel --bar --jobs ${THREADS} -a ${DATA_DIR}/trf/split/${RM_LIBRARY}_split.txt trf ${DATA_DIR}/trf/split/{} 2 7 7 80 10 50 500 -d -h -ngs ">" ${DATA_DIR}/trf/split/{}.trf
parallel --bar --jobs ${THREADS} -a ${DATA_DIR}/trf/split/${RM_LIBRARY}_split.txt python3 ${STRAIN_SCRIPTS}/trf_parser.py --trf ${DATA_DIR}/trf/split/{}.trf --out ${DATA_DIR}/trf/split/{}.trf.tsv
find ${DATA_DIR}/trf/split/ -type f -name "*trf.tsv" -exec cat {} + | cat > ${DATA_DIR}/trf/${RM_LIBRARY}.trf
# Run SA-SSR
echo "Running SA-SSR"
sa-ssr -e -l 20 -L 50000 -m 1 -M 5000 -t ${THREADS} ${DATA_DIR}/${RM_LIBRARY} ${DATA_DIR}/trf/${RM_LIBRARY}.sassr
# Run and compile mreps
echo "Running mreps"
parallel --bar --jobs ${THREADS} -a ${DATA_DIR}/trf/split/${RM_LIBRARY}_split.txt bash ${STRAIN_SCRIPTS}/mreps_parser.sh -i ${DATA_DIR}/trf/split/{} "2>"/dev/null
find ${DATA_DIR}/trf/split/ -type f -name "*mreps" -exec cat {} + | cat > ${DATA_DIR}/trf/${RM_LIBRARY}.mreps
# Interpret mreps, TRF and SA-SSR
echo "Trimming and sorting based on mreps, TRF, SA-SSR"
if [ ! -f ${DATA_DIR}/trf/${RM_LIBRARY}.sassr ]; then
   touch ${DATA_DIR}/trf/${RM_LIBRARY}.sassr
fi

$(which earlGrey | gsed 's/earlGrey/Rscript/g') ${STRAIN_SCRIPTS}/simple_repeat_filter_trim.R -i ${DATA_DIR}/${RM_LIBRARY} -d ${DATA_DIR}

# Delete temp files
echo "Removing temporary files"
rm -r ${DATA_DIR}/*/split/
find ${DATA_DIR}/ -mindepth 1 -name "run_*" -exec rm -r {} +
if [[ $RUNS -gt 0 ]]; then rm ${GENOME}.n*; fi

# Classify improved consensi using RepeatModeler's RepeatClassifier
echo "Reclassifying repeats"
mkdir -p ${DATA_DIR}/classify/
cp ${DATA_DIR}/trf/${RM_LIBRARY}.nonsatellite ${DATA_DIR}/classify/
cd ${DATA_DIR}/classify/
if [ -s ${RM_LIBRARY}.nonsatellite ]; then
    RepeatClassifier -threads ${THREADS} -consensi ${RM_LIBRARY}.nonsatellite
fi
# legacy command for older installations - DEPRECATED
# RepeatClassifier -pa ${THREADS} -consensi ${RM_LIBRARY}.nonsatellite
cd -
# Compile classified files
if [ -f ${DATA_DIR}/classify/${RM_LIBRARY}.nonsatellite.classified ]; then
    cp ${DATA_DIR}/classify/${RM_LIBRARY}.nonsatellite.classified ${DATA_DIR}/${RM_LIBRARY}.strained
else
    touch ${DATA_DIR}/${RM_LIBRARY}.strained
fi
echo "Compiling library"
if [ -f ${DATA_DIR}/trf/${RM_LIBRARY}.satellites ]; then
    cat ${DATA_DIR}/trf/${RM_LIBRARY}.satellites >> ${DATA_DIR}/${RM_LIBRARY}.strained
fi
```

Save the file with `CTRL+X` then press `Y` when asked to overwrite the file.

Make sure the updated file is executable:
```
chmod a+x $(which earlGrey | gsed 's|bin.*|share/earlgrey-5.0.0-0/scripts/TEstrainer/TEstrainer_for_earlGrey.sh|g')
```

Edit the script directory path in this file by running the following:
```
gsed -i "s|INSERT_FILENAME_HERE|$(which earlGrey | gsed 's:bin.*:share/earlgrey-5.0.0-0/scripts/TEstrainer/scripts/:g')|g" $(which earlGrey | gsed 's|bin.*|share/earlgrey-5.0.0-0/scripts/TEstrainer/TEstrainer_for_earlGrey.sh|g')
```

Edit famdb.py for use with our environment:
```
gsed -i 's/python3/python/g' $(which earlGrey | gsed 's|bin.*|share/RepeatMasker/famdb.py|g')
```

Edit LTR_FINDER_PARALLEL to be compatible with zsh
```
gsed -i "s|\`timeout $timeout|\`gtimeout $timeout|g" $(which earlGrey | gsed 's|bin.*|share/earlgrey-5.0.0-0/scripts/LTR_FINDER_parallel|g')
```

Install LTR_Finder from source
```
cd $(which earlGrey | gsed 's|bin.*|share/earlgrey-5.0.0-0/scripts/bin|g')
git clone https://github.com/xzhub/LTR_Finder
cd ./LTR_Finder/source
make
cp * ../../LTR_FINDER.x86_64-1.0.7/
```

Edit rcMergeRepeatsLoose:
```
gsed -i 's|sed|gsed|g' $(which earlGrey | gsed 's|bin.*|share/earlgrey-5.0.0-0/scripts/rcMergeRepeatsLoose|g')
var=$(which earlGrey | gsed "s/earlGrey/Rscript/g")
gsed -i "s|Rscript|${var}|g" $(which earlGrey | gsed 's|bin.*|share/earlgrey-5.0.0-0/scripts/rcMergeRepeatsLoose|g')
```

Edit main earlGrey script:
```
gsed -i "s|Rscript|${var}|g" $(which earlGrey | gsed 's|bin.*|share/earlgrey-5.0.0-0/earlGrey|g')
```

Add an important directory to PERL5LIB (for RepeatMasker)
```
echo "export PERL5LIB=$(which earlGrey | sed 's:bin.*:share/RepeatMasker:g')" >> ~/.start_miniconda3.sh
```

You are ready to go! Just remember to activate the _intel_ terminal, then the conda environment before running Earl Grey.

# A [biocontainer](https://biocontainers.pro/tools/earlgrey) has been generated from the Earl Grey Bioconda package.

## Get the preconfigured docker container from biocontainers
## This image will be configured with Dfam curated elements only
In this case, we need to bind a system directory to the docker container. In the line below, we are binding a directory call `host_data` that is found on our current path to `/data/` in the docker container. Please replace the file path before `:` to the directory you wish to bind to `/data/` in the container. This container must be run in interactive mode the first time you use it.

```
docker run -it -v `pwd`/host_data/:/data/ quay.io/biocontainers/earlgrey:5.0.0--h4ac6f70_0
```

## If you are running the container for the first time, you need to enable Earl Grey to configure the Dfam libraries correctly in interactive mode.
This is done by running Earl Grey for the first time. You will be prompted if Dfam libraries have not been configured properly. Answer `y` to configure Dfam.
Ensure the genome you wish to annotate is found in the directory you want to bind to the container. It will then be found in `/data` in the container. Alternatively, if you source a genome assembly from a database, download it to `/data/`. Always output results to the bound directory. In our case `/data/`. The results will then be found in your system directory after exiting the container.

```
earlGrey -g /data/genome.fasta -s test_genome -t 8 -o /data/
```

## If you need the container to run offline and/or without interactive mode
I try to keep an up-to-date container in docker hub, but this might not always be the case depending on if I have had time to build and upload a new image. Currently, there is an image with Dfam 3.7 curated elements only, and this is version 5.0.0. You can use this image by pulling the container:

```
# Interactive mode
docker run -it -v `pwd`:/data/ tobybaril/earlgrey_dfam3.7:latest

# Non interactive mode example:
docker run -v `pwd`:/data/ tobybaril/earlgrey_dfam3.7:latest earlGrey -g /data/NC_045808_EarlWorkshop.fasta -s nonInteractiveTest -o /data/ -t 8
```

# Dfam Workaround For Latest Release

With the latest release of RepeatMasker (v4.1.7), Dfam 3.8 has been reorganised into partitions containing both curated and uncurated sequences for specific taxonomic groups (see https://dfam.org/releases/Dfam_3.8/families/FamDB/README.txt). Consequently, it is challenging to provide a stable conda release for RepeatMasker 4.1.7. 

I have built a container for Earl Grey with RepeatMasker 4.1.7 and the root partition (partition 0) of Dfam version 3.8 preconfigured. This is particularly useful for those who require working with docker containers on their HPC infrastructure.

```
# to run in interactive mode with bound directories
docker run -it -v `pwd`/host_data/:/data/ tobybaril/earlgrey_dfam3.8:latest

# to run the container passing the required input file
docker run -w /data/ -v /path/to/system/directory/inputGenome.fasta:/data/inputGenome.fasta --name=earlGreyContainer tobybaril/earlgrey_dfam3.8:latest earlGrey -g /data/inputGenome.fasta -s input -o /data/ -t $threads

# to get the output files from the stopped container
docker cp earlGreyContainer:/data/input_EarlGrey /path/to/system/directory
```

If you require Dfam 3.8 for your studies, I have devised a workaround that remains functional. HOWEVER, I do not recommend attempting this unless you are comfortable with altering files within conda environments, and have a good level of experience in configuring tools. Undertake the below at your own risk!


First, locate your conda environment installation. It should be something like the below:
```
/home/user/anaconda3/envs/earlgrey/
```

Next, find the directory containing RepeatMasker and associated files. It will be inside the share directory
```
cd /home/user/anaconda3/envs/earlgrey/share/
```

Compress and back up RepeatMasker in case anything goes wrong, then delete the RepeatMasker Directory
```
tar -czvf RepeatMasker.bak.tar.gz RepeatMasker/ && rm -r ./RepeatMasker/
```

Download and unpack RepeatMasker 4.1.7
```
wget https://repeatmasker.org/RepeatMasker/RepeatMasker-4.1.7-p1.tar.gz
tar -zxvf RepeatMasker-4.1.7-p1.tar.gz
```

Go to the famdb directory and fetch at least the root partition. check the [README](https://dfam.org/releases/Dfam_3.8/families/FamDB/README.txt) for which partition contains which species, and download all the ones you want.
```
cd ./RepeatMasker/Libraries/famdb/
wget https://dfam.org/releases/Dfam_3.8/families/FamDB/dfam38_full.0.h5.gz
gunzip *.gz
```

Move back to the RepeatMasker directory and reconfigure RepeatMasker
```
cd ../../
perl ./configure

# when prompted, you will need the full path to trf. This is in the conda environment bin directory:
## e.g in this example it would be:
/home/user/anaconda3/envs/earlgrey/bin/trf

# NOTE: In some cases I ran into an issue where the configuration file needed modifying. Check the following:
nano ./RepeatMaskerConfig.pm

# Check that the DEFAULT_SEARCH_ENGINE block looks like this:
'DEFAULT_SEARCH_ENGINE' => {
                                       'command_line_override' => 'default_search_engine',
                                       'description' => 'The default search engine to use',
                                       'param_type' => 'value',
                                       'required' => 1,
                                       'value' => 'rmblast'
                                     },
```

# Alternative installation methods

Before using Earl Grey, please ensure RepeatMasker (>=version 4.1.4) and RepeatModeler (>=version 2.0.4) are installed and configured. If these are not, please follow the instructions below to install these before continuing with Earl Grey Installation. 
NOTE: These instructions are provided to install RepeatMasker, RepeatModeler and related programs with sudo privileges. If you are working on a shared cluster, please request installation of RepeatMasker and RepeatModeler by your sysadmin before working with Earl Grey. Earl Grey will function with RepeatMasker and RepeatModeler installed in the local path environment, or when the modules are loaded on a HPC cluster.


#==============================================================================================================================================================================#

## Earl Grey is available in a [Docker container](./Docker) preconfigured with Dfam curated version 3.7

To use this container, please make sure Docker is installed and configured on your system. All instructions relating to the docker installation are found within the [Docker directory](./Docker) in this repository. Please consult the README in the Docker directory for installation instructions.

#==============================================================================================================================================================================#

## Earl Grey is available for [Singularity](./Singularity) preconfigured with Dfam curated version 3.7

To use this container, please make sure Singularity is installed and configured on your system. All instructions relating to the singularity installation are found within the [Singularity directory](./Singularity) in this repository. Please consult the README in the Singularity directory for installation instructions.

#==============================================================================================================================================================================#

## Earl Grey Installation and Configuration (If you already have RepeatMasker and RepeatModeler)

If you do not currently have RepeatMasker and RepeatModeler installed, the instructions are provided further down this page. If you do have them installed, **ensure the executables are in your PATH environment, including the RepeatMasker/util/ directory!**

All of the scripts and associated modules are contained within this github repository. Earl Grey runs inside an anaconda environment to ensure all required packages are present and are the correct version. Therefore to run Earl Grey, you will require anaconda to be installed on your system. 

**If anaconda is NOT installed on your system, please install it following these instructions:**

```
# Change to /tmp directory as we won't need the script after running it
cd /tmp

# Download the anaconda installation script
curl https://repo.anaconda.com/archive/Anaconda3-2020.02-Linux-x86_64.sh --output anaconda.sh

# Run the script to install anaconda
bash anaconda.sh
# answer yes when asked, and install anaconda3 in the specified location (recommended) unless you want it to be installed elsewhere.
# When asked "Do you wish the installer to initialize Anaconda3 by running conda init?", answer "yes" for ease of use.
# Activate conda by refreshing terminal session

source ~/.bashrc

# If successful, you should now see (base) on the left of your username on the command line
```

**Now that Anaconda is installed, we can install Earl Grey**

Clone the Earl Grey github repo

```
# Clone into a home directory, or somewhere you want to install Earl Grey
git clone https://github.com/TobyBaril/EarlGrey
```

Enter the Earl Grey directory and configure the program

```
cd ./EarlGrey
chmod +x ./configure
./configure
```

Once this is complete, remember to activate the earlGrey conda environment before attempting to run the Earl Grey pipeline

```
conda activate earlGrey

earlGrey -g genome.fasta -s speciesName -o outputDirectory -t threads
```

For suggestions or questions, please use the discussion and issues functions in github.

Thank you for trying Earl Grey!

#==============================================================================================================================================================================#

## Earl Grey Installation and Configuration (If you DO NOT have RepeatMasker and RepeatModeler) - WITH SUDO PRIVILEGES

These instructions will guide you through configuring all required programs and scripts to run Earl Grey.

All of the scripts and associated modules are contained within this github repository. Earl Grey runs inside an anaconda environment to ensure all required packages are present and are the correct version. Therefore to run Earl Grey, you will require anaconda to be installed on your system. 

**If anaconda is NOT installed on your system, please install it following these instructions:**

```
# Change to /tmp directory as we won't need the script after running it
cd /tmp

# Download the anaconda installation script
curl https://repo.anaconda.com/archive/Anaconda3-2020.02-Linux-x86_64.sh --output anaconda.sh

# Run the script to install anaconda
bash anaconda.sh
# answer yes when asked, and install anaconda3 in the specified location (recommended) unless you want it to be installed elsewhere.
# When asked "Do you wish the installer to initialize Anaconda3 by running conda init?", answer "yes" for ease of use.
# Activate conda by refreshing terminal session

source ~/.bashrc

# If successful, you should now see (base) on the left of your username on the command line
```

**Now that Anaconda is installed, we can install Earl Grey**

### To install RepeatMasker

RepeatMasker can be downloaded from: http://www.repeatmasker.org/RepeatMasker/. Installation instructions can be found on the website. Alternatively, please use the code below:

You will need to download and install a couple of programs for RepeatMasker to work. 

Download rmblast

```wget https://www.repeatmasker.org/rmblast/rmblast-2.13.0+-x64-linux.tar.gz```

Extract the rmblast package

```tar -zxvf rmblast-2.13.0+-x64-linux.tar.gz```

Make a note of the full path to rmblast-2.13.0/bin/ as you will need this in the RepeatMasker configuration
If you are not certain of the full path, please run the following command

```realpath ./rmblast-2.13.0/bin/```

Download RepeatMasker (this will download it to the current directory). NOTE some extra steps are required for Dfam 3.8, please refer to repeatmasker.org/RepeatMasker.

```wget https://www.repeatmasker.org/RepeatMasker/RepeatMasker-4.1.6.tar.gz```

Copy the RepeatMasker package to /usr/local/, or somewhere that all users will be able to access the installation. 
Copying to /usr/local/ might require sudo privileges

```sudo cp RepeatMasker-4.1.6.tar.gz /usr/local/```

Change directory to /usr/local, and extract the RepeatMasker package. This might require sudo privileges.

```cd /usr/local/```

```sudo tar -zxvf RepeatMasker-4.1.6.tar.gz```

Install the required RepeatMasker libraries - Earl Grey has been tested with Dfam 3.4 onwards and RepBase. 
Unfortunately, RepBase is now behind a paywall, but to ensure Earl Grey remains open it does not rely on RepBase, although inclusion of RepBase can improve classification of repeats by RepeatModeler. If you have access to this database, please include it in your configuration of RepeatMasker. 
We recommend that you download Dfam 3.8 before using Earl Grey, and this is required for RepeatMasker 4.1.6. The Dfam library is large - this could take a while!

We recommend downloading Dfam into your home directory (~/) or a subdirectory of home

Change directory to home

```cd ~/```

Download lastest Dfam release - This may take a while
*** NOTE: There are many partitions for Dfam, each containing curated and uncurated TE consensi split by taxonomy. As a minimum, you need to configure RepeatMasker with the root partition, but you can also include others. Note: There may be erroneous annotations in the uncurated Dfam database, use at your own risk***

For whole Dfam library run all of the below, or run each part as required (minimum the root (0) partition):
```
wget https://dfam.org/releases/Dfam_3.8/families/FamDB/dfam38_full.0.h5.gz
wget https://dfam.org/releases/Dfam_3.8/families/FamDB/dfam38_full.1.h5.gz
wget https://dfam.org/releases/Dfam_3.8/families/FamDB/dfam38_full.2.h5.gz
wget https://dfam.org/releases/Dfam_3.8/families/FamDB/dfam38_full.3.h5.gz
wget https://dfam.org/releases/Dfam_3.8/families/FamDB/dfam38_full.4.h5.gz
wget https://dfam.org/releases/Dfam_3.8/families/FamDB/dfam38_full.5.h5.gz
wget https://dfam.org/releases/Dfam_3.8/families/FamDB/dfam38_full.6.h5.gz
wget https://dfam.org/releases/Dfam_3.8/families/FamDB/dfam38_full.7.h5.gz
wget https://dfam.org/releases/Dfam_3.8/families/FamDB/dfam38_full.8.h5.gz
```

Unzip the Dfam release - This may take a while with no indication that anything is happening, please be patient!

For whole Dfam library (including uncurated elements):
```gunzip *.h5.gz```

Move the Dfam library to the RepeatMasker library folder
NOTE, a warning might come up that this will overwrite the existing file, allow this by pressing "y" then Enter

```mv *.h5 /usr/local/RepeatMasker/Libraries/famdb/```

DO NOT configure RepeatMasker just yet...

Add RepeatMasker and util directory to your path environment

```echo 'export PATH=$PATH:/usr/local/RepeatMasker:/usr/local/RepeatMasker/util/' >> ~/.bashrc```

### To install RepeatModeler

RepeatModeler can be downloaded from: http://www.repeatmasker.org/RepeatModeler/. Installation instructions can be found on the website. Alternatively, please use the code below:

You will need to download and install a couple of programs for RepeatModeler to work. 

Download RECON and RepeatScout

```
wget http://www.repeatmasker.org/RepeatModeler/RECON-1.08.tar.gz
wget http://www.repeatmasker.org/RepeatScout-1.0.6.tar.gz
```

Extract the RECON and RepeatScout packages

```
tar -zxvf RECON-1.08.tar.gz
tar -zxvf RepeatScout-1.0.6.tar.gz
```

Enter the RECON directory and make from source

```
cd ./RECON-1.08/src/
make
make install
```

Enter the RepeatScout directory and make from source

```
cd ../../RepeatScout-1.0.6/
make
```

Make a note of the full path to RECON-1.08/bin/ and RepeatScout-1.0.6/ as you will need these in the RepeatModeler configuration
If you are not certain of the full paths, please run the following commands

```
realpath ./RECON-1.08/bin/
realpath ./RepeatScout-1.0.6/
```

Install UCSC TwoBit Tools. Make a note of the full path to the directory where these have been installed as you will need these in the RepeatModeler configuration

```
mkdir ~/ucscTwoBitTools/
cd ~/ucscTwoBitTools/
wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/faToTwoBit
wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/twoBitToFa
wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/twoBitInfo
wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/twoBitMask
wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/twoBitDup

# make scripts executable
chmod +x ~/ucscTwoBitTools/*
```

For Earl Grey, the LTR structural pipeline of RepeatModeler is not required, as this is run as part of Earl Grey's defragmentation process.

Download RepeatModeler

```wget http://www.repeatmasker.org/RepeatModeler/RepeatModeler-2.0.4.tar.gz```

Copy the RepeatModeler package to /usr/local/, or somewhere that all users will be able to access the installation. 
Copying to /usr/local/ might require sudo privileges

```sudo cp RepeatModeler-2.0.4.tar.gz /usr/local/```

Change directory to /usr/local, and extract the RepeatModeler package. This might require sudo privileges.

```cd /usr/local/```

```sudo tar -zxvf RepeatModeler-2.0.4.tar.gz```

DO NOT configure RepeatModeler just yet...

Add RepeatModeler directory to your path environment

```echo 'export PATH=$PATH:/usr/local/RepeatModeler-2.0.4/' >> ~/.bashrc```

### To Install Earl Grey and Configure all Programs

Clone the Earl Grey github repo

```
# Clone into a home directory, or somewhere you want to install Earl Grey
git clone https://github.com/TobyBaril/EarlGrey
```

Enter the Earl Grey directory and configure the program

```
cd ./EarlGrey
chmod +x ./configure
./configure
```

Activate the conda environment

```
conda activate earlGrey
```

For configuring RepeatMasker with sudo, h5py is required in the root path. The easiest way to do this is to run the following two lines:

```
sudo apt install python3-pip
sudo pip install h5py
```

The linux64 version of TRF is included in the modules directory of Earl Grey. Make a note of the full path to this file as this will be needed in the RepeatMasker configuration.

```realpath ./modules/trf409*```

NOTE: For the latest version of RepeatModeler, this must be called trf. We can make a symlinked version: ```ln -s ./modules/trf409* ./modules/trf```

Now to finally configure RepeatMasker! - Run these lines and then follow the on-screen prompts from RepeatMasker. Again, this might require sudo privileges.

```
cd /usr/local/RepeatMasker
# check which perl interpreter you should use 
which perl
# replace [perl] in the below command with the path printed from the command above and use this as your perl interpreter (this should be the anaconda perl interpreter)
sudo [perl] ./configure
```

Time to configure RepeatModeler! You need to enter the paths to lots of the directories of programs we have installed, please note them down before running the configuration script!! 


```
cd /usr/local/RepeatModeler-2.0.4/
# check which perl interpreter you should use - COPY THIS
which perl
# replace perl in the below command with the path printed from the command above and use this as your perl interpreter
sudo perl ./configure
```

Once this is complete, remember to activate the earlGrey conda environment before attempting to run the Earl Grey pipeline

```
conda activate earlGrey

earlGrey -g genome.fasta -s speciesName -o outputDirectory -t threads
```

#==============================================================================================================================================================================#

## Earl Grey Installation and Configuration (If you DO NOT have RepeatMasker and RepeatModeler) - WITHOUT SUDO PRIVILEGES

These instructions will guide you through configuring all required programs and scripts to run Earl Grey.

All of the scripts and associated modules are contained within this github repository. Earl Grey runs inside an anaconda environment to ensure all required packages are present and are the correct version. Therefore to run Earl Grey, you will require anaconda to be installed on your system. 

**If anaconda is NOT installed on your system, please install it following these instructions:**

```
# Change to /tmp directory as we won't need the script after running it
cd /tmp

# Download the anaconda installation script
curl https://repo.anaconda.com/archive/Anaconda3-2020.02-Linux-x86_64.sh --output anaconda.sh

# Run the script to install anaconda
bash anaconda.sh
# answer yes when asked, and install anaconda3 in the specified location (recommended) unless you want it to be installed elsewhere.
# When asked "Do you wish the installer to initialize Anaconda3 by running conda init?", answer "yes" for ease of use.
# Activate conda by refreshing terminal session

source ~/.bashrc

# If successful, you should now see (base) on the left of your username on the command line
```

**Now that Anaconda is installed, we can install Earl Grey**

### To install RepeatMasker

RepeatMasker can be downloaded from: http://www.repeatmasker.org/RepeatMasker/. Installation instructions can be found on the website. Alternatively, please use the code below:

You will need to download and install a couple of programs for RepeatMasker to work. 

Download rmblast

```wget https://www.repeatmasker.org/rmblast/rmblast-2.13.0+-x64-linux.tar.gz```

Extract the rmblast package

```tar -zxvf rmblast-2.13.0+-x64-linux.tar.gz```

Make a note of the full path to rmblast-2.13.0/bin/ as you will need this in the RepeatMasker configuration
If you are not certain of the full path, please run the following command

```realpath ./rmblast-2.13.0/bin/```

Download RepeatMasker (this will download it to the current directory). NOTE some extra steps are required if you wish to use Dfam 3.7, please refer to repeatmasker.org/RepeatMasker for these extra steps if required.

```wget https://www.repeatmasker.org/RepeatMasker/RepeatMasker-4.1.6.tar.gz```

Extract the RepeatMasker package.

```
tar -zxvf RepeatMasker-4.1.6.tar.gz
```

Install the required RepeatMasker libraries - Earl Grey has been tested with Dfam 3.3 and RepBase. 
Unfortunately, RepBase is now behind a paywall, but to ensure Earl Grey remains open it does not rely on RepBase, although inclusion of RepBase can improve classification of repeats by RepeatModeler. If you have access to this database, please include it in your configuration of RepeatMasker. 
We recommend that you download Dfam 3.8 before using Earl Grey, and this is required for RepeatMasker 4.1.6. The Dfam library is large - this could take a while!

Change directory to RepeatMasker Libraries

```cd ./RepeatMasker/Libraries/famdb/```

Download lastest Dfam release - This may take a while
*** NOTE: There are many partitions for Dfam, each containing curated and uncurated TE consensi split by taxonomy. As a minimum, you need to configure RepeatMasker with the root partition, but you can also include others. Note: There may be erroneous annotations in the uncurated Dfam database, use at your own risk***

For whole Dfam library run all of the below, or run each part as required (minimum the root (0) partition):
```
wget https://dfam.org/releases/Dfam_3.8/families/FamDB/dfam38_full.0.h5.gz
wget https://dfam.org/releases/Dfam_3.8/families/FamDB/dfam38_full.1.h5.gz
wget https://dfam.org/releases/Dfam_3.8/families/FamDB/dfam38_full.2.h5.gz
wget https://dfam.org/releases/Dfam_3.8/families/FamDB/dfam38_full.3.h5.gz
wget https://dfam.org/releases/Dfam_3.8/families/FamDB/dfam38_full.4.h5.gz
wget https://dfam.org/releases/Dfam_3.8/families/FamDB/dfam38_full.5.h5.gz
wget https://dfam.org/releases/Dfam_3.8/families/FamDB/dfam38_full.6.h5.gz
wget https://dfam.org/releases/Dfam_3.8/families/FamDB/dfam38_full.7.h5.gz
wget https://dfam.org/releases/Dfam_3.8/families/FamDB/dfam38_full.8.h5.gz
```

Unzip the Dfam release - This may take a while with no indication that anything is happening, please be patient!

For whole Dfam library (including uncurated elements):
```gunzip *.h5.gz```

DO NOT configure RepeatMasker just yet...

Add RepeatMasker and util directory to your path environment (Replace /path/to/ with the full path to your installation directory)

```echo 'export PATH=$PATH:/path/to/RepeatMasker/:/path/to/RepeatMasker/util/' >> ~/.bashrc```

### To install RepeatModeler

RepeatModeler can be downloaded from: http://www.repeatmasker.org/RepeatModeler/. Installation instructions can be found on the website. Alternatively, please use the code below:

You will need to download and install a couple of programs for RepeatModeler to work. 

Download RECON and RepeatScout (Make sure you are no longer inside your RepeatMasker directory!)

```
wget http://www.repeatmasker.org/RepeatModeler/RECON-1.08.tar.gz
wget http://www.repeatmasker.org/RepeatScout-1.0.6.tar.gz
```

Extract the RECON and RepeatScout packages

```
tar -zxvf RECON-1.08.tar.gz
tar -zxvf RepeatScout-1.0.6.tar.gz
```

Enter the RECON directory and make from source

```
cd ./RECON-1.08/src/
make
make install
```

Enter the RepeatScout directory and make from source

```
cd ../../RepeatScout-1.0.6/
make
```

Make a note of the full path to RECON-1.08/bin/ and RepeatScout-1.0.6/ as you will need these in the RepeatModeler configuration
If you are not certain of the full paths, please run the following commands

```
realpath ./RECON-1.08/bin/
realpath ./RepeatScout-1.0.6/
```

Install UCSC TwoBit Tools. Make a note of the full path to the directory where these have been installed as you will need these in the RepeatModeler configuration

```
mkdir ./ucscTwoBitTools/
cd ./ucscTwoBitTools/
wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/faToTwoBit
wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/twoBitToFa
wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/twoBitInfo
wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/twoBitMask
wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/twoBitDup

# make scripts executable
chmod +x ./ucscTwoBitTools/*
```

For Earl Grey, the LTR structural pipeline of RepeatModeler is not required, as this is run as part of Earl Grey's defragmentation process.

Download RepeatModeler

```wget http://www.repeatmasker.org/RepeatModeler/RepeatModeler-2.0.4.tar.gz```

Unpack RepeatModeler

```sudo tar -zxvf RepeatModeler-2.0.4.tar.gz```

DO NOT configure RepeatModeler just yet...

Add RepeatModeler directory to your path environment (Replace /path/to/ with the full path to your installation directory)


```echo 'export PATH=$PATH:/path/to/RepeatModeler-2.0.4/' >> ~/.bashrc```

### To Install Earl Grey and Configure all Programs

Clone the Earl Grey github repo

```
# Clone into a home directory, or somewhere you want to install Earl Grey
git clone https://github.com/TobyBaril/EarlGrey
```

Enter the Earl Grey directory and configure the program

```
cd ./EarlGrey
chmod +x ./configure
./configure
```

Activate the conda environment

```
conda activate earlGrey
```

The linux64 version of TRF is included in the modules directory of Earl Grey. Make a note of the full path to this file as this will be needed in the RepeatMasker configuration.

```realpath ./modules/trf409*```

NOTE: For the latest version of RepeatModeler, this must be called trf. We can make a symlinked version: ```ln -s ./modules/trf409* ./modules/trf```

Now to finally configure RepeatMasker! - Run these lines and then follow the on-screen prompts from RepeatMasker. Again, this might require sudo privileges.

```
cd /path/to/RepeatMasker
perl ./configure
```

Time to configure RepeatModeler! You need to enter the paths to lots of the directories of programs we have installed, please note them down before running the configuration script!! 


```
cd /path/to/RepeatModeler-2.0.4/
perl ./configure
```

Once this is complete, remember to activate the earlGrey conda environment before attempting to run the Earl Grey pipeline

```
conda activate earlGrey

earlGrey -g genome.fasta -s speciesName -o outputDirectory -t threads
```


