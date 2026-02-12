![earlGreyIcon](https://user-images.githubusercontent.com/46785187/136248346-21e980ee-1154-48c2-9398-70938bbe2404.png)

[![DOI](https://zenodo.org/badge/412126708.svg)](https://zenodo.org/badge/latestdoi/412126708) [![Anaconda-Server Badge](https://anaconda.org/bioconda/earlgrey/badges/version.svg)](https://anaconda.org/bioconda/earlgrey) [![Anaconda_downloads](https://anaconda.org/bioconda/earlgrey/badges/downloads.svg)](
https://anaconda.org/bioconda/earlgrey) [![Anaconda_platforms](https://anaconda.org/bioconda/earlgrey/badges/platforms.svg)](https://anaconda.org/bioconda/earlgrey) [![Anaconda-Server Badge](https://anaconda.org/bioconda/earlgrey/badges/latest_release_relative_date.svg)](https://anaconda.org/bioconda/earlgrey)

# Earl Grey

Earl Grey is a full-automated transposable element (TE) annotation pipeline, leveraging the most widely-used tools and combining these with a consensus elongation process to better define _de novo_ consensus sequences when annotating new genome assemblies.

# Contents

[Changes in Latest Release](#changes-in-latest-release)

[Example Run](#example)

[References and Acknowledgements](#references-and-acknowledgements)

[Usage Without Installation](#usage-without-installation)

[Recommended Installation](#recommended-installation-with-conda-or-mamba)

[Docker Container](#docker-container)

<!-- toc -->

# Important Considerations

Earl Grey version 6 uses Dfam 3.9. After installation, you MUST configure Dfam partitions as needed. Earl Grey will generate the script to do this and provide guidance when you run it for the first time. You need to specify which partitions of Dfam and/or RepBase to configure Earl Grey with. Choose partitions carefully as the combination will highly influence your results, especially if you want to pre-mask your input genome. Please make use of issues and discussions tabs if you have questions about this, we are always happy to help!

# Notes / Updates

We often get questions related to runtime. TE curation and annotation remains resource and time intensive. Fast is not necessarily better, and runtime is highly dependent on genome size, complexity, and repeat content. Runs will likely take longer than you might expect, and be very RAM-hungry. As some generic benchmarks, a 40Mb genome can take anywhere from a few hours to a day, 400Mb up to around 4-5 days, a 3Gb genome ~a week, and a 25Gb genome several weeks! Things will be running even if it doesn't look like they are. Each step checkpoints, so if you have server limits, you can resubmit the same script with the same parameters, and Earl Grey will skip completed steps. `TEstrainer` and the final `divergence calculator` use a lot of memory. Check carefully for OOM errors in the logs! As a rule of thumb, you need _at least_ 3GB of RAM _per thread_, with more being better. Therefore, 16 threads requires at least 48GB of RAM depending on repeat complexity of the input genome.

We have been made aware of some instability in repeat annotation percentages when high numbers of CPUs are employed in certain server environments. Please be sure to check logs carefully for instances of interruption. Known cases so far will show the following message:

```
OpenBLAS blas_thread_init: pthread_create failed for thread X of X: Resource temporarily unavailable
OpenBLAS blas_thread_init: ensure that your address space and process count limits are big enough (ulimit -a)
OpenBLAS blas_thread_init: or set a smaller OPENBLAS_NUM_THREADS to fit into what you have available
```

If you see this message, re-run analysis with less threads. Alternatively, you can modify your instance of the TEstrainer script `initial_mafft_setup.py` to add the following after `import os`:

```
os.environ['OPENBLAS_NUM_THREADS'] = '1'
```

# Changes in Latest Release

Earl Grey v7.0.2 adds RepeatLandscapes for _Penelope_-like elements and SINEs. Importantly, the `-norna` option in RepeatMasker is no longer invoked as default behaviour, which will improver the detection and masking of small tRNA-derived SINEs.

### Previous Changes

Earl Grey v7.0.1 patches the summary table generation, where LINEs and Penelopes were being counted in both categories for nested repeats only.

☕ Earl Grey v7.0.0 is here!

Some long-awaited changes in this release—thank you for your patience while I found the time to properly test and implement them.

🐞 RepeatCraft fixes
First, I’ve fixed an issue in RepeatCraft where distal TEs could be grouped erroneously. This occurred in rare edge cases where the internal counter failed to iterate when new distal copies of an existing TE were detected on the same contig. This should now be fully resolved.

🧬 Fully nested TE handling (major update)
I’ve also completely revamped the Earl Grey post-processing and summary steps to properly deal with Fully Nested TEs.
Partially overlapping TEs are handled as before, following the approach described in the original Earl Grey paper.
Fully nested TEs are now identified using a new iterative process:
The GFF is scanned for nested TEs
These are labelled and stored in a separate file
The nested TE is removed and the GFF is re-scanned to detect deeper (multi-level) nesting
This continues until no new nested TEs are found

📊 Changes to coverage and summaries (important)
Nested TEs are no longer included in the TE coverage calculation used to generate the pie chart in the `summaryFiles` directory.

Instead:
Summary tables now include new categories showing how many base pairs are comprised of nested TEs
These base pairs are not counted toward Total Interspersed Repeats, as doing so would result in double-counting genomic space.

⚠️ This represents a substantial change from previous versions, so please be aware of this difference when upgrading to v7.

The output table summary now has the following format:

```
|TE Classification                                 | Coverage (bp)|Copy Number   | % Genome Coverage| Genome Size| TE Family Count|
|:-------------------------------------------------|-------------:|:-------------|-----------------:|-----------:|---------------:|
|DNA                                               |         80886|326           |         0.7607795|    10631990|             326|
|DNA-nested                                        |          1449|29            |         0.0136287|    10631990|              29|
|Rolling Circle                                    |          8022|31            |         0.0754515|    10631990|              31|
|Rolling Circle-nested                             |             0|0             |         0.0000000|    10631990|              NA|
|Penelope                                          |        112822|812           |         1.0611560|    10631990|             812|
|Penelope-nested                                   |          1268|13            |         0.0119263|    10631990|              13|
|LINE                                              |        132133|426           |         1.2427871|    10631990|             426|
|LINE-nested                                       |          1268|13            |         0.0119263|    10631990|               4|
|SINE                                              |             0|0             |         0.0000000|    10631990|              NA|
|SINE-nested                                       |             0|0             |         0.0000000|    10631990|              NA|
|LTR                                               |             0|0             |         0.0000000|    10631990|              NA|
|LTR-nested                                        |             0|0             |         0.0000000|    10631990|              NA|
|Other (Simple Repeat, Microsatellite, RNA)        |        189999|4241          |         1.7870502|    10631990|            4241|
|Other (Simple Repeat, Microsatellite, RNA)-nested |          1225|34            |         0.0115218|    10631990|              34|
|Unclassified                                      |        274026|1291          |         2.5773726|    10631990|            1291|
|Unclassified-nested                               |             0|0             |         0.0000000|    10631990|              NA|
|Total Interspersed Repeat                         |        607889|2886          |         5.7175468|    10631990|              NA|
|Non-Repeat                                        |      10024101|notApplicable |        94.2824532|    10631990|              NA|
```

In the final GFF output, nested TEs are clearly labelled with a `nested=FULLY_NESTED` attribute in column 9, enabling quick identification and downstream filtering.

```
NC_045808.1	Earl_Grey	Simple_repeat	51840	51871	14	+	.	ID=(CGCA)N_33;Name=(CGCA)N;tstart=1;tend=31;shortte=F;nested=FULLY_NESTED-ROUND1
```

As always, thank you to the TE community for your enthusiasm in using Earl Grey, and for your invaluable feedback and bug reports. I’ll continue to incorporate improvements and fixes as quickly—and carefully—as possible.

Happy New Year! 🎉

Earl Grey v6.3.6 patches a RepeatCraft bug that arises extremely rarely in specific genomes, linked to dictionary initialisation. 

Earl Grey v6.3.5 patches the annotation only pipeline to use the correct divergence calculator when a custom library is used without a RepBase term.

Earl Grey v6.3.4 adds small patches to correct TE family count table, which was previously using ID rather than Name. Also, a new awk one-liner has been added to clean the final GFF to make all attributes compliant with genometools, geneious etc.

Earl Grey v6.3.3 adds small updates to improve user-friendliness when installing and running for the first time.

Earl Grey v6.3.0 makes some important changes. Firstly, RepeatMasker now runs without checking for IS element contamination `-no_is`. In final GFF files, each insertion is now given a unique ID with `ID=`. The TE family is designated with `Name=`. This should enable parsing with Geneious or other visualisers. Other script changes are related to this change in field designation.

Earl Grey v6.2.0 now includes better checkpoints for TEstrainer. In the event of an interrupted run, resubmit your `earlGrey` command with _exactly_ the same parameters and Earl Grey will skip successfully completed steps. With this new update, TEstrainer intermediate files will no longer be deleted and a new run started from scratch. Now, TEstrainer can be recovered mid-run to save time and resources.

Earl Grey v6.1.1 patches a bug where threads set to <4 caused TEstrainer to crash (only present in 6.1.0 and not earlier versions)

Earl Grey v6.1.0 reintroduces the `--curated` flag when known elements are used to pre-mask your input genome. As usual, unless a good quality TE library already exists for your species of interest, or a _very_ closely related one, we do not recommend pre-masking your input genome with known repeats. This reduces the amount of information available for _de novo_ annotation, and can lead to overestimations of TE divergence and lower-quality consensus sequences. If in doubt, leave this option out!

Earl Grey v6.0.3 reduces CPU usage for TEstrainer to reduce memory pressure.

Earl Grey v6.0.2 patches an issue where the use of existing libraries did not work with the new `famdb` formats. 

Earl Grey v6.0.1 contains small bug fixes to verify installed RepeatMasker Libraries correctly. There is now a Docker container for Earl Grey v6.0.1 that contains all partitions of Dfam (so it is BIG!).

*Earl Grey v6.0.0 is here!* 

There are some relatively large changes in this release, resulting in the jump to v6.0.0. 

Importantly, Earl Grey has been updated to use *Dfam version 3.9*, with RepeatMasker 4.1.8 and famdb 2.0.1. This means that there is some extra configuration required to get the pipeline running. Upon first installation and running of Earl Grey, the pipeline will check whether RepeatMasker has been configured with the correct Dfam database partitions. If not, it will warn you, generate a script that you can modify and run to configure RepeatMasker, and provide instructions to `stdout` if you want to do this yourself.

Please take care to configure Earl Grey v6 with ALL required partitions. More information on the partitioning can be found at [Dfam.org](https://dfam.org/releases/current/families/FamDB/README.txt).

Earl Grey v5.1.1 will continue to work for those who are happy with Dfam v3.7, but we recommend upgrading to v6.0.0 to keep up to date with the latest improvements to the database.

Earl Grey v5.1.1 contains very small patches to improve compatibility with publicly available genome sequencing data. In rare instances, strange characters in fasta headers were causing issues preventing the pipeline from running. These have been resolved in the preparation step. 

In addition, new _pretty_ tables are now generated in the `summaryFiles` directory and at the end of a successful run. These contain the same information as the `txt` tables, but in the familiar pipe format for readability. These can be added to markdown files if required. One is produced for the high level count as well as for the family level count. Below is an example of the table that is printed at the end of Earl Grey runs as of `v5.1.1`:

```
|TE Classification                          | Coverage (bp)| Copy Number| % Genome Coverage| Genome Size| TE Family Count|
|:------------------------------------------|-------------:|-----------:|-----------------:|-----------:|---------------:|
|DNA                                        |         69284|         190|         0.6516560|    10631990|               9|
|Rolling Circle                             |         71788|         316|         0.6752076|    10631990|              12|
|Penelope                                   |        139334|         869|         1.3105167|    10631990|               6|
|LINE                                       |        141841|         491|         1.3340964|    10631990|              22|
|SINE                                       |         44137|         191|         0.4151339|    10631990|               1|
|Other (Simple Repeat, Microsatellite, RNA) |        191645|        4235|         1.8025318|    10631990|             940|
|Unclassified                               |        155836|         905|         1.4657275|    10631990|              35|
|Non-Repeat                                 |       9818125|          NA|        92.3451301|    10631990|              NA|
```

Earl Grey v5.1.0 contains small changes that drastically improve memory usage in the divergence calculator. We have replaced the use of EMBOSS `water` with EMBOSS `matcher`, which reduces memory consumption on large alignments whilst remaining rigorous. For more information, please see the notes section in the [EMBOSS Manual](https://www.bioinformatics.nl/cgi-bin/emboss/help/matcher). This should prevent jobs running out of memory, particularly when using queuing systems and shared resources. 

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
		-l == Starting consensus library for an initial mask (in fasta format)
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

Earl Grey version 6 uses Dfam 3.9. After installation, you MUST configure Dfam partitions as needed. Earl Grey will generate the script to do this and provide guidance when you run it for the first time. You need to specify which partitions of Dfam and/or RepBase to configure Earl Grey with. Choose partitions carefully as the combination will highly influence your results, especially if you want to pre-mask your input genome.

Earl Grey version 7.0.1 (latest stable release) with all required and configured dependencies is found in the `biooconda` conda channel. To install, simply run the following depending on your installation:
```
# With conda
conda create -n earlgrey -c conda-forge -c bioconda earlgrey=7.0.1

# With mamba
mamba create -n earlgrey -c conda-forge -c bioconda earlgrey=7.0.1

# Then run
earlGrey

# a script will be output to stdout and generated in the current directory to aid in setup
```

# Recommended Installation on ARM-based Mac Systems (M chips) Using Docker

It is possible to install x86-based Docker environments on M-chip Mac systems using Rosetta. The below form a guide that should work, but please reach out if you have any trouble!

Please first follow the Docker on Mac installation [instructions here](https://docs.docker.com/desktop/setup/install/mac-install/). Ensure you have installed Rosetta2, as this is required to get Earl Grey to behave as expected.

Next, we will create aliases to switch between arm and rosetta:

We want a simple way to activate the x86 architecture. We can do this by adding the following to `~/.zshrc`:
```
alias arm="env /usr/bin/arch -arm64 /bin/zsh --login"
alias intel="env /usr/bin/arch -x86_64 /bin/zsh --login"
```

Then close the terminal and open a new one.

To activate the intel/rosetta environment, run the following in a new terminal:
```
intel
```

Next, get the Docker installation and run in an interactive terminal following the instructions in the next section below. You will only need to pull the container once, then can use it for all your Earl Grey needs.

After this, you are ready to go! Just remember to activate the _intel_ terminal before starting the interactive container and running Earl Grey. 

# Docker Container 

A Docker container has been generated with none of Dfam 3.9, but with script generation to source required partitions

I try to keep an up-to-date container in docker hub, but this might not always be the case depending on if I have had time to build and upload a new image. Currently, the recommended image ready for use is `-nodfam` version. Upon running the container interactively and running the command `earlGrey`, instructions will print to `stdout` and a script that you can use will be placed in your current working directory. After an initial setup and configuration in an interative version of the container, you can commit the changes (i.e. the Dfam configuration) using `docker commit [container_ID] yourdockerusername/earlgrey:version7.0.1-configured`. Then, you can run this container interactively, or non-interatively, to annotate focal genomes.

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

# to quit the container and leave it running so you can commit the configured changes
# ctrl + p; ctrl + q

# get container ID to commit (not necessarily required if you leave the container running. If you kill the container and start again, you would need to reconfigure the DFAM libraries, so this is recommended. However, the container will be BIG depending on which parts of DFAM you pull.
docker ps -a

# commit the modified container so you can use at will (replace yourdockerusername with your docker username)
docker commit [container_ID] yourdockerusername/earlgrey:version7.0.1-configured

# you can then run non-interatively if required:
docker run -v 'pwd':/data/ yourdockerusername/earlgrey:version7.0.1-configured earlGrey -g /data/GENOME.fasta -s nonInteractiveTest -o /data/ -t 8

# alternatively you can still run interactive sessions
docker run -it -v 'pwd':/data/ yourdockerusername/earlgrey:version7.0.1-configured
``` 

