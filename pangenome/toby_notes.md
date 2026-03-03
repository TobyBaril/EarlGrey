# Development notes for the pangenome pipeline of Earlgrey

### February 2026
We decided to get rid of the EarlGrey set up steps, such a the DFAM configuration etc, in the pangenome pipeline. Before using the pangenome pipeline, a user would have to run the normal installation and set up script. So there is no need to repeat the steps here. The only thing left to do in the pangenome pipeline is to carry out the checks. 

The pipeline expects a config.yaml file called `config/config.yaml`. This file contains the necessary information to run the pipeline such as the species and their genome files. It can also contain some additional options, but for some defaults will be applied (see the validate_parameters function for more details).

I am testing with snakemake version 9.9.0. 

`snakemake --cores 1 --dag ` will run the parameter checks and create the dag. 

Now copy all the text which is not part of the check messages into a `temp.txt` file and use the following command to create the DAG plot. 

`cat temp.txt | dot -Tsvg > dag.svg` 

If you want a diagram that represents a summary of the DAG, independently of the different samples, you can use `--rulegraph` instead. 

`snakemake --cores 1` to run the pipeline.

### March 2026 - Toby Testing and Development Notes

I will first build a mamba env with earlgrey 7.0.3 & snakemake

```bash
mamba create -c conda-forge -c bioconda -n earlgrey-pan-dev earlgrey=7.0.3 snakemake
```

Configure earl grey properly by running `earlGrey`

```bash
earlGrey 
# edited /data/toby/EarlGrey_pangenome$ bash configure_dfam39.sh to only be part 0
bash configure_dfam39.sh
```

Now, make a test directory with two test genomes:

```bash
mkdir -p /data/toby/EarlGrey_pangenome/test
cd /data/toby/EarlGrey_pangenome/test

# get some test genomes
cp /data/toby/tools/test.fasta /data/toby/EarlGrey_pangenome/test/genome1.fasta
cp /data/toby/tools/test.fasta /data/toby/EarlGrey_pangenome/test/genome2.fasta
```

Edit the config file to point to these two genomes:

```yaml
genome: 
  genome1: "/data/toby/EarlGrey_pangenome/test/genome1.fasta"
  genome2: "/data/toby/EarlGrey_pangenome/test/genome2.fasta"
species: ["genome1", "genome2"]
output_dir: "/data/toby/EarlGrey_pangenome/test/results"
threads: 8
repeatmasker_species: ""  # e.g. "arthropoda" or "" for none
custom_library: ""        # path to custom library or "" for none
iterations: 10
flank: 1000
max_consensus_seqs: 20
min_consensus_seqs: 3
consensus_library : "results/combined_all_species.fa"
script_dir: "/data/toby/miniforge3/envs/earlgrey-pan-dev/share/earlgrey-7.0.3-0/scripts"
run_heliano: False
```

I saved the original as: `/data/toby/EarlGrey_pangenome/pangenome/config/config.yaml.bak`

Try and run the pipeline:

```bash
cd /data/toby/EarlGrey_pangenome/pangenome
snakemake --cores 1 --dag 2>&1 | sed -n '/^digraph/,$p' > /data/toby/EarlGrey_pangenome/test/temp.txt
cat /data/toby/EarlGrey_pangenome/test/temp.txt | dot -Tsvg > /data/toby/EarlGrey_pangenome/test/dag.svg
```

This should create a DAG plot in the test directory.

Now, run the pipeline with a dry run to check for any errors:

```bash
snakemake --cores 1 --dry-run
```

Now, run the pipeline for real:

```bash
snakemake --cores 1
```

**Issue Fixed (March 3, 2026):** RepeatMasker was running even when `repeatmasker_species` and `custom_library` were not specified, causing errors because it was called with empty `-species` parameter.

**Solution:** Modified `rules/lib_construct.smk` to:
1. Added `get_masked_genome_input()` function that returns the appropriate input:
   - Masked genome if `repeatmasker_species` or `custom_library` is specified
   - Prep genome directly if neither is specified
2. Updated `build_db` rule to use this conditional input function
3. Added `ruleorder` directive to prioritize the correct repeatmasker rule when needed

Now the workflow correctly skips the initial RepeatMasker step when no masking parameters are provided, going directly from `prep_genome` to `build_db`.

Now, try another dry run to confirm the issue is resolved:

```bash
snakemake --cores 1 --dry-run
```

Now, run the pipeline again to ensure it executes without errors:

```bash
snakemake --cores 1
```

This failed because I think the pipeline is using the wrong earl grey installation. it should be using the version in the active conda env, not in the EarlGrey_pangenome directory.

**Issue Fixed (March 3, 2026):** The pipeline was using hardcoded script paths pointing to the local repository instead of the conda environment's EarlGrey installation.

**Solution:** Modified `rules/lib_construct.smk` to:
1. Added `script_dir` parameter to `prep_genome` rule and changed the headSwap.sh call from `../scripts/headSwap.sh` to `{params.script_dir}/headSwap.sh`
2. Added `script_dir` parameter to `testrainer` rule and changed the TEstrainer call from `scripts/TEstrainer/TEstrainer_for_earlGrey.sh` to `{params.script_dir}/TEstrainer/TEstrainer_for_earlGrey.sh`
3. Also updated `/data/toby/EarlGrey_pangenome/scripts/headSwap.sh` to dynamically determine its own directory instead of hardcoding it (though this won't be used when the conda env version is used)

Now the pipeline correctly uses the scripts from the conda environment specified in `config["script_dir"]`.

Clean previous outputs and try again:

```bash
snakemake --cores 1 --delete-all-output
snakemake --cores 1
```

**Issue Fixed (March 3, 2026):** BuildDatabase was creating files in the wrong directory and the output specification was incorrect.

**Problem:** BuildDatabase creates multiple files (`.nhr`, `.nin`, `.nsq`, `.nnd`, `.nni`, `.nog`, `.translation`, `.njs`) in the current working directory, not a single `.db` file.

**Solution:** Modified `rules/lib_construct.smk` to:
1. Updated `build_db` rule to output the actual files BuildDatabase creates (`.nhr`, `.nin`, `.nsq`) instead of a non-existent `.db` file
2. Added `cd {params.outdir}` to change to the output directory before running BuildDatabase
3. Updated `repeatmodeler` rule to:
   - Accept the database files as input instead of a single `.db` file
   - Change to the database directory before running RepeatModeler
   - Pass just the database name (not full path) to RepeatModeler

Now the workflow correctly creates database files in the expected output directories and RepeatModeler can find them.

The pipeline is now working up to the clustering point. So far, I have added `snakemake`.

At the moment, there are too many conflicts to add `mmseqs2` to the mamba env, but this would eventually be preferred. Thus, I will continue with clustering with cd-hit for now.

---

## March 3, 2026 - Major Pipeline Refinement Session

### Summary of Changes

Today I worked through a comprehensive debugging and refinement session to make the pangenome pipeline fully functional. Below are all the major changes made:

#### 1. Clustering Implementation
- **Added genome identifiers to sequences**: Modified `clustering.smk` to prefix each sequence header with the genome name (e.g., `>genome1_family-13`) to track sequence origins
- **Added library source prefixes**: RepeatMasker libraries get `REPMASKER_{species}_` prefix (e.g., `>REPMASKER_fungi_TE_name#DNA/TcMar-Tc1`), custom libraries get `CUSTOM_` prefix. This makes it easy to identify both the source and the specific RepeatMasker species/clade used.
- **Fixed library extraction**: RepeatMasker species/clade library is now extracted once to `{OUTDIR}/{REPSPEC}.RepeatMasker.lib` instead of per-genome
- **Conditional library inclusion**: Clustering now conditionally includes RepeatMasker species libraries and/or custom libraries when specified in config
- **Single cd-hit-est run**: Changed from two-pass clustering to single pass with parameters: `-d 0 -aS 0.8 -c 0.8 -G 0 -g 1 -b 500 -r 1`
- **Fixed output location**: Clustered library outputs to `{OUTDIR}/combinedLibraries/combined_all_species.clstrd.fa`
- **Added cleanup**: Automatically removes temporary concatenated file and `.clstr` intermediate files

#### 2. RepeatMasker Annotation Fixes
- **Fixed library path**: Updated `annotate_simple.smk` to use absolute path with `realpath` for library input
- **Fixed working directory**: RepeatMasker now runs with `cd` into output directory and `-dir` parameter to ensure outputs go to correct location
- **Library reference**: Changed to use clustered pangenome library instead of per-genome libraries

#### 3. Switch to Full Annotation Pipeline
- **Activated annotate.smk**: Changed `Snakefile` to include `rules/annotate.smk` instead of `annotate_simple.smk`
- **Updated library references**: Modified `repeatmasker_annotation`, `calculate_divergence`, and `sweep_up_files` rules to use pangenome clustered library
- **Commented out obsolete rules**: Disabled `create_repeatmasker_library` and `combine_libraries` rules that are redundant in pangenome approach
- **Fixed HELIANO variable**: Added proper handling for both `heliano` and `run_heliano` config keys with boolean/string conversion

#### 4. Summary Charts and Divergence
- **Fixed autoPie.sh inputs**: Changed from `.bed` file to `.summary` file as required by the script
- **Added `.summary` output**: Updated `merge_repeats` rule to declare `.filteredRepeats.summary` as an output
- **Added divergence to workflow**: Updated `rule all` to request `{species}_divergence_summary_table.tsv` outputs
- **Fixed divergence file handling**: Changed from `mv` to `cp` for withDivergence.gff to keep output file in expected location while updating main GFF
- **Added divergence dependency**: Made `sweep_up_files` depend on divergence summary to ensure proper execution order

#### 5. Output Organization
- **Removed duplicate libraries**: Eliminated per-genome copies of pangenome library in summary directories - now only exists in `combinedLibraries/`
- **Single pangenome library**: Final clustered library at `{OUTDIR}/combinedLibraries/combined_all_species.clstrd.fa` is the single source of truth

#### 6. Bug Fixes and Refinements
- **Fixed variable naming**: Changed `{output.count}` to `{output.highLevelCount}` in autoPie.sh call
- **Removed duplicate closing**: Fixed syntax error from duplicate `"""` in sweep_up_filesrule
- **Disabled make_directories**: Commented out `make_directories()` call in onstart that was creating incorrectly-named directory with literal list string
- **Fixed REPSPEC usage**: Changed from treating `{repspec}` as wildcard to using `REPSPEC` as constant throughout pipeline
- **Added missing HELI variable**: Added `HELI = config.get("run_heliano", False)` definition to Snakefile

#### 7. Workflow Execution
- **Updated rule all**: Now requests final summary files including divergence outputs, not just intermediate RepeatMasker files
- **Proper dependencies**: Ensured correct execution order through input/output dependencies

### Pipeline Architecture

The final pangenome pipeline flow:
1. **Prep genomes** → Create cleaned, indexed genome files
2. **Build databases** → Create BLAST databases for RepeatModeler
3. **RepeatModeler** → De novo repeat discovery with fallback for classification failures
4. **TEstrainer** → Refine consensus sequences
5. **Extract libraries** → Extract RepeatMasker species library if specified (once, not per-genome)
6. **Clustering** → Combine all strained sequences with genome prefixes, plus RepeatMasker/custom libraries with source prefixes, cluster with cd-hit-est
7. **Annotation** → Annotate each genome with pangenome clustered library
8. **Merge repeats** → Merge overlapping repeat annotations
9. **Summary charts** → Generate pie charts and counts
10. **Divergence calculation** → Calculate divergence and generate landscape plots
11. **Sweep up** → Copy final results to summary directories

### Key Files Modified
- `Snakefile` - Changed to use annotate.smk, updated rule all, disabled make_directories
- `rules/lib_construct.smk` - Added extract_repeatmasker_library rule, fixed script paths
- `rules/clustering.smk` - Complete rewrite with genome prefixes, library prefixes, single cd-hit run
- `rules/annotate.smk` - Updated for pangenome library, fixed all input/output paths, added .summary output
- `config/config.yaml` - Configured for test genomes with proper script_dir path

### Testing
Successfully ran full pipeline on two test genomes (genome1.fasta and genome2.fasta, both copies of test.fasta) with all stages completing:
- Library construction: 40 sequences per genome after TEstrainer
- Clustering: Combined 80 sequences → 35 unique sequences in pangenome library
- Annotation: Both genomes annotated with pangenome library
- Divergence: Landscape plots and summary tables generated
- Final outputs: All summary files, charts, and annotations created successfully
