import os

configfile: "config/config.yaml"

# Configuration variables
GENOME = config["genome"]
SPECIES_LIST = config["species"]
OUTDIR = config["output_dir"]
THREADS = config["threads"]
REPSPEC = config.get("repeatmasker_species", None)
CONSENSUS_LIB = config["consensus_library"]
MARGIN = config.get("margin", "no")
SOFTMASK = config.get("softmask", "no")
# Handle both heliano and run_heliano config keys, convert boolean to yes/no
_heliano_val = config.get("heliano", config.get("run_heliano", False))
HELIANO = "yes" if (_heliano_val is True or _heliano_val == "yes") else "no"
SCRIPT_DIR = config["script_dir"]  # Path to EarlGrey scripts directory

# rule all:
#     input:
#         expand(
#             "{outdir}/{species}_EarlGrey/{species}_summaryFiles/{species}.highLevelCount.txt",
#             outdir=OUTDIR,
#             species=SPECIES_LIST
#         ),
#         expand(
#             "{outdir}/{species}_EarlGrey/{species}_summaryFiles/{species}.summaryPie.pdf",
#             outdir=OUTDIR,
#             species=SPECIES_LIST
#         ),
#         # Optional softmasked genome
#         expand(
#             "{outdir}/{species}_EarlGrey/{species}_summaryFiles/{species}.softmasked.fasta",
#             outdir=OUTDIR,
#             species=SPECIES_LIST
#         ) if SOFTMASK == "yes" else []

# rule prep_genome:
#     input:
#         genome=lambda wildcards: GENOME[wildcards.species]
#     output:
#         prep="{outdir}/{species}_EarlGrey/{species}.prep",
#         dict="{outdir}/{species}_EarlGrey/{species}.dict",
#         backup="{outdir}/{species}_EarlGrey/{species}.bak.gz"
#     params:
#         script_dir=SCRIPT_DIR
#     shell:
#         """
#         # Create output directory
#         mkdir -p {wildcards.outdir}/{wildcards.species}_EarlGrey
        
#         # Copy and backup original genome
#         cp {input.genome} {output.prep}.orig
#         cp {input.genome} {output.prep}.orig.bak && gzip -f {output.prep}.orig.bak
#         mv {output.prep}.orig.bak.gz {output.backup}
        
#         # Clean genome headers and create prep file
#         sed '/>/ s/[[:space:]].*//g; /^$/d' {output.prep}.orig > {output.prep}.tmp
#         {params.script_dir}/headSwap.sh -i {output.prep}.tmp -o {output.prep}
#         rm {output.prep}.tmp {output.prep}.orig
        
#         # Move dictionary file to correct location
#         mv {output.prep}.tmp.dict {output.dict}
        
#         # Replace ambiguous nucleotides
#         sed -i '/^>/! s/[DVHBPE]/N/g' {output.prep}
#         """

# These rules are commented out for pangenome pipeline - library creation handled by clustering.smk
# rule create_repeatmasker_library:
#     output:
#         replib="{outdir}/{species}_EarlGrey/{species}_Curated_Library/{repspec}.RepeatMasker.lib"
#     params:
#         script_dir=SCRIPT_DIR,
#         libpath=lambda wildcards: "$(which RepeatMasker | sed 's|bin/RepeatMasker|share/RepeatMasker/Libraries/famdb/|')" if REPSPEC else ""
#     shell:
#         """
#         if [[ $(which RepeatMasker) == *"bin"* ]]; then
#             libpath="$(which RepeatMasker | sed 's|bin/RepeatMasker|share/RepeatMasker/Libraries/famdb/|')"
#             export PATH=$PATH:"$(which RepeatMasker | sed 's|bin/RepeatMasker|share/RepeatMasker/|g')"
#         else
#             libpath="$(which RepeatMasker | sed 's|/[^/]*$||g')/Libraries/famdb/"
#         fi
#         
#         mkdir -p {wildcards.outdir}/{wildcards.species}_EarlGrey/{wildcards.species}_Curated_Library
#         famdb.py -i $libpath families -f fasta_name --include-class-in-name -a -d --curated {wildcards.repspec} > {output.replib}
#         """

# rule combine_libraries:
#     input:
#         consensus=CONSENSUS_LIB,
#         replib=lambda wildcards: f"{wildcards.outdir}/{wildcards.species}_EarlGrey/{wildcards.species}_Curated_Library/{REPSPEC}.RepeatMasker.lib" if REPSPEC else []
#     output:
#         combined="{outdir}/{species}_EarlGrey/{species}_Curated_Library/{species}_combined_library.fasta"
#     shell:
#         """
#         mkdir -p {wildcards.outdir}/{wildcards.species}_EarlGrey/{wildcards.species}_Curated_Library
#         if [ -f "{input.replib}" ]; then
#             cat {input.consensus} {input.replib} > {output.combined}
#         else
#             cp {input.consensus} {output.combined}
#         fi
#         """

rule repeatmasker_annotation:
    input:
        genome="{outdir}/{species}_EarlGrey/{species}.prep",
        library=f"{OUTDIR}/combinedLibraries/combined_all_species.clstrd.fa"
    output:
        masked="{outdir}/{species}_EarlGrey/{species}_RepeatMasker_Against_Custom_Library/{species}.prep.masked",
        out="{outdir}/{species}_EarlGrey/{species}_RepeatMasker_Against_Custom_Library/{species}.prep.out",
        tbl="{outdir}/{species}_EarlGrey/{species}_RepeatMasker_Against_Custom_Library/{species}.prep.tbl"
    params:
        outdir="{outdir}/{species}_EarlGrey/{species}_RepeatMasker_Against_Custom_Library",
        threads=lambda wildcards: int(THREADS/4)
    shell:
        """
        mkdir -p {params.outdir}
        cd {params.outdir}
        RepeatMasker -lib $(realpath {input.library}) -norna -no_is -lcambig -s -a -pa {params.threads} \
                     -dir {params.outdir} $(realpath {input.genome})
        """

rule heliano_detection:
    input:
        genome="{outdir}/{species}_EarlGrey/{species}.prep"
    output:
        helitron_gff="{outdir}/{species}_EarlGrey/{species}_heliano/RC.representative.gff"
    params:
        heliano_dir="{outdir}/{species}_EarlGrey/{species}_heliano",
        threads=THREADS
    shell:
        """
        if [ "{HELIANO}" == "yes" ]; then
            mkdir -p {params.heliano_dir}
            cd {params.heliano_dir}
            timestamp=$(date +"%Y%m%d_%H%M")
            heliano -g {input.genome} --nearest -dn 6000 -flank_sim 0.5 \
                    -o {params.heliano_dir}/HEL_$timestamp -w 10000 -n {params.threads}
            awk '{{OFS="\t"}}{{print $1, "HELIANO", "RC/Helitron", $2+1, $3, $5, $6, ".", "ID="$9"_"$11";shortTE=F"}}' \
                {params.heliano_dir}/HEL_$timestamp/RC.representative.bed > {output.helitron_gff}
        else
            mkdir -p {params.heliano_dir}
            touch {output.helitron_gff}
        fi
        """

rule merge_repeats:
    input:
        genome="{outdir}/{species}_EarlGrey/{species}.prep",
        dict="{outdir}/{species}_EarlGrey/{species}.dict",
        out="{outdir}/{species}_EarlGrey/{species}_RepeatMasker_Against_Custom_Library/{species}.prep.out",
        tbl="{outdir}/{species}_EarlGrey/{species}_RepeatMasker_Against_Custom_Library/{species}.prep.tbl",
        helitron_gff="{outdir}/{species}_EarlGrey/{species}_heliano/RC.representative.gff" if HELIANO == "yes" else []
    output:
        bed="{outdir}/{species}_EarlGrey/{species}_mergedRepeats/looseMerge/{species}.filteredRepeats.bed",
        gff="{outdir}/{species}_EarlGrey/{species}_mergedRepeats/looseMerge/{species}.filteredRepeats.gff",
        summary="{outdir}/{species}_EarlGrey/{species}_mergedRepeats/looseMerge/{species}.filteredRepeats.summary"
    params:
        script_dir=SCRIPT_DIR,
        outdir="{outdir}/{species}_EarlGrey/{species}_mergedRepeats/looseMerge",
        threads=THREADS,
        margin=MARGIN,
        helitron_param=lambda wildcards, input: f"-e {input.helitron_gff}" if HELIANO == "yes" and os.path.getsize(input.helitron_gff if HELIANO == "yes" else "/dev/null") > 0 else ""
    shell:
        """
        mkdir -p {params.outdir}
        
        # Try loose merge first
        {params.script_dir}/rcMergeRepeatsLoose -f {input.genome} -s {wildcards.species} \
            -d {params.outdir} -u {input.out} -q {input.tbl} -t {params.threads} \
            -b {input.dict} -m {params.margin} {params.helitron_param}
        
        # Fix GFF formatting
        if [ -f "{output.gff}" ]; then
            awk '{{OFS="\t"}}{{print $1, $2, $3, $4, $5, $6, $7, $8, toupper($9)}}' {output.gff} > {output.gff}.tmp
            mv {output.gff}.tmp {output.gff}
        fi
        
        # If loose merge failed, try strict merge
        if [ ! -f "{output.bed}" ]; then
            echo "Loose merge failed, trying strict merge..."
            {params.script_dir}/rcMergeRepeats -f {input.genome} -s {wildcards.species} \
                -d {wildcards.outdir}/{wildcards.species}_EarlGrey/{wildcards.species}_mergedRepeats \
                -u {input.out} -q {input.tbl} -t {params.threads} \
                -b {input.dict} -m {params.margin} {params.helitron_param}
            
            # Move strict merge results to expected location if loose merge failed
            if [ -f "{wildcards.outdir}/{wildcards.species}_EarlGrey/{wildcards.species}_mergedRepeats/{wildcards.species}.filteredRepeats.bed" ]; then
                mkdir -p {params.outdir}
                mv {wildcards.outdir}/{wildcards.species}_EarlGrey/{wildcards.species}_mergedRepeats/{wildcards.species}.filteredRepeats.* {params.outdir}/
            fi
        fi
        """

rule generate_summary_charts:
    input:
        summary="{outdir}/{species}_EarlGrey/{species}_mergedRepeats/looseMerge/{species}.filteredRepeats.summary",
        tbl="{outdir}/{species}_EarlGrey/{species}_RepeatMasker_Against_Custom_Library/{species}.prep.tbl"
    output:
        pie="{outdir}/{species}_EarlGrey/{species}_summaryFiles/{species}.summaryPie.pdf",
        highLevelCount="{outdir}/{species}_EarlGrey/{species}_summaryFiles/{species}.highLevelCount.txt"
    params:
        script_dir=SCRIPT_DIR,
        outdir="{outdir}/{species}_EarlGrey/{species}_summaryFiles"
    shell:
        """
        mkdir -p {params.outdir}
        cd {params.outdir}
        {params.script_dir}/autoPie.sh -i {input.summary} -t {input.tbl} \
                                       -p {output.pie} -o {output.highLevelCount}
        """

rule calculate_divergence:
    input:
        library=f"{OUTDIR}/combinedLibraries/combined_all_species.clstrd.fa",
        genome_orig=lambda wildcards: GENOME[wildcards.species],
        gff="{outdir}/{species}_EarlGrey/{species}_mergedRepeats/looseMerge/{species}.filteredRepeats.gff"
    output:
        div_gff="{outdir}/{species}_EarlGrey/{species}_RepeatLandscape/{species}.filteredRepeats.withDivergence.gff",
        div_summary="{outdir}/{species}_EarlGrey/{species}_summaryFiles/{species}_divergence_summary_table.tsv"
    params:
        script_dir=SCRIPT_DIR,
        landscape_dir="{outdir}/{species}_EarlGrey/{species}_RepeatLandscape",
        summary_dir="{outdir}/{species}_EarlGrey/{species}_summaryFiles",
        threads=THREADS
    shell:
        """
        mkdir -p {params.landscape_dir}
        cd {params.landscape_dir}
        
        # Calculate divergence
        python {params.script_dir}/divergenceCalc/divergence_calc.py \
            -l {input.library} -g {input.genome_orig} -i {input.gff} \
            -o {output.div_gff} -t {params.threads}
        
        # Generate divergence plots
        Rscript {params.script_dir}/divergenceCalc/divergence_plot.R \
            -s {wildcards.species} -g {output.div_gff} -o {params.landscape_dir}
        
        # Copy results to summary directory
        mkdir -p {params.summary_dir}
        cp {params.landscape_dir}/*.pdf {params.summary_dir}/ || true
        cp {params.landscape_dir}/*_summary_table.tsv {output.div_summary} || true
        
        # Update main GFF with divergence info (copy instead of move to keep output file)
        cp {output.div_gff} {input.gff}
        
        # Cleanup
        rm -rf {params.landscape_dir}/tmp/ || true
        """

rule sweep_up_files:
    input:
        bed="{outdir}/{species}_EarlGrey/{species}_mergedRepeats/looseMerge/{species}.filteredRepeats.bed",
        gff="{outdir}/{species}_EarlGrey/{species}_mergedRepeats/looseMerge/{species}.filteredRepeats.gff",
        highLevelCount="{outdir}/{species}_EarlGrey/{species}_summaryFiles/{species}.highLevelCount.txt",
        pie="{outdir}/{species}_EarlGrey/{species}_summaryFiles/{species}.summaryPie.pdf",
        div_summary="{outdir}/{species}_EarlGrey/{species}_summaryFiles/{species}_divergence_summary_table.tsv"
    output:
        summary_bed="{outdir}/{species}_EarlGrey/{species}_summaryFiles/{species}.filteredRepeats.bed",
        summary_gff="{outdir}/{species}_EarlGrey/{species}_summaryFiles/{species}.filteredRepeats.gff"
    params:
        summary_dir="{outdir}/{species}_EarlGrey/{species}_summaryFiles"
    shell:
        """
        # Copy final results to summary directory
        cp {input.bed} {output.summary_bed}
        cp {input.gff} {output.summary_gff}
        """

rule generate_softmasked_genome:
    input:
        bed="{outdir}/{species}_EarlGrey/{species}_summaryFiles/{species}.filteredRepeats.bed",
        backup="{outdir}/{species}_EarlGrey/{species}.bak.gz"
    output:
        softmasked="{outdir}/{species}_EarlGrey/{species}_summaryFiles/{species}.softmasked.fasta"
    shell:
        """
        if [ "{SOFTMASK}" == "yes" ]; then
            gunzip -c {input.backup} > {input.backup}.tmp
            bedtools maskfasta -fi {input.backup}.tmp -bed {input.bed} \
                              -fo {output.softmasked} -soft
            rm {input.backup}.tmp
        else
            touch {output.softmasked}
        fi
        """

# Commented out for pangenome pipeline
# # Handle optional rules based on configuration
# if REPSPEC:
#     # If RepeatMasker species is specified, create the library first
#     ruleorder: create_repeatmasker_library > combine_libraries
# else:
#     # If no RepeatMasker species, just use consensus library
#     rule combine_libraries_no_repspec:
#         input:
#             consensus=CONSENSUS_LIB
#         output:
#             combined="{outdir}/{species}_EarlGrey/{species}_Curated_Library/{species}_combined_library.fasta"
#         shell:
#             """
#             mkdir -p {wildcards.outdir}/{wildcards.species}_EarlGrey/{wildcards.species}_Curated_Library
#             cp {input.consensus} {output.combined}
#             """