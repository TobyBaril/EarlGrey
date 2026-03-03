import os

GENOME = config["genome"]
SPECIES_LIST = config["species"]
OUTDIR = config["output_dir"] #os.path.join(config["output_dir"], f"{SPECIES}_EarlGrey")
THREADS = config["threads"]
REPSPEC = config["repeatmasker_species"]
CUSTOM_LIB = config["custom_library"]
ITER = config["iterations"]
FLANK = config["flank"]
MAX_SEQ = config["max_consensus_seqs"]
MIN_SEQ = config["min_consensus_seqs"]
HELI = config.get("run_heliano", False)
SCRIPT_DIR = config["script_dir"]

# Rule priority: if CUSTOM_LIB is set, use repeatmasker_custom; otherwise use repeatmasker
if CUSTOM_LIB:
    ruleorder: repeatmasker_custom > repeatmasker
elif REPSPEC:
    ruleorder: repeatmasker > repeatmasker_custom

def get_masked_genome_input(wildcards):
    """Return appropriate input based on config settings"""
    if REPSPEC:
        return f"{wildcards.outdir}/{wildcards.species}_EarlGrey/{wildcards.species}_RepeatMasker/{wildcards.species}.masked"
    elif CUSTOM_LIB:
        return f"{wildcards.outdir}/{wildcards.species}_EarlGrey/{wildcards.species}_RepeatMasker/{wildcards.species}.masked"
    else:
        # No masking needed, use prep genome directly
        return f"{wildcards.outdir}/{wildcards.species}_EarlGrey/{wildcards.species}.prep"

rule prep_genome:
    input:
        genome=lambda wildcards: GENOME[wildcards.species]
    output:
        gen_prep="{OUTDIR}/{species}_EarlGrey/{species}.prep",
        gen_dict="{OUTDIR}/{species}_EarlGrey/{species}.dict"
    params:
        script_dir=SCRIPT_DIR
    shell:
        """
        cp {input.genome} {input.genome}.bak && gzip -f {input.genome}.bak
        sed '/>/ s/[[:space:]].*//g; /^$/d' {input.genome} > {input.genome}.tmp
        {params.script_dir}/headSwap.sh -i {input.genome}.tmp -o {output.gen_prep} && rm {input.genome}.tmp
        mv {input.genome}.tmp.dict {output.gen_dict}
        sed -i.bak '/^>/! s/[DVHBPE]/N/g' {output.gen_prep}

        """

def get_masked_genome_input(wildcards):
    """Return appropriate input based on config settings"""
    if REPSPEC:
        return f"{wildcards.outdir}/{wildcards.species}_EarlGrey/{wildcards.species}_RepeatMasker/{wildcards.species}.masked"
    elif CUSTOM_LIB:
        return f"{wildcards.outdir}/{wildcards.species}_EarlGrey/{wildcards.species}_RepeatMasker/{wildcards.species}.masked"
    else:
        # No masking needed, use prep genome directly
        return f"{wildcards.outdir}/{wildcards.species}_EarlGrey/{wildcards.species}.prep"

rule repeatmasker:
    input:
        genome="{outdir}/{species}_EarlGrey/{species}.prep"
    output:
        masked="{outdir}/{species}_EarlGrey/{species}_RepeatMasker/{species}.masked"
    params:
        outdir=OUTDIR,
        rep_spec=REPSPEC,
        threads=int(THREADS/4)
    shell:
        """
        RepeatMasker \
           -species {params.rep_spec} \
           -norna -no_is -lcambig -s -a -pa {params.threads} \
           -dir {params.outdir}/{wildcards.species}_EarlGrey/{wildcards.species}_RepeatMasker \
           {input.genome}
        """

rule repeatmasker_custom:
    input:
        genome="{outdir}/{species}_EarlGrey/{species}.prep",
        lib=CUSTOM_LIB
    output:
        masked="{outdir}/{species}_EarlGrey/{species}_RepeatMasker/{species}.masked"
    params:
        threads=int(THREADS/4),
        outdir=OUTDIR
    shell:
        """
        RepeatMasker \
            -lib {input.lib} \
            -norna -no_is -lcambig -s -a -pa {params.threads} \
            -dir {params.outdir}/{wildcards.species}_RepeatMasker \
            {input.genome}
        """

rule extract_repeatmasker_library:
    output:
        replib=f"{{outdir}}/{REPSPEC}.RepeatMasker.lib"
    params:
        repspec=REPSPEC,
        outdir=OUTDIR
    shell:
        """
        # Determine RepeatMasker library path
        if [[ $(which RepeatMasker) == *"bin"* ]]; then
            libpath="$(which RepeatMasker | sed 's|bin/RepeatMasker|share/RepeatMasker/Libraries/famdb/|')"
            export PATH=$PATH:"$(which RepeatMasker | sed 's|bin/RepeatMasker|share/RepeatMasker/|g')"
        else
            libpath="$(which RepeatMasker | sed 's|/[^/]*$||g')/Libraries/famdb/"
        fi
        
        # Create output directory
        mkdir -p {params.outdir}
        
        # Extract RepeatMasker library for specified species/clade
        famdb.py -i $libpath families -f fasta_name --include-class-in-name -a -d --curated {params.repspec} > {output.replib}
        """

rule build_db:
    input:
        masked=get_masked_genome_input
    output:
        db="{outdir}/{species}_EarlGrey/{species}_Database/{species}.nhr",
        nin="{outdir}/{species}_EarlGrey/{species}_Database/{species}.nin",
        nsq="{outdir}/{species}_EarlGrey/{species}_Database/{species}.nsq"
    params:
        outdir="{outdir}/{species}_EarlGrey/{species}_Database",
        name="{species}"
    shell:
        """
        mkdir -p {params.outdir}
        cd {params.outdir}
        BuildDatabase -name {params.name} {input.masked}
        """

rule repeatmodeler:
    input:
        db="{outdir}/{species}_EarlGrey/{species}_Database/{species}.nhr",
        nin="{outdir}/{species}_EarlGrey/{species}_Database/{species}.nin",
        nsq="{outdir}/{species}_EarlGrey/{species}_Database/{species}.nsq"
    output:
        families="{outdir}/{species}_EarlGrey/{species}_Database/{species}-families.fa"
    params:
        db_dir="{outdir}/{species}_EarlGrey/{species}_Database",
        db_name="{species}",
        threads=THREADS
    shell:
        """
        cd {params.db_dir}
        RepeatModeler -threads {params.threads} -database {params.db_name}
        """

rule testrainer:
    input:
        genome="{outdir}/{species}_EarlGrey/{species}.prep",
        families="{outdir}/{species}_EarlGrey/{species}_Database/{species}-families.fa"
    output:
        strained="{outdir}/{species}_EarlGrey/{species}_strainer/{species}-families.fa.strained",
        summary="{outdir}/{species}_EarlGrey/{species}_summaryFiles/{species}-families.fa.strained"
    params:
        outdir=OUTDIR,
        threads=int(THREADS/4),
        flank=FLANK,
        iter=ITER,
        max_seq=MAX_SEQ,
        min_seq=MIN_SEQ,
        script_dir=SCRIPT_DIR,
        strainer_dir="{outdir}/{species}_EarlGrey/{species}_strainer"
    shell:
        """
        mkdir -p {params.strainer_dir}
        cd {params.strainer_dir}
        {params.script_dir}/TEstrainer/TEstrainer_for_earlGrey.sh \
           -g {input.genome} -l {input.families} \
           -t {params.threads} -f {params.flank} \
           -r {params.iter} -n {params.max_seq} \
           -m {params.min_seq}

        # Find and copy the latest TEstrainer output from subdirectory
        latestDir=$(ls -td {params.strainer_dir}/*/ 2>/dev/null | head -n 1)
        if [ -n "$latestDir" ]; then
            latestFile="${{latestDir}}{wildcards.species}-families.fa.strained"
            if [ -f "$latestFile" ]; then
                cp "$latestFile" {output.strained}
            fi
        fi

        # Add species name to fasta headers
        sed -i.bak "s/>/>{wildcards.species}_/g" {output.strained}

        mv {output.strained}.bak {output.summary}
        """

