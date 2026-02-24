import os

configfile: "config/config.yaml"

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

rule check_and_configure_env:
    params:
        outdir = OUTDIR,
        species = SPECIES_LIST,
        repspec = REPSPEC,
        startCust = CUSTOM_LIB,
        heli = HELI,
        script_dir = SCRIPT_DIR
    output:
        envchecked = touch("{outdir}/.envchecked")
    run:
        # Validate parameters
        validated_config = validate_parameters(config)
        
        # Check directories
        check_script_directories(script_dir)

        # Check biocontainer
        check_biocontainer()
        
        # Check Dfam 3.9
        check_dfam39()
        
        print("All startup checks completed successfully!")

        print("Making output directories...")
        make_directories(outdir, species, RepSpec=repspec, startCust=startCust, heli=heli)

rule prep_genome:
    input:
        genome=lambda wildcards: GENOME[wildcards.species],
        dfam_ready = r"{outdir}/.envchecked" 
    params:
        outdir = OUTDIR
    output:
        gen_prep="{outdir}/{species}_EarlGrey/{species}.prep",
        gen_dict="{outdir}/{species}_EarlGrey/{species}.dict"
    shell:
        """
        cp {input.genome} {input.genome}.bak && gzip -f {input.genome}.bak
        sed '/>/ s/[[:space:]].*//g; /^$/d' {input.genome} > {input.genome}.tmp
        scripts/headSwap.sh -i {input.genome}.tmp -o {output.gen_prep} && rm {input.genome}.tmp
        mv {input.genome}.tmp.dict {output.gen_dict}
        sed -i.bak '/^>/! s/[DVHBPE]/N/g' {output.gen_prep}

        """

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
           -dir {params.outdir}/{species}_EarlGrey/{wildcards.species}_RepeatMasker \
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

rule build_db:
    input:
        masked="{outdir}/{species}_EarlGrey/{species}_RepeatMasker/{species}.masked"
    output:
        db="{outdir}/{species}_EarlGrey/{species}_Database/{species}.db"
    params:
        outdir=OUTDIR
    shell:
        """
        BuildDatabase -name {wildcards.species} {input.masked}
        """

rule repeatmodeler:
    input:
        db="{outdir}/{species}_EarlGrey/{species}_Database/{species}.db"
    output:
        families="{outdir}/{species}_EarlGrey/{species}_Database/{species}-families.fa"
    params:
        outdir=OUTDIR,
        threads=THREADS
    shell:
        """
        RepeatModeler -threads {params.threads} -database {input.db}
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
        min_seq=MIN_SEQ
    shell:
        """
        scripts/TEstrainer/TEstrainer_for_earlGrey.sh \
           -g {input.genome} -l {input.families} \
           -t {params.threads} -f {params.flank} \
           -r {params.iter} -n {params.max_seq} \
           -m {params.min_seq}

        # Add species name to fasta headers
        sed -i.bak "s/>/>{wildcards.species}_/g" {output.strained}

        mv {output.strained}.bak {output.summary}
        """

