rule cluster_all_species:
    input:
        expand("{outdir}/{species}_EarlGrey/{species}_summaryFiles/{species}-families.fa.strained", outdir=OUTDIR, species=SPECIES_LIST)
    output:
        combined="results/combined_all_species.fa"
    params:
        outdir=OUTDIR
    shell:
        """
        cat {input} > {output.combined}
        cd-hit-est -i {output.combined} -o {output.combined}.clstr -c 0.8 -n 5 -M 16000 -T {THREADS}
        cd-hit-est \
           -d 0 -aS 0.8 -c 0.8 -G 0 -g 1 -b 500 -r 1 \
           -i {output.combined} -o ${output.combined}.clstrd.fa
        """

# TODO: in -r and -L options, add previous consensus to clustering so that the new library is built on both new and old consensus sequences