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