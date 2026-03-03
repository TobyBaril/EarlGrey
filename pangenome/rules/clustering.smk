rule cluster_all_species:
    input:
        strained=expand("{outdir}/{species}_EarlGrey/{species}_summaryFiles/{species}-families.fa.strained", 
                       outdir=OUTDIR, species=SPECIES_LIST),
        replib=f"{OUTDIR}/{REPSPEC}.RepeatMasker.lib" if REPSPEC else [],
        custom=CUSTOM_LIB if CUSTOM_LIB else []
    output:
        combined="{outdir}/combinedLibraries/combined_all_species.clstrd.fa"
    params:
        threads=THREADS,
        outdir=OUTDIR,
        species=SPECIES_LIST,
        repspec=REPSPEC
    run:
        import os
        shell("mkdir -p {params.outdir}/combinedLibraries")
        
        # Define temporary input file path
        tmp_file = f"{params.outdir}/combinedLibraries/combined_all_species.fa"
        
        # Create temporary file with genome-prefixed sequences
        with open(tmp_file, 'w') as outf:
            # Add strained sequences with genome prefix
            for i, species in enumerate(params.species):
                with open(input.strained[i], 'r') as inf:
                    for line in inf:
                        if line.startswith('>'):
                            outf.write(f">{species}_{line[1:]}")
                        else:
                            outf.write(line)
            
            # Add RepeatMasker library if present with REPMASKER_{species} prefix
            if input.replib and os.path.exists(input.replib):
                with open(input.replib, 'r') as inf:
                    for line in inf:
                        if line.startswith('>'):
                            outf.write(f">REPMASKER_{params.repspec}_{line[1:]}")
                        else:
                            outf.write(line)
            
            # Add custom library if present with CUSTOM prefix
            if input.custom and os.path.exists(input.custom):
                with open(input.custom, 'r') as inf:
                    for line in inf:
                        if line.startswith('>'):
                            outf.write(f">CUSTOM_{line[1:]}")
                        else:
                            outf.write(line)
        
        # Run single cd-hit-est clustering
        shell(f"cd-hit-est -d 0 -aS 0.8 -c 0.8 -G 0 -g 1 -b 500 -r 1 "
              f"-i {tmp_file} -o {{output.combined}} "
              f"-M 16000 -T {{params.threads}}")
        
        # Clean up intermediate files
        shell(f"rm -f {tmp_file}")