import os
from os.path import exists, getsize
import sys
import argparse
import pandas as pd
import multiprocessing
import pybedtools
import subprocess
import shlex
import shutil
from Bio import AlignIO, SeqIO
from math import log, sqrt
from functools import partial
from time import time
from re import sub

parser = argparse.ArgumentParser()
parser.add_argument('-l', '--repeat_library', type=str, required=True,
                    help='repeat_library')
parser.add_argument('-i', '--in_gff', type=str, required=True,
                    help='Path to gff')
parser.add_argument('-g', '--genome', type=str, required=True,
                    help='Path to genome')
parser.add_argument('-o', '--out_gff', type=str, required=True,
                    help='Output gff')
parser.add_argument('-tmp', '--temp_dir', type=str, default='tmp/',
                    help='Temporary directory')
parser.add_argument('-t', '--cores', type=int, default=4,
                    help='Number of cores')
parser.add_argument('-k', '--timeout', type=int, default=30,
                    help='Seconds after which matcher will be cancelled and repeat treated as unalignable')

args = parser.parse_args()

def file_check(repeat_library, in_gff, genome, out_gff, temp_dir):
    if(exists(repeat_library) == False or exists(in_gff) == False or exists(genome) == False):
        sys.exit('Files not found. Requires the repeat library, path to the genome, and path to gff containing coordinates and corresponding repeat files')
    if(exists(genome+".fai") == False):
        print("Indexing genome")
        subprocess.run(["samtools","faidx",genome], stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)
    if(exists(temp_dir) == False):
        os.mkdir(temp_dir)
    if(exists(temp_dir+"/qseqs") == False):
        os.mkdir(temp_dir+"/qseqs")
    if(exists(temp_dir+"/split_library/") == False):
        os.mkdir(temp_dir+"/split_library/")
        

def splitter(in_seq, temp_dir):
    with open(in_seq, 'r') as handle:
        for record in SeqIO.parse(handle, "fasta"):
            repeat_name = record.name.split(sep="#")[0]
            repeat_name = repeat_name.lower()
            file_name = (temp_dir+"/split_library/"+repeat_name+".fasta")
            SeqIO.write(record, file_name, "fasta-2line")

def parse_gff(in_gff):
    gff = pd.read_table(in_gff, header = None, names=['seqnames', 'tool', 'repeat_class', 'start', 'end', 'score', 'strand', 'phase', 'metadata'])
    simple_gff = gff[gff['repeat_class'].str.contains('Simple_repeat|Satellite|Low_complexity')].reset_index()
    gff = gff[~gff['repeat_class'].str.contains('Simple_repeat|Satellite|Low_complexity')].reset_index()
    other_gff = gff[~gff['tool'].str.contains('Earl_Grey|RepeatMasker')].reset_index()
    gff = gff[gff['tool'].str.contains('Earl_Grey|RepeatMasker')].reset_index()
    gff['metadata_tmp'] = gff['metadata'].str.replace(';TSTART.*', '', regex=True)
    gff[['repeat_id', 'repeat_family']] = gff['metadata_tmp'].str.split(';', n=2, expand=True)
    gff = gff.drop(columns = ['metadata_tmp', 'repeat_id'])
    gff['repeat_family'] = gff['repeat_family'].str.replace('NAME=', '', regex=True)
    gff['repeat_family'] = gff['repeat_family'].str.lower()
    return(gff, simple_gff, other_gff)

def file_name_generator():
    import random
    import string
    file_name = ''.join(random.sample(string.ascii_letters, 12))+'.tmp'
    return(file_name)

def Kimura80(qseq, sseq):
    """
    Calculations adapted from https://github.com/kgori/python_tools_on_github/blob/master/pairwise_distances.py
    """
    # define transitions, transversions, matches
    transitions = [ "AG", "GA", "CT", "TC"]
    transversions = [ "AC", "CA", "AT", "TA",
                    "GC", "CG", "GT", "TG" ]
    matches = [ "AA", "GG", "CC", "TT"]
    # set counters to 0
    m,ts,tv=0,0,0
    # count transitions, transversions, matches
    for i, j in zip(qseq, sseq):
        if i+j in matches: m+=1
        if i+j in transitions: ts+=1
        if i+j in transversions: tv+=1
    # count number of bp which align (excludes gaps, Ns)
    aln_len = m + ts + tv
    
    if aln_len != 0:
        # calculate p and q 
        p = ts/aln_len
        q = tv/aln_len
    
        # calculate Kimura distance
        Kimura_dist = -0.5 * log((1 - 2*p - q) * sqrt( 1 - 2*q ))
    else:
        Kimura_dist = "NA"
    
    return(Kimura_dist)

def outer_func(genome_path, temp_dir, timeoutSeconds, gff):
    generated_name = file_name_generator()
    holder_file_name = temp_dir+generated_name
    failed_file_name = temp_dir+"failed_"+generated_name
    with open(holder_file_name, 'w') as tmp_out:
        header = list(gff.columns.values)[1:] + ["Kimura"]
        header = "\t".join(header)+"\n"
        tmp_out.write(header)
        for row in gff.iterrows():
            # Set index
            idx = row[0]
            # Set scaffold, coordinates, strand, repeat family
            seqnames, start, end, strand, repeat_family = row[1]['seqnames'], str(row[1]['start'] - 1), str(row[1]['end']), row[1]['strand'], row[1]['repeat_family']
            # Create BED string for BEDtools
            bed_str = " ".join([seqnames, start, end, ".", ".", strand])
            # Set path for query sequence
            query_path = temp_dir+"/qseqs/"+str(idx)
            # Create bedtools command and getfasta
            a=pybedtools.BedTool(bed_str, from_string=True)
            try:
                a = a.sequence(fi=genome_path, fo=query_path, s=True)
            except: # Occasionally a samtools error occurs, this overcomes this
                with open(failed_file_name, "a") as failed_file:
                    failed_file.write(seqnames+":"+start+"-"+end+"_"+strand+"_"+repeat_family+"\n")
            if exists(query_path) is True and getsize(query_path) > 0:
                # Set path to subject sequence
                subject_path=temp_dir+"/split_library/"+repeat_family+".fasta"
                # Run matcher, with timeout exception
                test_command = shlex.split("matcher "+query_path+" "+subject_path+" -outfile "+query_path+".matcher -aformat fasta")
                # Run test and kill if it takes more than 10 seconds
                alignment_p = subprocess.Popen(test_command, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)
                try:
                    alignment_p.wait(timeoutSeconds)
                except subprocess.TimeoutExpired:
                    # if matcher fails to complete before timeout, kill and move on
                    with open(failed_file_name, "a") as failed_file:
                        failed_file.write(seqnames+":"+start+"-"+end+"_"+strand+"_"+repeat_family+"\n")
                    alignment_p.kill()

            if exists(query_path+".matcher") is False or getsize(query_path+".matcher") == 0:
                # If no alignment is possible, set distances to NA and alignment length to 0
                Kdist = "NA"
                os.remove(query_path)
                if exists(query_path+".matcher") is True:
                    os.remove(query_path+".matcher")
            else:
                # Read in alignments
                aln = list(SeqIO.parse(query_path+".matcher", 'fasta'))
                ref_seq, gen_seq = str(aln[0].seq).upper(), str(aln[1].seq).upper()
                # Check ref and genome sequence are same length, set Kdist to NA if not
                if len(ref_seq) == len(gen_seq):
                    # Calculate distances based on model
                    Kdist = Kimura80(ref_seq, gen_seq)
                    # Convert numbers to strings
                    if Kdist != "NA":
                        Kdist = str(round(Kdist, 4))
                else:
                    Kdist = "NA"
                # Delete temporary files
                os.remove(query_path+".matcher")
                os.remove(query_path)
            # Make line for temporary file and write to file
            tmp_holder = row[1].to_list()[1:]
            tmp_holder = "\t".join(str(x) for x in tmp_holder)+"\t"+Kdist+"\n"
            tmp_out.write(tmp_holder)

    return(holder_file_name)

def tmp_out_parser(file_list, simple_gff, other_gff):
    # Loop through results 
    gff=pd.DataFrame()
    for file in file_list:
        # read in gff
        in_gff = pd.read_csv(file, sep = "\t")
        # concatenate gff
        gff = pd.concat([gff, in_gff], ignore_index=True)
    # Convert numbers to strings for concatenation
    gff['Kimura'] = gff['Kimura'].astype(str)
    # Convert new data onto metadata
    gff['metadata'] = gff['metadata'] + ";KIMURA80=" + gff['Kimura']
    # Remove unnecessary rows
    gff = gff.drop(columns = ['Kimura', 'repeat_family'])
    # Combine columns, sort and drop unneccessary columns
    gff = pd.concat([gff, simple_gff, other_gff], ignore_index=True)
    gff = gff.drop(columns = ['level_0', 'index'])
    gff = gff.sort_values(by=['seqnames', 'start'])
    gff = gff.reset_index()

    return(gff)

if __name__ == "__main__":
    
    start_time = time()

    # check files exist
    file_check(args.repeat_library, args.in_gff, args.genome, args.out_gff, args.temp_dir)
    
    # split library file
    print("Splitting repeat library")
    splitter(args.repeat_library, args.temp_dir)

    # read in gff and take head
    print("Reading in gff")
    in_gff, simple_gff, other_gff = parse_gff(args.in_gff)
    
    # create as many processes as instructed cores
    num_processes = args.cores

    # calculate the chunk size as an integer
    chunk_size = int(in_gff.shape[0]/num_processes)

    # break into chunks
    chunks = [in_gff.iloc[in_gff.index[i:i + chunk_size]] for i in range(0, in_gff.shape[0], chunk_size)]

    # set pybedtools temp path
    try:
        os.mkdir(args.temp_dir+"/pybedtools/")
    except FileExistsError:
        pass
    pybedtools.set_tempdir(args.temp_dir+'/pybedtools')

    print("Starting calculations") 
    # Peform calulations in parallel
    func = partial(outer_func, args.genome, args.temp_dir, args.timeout)
    pool = multiprocessing.Pool(processes=num_processes)
    results = pool.map(func, chunks)
    pool.close()
    pool.join()
    print("Finished calculations") 

    # Free up memory (necessary with very large gffs and low memory machines)
    del chunks
    del in_gff
  
    # Read in temp files, fix metadata, add simple repeats back, and sort
    calc_gff = tmp_out_parser(results, simple_gff, other_gff)
        
    # remove first column and write to file
    calc_gff.drop(columns = ['index']).to_csv(args.out_gff, sep = "\t", header = False, index=False)

    # print run time for number of rows
    run_time = time() - start_time
    print("Total run time for ", len(calc_gff), " rows was ", run_time, " seconds")

    # Delete folder of split library
    shutil.rmtree(args.temp_dir+"/split_library/", ignore_errors=True)
