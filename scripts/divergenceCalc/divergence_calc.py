import os
from os.path import exists, getsize
import sys
import argparse
import numpy as np
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
parser.add_argument('-c', '--chunk_factor', type=int, default=4,
                    help='Number of work chunks to create per core for load balancing (default 4). '
                         'Higher values give finer-grained scheduling when per-repeat alignment '
                         'times are uneven, at the cost of slightly more worker restarts.')

args = parser.parse_args()

def file_check(repeat_library, in_gff, genome, out_gff, temp_dir):
    if(exists(repeat_library) == False or exists(in_gff) == False or exists(genome) == False):
        sys.exit('Files not found. Requires the repeat library, path to the genome, and path to gff containing coordinates and corresponding repeat files')
    if(exists(genome+".fai") == False):
        print("Indexing genome")
        subprocess.run(["samtools","faidx",genome], stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)
    if(exists(temp_dir) == False):
        os.mkdir(temp_dir)
    if(exists(os.path.join(temp_dir, "qseqs")) == False):
        os.mkdir(os.path.join(temp_dir, "qseqs"))
    if(exists(os.path.join(temp_dir, "split_library")) == False):
        os.mkdir(os.path.join(temp_dir, "split_library"))
        

def splitter(in_seq, temp_dir):
    with open(in_seq, 'r') as handle:
        for record in SeqIO.parse(handle, "fasta"):
            repeat_name = record.name.split(sep="#")[0]
            repeat_name = repeat_name.lower()
            file_name = os.path.join(temp_dir, "split_library", repeat_name + ".fasta")
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

# Map each ordered base pair to its category once, at import time, so Kimura80
# does a single O(1) dict lookup per aligned base instead of three list scans.
_BASE_PAIR_CATEGORY = {pair: "MATCH" for pair in ("AA", "GG", "CC", "TT")}
_BASE_PAIR_CATEGORY.update({pair: "TRANSITION" for pair in ("AG", "GA", "CT", "TC")})
_BASE_PAIR_CATEGORY.update({pair: "TRANSVERSION" for pair in (
    "AC", "CA", "AT", "TA", "GC", "CG", "GT", "TG")})


def Kimura80(qseq, sseq):
    """
    Calculations adapted from https://github.com/kgori/python_tools_on_github/blob/master/pairwise_distances.py
    """
    # set counters to 0
    m,ts,tv=0,0,0
    # count transitions, transversions, matches via a single dict lookup per base
    # (pairs not in the table -- gaps, Ns, mismatching ambiguity codes -- score 0)
    for i, j in zip(qseq, sseq):
        category = _BASE_PAIR_CATEGORY.get(i+j)
        if category == "MATCH": m+=1
        elif category == "TRANSITION": ts+=1
        elif category == "TRANSVERSION": tv+=1
    # count number of bp which align (excludes gaps, Ns)
    aln_len = m + ts + tv
    
    if aln_len != 0:
        # calculate p and q
        p = ts/aln_len
        q = tv/aln_len

        # The Kimura-2-parameter distance is only defined while both inner terms
        # stay positive. Saturated (highly divergent) alignments drive them to or
        # below zero, where the distance diverges -- guard against this so log()/
        # sqrt() never raise a math domain error and crash the worker. Such cases
        # are reported as NA, consistent with the unalignable case below.
        term1 = 1 - 2*p - q
        term2 = 1 - 2*q
        if term1 > 0 and term2 > 0:
            Kimura_dist = -0.5 * log(term1 * sqrt(term2))
        else:
            Kimura_dist = "NA"
    else:
        Kimura_dist = "NA"

    return(Kimura_dist)

def batch_getfasta(gff, genome_path):
    """Extract every query sequence for a chunk in a single bedtools getfasta call.

    Returns a dict mapping the GFF row index (as a string) to its extracted
    sequence. Extraction is strand-aware (-s), so minus-strand features are
    reverse-complemented exactly as the previous per-row code did. The row index
    is carried in the BED name column and recovered from the FASTA header, so the
    mapping is robust to bedtools reordering, header-format differences between
    bedtools versions, or any record it declines to emit.
    """
    bed_lines = []
    for idx, r in gff.iterrows():
        # GFF is 1-based inclusive; BED is 0-based half-open -> start - 1
        bed_lines.append("\t".join([
            str(r['seqnames']), str(int(r['start']) - 1), str(int(r['end'])),
            str(idx), ".", str(r['strand'])
        ]))
    if not bed_lines:
        return {}
    seqs = {}
    try:
        bed = pybedtools.BedTool("\n".join(bed_lines) + "\n", from_string=True)
        fasta = bed.sequence(fi=genome_path, s=True, name=True)
        for record in SeqIO.parse(fasta.seqfn, "fasta"):
            # Header is one of 'idx', 'idx::chr:start-end(strand)', 'idx(strand)'
            # depending on bedtools version; strip the coordinate/strand suffixes.
            key = record.id.split("::")[0].split("(")[0]
            seqs[key] = str(record.seq)
    except Exception:
        # A whole-batch failure (rare samtools/bedtools error) returns an empty
        # map; the caller then marks each affected row as failed -> NA, matching
        # the original per-row error handling.
        return {}
    return seqs


def outer_func(genome_path, temp_dir, timeoutSeconds, chunk_path):
    # Set pybedtools temp directory within this worker (required when using forkserver).
    pybedtools.set_tempdir(os.path.join(temp_dir, 'pybedtools'))
    # Load the chunk from its temp file and delete it immediately to reclaim disk space.
    # index_col=0 restores the original GFF row indices, which are used as unique
    # per-row filenames in qseqs/. Without this, every chunk restarts at 0 and
    # concurrent workers collide on the same paths.
    gff = pd.read_csv(chunk_path, sep="\t", index_col=0)
    os.remove(chunk_path)
    generated_name = file_name_generator()
    holder_file_name = os.path.join(temp_dir, generated_name)
    failed_file_name = os.path.join(temp_dir, "failed_" + generated_name)
    # Extract all query sequences for this chunk in a single bedtools getfasta
    # call (one genome open per chunk instead of one per row).
    query_seqs = batch_getfasta(gff, genome_path)
    # Release the single getfasta temp file now that it has been read into memory.
    pybedtools.cleanup(remove_all=False)
    with open(holder_file_name, 'w') as tmp_out:
        header = list(gff.columns.values)[1:] + ["Kimura"]
        header = "\t".join(header)+"\n"
        tmp_out.write(header)
        for row in gff.iterrows():
            # Set index
            idx = row[0]
            # Set scaffold, coordinates, strand, repeat family
            seqnames, start, end, strand, repeat_family = str(row[1]['seqnames']), str(row[1]['start'] - 1), str(row[1]['end']), str(row[1]['strand']), str(row[1]['repeat_family'])
            # Set path for query sequence
            query_path = os.path.join(temp_dir, "qseqs", str(idx))
            # Write this row's query sequence (from the batched extraction) to file
            # for matcher; a missing sequence is treated as a getfasta failure.
            seq = query_seqs.get(str(idx))
            if seq:
                with open(query_path, "w") as qf:
                    qf.write(">" + str(idx) + "\n" + seq + "\n")
            else:
                with open(failed_file_name, "a") as failed_file:
                    failed_file.write(seqnames+":"+start+"-"+end+"_"+strand+"_"+repeat_family+"\n")
            if exists(query_path) is True and getsize(query_path) > 0:
                # Set path to subject sequence
                subject_path = os.path.join(temp_dir, "split_library", repeat_family + ".fasta")
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
                if exists(query_path):
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
                if exists(query_path):
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
    # Use forkserver so workers start from a clean process rather than inheriting
    # a fork-copy of the parent address space (which holds the full GFF DataFrame).
    # This is the primary mechanism for reducing peak RAM at high thread counts.
    if multiprocessing.get_start_method(allow_none=True) is None:
        try:
            multiprocessing.set_start_method('forkserver')
        except (RuntimeError, ValueError):
            # Fall back to the existing/default start method if forkserver is unavailable
            # or if a start method has already been set by the environment.
            pass

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

    # Split into more chunks than workers (chunk_factor per core) so imap_unordered
    # can rebalance dynamically: a worker that finishes a light chunk immediately
    # picks up another instead of idling while one worker grinds through a chunk
    # full of slow or timeout-prone alignments. Cap so there is at most one chunk
    # per row (and at least one chunk overall).
    n_chunks = min(max(num_processes * args.chunk_factor, 1), max(len(in_gff), 1))
    chunks = [in_gff.iloc[idx] for idx in np.array_split(range(len(in_gff)), n_chunks)]

    # Write chunks to temp TSV files so the parent DataFrame can be freed before workers
    # are created. Workers read from disk rather than receiving pickled DataFrames via IPC.
    print("Writing chunks to disk")
    chunk_files = []
    for i, chunk in enumerate(chunks):
        chunk_path = os.path.join(args.temp_dir, f"chunk_{i}.tsv")
        chunk.to_csv(chunk_path, sep="\t", index=True)
        chunk_files.append(chunk_path)

    # Free the main GFF DataFrame and chunk list from the parent process before the pool
    # is created. With forkserver this is already avoided, but freeing here also reduces
    # parent RSS during the pool run, which matters on memory-constrained machines.
    del chunks
    del in_gff

    # set pybedtools temp path (also set per-worker inside outer_func for forkserver)
    try:
        os.mkdir(os.path.join(args.temp_dir, "pybedtools"))
    except FileExistsError:
        pass
    pybedtools.set_tempdir(os.path.join(args.temp_dir, 'pybedtools'))

    print("Starting calculations") 
    # Perform calculations in parallel. maxtasksperchild=1 restarts each worker after
    # processing one chunk, releasing any lingering pybedtools handles or cached objects.
    func = partial(outer_func, args.genome, args.temp_dir, args.timeout)
    pool = multiprocessing.Pool(processes=num_processes, maxtasksperchild=1)
    results = list(pool.imap_unordered(func, chunk_files))
    pool.close()
    pool.join()
    print("Finished calculations") 

    # Read in temp files, fix metadata, add simple repeats back, and sort
    calc_gff = tmp_out_parser(results, simple_gff, other_gff)
        
    # remove first column and write to file
    calc_gff.drop(columns = ['index']).to_csv(args.out_gff, sep = "\t", header = False, index=False)

    # print run time for number of rows
    run_time = time() - start_time
    print("Total run time for ", len(calc_gff), " rows was ", run_time, " seconds")

    # Delete folder of split library
    shutil.rmtree(os.path.join(args.temp_dir, "split_library"), ignore_errors=True)
