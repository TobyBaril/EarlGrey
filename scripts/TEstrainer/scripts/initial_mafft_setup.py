#!/usr/bin/env python

import argparse
import sys
import os
from os.path import exists

parser = argparse.ArgumentParser()
parser.add_argument('-d', '--directory', type=str, required=True,
                    help='Directory housing TEstrainer run')
parser.add_argument('-r', '--iteration', type=str, required=True,
                    help='Interation number of TEstrainer curation')
parser.add_argument('-s', '--seq_name', type=str, required=True,
                    help='Sequence being prepared for alignment')
parser.add_argument('-g', '--genome', type=str, required=True,
                    help='Path to genome sequence')
parser.add_argument('-f', '--flank', type=int, default=1000,
                    help='Length of flank to be extended')
parser.add_argument('-n', '--no_seq', type=int, default=20,
                    help='Number of sequences to use for alignment')
parser.add_argument('-D', '--debug', action='store_true',
                    help='Set for full messaging')
args = parser.parse_args()

# function to check files exist
def file_check(file_name, debug):
  if(exists(file_name) == False):
    if(debug == False):
      exit()
    else:
      sys.exit((file_name+" not found"))

# check input files/folders exist
file_check((args.directory+"/run_"+args.iteration+"/raw/"+args.seq_name), args.debug) # current consensus
file_check((args.directory+"/run_0/og/"+args.seq_name), args.debug) # og consensus
file_check(args.genome, args.debug) # genome sequence
file_check((args.directory+'/run_'+args.iteration+'/self_search/'), args.debug) # self search folder
file_check((args.directory+'/run_'+args.iteration+'/to_align/'), args.debug) # final output folder
file_check((args.directory+'/run_'+args.iteration+'/TEtrim_complete/'), args.debug) # alternate final output folder

import string
from os import system
import statistics
import re
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import pandas as pd
import numpy as np
import pyranges as pr

# function to output current consensus and exit if less than a certain number of sequences are found
def size_check(var_name, size):
  if(len(var_name) < size):
      with open((args.directory+'/run_'+args.iteration+'/TEtrim_complete/'+args.seq_name), "w") as o:
        SeqIO.write(start_seq, o, "fasta-2line")
      exit()

# Function for writing sequence from df to file    
def df_to_fasta(df, path, name, overwrite):
  df = df.astype('str')
  if('Chromosome' not in df.columns) or ('Start' not in df.columns) or ('End' not in df.columns) or ('Seq' not in df.columns):
    sys.exit(status="Necessary sequence columns not found in df")
  # Name sequences
  if('Name' not in df.columns) or (name is True):
    if('Strand' in df.columns):
      df['Name']=(df.Chromosome+":"+df.Start+"-"+df.End+"("+df.Strand+")")
    else:
      df['Name']=(df.Chromosome+":"+df.Start+"-"+df.End)
  # Write sequences
  if(overwrite is True):
    with open((path), "w") as o:
      for x in range(len(df)):
        SeqIO.write(SeqRecord(seq=Seq(df.Seq[x]),id=df.Name[x],name=df.Name[x], description=""), o, "fasta-2line")
  else:
    with open((path), "a") as o:
      for x in range(len(df)):
        SeqIO.write(SeqRecord(seq=Seq(df.Seq[x]),id=df.Name[x],name=df.Name[x], description=""), o, "fasta-2line")

# Function for adjusting blast output for coordinates direction, convert to base 0
def blast_to_bed(df):
  rev = df['Start'] > df['End']
  fwd = df['Start'] < blast_df['End']
  df.loc[rev, ['Start', 'End']] = (df.loc[rev, ['End', 'Start']].values)
  df.loc[fwd, ['Start']] = df.loc[fwd, ['Start']] - 1
  df.loc[rev, ['End']] = df.loc[rev, ['End']] +1
  return(df)

# read in starting seq
start_seq = SeqIO.read((args.directory+"/run_"+args.iteration+"/raw/"+args.seq_name), "fasta")

# check if sequence is largely tandem repeat
# make pd.df to use as holder for trf coordinates
te_len=len(start_seq.seq)
trf_df = pd.DataFrame(columns=['Chromosome', 'Start', 'End'])
with open((args.directory+'/run_'+args.iteration+'/raw/'+args.seq_name+'.trf'), 'r') as trf:
  for line in trf:
    if line.startswith('@'):
      seqnames = re.sub('@', '', line.split()[0])
    elif float(line.split()[3]) > 5:
      trf_df = pd.concat([trf_df, pd.DataFrame({'Chromosome':[seqnames], 'Start':[int(line.split()[0])], 'End':[int(line.split()[1])]})])
# Merge trf pr
trf_pr=pr.PyRanges(df=trf_df).merge()
# Calculate total length of satellite sequences
# TEs >90% consider satellites
if trf_pr.length/te_len > 0.9:
  with open((args.directory+'/run_'+args.iteration+'/TEtrim_complete/'+args.seq_name), "w") as o:
    if args.debug is True:
      print('Hold your horses, '+start_seq.name+' is likely a tandem repeat')
    SeqIO.write(start_seq, o, "fasta-2line")
  exit()
# TEs >50% consider TEs needing trimming
elif (trf_pr.length/te_len > 0.5) & (len(start_seq.seq) > 500):
  if args.debug is True:
      print('Whoa boy, '+start_seq.name+' is inside tandem repeat')
  # make pr object of whole repeat
  te_pr=pr.from_dict({'Chromosome': [seqnames], 'Start': [0], 'End': [te_len]})
  # get inverse or trf, i.e. non tandem repeat section
  trf_inverse_pr=te_pr.subtract(trf_pr)
  # convert to df
  trf_inverse_df=trf_inverse_pr.as_df()
  # calculate width
  trf_inverse_df['Width']=trf_inverse_df.End-trf_inverse_df.Start+1
  # select longest non TR section and convert to ranges
  non_trf_longest_df=trf_inverse_df[trf_inverse_df.Width==max(trf_inverse_df.Width)].head(n=1)
  non_trf_longest_pr=pr.PyRanges(df=non_trf_longest_df)
  # get sequence of non TR section, write to file
  non_trf_seq=SeqRecord(seq=Seq(str(start_seq.seq)[int(non_trf_longest_df['Start']):int(non_trf_longest_df['End'])]), id=start_seq.id, name=start_seq.name)
  with open((args.directory+'/run_'+args.iteration+'/TEtrim_complete/'+args.seq_name), "w") as o:
    SeqIO.write(non_trf_seq, o, "fasta-2line")
  exit()

# perform initial blast
system("blastn -task dc-megablast -query "+args.directory+"/run_"+args.iteration+"/raw/"+args.seq_name+" -db "+args.genome+" -evalue 1e-5 -outfmt \"6 qseqid sseqid pident length qstart qend qlen sstart send slen evalue bitscore qcovs\" -out "+args.directory+"/run_"+args.iteration+"/initial_blast/"+args.seq_name+".out -num_threads 1")

# check if any hits found, if not write to file and exit
if os.path.getsize(args.directory+'/run_'+args.iteration+'/initial_blast/'+args.seq_name+'.out') == 0:
  with open((args.directory+'/run_'+args.iteration+'/TEtrim_complete/'+args.seq_name), "w") as o:
        SeqIO.write(start_seq, o, "fasta-2line")
  sys.exit()

# read in blast table and filter
blast_df = pd.read_table((args.directory+'/run_'+args.iteration+'/initial_blast/'+args.seq_name+'.out'), names=['qseqid', 'Chromosome', 'pident', 'length', 'qstart', 'qend', 'qlen', 'Start', 'End', 'slen', 'evalue', 'bitscore', 'qcovs'])
blast_df = blast_df.query('pident >= 70 & qcovs >= 50').copy()
size_check(blast_df, 3)

# adjust for coordinates direction, convert to base 0
blast_to_bed(blast_df)

# select 50 best, add flanks and convert to ranges
blast_df = blast_df.sort_values(by = 'bitscore', ascending=False)
blast_df = blast_df.iloc[:args.no_seq]
blast_df['Start'] = blast_df['Start'] - args.flank
blast_df['End'] = blast_df['End'] + args.flank
blast_df.loc[blast_df['Start'] < 1, 'Start'] = 0
blast_df.loc[blast_df['End'] > blast_df['slen'], 'End'] = blast_df['slen']
blast_df = blast_df.sort_values(by = ['Chromosome', 'Start'])
blast_gr = pr.from_dict({"Chromosome": blast_df.Chromosome, "Start": blast_df.Start, "End": blast_df.End, "Bitscore": blast_df.bitscore, "slen" : blast_df.slen})

# reduce/merge ranges
best_hits_df = blast_gr.cluster(strand=False).df.groupby(['Cluster']).agg({'Chromosome':'first', 'Start':'min', 'End':'max', 'Bitscore':'max'})[['Chromosome','Start','End','Bitscore']].reset_index()
best_hits_py = pr.from_dict({"Chromosome": best_hits_df.Chromosome, "Start": best_hits_df.Start, "End": best_hits_df.End, "Bitscore": best_hits_df.Bitscore})

# get sequence
best_hits_seq = pr.get_fasta(best_hits_py, path=args.genome)

# convert df to strings
best_hits_df = best_hits_df.astype('str')
best_hits_df['Name'] = (best_hits_df.Chromosome+":"+best_hits_df.Start+"-"+best_hits_df.End+"#"+best_hits_df.Bitscore)
best_hits_df['Seq'] = best_hits_seq
best_hits_outpath = (args.directory+'/run_'+args.iteration+'/self_search/'+args.seq_name+'_check_1')

# write sequences to file
df_to_fasta(best_hits_df, best_hits_outpath, False, True)

# blast against og sequence,read in and filter
system("blastn -task dc-megablast -query "+args.directory+'/run_0/og/'+args.seq_name+" -subject "+args.directory+"/run_"+args.iteration+"/self_search/"+args.seq_name+"_check_1 -outfmt \"6 sseqid pident qcovs\" -out "+args.directory+"/run_"+args.iteration+"/self_search/"+args.seq_name+"_check_1.out")
check_df = pd.read_table((args.directory+"/run_"+args.iteration+"/self_search/"+args.seq_name+"_check_1.out"), names=['Chromosome', 'pident', 'qcovs'])
check_df = check_df.query('pident >= 70 & qcovs >= 50')
size_check(check_df, 3)

# select one instance of each sequence
check_df = check_df.groupby(['Chromosome']).head(n=1).reset_index()
check_df[['Chromosome', 'bitscore']] = check_df['Chromosome'].str.split('#', n=1, expand = True)
best_hits_df[['Name', 'bitscore']] = best_hits_df['Name'].str.split('#', n=1, expand = True)
size_check(check_df, 3)

# sort by bitscore and select no_seq (default 20) highest bitscores
check_df = check_df.astype({'bitscore': 'float'})
check_df = check_df.sort_values(by = 'bitscore', ascending=False)

# select accurate and write ready for self blast
correct_df = best_hits_df.loc[best_hits_df.Name.isin(check_df.Chromosome)].reset_index()
correct_df_path = (args.directory+'/run_'+args.iteration+'/self_search/'+args.seq_name+'_check_2')
df_to_fasta(correct_df, correct_df_path, False, True)

# self blast
system("blastn -task dc-megablast -query "+correct_df_path+" -subject "+correct_df_path+" -outfmt \"6 qseqid sseqid length pident qstart qend bitscore\" -out "+correct_df_path+".out")
self_blast_df = pd.read_table((args.directory+"/run_"+args.iteration+"/self_search/"+args.seq_name+"_check_2.out"), names=['Chromosome', 'sseqid', 'length', 'pident', 'Start', 'End', 'bitscore'])
size_check(self_blast_df, 3)

# fix coordinates
self_blast_df.Start = self_blast_df.Start - 1

# filter self hits, small hits
self_blast_df = self_blast_df[self_blast_df.Chromosome != self_blast_df.sseqid].copy()
self_blast_df = self_blast_df[self_blast_df.length >= (0.5*len(start_seq))].copy()
size_check(self_blast_df, 3)

# calculate and filter quantiles
self_blast_df['q1'] = self_blast_df.groupby('Chromosome').Start.transform(lambda x: x.quantile(0.1))
self_blast_df['q9'] = self_blast_df.groupby('Chromosome').End.transform(lambda x: x.quantile(0.9))
self_blast_df = self_blast_df[(self_blast_df.Start.astype(float) >= self_blast_df.q1) & (self_blast_df.End.astype(float) <= self_blast_df.q9)].copy()
size_check(self_blast_df, 3)

# create and reduce/merge ranges
self_blast_gr = pr.from_dict({"Chromosome": self_blast_df.Chromosome, "Start": self_blast_df.Start, "End": self_blast_df.End})
self_hits_trimmed_df = self_blast_gr.cluster(strand=False).df.groupby(['Cluster']).agg({'Chromosome':'first', 'Start':'min', 'End':'max'})[['Chromosome','Start','End']].reset_index()
self_hits_trimmed_py = pr.from_dict({"Chromosome": self_hits_trimmed_df.Chromosome, "Start": self_hits_trimmed_df.Start, "End": self_hits_trimmed_df.End})
size_check(self_hits_trimmed_py, 3)
 
# get trimmed sequence
self_hits_trimmed_seq = pr.get_fasta(self_hits_trimmed_py, path=(args.directory+'/run_'+args.iteration+'/self_search/'+args.seq_name+'_check_2'))
self_hits_trimmed_df['Seq'] = self_hits_trimmed_seq
self_hits_trimmed_outpath = (args.directory+'/run_'+args.iteration+'/to_align/'+args.seq_name)

# write current consensus fasta to file
with open(self_hits_trimmed_outpath, "w") as o:
  SeqIO.write(start_seq, o, "fasta-2line")

# write trimmed fasta to file
df_to_fasta(self_hits_trimmed_df, self_hits_trimmed_outpath, True, False)
