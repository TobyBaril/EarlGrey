#!/usr/bin/env python

import os
from re import sub
from os.path import exists
import sys
import argparse
from Bio import SeqIO

parser = argparse.ArgumentParser()
parser.add_argument('-l', '--library', type=str, required=True,
                    help='Input library fasta')
parser.add_argument('-d', '--directory', type=str, required=True,
                    help='Path to working data directory')                    
args = parser.parse_args()

# check file/folder exists
if(exists(args.library) == False):
  sys.exit('Input fasta not found')
if(exists(args.directory) == False):
  os.mkdir(args.out_dir)

# name of out file
out_library=args.directory+'/'+sub('.*/', '', args.library)
nondfam_library=args.directory+'/non-dfam.'+sub('.*/', '', args.library)

# write Dfam pseudo-curated sequences to out_library
# sequence headers beginning with DR and having digits between "DR" and "#" treated as as Dfam
with open(out_library, 'w') as out_handle:
  with open(args.library, 'r') as in_handle:
    for record in SeqIO.parse(in_handle, "fasta"):
      if record.name.startswith("DR"):
        if sub('DR', '', sub('#.*', '', record.name)).isdigit():
          SeqIO.write(record, out_handle, "fasta")

# write non-Dfam sequences to out_library
# sequence headers beginning with DR and having digits between "DR" and "#" treated as as Dfam
with open(nondfam_library, 'w') as out_handle:
  with open(args.library, 'r') as in_handle:
    for record in SeqIO.parse(in_handle, "fasta"):
      if record.name.startswith("DR") is False:
        SeqIO.write(record, out_handle, "fasta")
      elif sub('DR', '', sub('#.*', '', record.name)).isdigit() is False:
        SeqIO.write(record, out_handle, "fasta")
