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
  os.mkdir(args.directory)

# name of out file
out_library=args.directory+'/'+sub('.*/', '', args.library)
nondfam_library=args.directory+'/non-dfam.'+sub('.*/', '', args.library)

# Split the library into Dfam pseudo-curated and non-Dfam sequences in a single pass.
# A sequence is treated as Dfam when its header begins with "DR" and the text
# between "DR" and "#" is all digits; everything else is non-Dfam.
with open(out_library, 'w') as dfam_handle, open(nondfam_library, 'w') as nondfam_handle:
  with open(args.library, 'r') as in_handle:
    for record in SeqIO.parse(in_handle, "fasta"):
      is_dfam = record.name.startswith("DR") and sub('DR', '', sub('#.*', '', record.name)).isdigit()
      SeqIO.write(record, dfam_handle if is_dfam else nondfam_handle, "fasta")
