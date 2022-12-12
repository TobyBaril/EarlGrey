#!/usr/bin/env python

import os
import re
from os.path import exists
import sys
import argparse
from Bio import SeqIO

parser = argparse.ArgumentParser()
parser.add_argument('-i', '--in_seq', type=str, required=True,
                    help='Input multi-fasta to be split')
parser.add_argument('-o', '--out_dir', type=str, required=True,
                    help='Output directory')                    
args = parser.parse_args()
if(exists(args.in_seq) == False):
  sys.exit('File not found')
if(exists(args.out_dir) == False):
  os.mkdir(args.out_dir)
file_list=[]

# split fasta file
with open(args.in_seq, 'r') as handle:
    for record in SeqIO.parse(handle, "fasta"):
        file_name = (args.out_dir+"/"+record.name.split(sep="#")[0]+".fasta")
        file_list.append(record.name.split(sep="#")[0]+".fasta")
        SeqIO.write(record, file_name, "fasta-2line")
# write file list
with open((args.out_dir+"/"+re.sub('.*/', '', args.in_seq)+"_split.txt"), 'w') as fp:
    for item in file_list:
        # write each item on a new line
        fp.write("%s\n" % item)
