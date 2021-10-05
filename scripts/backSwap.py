#!/home/toby/anaconda3/bin/python

import pandas as pd
import sys
import csv

args = sys.argv

input = args[1]
dictionary = args[2]
output = args[3]

orig = dict(csv.reader(open(dictionary), delimiter="\t"))
dict = {value:key for key, value in orig.items()}

table = pd.read_csv(input, names = ['scaf', 'start', 'end', 'repeat', 'score', 'strand'], delim_whitespace = True, header = None)
table['scaf'] = table['scaf'].astype(str)
table['scaf'].replace(dict, inplace = True)

table.to_csv(output, index = False, header = False, sep = '\t') 
