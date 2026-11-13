import pandas as pd
import sys
import csv

args = sys.argv

input = args[1]
dictionary = args[2]
output = args[3]

orig = dict(csv.reader(open(dictionary), delimiter="\t"))
dict = {value:key for key, value in orig.items()}

table = pd.read_csv(input, names = ['scaf', 'two', 'three', 'four', 'five', 'six', 'seven', 'eight', 'nine'], sep='\t', header = None)
table['scaf'] = table['scaf'].astype(str)
# table['scaf'].replace(dict, inplace = True)
table['scaf'] = table['scaf'].replace(dict)

table.to_csv(output, index = False, header = False, sep = '\t') 
