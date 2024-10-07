#!/usr/bin/env python3

import sys
import csv

args = sys.argv

dictionary = args[1]
fastaIn = args[2]

a=dict(csv.reader(open(dictionary),delimiter="\t"))
lines=open(fastaIn).read().split("\n")
lines=map(lambda x: ">"+a[x[1:]] if x and x[0]==">" else x,lines);print("\n".join(lines))
