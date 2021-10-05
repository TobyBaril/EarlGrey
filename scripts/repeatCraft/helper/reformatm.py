import sys
import re

def reformat(rmgff,rmout,outfile):

	# print track
	stdout = sys.stdout
	sys.stdout = open(outfile, 'w')

	classD = {}
	with open(rmout, "r") as f:
		for i in range(3):  # Skip header	
			next(f)
		for line in f:
			[_, _, _, _, _, _, _, _, _, repeat, repeatClass, _, _, _] = line.rstrip().split()[0:14]
			classD[repeat] = repeatClass

	# Rewrite the attr in repeatmasker gff
	print("## gff-version 3")
	with open(rmgff, "r") as f:
		for line in f:
			if line.startswith("#"):  # skip header
				next(f)
			else:
				[seqid, source, T, start, end, score, strand, phase, remark] = line.rstrip().split("\t")
				if re.search(".*Motif:.*",line):
					family = re.findall("Motif:(.*)\"", remark)[0]
					s,e = remark.split()[2:]
				if re.search("Target=.*",line):
					attr =  re.findall("Target=(.*)$",line)[0]
					family = attr.split()[0]
					s = attr.split()[1]
					e = attr.split()[2]

				c = classD[family]
				nremark = "Tstart=" + s + ";Tend=" + e + ";ID=" + family
				print(*[seqid, source, c, start, end, score, strand, phase, nremark], sep="\t")

	sys.stdout.close()
	sys.stdout = stdout

