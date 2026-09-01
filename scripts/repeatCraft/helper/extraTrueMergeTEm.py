import sys
import re
from collections import defaultdict
nested_dict = lambda: defaultdict(nested_dict)

def extratruemergete(gffp,outfile):


	gff = gffp

	# print track
	stdout = sys.stdout
	sys.stdout = open(outfile, 'w')

	d = nested_dict()
	lastchrom = ""


	# Check number of row of header
	cnt = 0
	with open(gff, "r") as f:
		for line in f:
			cnt += 1
			if not line.startswith("#"):
				cnt -= 1
				break

	print("##gff-version 3")

	with open(gff, "r") as f:
		for i in range(cnt):
			next(f)
		for line in f:
			col = line.rstrip().split("\t")

			# if changing to the last column, need to print the what havn't print out (all groups in the last chrom)
			if col[0] != lastchrom:
				if lastchrom != "":
					# print("new chrom, print remaining lastcol") # debug
					for family in d[lastchrom]:
						for group in d[lastchrom][family]:
							col2print = d[lastchrom][family][group]["firstcol"]
							col2print[4] = d[lastchrom][family][group]["end"]
							col2print[6] = d[lastchrom][family][group]["strand"]
							col2print[8] = re.sub("Tstart=.*?;", "Tstart=" + str(d[lastchrom][family][group]["Tstart"]) + ";", col2print[8])
							col2print[8] = re.sub("Tend=.*?;", "Tend=" + str(d[lastchrom][family][group]["Tend"]) + ";", col2print[8])
							print(*col2print,sep="\t")

			lastchrom = col[0]

			# Extract attribute
			cattrD = {}
			cattr = col[8].split(";")
			for i in cattr:
				k, v = i.split("=")
				cattrD[k] = v
			cattrD["Tstart"] = int(cattrD["Tstart"])
			cattrD["Tend"] = int(cattrD["Tend"])

			tegroup = re.findall(r"TEgroup=(.*)$", col[8])
			if len(tegroup) > 0:
				tegroup = tegroup[0]  # to string
			else:
				tegroup = False

			# if with tag, store it
			if tegroup:
				tchrom, tfamily, tnumber = tegroup.split("|")
				if d[tchrom][tfamily][tnumber]:
					# update group information
					d[tchrom][tfamily][tnumber]["end"] = int(col[4])
					if cattrD["Tstart"] < d[tchrom][tfamily][tnumber]["Tstart"]:
						d[tchrom][tfamily][tnumber]["Tstart"] = cattrD["Tstart"]
					if cattrD["Tend"] > d[tchrom][tfamily][tnumber]["Tend"]:
						d[tchrom][tfamily][tnumber]["Tend"] = cattrD["Tend"]
					if cattrD["Tend"] - cattrD["Tstart"] > d[tchrom][tfamily][tnumber]["maxsize"]:
						d[tchrom][tfamily][tnumber]["maxsize"] = cattrD["Tend"] - cattrD["Tstart"]
						d[tchrom][tfamily][tnumber]["strand"] = col[6]
				else: # first member of the group
					d[tchrom][tfamily][tnumber]["firstcol"] = col
					d[tchrom][tfamily][tnumber]["end"] = int(col[4])
					d[tchrom][tfamily][tnumber]["Tstart"] = cattrD["Tstart"]
					d[tchrom][tfamily][tnumber]["Tend"] = cattrD["Tend"]
					d[tchrom][tfamily][tnumber]["maxsize"] = cattrD["Tend"] - cattrD["Tstart"]
					d[tchrom][tfamily][tnumber]["strand"] = col[6]

			else:
				print(*col,sep="\t") # if no tag just print it

		# print group in last chrom
		for family in d[lastchrom]:
			for group in d[lastchrom][family]:
				col2print = d[lastchrom][family][group]["firstcol"]
				col2print[4] = d[lastchrom][family][group]["end"]
				col2print[6] = d[lastchrom][family][group]["strand"]
				col2print[8] = re.sub("Tstart=.*?;", "Tstart=" + str(d[lastchrom][family][group]["Tstart"]) + ";", col2print[8])
				col2print[8] = re.sub("Tend=.*?;", "Tend=" + str(d[lastchrom][family][group]["Tend"]) + ";", col2print[8])
				print(*col2print, sep="\t")

	sys.stdout.close()
	sys.stdout = stdout
