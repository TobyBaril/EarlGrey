import sys
import re

def truemergete(rmgff,outfile):

	gff = rmgff

	# print track
	stdout = sys.stdout
	sys.stdout = open(outfile, 'w')

	tag = False


	def attr2dict(attrcol):
		attr = attrcol.split(";")
		attrD = {}
		for i in attr:
			k, v = i.split("=")
			attrD[k] = v
		return (attrD)


	d = {
		"col": [],
		"start": 0,
		"end": 0,
		"Tstart": 0,
		"Tend": 0,
		"csize": 0,
		"strand": "",
		"TEgroup": ""
	}

	# Number of row of header
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
			tegroup = re.findall(r"TEgroup=(.*)$", col[8])
			if len(tegroup) > 0:
				tegroup = tegroup[0]  # to string
			else:
				tegroup = False

			if tegroup:  # have TEgroup tag
				# is the last row also has a tag?
				if tag:
					# is it the same tag?
					if tegroup == d["TEgroup"]:
						d["end"] = col[4]  # update the end
						# check if need to update Tstart, Tend and strand
						tmpattr = attr2dict(col[8])
						if tmpattr["Tend"] > d["Tend"]:
							d["Tend"] = tmpattr["Tend"]
						if (int(tmpattr["Tend"]) - int(tmpattr["Tstart"])) > d["csize"]:
							d["csize"] = int(tmpattr["Tend"]) - int(tmpattr["Tstart"])
							d["strand"] = col[6]
					else:
						# this is a new group
						# print the save row
						col2p = d["col"]
						col2p[3] = d["start"]
						col2p[4] = d["end"]
						col2p[6] = d["strand"]
						col2p[8] = re.sub("Tstart=.*?;", "Tstart=" + d["Tstart"] + ";", col2p[8])
						col2p[8] = re.sub("Tend=.*?;", "Tend=" + d["Tend"] + ";", col2p[8])
						# col2p[8] = re.sub("TEgroup=.*?$", "", col2p[8])
						print(*col2p, sep="\t")

						# update the d with current row
						d["col"] = col
						d["start"] = col[3]
						d["end"] = col[4]
						d["strand"] = col[6]
						tmpattr = attr2dict(col[8])
						d["Tstart"] = tmpattr["Tstart"]
						d["Tend"] = tmpattr["Tend"]
						d["TEgroup"] = tmpattr["TEgroup"]
						d["csize"] = int(tmpattr["Tend"]) - int(tmpattr["Tstart"])
						tag = True  # useless
				else:  # new group
					# update the d with current row
					d["col"] = col
					d["start"] = col[3]
					d["end"] = col[4]
					d["strand"] = col[6]
					tmpattr = attr2dict(col[8])
					d["Tstart"] = tmpattr["Tstart"]
					d["Tend"] = tmpattr["Tend"]
					d["TEgroup"] = tmpattr["TEgroup"]
					d["csize"] = int(tmpattr["Tend"]) - int(tmpattr["Tstart"])
					tag = True  # useless
			else:  # don't have TEgroup tag
				if tag:  # Last row has TEgroup, it is the end of the cluster, print last row
					col2p = d["col"]
					col2p[3] = d["start"]
					col2p[4] = d["end"]
					col2p[6] = d["strand"]
					col2p[8] = re.sub("Tstart=.*?;", "Tstart=" + d["Tstart"] + ";", col2p[8])
					col2p[8] = re.sub("Tend=.*?;", "Tend=" + d["Tend"] + ";", col2p[8])
					# col2p[8] = re.sub("TEgroup=.*?$", "", col2p[8])
					print(*col2p, sep="\t")

					# clean the d (for safe)
					d["col"] = []
					d["start"] = 0
					d["end"] = 0
					d["Tstart"] = 0
					d["Tend"] = 0
					d["csize"] = 0
					d["strand"] = ""
					d["TEgroup"] = ""
					tag = False

					print(*col, sep="\t") # print current row
				else:  # Last row has no tag
					print(*col, sep="\t")

	# print the last row
	if len(d["col"]) > 0:
		print(*d["col"],sep="\t")
	sys.stdout.close()
	sys.stdout = stdout
