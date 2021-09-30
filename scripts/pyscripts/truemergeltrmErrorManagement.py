import sys
import re

def trumergeLTR(rmgff,outfile):

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
		"strand": "",
		"LTRgroup": ""
	}

	# Number of row of header
	cnt = 0
	with open(gff, "r") as f:
		for line in f:
			cnt += 1
			if not line.startswith("#"):
				cnt -= 1
				break


	with open(gff, "r") as f:
		for i in range(cnt):
			next(f)
		for line in f:
			col = line.rstrip().split("\t")
			try: 
                ltrgroup = re.findall(r"LTRgroup=(.*)$", col[8])
                    if len(ltrgroup) > 0:
                    ltrgroup = ltrgroup[0]  # to string
                    else:
                        ltrgroup = False
            except:
                ltrgroup = False

			if ltrgroup:  # have TEgroup tag
				# is the last row also has a tag?
				if tag:
					# is it the same tag?
					if ltrgroup == d["LTRgroup"]:
						d["end"] = col[4]  # update the end
						# if size larger than the previous one, update strand
						if int(col[4]) - int(col[3]) > int(d["end"]) - int(d["start"]):
							d["strand"] = col[6]

					else:
						# this is a new group
						# print the save row
						col2p = d["col"]
						col2p[3] = d["start"]
						col2p[4] = d["end"]
						col2p[6] = d["strand"]
						# col2p[8] = re.sub("TEgroup=.*?$", "", col2p[8])
						print(*col2p, sep="\t")

						# update the d with current row
						d["col"] = col
						d["start"] = col[3]
						d["end"] = col[4]
						d["strand"] = col[6]
						tmpattr = attr2dict(col[8])
						d["LTRgroup"] = tmpattr["LTRgroup"]
						tag = True  # useless
				else:  # new group
					# update the d with current row
					d["col"] = col
					d["start"] = col[3]
					d["end"] = col[4]
					d["strand"] = col[6]
					tmpattr = attr2dict(col[8])
					d["LTRgroup"] = tmpattr["LTRgroup"]
					tag = True  # useless
			else:  # don't have TEgroup tag
				if tag:  # Last row has TEgroup, it is the end of the cluster, print last row
					col2p = d["col"]
					col2p[3] = d["start"]
					col2p[4] = d["end"]
					col2p[6] = d["strand"]
					# col2p[8] = re.sub("TEgroup=.*?$", "", col2p[8])
					print(*col2p, sep="\t")

					# clean the d (for safe)
					d["col"] = []
					d["start"] = 0
					d["end"] = 0
					d["strand"] = ""
					d["LTRgroup"] = ""
					tag = False

					print(*col, sep="\t") # print current row
				else:  # Last row has no tag
					print(*col, sep="\t")

	# print the last row
	#print(*d["col"],sep="\t")
	sys.stdout.close()
	sys.stdout = stdout
