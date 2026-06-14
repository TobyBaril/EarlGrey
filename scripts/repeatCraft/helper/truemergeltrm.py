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
			k, v = i.split("=", 1)
			attrD[k] = v
		return (attrD)


	d = {
		"col": [],
		"start": 0,
		"end": 0,
		"size": 0,
		"strand": "",
		"LTRgroup": ""
	}

	# Number of row of header
	#cnt = 0
	#with open(gff, "r") as f:
	#	for line in f:
	#		cnt += 1
	#		if not line.startswith("#"):
	#			cnt -= 1
	#			break
	print("##gff-version 3")
	with open(gff, "r") as f:
		for line in f:
			if line.startswith("#") or not line.strip():
				continue
			col = line.rstrip().split("\t")
			if len(col) < 9:
				print(
					f"Error parsing intermediate file {gff} — too few columns:",
					file=sys.stderr,
				)
				print(line, file=sys.stderr)
				# Prompt on the real terminal: sys.stdout is redirected to the
				# output GFF here, so write the prompt to stderr instead. When
				# running non-interactively (e.g. in a pipeline) stdin is closed
				# and input() raises EOFError -- default to skipping the line.
				print("Skip this line? (Y/N) ", end="", file=sys.stderr, flush=True)
				try:
					skip_line = input()
				except EOFError:
					skip_line = "Y"
				if skip_line.strip().upper() == "N":
					sys.exit("Error in merging fragment.")
				else:
					continue
			try:
				ltrgroup = re.findall(r"LTRgroup=(.*)$", col[8])
			except Exception as e:
				print(f"Error parsing col[8] in {gff}: {e}", file=sys.stderr)
				print(line, file=sys.stderr)
				continue

			if len(ltrgroup) > 0:
				ltrgroup = ltrgroup[0]  # to string
			else:
				ltrgroup = False

			if ltrgroup:  # have TEgroup tag
				# is the last row also has a tag?
				if tag:
					# is it the same tag?
					if ltrgroup == d["LTRgroup"]:
						d["end"] = col[4]  # update the end
						# adopt the strand of the largest fragment seen so far
						fragsize = int(col[4]) - int(col[3])
						if fragsize > d["size"]:
							d["size"] = fragsize
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
						d["size"] = int(col[4]) - int(col[3])
						d["strand"] = col[6]
						tmpattr = attr2dict(col[8])
						d["LTRgroup"] = tmpattr["LTRgroup"]
						tag = True
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
					d["size"] = 0
					d["strand"] = ""
					d["LTRgroup"] = ""
					tag = False

					print(*col, sep="\t") # print current row
				else:  # Last row has no tag
					print(*col, sep="\t")

	# print the last row (flush the final accumulated cluster, if any)
	if len(d["col"]) > 0:
		col2p = d["col"]
		col2p[3] = d["start"]
		col2p[4] = d["end"]
		col2p[6] = d["strand"]
		print(*col2p, sep="\t")
	sys.stdout.close()
	sys.stdout = stdout
