#!/usr/bin/env python3

import sys
import subprocess
import re
import os

from collections import defaultdict
nested_dict = lambda: defaultdict(nested_dict)

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
	print("##gff-version 3")
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

def filtershortTE(rmgff,m,tesize,mapfile,outfile):

	gff = rmgff

	# print track
	stdout = sys.stdout
	sys.stdout = open(outfile, 'w')

	if m == False:
		s = tesize
		sys.stderr.write("Missing mapfile, use unite size for all TEs except simple repeat, low complexity and satellite.\n")
	else:
		s = None

	# if mapfile.tsv is not available, use defalut map
	if s is not None:
		m = {
			"Unknown":s,
			"SINE":s,
			"SINE?":s,
			"LINE":s,
			"LTR" :s,
			"DNA" :s,
			"RC"  :s,
			"rRNA":s,
			"Simple_repeat":0,
			"Low_complexity":0,
			"Satellite":0,
			"snRNA": 0
		}
	else:
		m = {}
		with open(mapfile,"r") as f:
			for line in f:
				c,t = line.rstrip().split("\t")
				m[c] = int(t)
		# Assume SINE? = SINE
		m["SINE?"] = m["SINE"]

	# Find how many line to skip in gff (faster for long file)
	cnt = 0
	with open(gff,"r") as f:
		for line in f:
				cnt +=1
				if not line.startswith("#"):
					cnt -=1
					break

	# Check if size of TEs samller than lower threshold
	with open(gff,"r") as f:
		for i in range(cnt):
			next(f)
		for line in f:
			col = line.rstrip().split("\t")
			generalClass = col[2].split("/")[0]
			teSize = int(col[4]) - int(col[3])
			try:
				m[generalClass]
			except KeyError:
				m[generalClass] = 0
				sys.stderr.write("'" + generalClass + "' is not in the mapfile. Skip filtering " + generalClass)
			if teSize < m[generalClass]:
				col[8] =  col[8] + ";shortTE=T"
			else:
				col[8] = col[8] + ";shortTE=F"
			print(*col,sep = "\t")

	sys.stdout.close()
	sys.stdout = stdout

def fuseltr(rmgff,ltrgff_p,ltr_maxlength,ltr_flank,outfile):

	# print track
	stdout = sys.stdout
	sys.stdout = open(outfile, 'w')

	nested_dict = lambda: defaultdict(nested_dict)

	# Mark the position and group number of LTR identified by LTR_FINDER
	# Add new attribute (LTRgroup), export to new .gff
	ltrD = nested_dict()
	lastChrom = ""
	gpcnt = 1

	bname = os.path.basename(ltrgff_p)
	ltrout =  open(bname.replace(".gff","_label.gff"),"w")
	with open(ltrgff_p,"r") as f:
		for line in f:
			col = line.rstrip().split("\t")
			if int(col[4]) - int(col[3]) > ltr_maxlength:
				col = [str(i) for i in col]
				ltrout.write("\t".join(col) + "\n")
				continue
			else:
				if col[0] != lastChrom:
					gpcnt = 1
				for i in range(int(col[3])-ltr_flank,int(col[4])+ltr_flank+1):
					ltrD[col[0]][i] = gpcnt
				remark = col[0] + "_g" + str(gpcnt)
				col[8] = col[8] + ";LTRgroup=" + remark
				col = [str(i) for i in col]
				ltrout.write("\t".join(col) + "\n")
				lastChrom = col[0]
				gpcnt += 1
	ltrout.close()
	stder1 = "Updated LTR.gff with LTRgroup attribute to:" + bname
	sys.stderr.write(stder1)

	# repeatmasker gff
	print("##gff-version 3")
	with open(rmgff,"r") as f:
		for line in f:
			col = line.rstrip().split("\t")
			# skip non LTR rows
			generalClass = col[2].split("/")[0]
			if generalClass != "LTR":
				print(*col, sep="\t")
				continue
			# skip short sequence
			if re.search(r"shortTE:T",col[8]):
				print(*col,sep="\t")
				continue

			# Check if LTR:repeatmasker is complete covered by LTR:LTR_FINDER + flank
			n = {}
			outrange = False
			for i in range(int(col[3]),int(col[4]) + 1):
				if not isinstance(ltrD[col[0]][i],int):
					outrange = True
					break
				else:
					n[ltrD[col[0]][i]] = "T"

			if outrange == True:
				print(*col,sep="\t")
			else:
				if len(n.keys()) == 1:
					col[8] =  col[8] + ";LTRgroup=" + col[0] + "_g" + str(list(n.keys())[0])
				else:
					for i in range(len(list(n.keys()))):
						if i ==  0:
							col[8] = col[8] + ";LTRgroup=" + col[0] + "_g" + str(list(n.keys())[i])
						else:
							col[8] = col[8] + "," + col[0] + "_g" + str(list(n.keys())[i])
				print(*col,sep="\t")

	sys.stdout.close()
	sys.stdout = stdout

def fusete(gffp,outfile,gapsize=150):

	gff = gffp
	gapSize = gapsize

	# print track
	stdout = sys.stdout
	sys.stdout = open(outfile, 'w')

	# Main
	# Initial parameters dict
	P = {
		"pcol": [],
		"pchrom" : "",
		"pEnd": 0,
		"pFamily": "",
		"pCstart": 0,
		"pCend": 0,
		"plabel": ""
	}

	D = nested_dict()
	C = nested_dict() # consensus coverage


	def update_pcol(c, label=""):
		attr = c[8].split(";")
		attrD = {}
		for i in attr:
			k, v = i.split("=")
			attrD[k] = v

		P["pcol"] = c
		P["pchrom"] = c[0]
		P["pEnd"] = int(c[4])
		P["pFamily"] = attrD["ID"]
		P["pCstart"] = int(attrD["Tstart"])
		P["pCend"] = int(attrD["Tend"])
		P["plabel"] = label

	# quick check number of line of the file
	sh = subprocess.run(['wc', '-l', gff], stdout=subprocess.PIPE)
	totalline = str(int(sh.stdout.split()[0]))
	dcnt = 0

	# Number of row of header
	cnt = 0
	with open(gff,"r") as f:
		for line in f:
				cnt +=1
				if not line.startswith("#"):
					cnt -=1
					break

	print("##gff-version 3")
	with open(gff, "r") as f:
		for i in range(cnt):
			next(f)
		for line in f:
			# If current TE is far away from the last TE
			col = line.rstrip().split("\t")

			dcnt += 1
			sys.stderr.write("\rProgress:" + str(dcnt) + "/" + totalline + "...")

			if (int(col[3]) - P["pEnd"] > gapSize) or (col[0] != P["pchrom"]): # Also make sure they are in different chrom
				if P["pcol"]:  # Make sure not the first line
					print(*P["pcol"], sep="\t")  # print last row
				update_pcol(c=col, label="")
				C = nested_dict() # reset consensus dict
			else:
				# Extract attribute
				cattrD = {}
				cattr = col[8].split(";")
				for i in cattr:
					k, v = i.split("=")
					cattrD[k] = v

				# If current TE close to last TE, but not the same family
				if P["pFamily"] != cattrD["ID"]:
					if P["pcol"]:  # Make sure not the first line
						print(*P["pcol"], sep="\t")
					update_pcol(c=col, label="")
					C = nested_dict() # reset consensus dict
				# If current TE close to last TE and belong to same family
				# Since has compared the family, pcol must not be empty in the following block
				else:
					o = min(int(P["pCend"]), int(cattrD["Tend"])) - max(int(P["pCstart"]), int(cattrD["Tstart"]))
					if o > 0:  # Consensus position overlap
						print(*P["pcol"], sep="\t")
						update_pcol(c=col, label="")
						C = nested_dict()  # reset consensus dict
					else:
						# Group the two TE by adding group attr
						# If the previous TE is already in a group
						# further check if overlap with consensus in the group
						allo = False
						for i in range(int(cattrD["Tstart"]), int(cattrD["Tend"])):
							if i in list(C):
								allo = True
								break
						if allo:
							print(*P["pcol"], sep="\t")
							update_pcol(c=col, label="")
							C = nested_dict()  # reset consensus dict
						else:
							if P["plabel"]:

								# update consensus coverage
								for i in range(int(cattrD["Tstart"]), int(cattrD["Tend"])):
									C[i] = 1

								print(*P["pcol"], sep="\t")
								update_pcol(c=col, label=P["plabel"])  # keep using the same label
								P["pcol"][8] = P["pcol"][8] + ";" + P["plabel"]  # Add attribute
							# Previous TE is the first element in the group, add the attribute before print
							# Previous TE has no label yet
							else:
								# Check label
								if D[col[0]][col[2]][cattrD["ID"]]:  # Use a new label for the TE pair
									D[col[0]][col[2]][cattrD["ID"]] += 1
									groupID = D[col[0]][col[2]][cattrD["ID"]]
								# D[col[0]][col[2]][cattrD["ID"]] = groupID
								else:  # If this is the first family in this chrom
									D[col[0]][col[2]][cattrD["ID"]] = 1
									groupID = 1

								# Add the coverage to consensus dict
								for i in range(P["pCstart"], P["pCend"]):
									C[i] = 1
								for i in range(int(cattrD["Tstart"]), int(cattrD["Tend"])):
									C[i] = 1

								grouplabel = "TEgroup=" + col[0] + "|" + col[2] + "|" + cattrD["ID"] + "|" + str(groupID)
								P["pcol"][8] = P["pcol"][8] + ";" + grouplabel
								print(*P["pcol"], sep="\t")
								update_pcol(c=col, label=grouplabel)
								P["pcol"][8] = P["pcol"][8] + ";" + grouplabel  # Add attribute

	# Print the last row
	print(*P["pcol"], sep="\t")
	sys.stdout.close()
	sys.stdout = stdout

def combineGff(ltrgff,column,outfile):
	# the gff need to be sorted by contig and start position
	# `sort -k1,1 -k4,4n -k5,5n .gff` before running this script
	# ignore strand (all "+")
	# print track
	stdout = sys.stdout
	sys.stdout = open(outfile, 'w')

	with open(ltrgff, "r") as f:
		startpos = 0
		endpos = 0
		lastContig = "null"
		for line in f:
			if column == 9:
				[contig, prog, type, start, end, value, c, coding, attr] = line.rstrip().split("\t")
			else:
				[contig, prog, type, start, end, value, c, coding] = line.rstrip().split("\t")
			start = int(start)
			end = int(end)
			# debug
			# print("{}\t{}\t{}".format(contig,start,end))
			if contig != lastContig:  # if meet a new contig
				# print("!=")
				if startpos != 0:  # not the first
					print("{0}\t{1}\t{2}\t{3}\t{4}\t0\t+\t.\tNA".format(lastContig, prog, type, startpos,endpos))  # print the last contig
				startpos = start
				endpos = end
				lastContig = contig
				continue
			else:  # still on the same contig
				# print("same contig")
				if endpos >= start:  # link consecutive structure together, update the end position
					endpos = end
				else:  # the new structure is not connected with the previous one, print it out
					print("{0}\t{1}\t{2}\t{3}\t{4}\t0\t+\t.\tNA".format(contig, prog, type, startpos, endpos))
					startpos = start
					endpos = end

	sys.stdout.close()
	sys.stdout = stdout

def truefusete(gffp,gapsize,outfile):

	gff = gffp

	# print track
	stdout = sys.stdout
	sys.stdout = open(outfile, 'w')

	# quick check number of line of the file
	sh = subprocess.run(['wc', '-l',gff], stdout=subprocess.PIPE)
	totalline = str(int(sh.stdout.split()[0]))

	cnt = 0
	d = nested_dict()
	lastchrom = ""

	# Progress
	dcnt = 0

	# Check number of row of header
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
			# Progress
			dcnt += 1
			sys.stderr.write("\rProgress:" + str(dcnt) + "/"+ totalline+ "...")

			col = line.rstrip().split("\t")

			# Extract attribute
			cattrD = {}
			cattr = col[8].split(";")
			for i in cattr:
				k, v = i.split("=")
				cattrD[k] = v
			cattrD["Tstart"] = int(cattrD["Tstart"])
			cattrD["Tend"] = int(cattrD["Tend"])

			# if changing to the last column, need to print the what havn't print out (lastcol for all families in last chrom)
			if  col[0] != lastchrom:
				if lastchrom != "":
					#print("new chrom, print remaining lastcol") # debug
					for family in d[lastchrom]:
						print(*d[lastchrom][family]["lastcol"],sep="\t")

			lastchrom = col[0] # Update lastcol

			if d[col[0]][cattrD["ID"]]:  # not the first family on this chrom
				#print("not first family") # debug
				if (int(col[3]) - d[col[0]][cattrD["ID"]]["lastend"]) > gapsize or col[0] != d[col[0]][cattrD["ID"]]["lastcol"][
					0]:
					#print("larger than 150") # debug
					# don't need to group the two records
					# print the lastest record of the lastest family group without adding new label
					col2print = d[col[0]][cattrD["ID"]]["lastcol"]
					print(*col2print, sep = "\t")
					# update the dictionary
					d[col[0]][cattrD["ID"]]["lastcol"] = col
					d[col[0]][cattrD["ID"]]["lastend"] = int(col[4])
					d[col[0]][cattrD["ID"]]["Tstart"] = cattrD["Tstart"]
					d[col[0]][cattrD["ID"]]["Tend"] = cattrD["Tend"]
					d[col[0]][cattrD["ID"]]["lastTElabel"] = False
				else:
					#print("less than 150") # debug
					# Is the lastcol carrying a TEgroup label?
					if d[col[0]][cattrD["ID"]]["lastTElabel"]:
						#print("last one have label") # debug
						# check consensus information (all in last group)
						groupnumber = d[col[0]][cattrD["ID"]]["groupcnt"]
						o = False
						for i in range(cattrD["Tstart"],cattrD["Tend"]):
							if i in list(d[col[0]][cattrD["ID"]][groupnumber]):
								o = True
								break
						if o: # overlap with some copies in the group, break the grouping
							#print("consensus overlap") # debug
							print(*d[col[0]][cattrD["ID"]]["lastcol"], sep="\t")
							d[col[0]][cattrD["ID"]]["lastcol"] = col
							d[col[0]][cattrD["ID"]]["lastend"] = int(col[4])
							d[col[0]][cattrD["ID"]]["Tstart"] = cattrD["Tstart"]
							d[col[0]][cattrD["ID"]]["Tend"] = cattrD["Tend"]
							d[col[0]][cattrD["ID"]]["lastTElabel"] = True
						else:
							#print("consensus pass") # debug
							# print the lastcol directly
							print(*d[col[0]][cattrD["ID"]]["lastcol"], sep="\t")
							# Update consensus coverage
							for i in range(cattrD["Tstart"],cattrD["Tend"]):
								d[col[0]][cattrD["ID"]][groupnumber][i] = 1

							# Update last col using label from last label
							attr = ";TEgroup=" + col[0] + "|" + cattrD["ID"] + "|" + str(d[col[0]][cattrD["ID"]]["groupcnt"])
							col[8] = col[8] + attr
							d[col[0]][cattrD["ID"]]["lastcol"] = col
							d[col[0]][cattrD["ID"]]["lastend"] = int(col[4])
							d[col[0]][cattrD["ID"]]["Tstart"] = cattrD["Tstart"]
							d[col[0]][cattrD["ID"]]["Tend"] = cattrD["Tend"]
							d[col[0]][cattrD["ID"]]["lastTElabel"] = True


					else: # the lastcol is the first element is this group, just need to check if last and current copies overlap
						#print("last copy no label") # debug
						o = min(d[col[0]][cattrD["ID"]]["Tend"], cattrD["Tstart"]) - max(d[col[0]][cattrD["ID"]]["Tstart"], cattrD["Tend"])
						if o > 0:  # Consensus position overlap, don't need to group them just print
							#print("consensus overlap") # debug
							print(*d[col[0]][cattrD["ID"]]["lastcol"], sep="\t")
							d[col[0]][cattrD["ID"]]["lastcol"] = col
							d[col[0]][cattrD["ID"]]["lastend"] = int(col[4])
							d[col[0]][cattrD["ID"]]["Tstart"] = cattrD["Tstart"]
							d[col[0]][cattrD["ID"]]["Tend"] = cattrD["Tend"]
							d[col[0]][cattrD["ID"]]["lastTElabel"] = True

						else: # can open a new group now and update the attr of the last and current copies
							#print("consensus pass") # debug
							# Make a new label for lastcol and current col
							if d[col[0]][cattrD["ID"]]["groupcnt"]: # Is there a previous family group in the same chrom?
								d[col[0]][cattrD["ID"]]["groupcnt"] = d[col[0]][cattrD["ID"]]["groupcnt"] + 1
							else:
								d[col[0]][cattrD["ID"]]["groupcnt"] = 1

							# Mark down the consensus coverage
							groupnumber = d[col[0]][cattrD["ID"]]["groupcnt"]
							for i in range(d[col[0]][cattrD["ID"]]["Tstart"],d[col[0]][cattrD["ID"]]["Tend"]):
								d[col[0]][cattrD["ID"]][groupnumber][i] = 1
							for i in range(cattrD["Tstart"],cattrD["Tend"]):
								d[col[0]][cattrD["ID"]][groupnumber][i] = 1

							# Print lastcol
							lastcol2print = d[col[0]][cattrD["ID"]]["lastcol"]
							attr = ";TEgroup=" + col[0] + "|" + cattrD["ID"] + "|" + str(d[col[0]][cattrD["ID"]]["groupcnt"])
							lastcol2print[8] = lastcol2print[8] + attr
							print(*lastcol2print,sep="\t")

							# Update lastcol
							col[8] = col[8] + attr
							d[col[0]][cattrD["ID"]]["lastcol"] = col
							d[col[0]][cattrD["ID"]]["lastend"] = int(col[4])
							d[col[0]][cattrD["ID"]]["Tstart"] = cattrD["Tstart"]
							d[col[0]][cattrD["ID"]]["Tend"] = cattrD["Tend"]
							d[col[0]][cattrD["ID"]]["lastTElabel"] = True

			else: # first family on this chrom
				#print("first element on this chrom") # debug
				d[col[0]][cattrD["ID"]]["lastcol"] = col
				d[col[0]][cattrD["ID"]]["lastend"] = int(col[4])
				d[col[0]][cattrD["ID"]]["Tstart"] = cattrD["Tstart"]
				d[col[0]][cattrD["ID"]]["Tend"] = cattrD["Tend"]
				d[col[0]][cattrD["ID"]]["lastTElabel"] = False

		# print the last record for all families from the last chrom
		#print("print remaining lastcol") # debug
		for family in d[lastchrom]:
			print(*d[lastchrom][family]["lastcol"],sep="\t")

	sys.stdout.close()
	sys.stdout = stdout

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


	print("##gff-version 3")
	with open(gff, "r") as f:
		for i in range(cnt):
			next(f)
		for line in f:
			col = line.rstrip().split("\t")
			ltrgroup = re.findall(r"LTRgroup=(.*)$", col[8])
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

def rcstat(rclabelp,rmergep,outfile,ltrgroup=True):

	rlabel =  rclabelp
	rmerge = rmergep

	# print track
	stdout = sys.stdout
	sys.stdout = open(outfile, 'w')

	# Read rlabel
	# Stat variables
	telabel = 0
	ltrlabel = 0
	teltrlabel = 0
	teD = {}
	ltrD ={}
	teltrD = {}

	rowRaw = {}
	rowMerge = {}

	# flag
	teflag = False
	ltrflag = False


	# Check number of row of header
	cnt = 0
	with open(rlabel, "r") as f:
		for line in f:
			cnt += 1
			if not line.startswith("#"):
				cnt -= 1
				break

	with open(rlabel,"r") as f:
		for i in range(cnt):
			next(f)
		for line in f:
			col = line.rstrip().split("\t")

			'''
			# Extract attribute
			cattrD = {}
			cattr = col[8].split(";")
			for i in cattr:
				k, v = i.split("=")
				cattrD[k] = v
			'''

			if rowRaw.get(col[2]):
				rowRaw[col[2]] += 1
			else:
				rowRaw[col[2]] = 1

			if re.search("TEgroup=",col[8]):
				teflag = True
				telabel += 1

				if teD.get(col[2]):
					teD[col[2]] += 1
				else:
					teD[col[2]] = 1

			if re.search("LTRgroup=",col[8]):
				ltrflag = True
				ltrlabel += 1

				if ltrD.get(col[2]):
					ltrD[col[2]] += 1
				else:
					ltrD[col[2]] = 1

			if teflag and ltrflag:
				teltrlabel += 1
				if teltrD.get(col[2]):
					teltrD[col[2]] += 1
				else:
					teltrD[col[2]] = 1

			teflag = False
			ltrflag = False

	# Read rmerge
	# Check number of row of header
	cnt = 0
	with open(rmerge, "r") as f:
		for line in f:
			cnt += 1
			if not line.startswith("#"):
				cnt -= 1
				break

	with open(rmerge,"r") as f:
		for i in range(cnt):
			next(f)
		for line in f:
			col = line.rstrip().split("\t")

			if rowMerge.get(col[2]):
				rowMerge[col[2]] += 1
			else:
				rowMerge[col[2]] = 1

	print("#1. Number of repeats (by class) before and after merge")
	print("=============================================================")
	print(*["repeat class","no. before merge","no. after merge"], sep="\t")
	for c in list(rowRaw.keys()):
		print(*[c,rowRaw[c],rowMerge[c]],sep="\t")
	print("\n")
	print("#2. Number of repeats (by class) merged by TEgroup and LTRgroup")
	print("=============================================================")
	for c in list(teD.keys()):
		if ltrgroup:
			if re.search("LTR",c):
				try:
					print(*[c,teD[c],ltrD[c]],sep="\t")
				except:
					print(*[c,teD[c],""],sep="\t")
			else:
				print(*[c,teD[c],""],sep="\t")
		else:
			print(*[c, teD[c], ""], sep="\t")

	sys.stdout.close()
	sys.stdout = stdout
