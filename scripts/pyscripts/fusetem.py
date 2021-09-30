import sys
import subprocess
from collections import defaultdict

nested_dict = lambda: defaultdict(nested_dict)

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

