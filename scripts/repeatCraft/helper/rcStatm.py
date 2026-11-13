import sys
import re

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
		
	# delete shitty file headers
	with open(rlabel, "r") as f:
		lines = f.readlines()
	with open(rlabel, "w") as f:
		for line in lines:
			if line.strip("\n") != "##gff-version 3":
				f.write(line)

	with open(rmerge, "r") as f:
		lines = f.readlines()
	with open(rmerge, "w") as f:
		for line in lines:
			if line.strip("\n") != "##gff-version 3":
				f.write(line)

	with open(rlabel,"r") as f:
		for i in range(cnt):
			next(f)
		for line in f:
			col = line.rstrip().split("\t")

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

