import re
import sys
import os
from collections import defaultdict


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














