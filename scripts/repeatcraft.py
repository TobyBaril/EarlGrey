#!/usr/bin/env python3

import argparse
import sys
import re
import os
import subprocess
from subprocess import call

sys.path.append(os.path.join(sys.path[0],"helper"))

import reformatm
import combineGFFoverlapm
import filtershortm
import fuseltr
import fusetem
import truemergeltrm
import truemergetem
import extraFuseTEm
import extraTrueMergeTEm
import rcStatm

# Arguments
parser = argparse.ArgumentParser(description="RepeatCraft pipeline for improving repeat elements annotation by defragments closely spanced repeat elements,based on sequence similarity and structural features from different annotators")
parser.add_argument("-r", "--rmgff", help="RepeatMasker GFF", type=str)
parser.add_argument("-u", "--rmout", help="RepeatMasker OUT", type=str)
parser.add_argument("-c", "--config", help="Configuration file", type=str)
parser.add_argument("-o", "--output", help="Output file name", default="repeatcraft.gff", type=str)
parser.add_argument("-m", "--mode",
                    help="Merge mode. strict or loose. Default = loose",
                    default="loose", type=str)
args = parser.parse_args()

if len(sys.argv) <= 1:
	parser.print_help()
	sys.exit()

if args.rmgff is None:
	parser.print_help()
	sys.exit("Missing Repeatmasker GFF")

if args.rmout is None:
	parser.print_help()
	sys.exit("Missing Repeatmasker OUT")

if args.config is None:
	parser.print_help()
	sys.exit("Missing Configuration file")

rmgffp = args.rmgff
rmoutp = args.rmout
configp = args.config
outputname = args.output
mergemode = args.mode

param = {
	"shortTEsize": 0,
	"mapfile": "None",
	"ltrgff": "None",
	"maxltrsize": 0,
	"ltrflanksize": 0,
	"tegap": 150
}

# Parse config file
with open(configp, "r") as f:
	for line in f:
		if re.search("shortTE_size:.*", line):
			param["shortTEsize"] = int(re.findall("shortTE_size:(.*)", line.rstrip())[0].lstrip())
		if re.search("mapfile:.*", line):
			param["mapfile"] = re.findall("mapfile:(.*)", line.rstrip())[0].lstrip()
		if re.search("ltr_finder_gff:.*", line):
			param["ltrgff"] = re.findall("ltr_finder_gff:(.*)", line.rstrip())[0].lstrip()
		if re.search("max_LTR_size:.*", line):
			param["maxltrsize"] = int(re.findall("max_LTR_size:(.*)", line.rstrip())[0].lstrip())
		if re.search("LTR_flanking_size:.*", line):
			param["ltrflanksize"] = int(re.findall("LTR_flanking_size:(.*)", line.rstrip())[0].lstrip())
		if re.search("gap_size:.*", line):
			param["tegap"] = int(re.findall("gap_size:(.*)", line.rstrip())[0].lstrip())

# Check inputs
mfile = False
checkltr = False

if param["mapfile"] != "None":
	mfile = True
if param["ltrgff"] != "None":
	checkltr = True

# todo: print config param
# todo: check if any tmp file in the directory, ask if need to remove?

# define the tmpfile before try catch
outputnamelabel_tobesort = outputname + ".rclabel.gff.tmp"
outputnamemerge_tobesort = outputname + ".rmerge.gff.tmp"


try:
	# Reformat GFF
	sys.stderr.write("Step 1: Reformating GFF...\n")
	reformatm.reformat(rmgff=rmgffp, rmout=rmoutp, outfile="tmp01.unsort.gff")
	sortrgff = "sort -k1,1 -k4,4n -k5,5n tmp01.unsort.gff > tmp01.gff"
	subprocess.run(sortrgff,shell=True)

	# Reformat and fuse LTR_FINDER GFF
	if checkltr:
		sys.stderr.write("       Parsing LTR_FINDER GFF...\n")
		filterrow = "grep 'LTR_retrotransposon' " + param["ltrgff"] + " > parseltrfinder.tmp"
		cutcolumn = "cut -f1-8 " + "parseltrfinder.tmp" + "|sort|uniq > parseltrfinder.tmp2"
		sortgff = "sort -k1,1 -k4,4n -k5,5n parseltrfinder.tmp2 > parseltrfinder.tmp"
		subprocess.run(filterrow,shell=True)
		subprocess.run(cutcolumn,shell=True)
		subprocess.run(sortgff,shell=True)
		combineGFFoverlapm.combineGff(ltrgff="parseltrfinder.tmp",column=8,outfile="ltrfinder_reformat.gff")
		param["ltrgff"] =  "ltrfinder_reformat.gff"

	# Label short TEs
	sys.stderr.write("Step 2: Labelling short TEs...\n")
	if mfile:
		filtershortm.filtershortTE(rmgff="tmp01.gff", m=mfile, tesize=param["shortTEsize"], mapfile=param["mapfile"],outfile="tmp02.gff")
	else:
		filtershortm.filtershortTE(rmgff="tmp01.gff", m=False, tesize=param["shortTEsize"],outfile="tmp02.gff",mapfile=False)

	# Label LTR group and TE group
	if mergemode == "strict":
		outputnamelabel = outputname + ".rclabel.gff"
		if checkltr:
			sys.stderr.write("Step 3: Labelling LTR groups...\n")
			fuseltr.fuseltr(rmgff="tmp02.gff", ltrgff_p=param["ltrgff"], ltr_maxlength=param["maxltrsize"],ltr_flank=param["ltrflanksize"], outfile="tmp03.gff")
			sys.stderr.write("Step 4: Labelling TE groups...(strict mode)\n")
			fusetem.fusete(gffp="tmp03.gff", outfile=outputnamelabel, gapsize=param["tegap"])
		else:
			sys.stderr.write("Missing LTR_FINDER GFF, skip adding LTRgroup attribute label.\n")
			sys.stderr.write("Step 4: Labelling TE groups...(strict mode)\n")
			fusetem.fusete(gffp="tmp02.gff", outfile=outputnamelabel, gapsize=param["tegap"])
	else:
		outputnamelabel = outputname + ".rclabel.gff"
		outputnamelabel_tobesort = outputname + ".rclabel.gff.tmp"
		if checkltr:
			sys.stderr.write("Step 3: Labelling LTR groups...\n")
			fuseltr.fuseltr(rmgff="tmp02.gff", ltrgff_p=param["ltrgff"], ltr_maxlength=param["maxltrsize"],
			                ltr_flank=param["ltrflanksize"], outfile="tmp03.gff")
			sys.stderr.write("Step 4: Labelling TE groups...(loose mode)\n")
			extraFuseTEm.truefusete(gffp="tmp03.gff",outfile=outputnamelabel_tobesort,gapsize=param["tegap"])
			c = "sort -k1,1 -k4,4n -k5,5n " + outputnamelabel_tobesort + " >" + outputnamelabel
			subprocess.run(c, shell=True)
		else:
			sys.stderr.write("Missing LTR_FINDER GFF, skip adding LTRgroup attribute label.\n")
			sys.stderr.write("Step 4: Labelling TE groups...(loose mode)\n")
			extraFuseTEm.truefusete(gffp="tmp02.gff", outfile=outputnamelabel_tobesort, gapsize=param["tegap"])
			c = "sort -k1,1 -k4,4n -k5,5n " + outputnamelabel_tobesort + " >" + outputnamelabel
			subprocess.run(c, shell=True)


	# True merge
	sys.stderr.write("Step 5: Merging GFF records by labels...\n")
	strictoutputnamemerge = outputname + ".rmerge.gff"
	outputnamemerge = outputname + ".rmerge.gff"
	outputnamemerge_tobesort = outputname + ".rmerge.gff.tmp"
	if checkltr:
		truemergeltrm.trumergeLTR(rmgff=outputnamelabel, outfile="ltrmerge.tmp.gff")
		if mergemode == "strict":
			truemergetem.truemergete(rmgff="ltrmerge.tmp.gff", outfile=strictoutputnamemerge)
		else:
			extraTrueMergeTEm.extratruemergete(gffp="ltrmerge.tmp.gff",outfile=outputnamemerge_tobesort )
			c = "sort -k1,1 -k4,4n -k5,5n " + outputnamemerge_tobesort + " >" +outputnamemerge
			subprocess.run(c, shell=True)

	else:
		if mergemode == "strict":
			truemergetem.truemergete(rmgff=outputnamelabel, outfile=strictoutputnamemerge)
		else:
			extraTrueMergeTEm.extratruemergete(gffp=outputnamelabel, outfile=outputnamemerge_tobesort )
			c = "sort -k1,1 -k4,4n -k5,5n " + outputnamemerge_tobesort + " >" + outputnamemerge
			subprocess.run(c, shell=True)


	sys.stderr.write("Step 6: Writing stat file..")
	statfname =  outputname + ".summary.txt"
	if mergemode == "strict":
		if checkltr:
			rcStatm.rcstat(rclabelp=outputnamelabel,rmergep=strictoutputnamemerge,outfile=statfname, ltrgroup = True)
		else:
			rcStatm.rcstat(rclabelp=outputnamelabel,rmergep=strictoutputnamemerge,outfile=statfname, ltrgroup = False)
	else:
		if checkltr:
			rcStatm.rcstat(rclabelp=outputnamelabel,rmergep=outputnamemerge,outfile= statfname, ltrgroup = True)
		else:
			rcStatm.rcstat(rclabelp=outputnamelabel,rmergep=outputnamemerge, outfile=statfname, ltrgroup = False)
finally:
	# Remove tmp files
	sys.stderr.write("Removing tmp files...\n")
	os.remove("tmp01.unsort.gff")
	os.remove("tmp01.gff")
	os.remove("tmp02.gff")
	if mergemode != "strict":
		os.remove(outputnamelabel_tobesort)
		os.remove(outputnamemerge_tobesort)
	if checkltr:
		os.remove("tmp03.gff")
		os.remove("ltrmerge.tmp.gff")
		os.remove("parseltrfinder.tmp")
		os.remove("parseltrfinder.tmp2")
		os.remove("ltrfinder_reformat.gff")
	sys.stderr.write("Done\n")
