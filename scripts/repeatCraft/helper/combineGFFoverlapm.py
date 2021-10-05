import sys

# the gff need to be sorted by contig and start position
# `sort -k1,1 -k4,4n -k5,5n .gff` before running this script
# ignore strand (all "+")

def combineGff(ltrgff,column,outfile):

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

