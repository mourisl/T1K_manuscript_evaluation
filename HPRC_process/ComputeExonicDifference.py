#!/usr/bin/env python3

import sys

if (len(sys.argv) <= 1):
	print("usage: a.py exon_coord_file.fa vcf > exonic_difference.txt")
	sys.exit(0)

alleleExonCoords = {}
fp = open(sys.argv[1])
for line in fp:
	header = line.rstrip()
	cols = header[1:].split()
	seq = fp.readline()
	
	allele = cols[0]
	alleleExonCoords[allele] = []
	for i in range(2, len(cols), 2):
		alleleExonCoords[allele].append([int(cols[i]), int(cols[i + 1])])
fp.close()

def RangeOverlap(start, end, exons):
	for e in exons:
		if (not(end <= e[0] or start >= e[1])):
			x = max(start, e[0])
			y = min(end, e[1])
			return y - x
	return 0

alleleExonDiff = {}
fp = open(sys.argv[2])
for line in fp:
	cols = line.rstrip().split()
	allele = cols[0] 
	if (allele not in alleleExonDiff):
		alleleExonDiff[allele] = 0
	start = int(cols[1])
	if (start == -1):
		continue 

	# the intervals are left close right open here
	end = start + max([len(cols[3]), len(cols[4])])
	alleleExonDiff[allele] += RangeOverlap(start, end, alleleExonCoords[allele])
fp.close()

for allele in alleleExonDiff:
	print(allele + "\t" + str(alleleExonDiff[allele]))
