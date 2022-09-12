#!/usr/bin/env python3

import sys

if (len(sys.argv) <= 1):
	print("usage: a.py exon_coord_file.fa vcf allele_list... > alleleExonDiff.vcf")
	sys.exit(0)

selectedAlleles = {}
for i in range(3, len(sys.argv)):
	selectedAlleles[sys.argv[i]] = 1
if (len(selectedAlleles) == 0):
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
		alleleExonCoords[allele].append([int(cols[i]), int(cols[i + 1]) + 1])
fp.close()

def RangeOverlap(start, end, exons):
	for e in exons:
		if (not(end <= e[0] or start >= e[1])):
			x = max(start, e[0])
			y = min(end, e[1])
			return y - x
	return 0

def ConvertToExonCoord(pos, exons):
	ret = 0 
	for e in exons:
		if (pos >= e[1]):
			ret += e[1] - e[0] 
		else:
			ret += pos - e[0]
			break
	return ret

alleleExonVcf = {}
fp = open(sys.argv[2])
for line in fp:
	cols = line.rstrip().split()
	allele = cols[0] 
	start = int(cols[1])
	if (start == -1 or allele not in selectedAlleles):
		continue 

	if (allele not in alleleExonVcf):
		alleleExonVcf[allele] = []

	# the intervals are left close right open here
	end = start + max([len(cols[3]), len(cols[4])])
	#alleleExonDiff[allele] += RangeOverlap(start, end, alleleExonCoords[allele])
	if (RangeOverlap(start, end, alleleExonCoords[allele]) > 0):
		cols[1] = str(ConvertToExonCoord(start, alleleExonCoords[allele]) + 1)
		alleleExonVcf[allele].append("\t".join(cols))
fp.close()

for allele in alleleExonVcf:
	for vcf in alleleExonVcf[allele]:
		print(vcf) # this vcf file is 1-based

