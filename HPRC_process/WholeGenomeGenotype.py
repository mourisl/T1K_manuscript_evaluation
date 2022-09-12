#!/usr/bin/env python3

import sys

if (len(sys.argv) <= 1):
	print("Usage: python3 a.py minimap2.paf")
	sys.exit(1) 

def RangeOverlap(r1, r2):
	overlapSize = 0
	if (r2[1] <= r1[0] or r2[0] >= r1[1]):
		return False
	elif (r1[0] <= r2[0] and r1[1] >= r2[1]):
		overlapSize = r2[1] - r2[0] 
	elif (r2[0] <= r1[0] and r2[1] >= r1[1]):
		overlapSize = r1[1] - r1[0] 
	elif (r1[0] >= r2[0]):
		overlapSize = r2[1] - r1[0]
	elif (r1[0] < r2[0]):
		overlapSize = r1[1] - r2[0] 
	#print(r1, r2, overlapSize)
	# Similarity >= 0.95
	if (overlapSize > 0.95 * (r1[1] - r1[0])
			or overlapSize > 0.95 * (r2[1] - r2[0])):
		return True 
	return False

ranges = []
fp = open(sys.argv[1])
for line in fp:
	cols = line.rstrip().split()
	allele = cols[0]
	# Coverage >= 99%
	if (int(cols[3]) - int(cols[2]) < 0.99 * int(cols[1])):
		continue
	
	grange = [int(cols[7]), int(cols[8]), int(cols[9]), float(cols[9]) / float(cols[10]), allele]
	addNew = True
	for i, r in enumerate(ranges):
		if (RangeOverlap(grange, r)):
			if (grange[2] > r[2] or (grange[2] == r[2] and grange[3] > r[3])):
				ranges[i][:] = grange[:]
			addNew = False
			break
	if (addNew):
		ranges.append(grange)
fp.close()

for i, r in enumerate(ranges):
	print(r[-1])
