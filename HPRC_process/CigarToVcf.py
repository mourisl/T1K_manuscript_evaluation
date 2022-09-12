#!/usr/bin/env python3

# Find the difference on the query

import sys

if (len(sys.argv) <= 1):
	print("usage: a.py ref.fa query.fa align.sam")
	sys.exit(0)

coverageThreshold = 0.9

def CigarToBlocks(cigar):
	blocks = []
	n = 0
	for c in cigar:
		tmp = ord(c) - ord("0")
		if (tmp >= 0 and tmp <= 9):
			n = n * 10 + tmp
		else:
			blocks.append([n, c]) 
			n = 0
	return blocks

def ReverseComplement(s):
	rc = []
	rcMap = {"A":"T", "C":"G", "G":"C", "T":"A"}
	for i in range(len(s)):
		rc.append(rcMap[s[i]])
	return "".join(rc[::-1])

refSeq = {}
querySeq = {}

fp = open(sys.argv[1]) 
for line in fp:
	header = line.rstrip()
	seq = fp.readline().rstrip()
	refSeq[header[1:]] = seq
fp.close()

fp = open(sys.argv[2])
for line in fp:
	header = line.rstrip()
	seq = fp.readline().rstrip()
	querySeq[header.split()[0][1:]] = seq
fp.close()

fp = open(sys.argv[3])
for line in fp:
	if (line[0] == "@"):
		continue
	cols = line.rstrip().split("\t")
	blocks = CigarToBlocks(cols[5])
	start = int(cols[3])
	flag = int(cols[1])
	
	query = querySeq[cols[0]] 
	if (flag & 0x10):
		query = ReverseComplement(query)
	ref = refSeq[cols[2]]

	matchableLength = 0
	for b in blocks:
		if (b[1] == "M"):
			matchableLength += b[0]
	if (matchableLength < coverageThreshold * len(query)):
		continue

	refPos = start - 1
	queryPos = 0 
	vcf = [] # position (0-based) on query, ref, query
	#print(refPos, queryPos, blocks)
	for b in blocks:
		if (b[1] == "S"):
			queryPos += b[0]
			#refPos += b[0]
		elif (b[1] == "M"):
			# search for mismatches
			for i in range(b[0]):
				if (ref[refPos + i] != query[queryPos + i]):
					vcf.append([queryPos + i, query[queryPos + i], ref[refPos + i]])
			queryPos += b[0]
			refPos += b[0]
		elif (b[1] == "I" or b[1] == "H"):
			vcf.append([queryPos, query[queryPos:(queryPos + b[0])], "."])
			queryPos += b[0]
		elif (b[1] == "D"):
			vcf.append([queryPos, ".", ref[refPos:(refPos + b[0])]])
			refPos += b[0]
		elif (b[1] == "N"):
			refPos += b[0]
		else:
			print("Unknown CIGAR op " + b[1])
			sys.exit(1)
	if (len(vcf) == 0): # for the perfect alignment
		vcf.append([-1, ".", "."]) # -1 is definitely outside any interesting region
	# Reverse complement
	queryLen = len(query)
	strand = "+"
	if (flag & 0x10):
		for i in range(len(vcf)):
			vcf[i][0] = queryLen - 1 - vcf[i][0]
			if (vcf[i][1] != "."):
				vcf[i][1] = ReverseComplement(vcf[i][1])
			if (vcf[i][2] != "."):
				vcf[i][2] = ReverseComplement(vcf[i][2])
		vcf = vcf[::-1]
		strand = "-"
	
	# output
	for i in range(len(vcf)):
		print("\t".join([cols[0], str(vcf[i][0]), ".", vcf[i][1], vcf[i][2], strand]))
fp.close()
