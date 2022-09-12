#!/usr/bin/env python3

import sys
import re

if (len(sys.argv) <= 1):
	# --formata STRING: format for file a - "tonek", "list", "arcasHLA"
	# --formatb STRING: format for file b - "tonek", "list", "arcasHLA"
	# --eval INT: bits represent to evalue a or b: first bit a, second bit b.
	# --type STRING: KIR or HLA
	print("usage: a.py a_genotype b_genotype [--formata tonek --formatb tonek --eval 3 --type KIR]")
	sys.exit(1)



qualityThreshold = 0
genotypeA = {}
genotypeB = {}
genotypeACluster = {} # most useful for list, each line will corresponds to a cluster of alleles
genotypeBCluster = {} 
alist = False
formatA = "tonek"
formatB = "tonek"
geneType = "KIR"
evalSamples = 3

for i in range(3, len(sys.argv)):
	if (sys.argv[i][0] != '-'):
		continue
	if (sys.argv[i] == "--formata"):
		formatA = sys.argv[i + 1]
	if (sys.argv[i] == "--formatb"):
		formatB = sys.argv[i + 1]
	elif (sys.argv[i] == "--eval"):
		evalSamples = int(sys.argv[i + 1])
	elif (sys.argv[i] == "--type"):
		geneType = sys.argv[i + 1]

def ParseAlleleName(n):
	cols = n.split("*")
	if (geneType == "HLA"):
		subcols = cols[1].split(":")
		return cols[0], cols[0] + "*" + ":".join(subcols[0:3]), n
	else:
		return cols[0], cols[0] + "*" + cols[1][0:3], n

allGenes = {}
for k in range(2):
	fp = open(sys.argv[k + 1])
	genotype = genotypeA
	genotypeCluster = genotypeACluster
	fmt = formatA 
	if (k == 1):
		genotype = genotypeB 
		genotypeCluster = genotypeBCluster
		fmt = formatB
	
	lineCnt = 0
	for line in fp:
		cols = line.rstrip().split()
		if (fmt == "tonek"):
			gene = cols[0]
			if (gene not in genotype):
				genotype[gene] = {}
				genotypeCluster[gene] = {}
			if (gene not in allGenes):
				allGenes[gene] = 1
			if (int(cols[1]) >= 1):
				for field in cols[2].split(","):
					gene, majorAllele, allele = ParseAlleleName(field)
					genotype[gene][majorAllele] = float(cols[4])
					genotypeCluster[gene][majorAllele] = [0]
			if (int(cols[1]) >= 2):
				for field in cols[5].split(","):
					gene, majorAllele, allele = ParseAlleleName(field)
					genotype[gene][majorAllele] = float(cols[7])
					genotypeCluster[gene][majorAllele] = [1]
		elif (fmt == "list"):
			#gene, majorAllele, allele = ParseAlleleName(line.rstrip())
			cols = line.rstrip().split(",")
			
			for a in cols:
				if (len(a) == 0):
					continue
				gene, majorAllele, allele = ParseAlleleName(a)
				if (gene not in genotype):
					genotype[gene] = {}
					genotypeCluster[gene] = {}
				if (gene not in allGenes):
					allGenes[gene] = 1
				genotype[gene][majorAllele] = 60
				if (majorAllele not in genotypeCluster[gene]):
					genotypeCluster[gene][majorAllele] = []
				genotypeCluster[gene][majorAllele].append(lineCnt)

		elif (fmt == "arcasHLA"):
			patterns = re.findall('\".*?\"', line)	
			i = 0
			for p in patterns:
				if "*" not in p:
					continue
				allele = p[1:-1] # remove the " symbols
				gene, majorAllele, allele = ParseAlleleName(geneType + "-" + allele)
				if (gene not in genotype):
					genotype[gene] = {}
					genotypeCluster[gene] = {}
				if (gene not in allGenes):
					allGenes[gene] = 1
				genotype[gene][majorAllele] = 60
				genotypeCluster[gene][majorAllele] = [i]
				i += 1
		else:
			print("Unknown format", fmt)
			sys.exit(1)
		lineCnt += 1
	fp.close()


for g in allGenes:
	if (g not in genotypeA):
		genotypeA[g] = {}
		genotypeACluster[g] = {}
	if (g not in genotypeB):
		genotypeB[g] = {}
		genotypeBCluster[g] = {}

def AlleleMatch(allele, alleleScore, otherGenotype, otherGenotypeCluster, used):
	if (alleleScore <= qualityThreshold):
		return 3
	if (allele not in otherGenotype):
		return 2
	
	# check whether it is aligned.
	filterFlag = True
	for k in range(len(otherGenotypeCluster[allele])):
		if (otherGenotypeCluster[allele][k] not in used):
			used[otherGenotypeCluster[allele][k]] = 1
			filterFlag = False
			break

	if (filterFlag):
		return 2

	if (otherGenotype[allele] > qualityThreshold):
		return 0
	else:
		return 1

def EvaluateGeneInOne(gA, gACluster, gB, gBCluster, matchCnt, allowRepetitiveMatch=False):
	matchResult = []
	used = {}
	
	if (len(gA) == 0):
		return

	aClusterCnt = max([max(v) for v in gACluster.values()]) + 1
	
	aClusterAlleles = []
	for i in range(aClusterCnt):
		aClusterAlleles.append([])
	for allele in gA:
		for cluster in gACluster[allele]:
			aClusterAlleles[cluster].append(allele)
	
	for cluster in range(aClusterCnt): 
		if (len(aClusterAlleles[cluster]) == 0):
			continue
		if (allowRepetitiveMatch):
			used = {}
		matchResult.append(
				min([AlleleMatch(allele, gA[allele], gB, gBCluster, used) 
					for allele in aClusterAlleles[cluster]])) 
	used = {}
	matchResultRv = []
	for cluster in reversed(range(aClusterCnt)):
		if (len(aClusterAlleles[cluster]) == 0):
			continue
		if (allowRepetitiveMatch):
			used = {}
		matchResultRv.append(
				min([AlleleMatch(allele, gA[allele], gB, gBCluster, used) 
					for allele in aClusterAlleles[cluster]])) 
	
	# it is not correct when there is more than 2 group clusters in A, 
	# but we should only consider at most 2 alleles per gene for now.
	if (sum(matchResult) <= sum(matchResultRv)):
		for r in matchResult:
			matchCnt[r] += 1
	else:
		for r in matchResultRv:
			matchCnt[r] += 1


totalMatchResult = [0, 0, 0, 0]
for gene in sorted(allGenes.keys()):
	matchCnt = [0, 0, 0, 0] # 0-match, 1-weak match, 2-mismatch, 3-no use	
	if (evalSamples & 1 == 1):
		EvaluateGeneInOne(genotypeA[gene], genotypeACluster[gene], genotypeB[gene], genotypeBCluster[gene], matchCnt, formatA == "list")
	if (evalSamples & 2 == 2):
		EvaluateGeneInOne(genotypeB[gene], genotypeBCluster[gene], genotypeA[gene], genotypeACluster[gene], matchCnt, formatB == "list")
	
	accuracy = 1
	weakAccuracy = 1 
	total = sum(matchCnt[0:3])
	if (total > 0):
		accuracy = matchCnt[0] / total
		weakAccuracy = (matchCnt[0] + matchCnt[1]) / total
		
	print("\t".join([gene] + ["%.4f"%accuracy, "%.4f"%weakAccuracy] + [str(x) for x in matchCnt]))
	for i in range(len(matchCnt)):
		totalMatchResult[i] += matchCnt[i]


accuracy = 1
weakAccuracy = 1 
total = sum(totalMatchResult[0:3])
if (total > 0):
	accuracy = totalMatchResult[0] / total
	weakAccuracy = (totalMatchResult[0] + totalMatchResult[1]) / total
print("\t".join(["Overall"] + ["%.4f"%accuracy, "%.4f"%weakAccuracy] + [str(x) for x in totalMatchResult]))
