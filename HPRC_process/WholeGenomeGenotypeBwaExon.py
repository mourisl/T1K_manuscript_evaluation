#!/usr/bin/env python3
import sys
import subprocess

if (len(sys.argv) <= 1):
	print("Usage: python3 a.py minimap2.paf genome.fa kir_ref.fa [-n diff -o prefix --saveall]")
	sys.exit(1) 

def ParseAlleleName(n):
	cols = n.split("*")
	return cols[0], cols[0] + "*" + cols[1][0:3], n

def RangeOverlap(r1, r2):
	overlapSize = 0
	if (r1[0] != r2[0]):
		return False

	if (r2[2] <= r1[1] or r2[1] >= r1[2]):
		return False
	elif (r1[1] <= r2[1] and r1[2] >= r2[2]):
		overlapSize = r2[2] - r2[1] 
	elif (r2[1] <= r1[1] and r2[2] >= r1[2]):
		overlapSize = r1[2] - r1[1] 
	elif (r1[1] >= r2[1]):
		overlapSize = r2[2] - r1[1]
	elif (r1[1] < r2[1]):
		overlapSize = r1[2] - r2[1] 
	#print(r1, r2, overlapSize)
	# Similarity >= 0.95
	if (overlapSize > 0.95 * (r1[2] - r1[1])
			or overlapSize > 0.95 * (r2[2] - r2[1])):
		return True 
	return False

secondaryHitDiff = 0
outputPrefix = ""
saveAll = False
for i in range(4, len(sys.argv)):
	if (sys.argv[i] == "-n"):
		secondaryHitDiff = int(sys.argv[i+1])
	elif (sys.argv[i] == "-o"):
		outputPrefix = sys.argv[i+1]
	elif (sys.argv[i] == "--saveall"):
		saveAll = True

ranges = []
fp = open(sys.argv[1])
for line in fp:
	cols = line.rstrip().split()
	allele = cols[0]
	# Coverage >= 99%
	if (int(cols[3]) - int(cols[2]) < 0.99 * int(cols[1])):
		continue
	
#KIR2DL1*0010101	14738	12	14728	-	HG002#1#JAHKSE010000050.1	19771003	958929	973648	13726	14719	60	tp:A:P	cm:i:1329	s1:i:13725	s2:i:10386	dv:f:0.0017	rl:i:944
	grange = [cols[5], int(cols[7]), int(cols[8]), int(cols[9]), float(cols[9]) / float(cols[10]), allele]
	addNew = True
	for i, r in enumerate(ranges):
		if (RangeOverlap(grange, r)):
			if (grange[3] > r[3] or (grange[3] == r[3] and grange[4] > r[4])):
				ranges[i][:] = grange[:]
			addNew = False
			break
	if (addNew):
		ranges.append(grange)
fp.close()

#for i, r in enumerate(ranges):
#	print(r[-1])
genome = {}
fp = open(sys.argv[2])
seq = ""
chrom = ""
for line in fp:
	line = line.rstrip()
	if (line[0] == '>'):
		if (chrom != ""):
			genome[chrom] = seq
		seq = ""
		chrom = line.split()[0][1:]
	else:
		seq += line
fp.close()
if (chrom != ""):
	genome[chrom] = seq

fp = open(sys.argv[3])
alleleSeq = {}
geneToAlleles = {}
for line in fp:
	header = line.rstrip()
	allele = header.split()[0][1:]
	seq = fp.readline().rstrip()
	alleleSeq[allele] = seq
	gene = ParseAlleleName(allele)[0]
	if gene not in geneToAlleles:
		geneToAlleles[gene] = {}
	geneToAlleles[gene][allele] = 1 
fp.close()


fpBwaQuery = open("bwa_query.fa", "w")
for allele in alleleSeq:
	fpBwaQuery.write(">" + allele + "\n" + alleleSeq[allele] + "\n")
fpBwaQuery.close()

# For each interval and gene, validate with blastn
calledAlleles = []
vcfFile = outputPrefix + "_exon.vcf" ;
subprocess.run("echo -n > " + vcfFile, shell=True)
for i, r in enumerate(ranges):
	calledAlleles.append([])
	fpBwaRef = open("bwa_ref.fa", "w")
	start = r[1] - 100
	if (start < 0):
		start = 0
	end = r[2] + 100
	fpBwaRef.write(">ref\n" + genome[r[0]][start:end])
	gene = ParseAlleleName(r[-1])[0]
	fpBwaRef.close()
	subprocess.run("bwa index -p bwa_idx bwa_ref.fa > /dev/null", shell=True)
	subprocess.run("bwa mem -t 8 bwa_idx bwa_query.fa > bwa.sam", shell=True)
	subprocess.run("python3 /liulab/lsong/projects/kir/kir/HPRC/CigarToVcf.py bwa_ref.fa bwa_query.fa bwa.sam > bwa.vcf", shell=True)
	subprocess.run("python3 /liulab/lsong/projects/kir/kir/HPRC/ComputeExonicDifference.py "+sys.argv[3]+" bwa.vcf > bwa_exondiff.out", shell=True)
	fp = open("bwa_exondiff.out")
	exonDiffs = []
	for line in fp:
		cols = line.rstrip().split()
		exonDiffs.append([cols[0], int(cols[1])])
	fp.close()
	if (len(exonDiffs) == 0):
		continue
	exonDiffs.sort(key = lambda x:x[1])
	calledAlleles[i].append(exonDiffs[0][0])
	for j in range(1, len(exonDiffs)):
		if (exonDiffs[j][1] <= exonDiffs[0][1] + secondaryHitDiff):
			calledAlleles[i].append(exonDiffs[j][0])
		else:
			break
	subprocess.run("python3 /liulab/lsong/projects/kir/kir/HPRC/FindExonVcf.py "+sys.argv[3]+" bwa.vcf " + " ".join(calledAlleles[i]) + ">>" + vcfFile, shell=True)
#for allele in sorted(calledAlleles.keys()):
#	print(allele)
for i, alleles in enumerate(calledAlleles):
	print(",".join(alleles))

