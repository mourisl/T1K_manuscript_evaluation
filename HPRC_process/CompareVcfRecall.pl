#!/bin/perl

use strict ;
use warnings ;

die "usage: a.pl sample_list.out [genome_path]\n" if (@ARGV == 0) ;

my $genomePath = "Genome" ;
if (scalar(@ARGV) > 1)
{
	$genomePath = $ARGV[1] ;
}

open FP, $ARGV[0] ;
my %sample ;
while (<FP>)
{
	chomp ; 
	$sample{$_} = 1 ;
}
close FP ;

my $path="./" ;
my $overallMatchCnt = 0 ;
my $overallCnt = 0 ;
foreach my $p (keys %sample)
{
	my $a1 = "$path/$genomePath/${p}.1_exon.vcf" ;
	my $l1 = "$path/$genomePath/${p}.1_bwa_exon_genotype.out" ;
	my $a2 = "$path/$genomePath/${p}.2_exon.vcf" ;
	my $l2 = "$path/$genomePath/${p}.2_bwa_exon_genotype.out" ;
	my $b = "$path/Illumina/${p}_allele.vcf" ;
	my $lb = "$path/Illumina/${p}_allele.tsv" ;

	#my $result = `perl /liulab/lsong/projects/kir/kir/PairSampleCompare.pl $path/$a/${a}_genotype.tsv $path/$b/${b}_genotype.tsv | grep -v ^K` ;
	`cat $a1 $a2 > tmp.vcf` ;
	`grep PASS $b > tmp_truth.vcf` ;
	#my $result = `python3 $path/ValidateVcf.py $b $lb tmp.vcf` ;
	my $result = `python3 $path/ValidateVcf.py tmp_truth.vcf $lb tmp.vcf` ;
	chomp $result ;
	my @cols = split /\s+/, $result ;
	print($p, "\t", $cols[2], "\t", $cols[3], "\t", $cols[5], "\n") ;
	$overallMatchCnt += $cols[2] ;
	$overallCnt += $cols[3] ;
}
print("Overall\t", $overallMatchCnt, "\t", $overallCnt, "\t", $overallMatchCnt/$overallCnt, "\n") ;
