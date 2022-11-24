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

	#my $result = `perl /liulab/lsong/projects/kir/kir/PairSampleCompare.pl $path/$a/${a}_genotype.tsv $path/$b/${b}_genotype.tsv | grep -v ^K` ;
	`cat $a1 $a2 > tmp.vcf` ;
	`cat $l1 $l2 > tmp.list` ;
	my $result = `python3 $path/ValidateVcf.py tmp.vcf tmp.list $b` ;
	chomp $result ;
	my @cols = split /\s+/, $result ;
	print($p, "\t", $cols[0], "\t", $cols[1], "\t", $cols[4], "\n") ;
	$overallMatchCnt += $cols[0] ;
	$overallCnt += $cols[1] ;
}
print("Overall\t", $overallMatchCnt, "\t", $overallCnt, "\t", $overallMatchCnt/$overallCnt, "\n") ;
