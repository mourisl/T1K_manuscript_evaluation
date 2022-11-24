#!/bin/perl

use strict ;
use warnings ;

die "usage: a.pl sample_list.out\n" if (@ARGV == 0) ;

open FP, $ARGV[0] ;
my %sample ;
while (<FP>)
{
	chomp ; 
	$sample{$_} = 1 ;
}
close FP ;


my %geneResult ;
my $path="./" ;
foreach my $p (keys %sample)
{
	my $a1 = "$path/Genome/${p}.1_bwa_exon_genotype.out" ;
	my $a2 = "$path/Genome/${p}.2_bwa_exon_genotype.out" ;
	my $b = "$path/Illumina/${p}_genotype.tsv" ;
	
#my $result = `perl /liulab/lsong/projects/kir/kir/PairSampleCompare.pl $path/$a/${a}_genotype.tsv $path/$b/${b}_genotype.tsv | grep -v ^K` ;
	#`cat $a1 $a2 > tmp.out` ;
	my $result = `cat $a1 $a2 | python3 PairSampleCompare.py /dev/stdin $b --formata list --eval 2` ;
	my @lines = split /\n/, $result ;
	foreach my $line (@lines)
	{
		my @cols = split /\s+/, $line ;
		my $gene = $cols[0] ;
		if (!defined $geneResult{$gene}) 
		{
			@{$geneResult{$gene}} = (0, 0, 0) ;
		}
		${$geneResult{$gene}}[0] += $cols[3] ;
		${$geneResult{$gene}}[1] += $cols[4] ;
		${$geneResult{$gene}}[2] += $cols[5] ;
	}
	#my $nomatchCnt = $cols[0] ;
	#my $wrongfilterCnt = $cols[1] ;
	#my $totalCnt = $cols[2] ;
}
foreach my $gene (sort(keys %geneResult))
{
	my $accuracy = 1 ;
	my $weakAccuracy = 1 ;
	my $total = 0 ;
	for (my $i = 0 ; $i < 3 ; ++$i)
	{
		$total += ${$geneResult{$gene}}[$i] ;
	}

	if ($total > 0)
	{
		$accuracy = ${$geneResult{$gene}}[0] / $total ;
		$weakAccuracy = (${$geneResult{$gene}}[0] + ${$geneResult{$gene}}[1]) / $total ;
	}
	print(join("\t", ($gene, $accuracy, $weakAccuracy, @{$geneResult{$gene}})), "\n") ;
}
