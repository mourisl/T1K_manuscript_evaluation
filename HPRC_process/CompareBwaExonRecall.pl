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

my $path="/liulab/lsong/projects/kir/kir/HPRC" ;
foreach my $p (keys %sample)
{
	my $a1 = "$path/$genomePath/${p}.1_bwa_exon_genotype.out" ;
	my $a2 = "$path/$genomePath/${p}.2_bwa_exon_genotype.out" ;
	my $b = "$path/Illumina/${p}/${p}_genotype.tsv" ;

	#my $result = `perl /liulab/lsong/projects/kir/kir/PairSampleCompare.pl $path/$a/${a}_genotype.tsv $path/$b/${b}_genotype.tsv | grep -v ^K` ;
	my $result = `cat $a1 $a2 | python3 /liulab/lsong/projects/kir/kir/PairSampleCompare.py /dev/stdin $b --formata list --eval 1` ;
	foreach my $line (split /\n/, $result)
	{
		my @cols = split /\s+/, $line ;
		#my $nomatchCnt = $cols[0] ;
		#my $wrongfilterCnt = $cols[1] ;
		#my $totalCnt = $cols[2] ;
		print($p, "\t", $cols[0], "\t", $cols[1], "\t", $cols[2], "\n") ;
	}
}
