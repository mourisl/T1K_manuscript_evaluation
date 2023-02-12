#!/bin/perl

use strict ;
use warnings ;


die "usage: a.pl hla.fa hla_freq.tsv > filtered.fa\n" if (@ARGV == 0) ;

sub GetMajorAllele
{
	my @cols = split /:/, $_[0] ;
	if (scalar(@cols) <= 1) 
	{
		return $_ ;
	}
	else
	{
		return join(":", @cols[0..1]) ;
	}
}

open FP, $ARGV[1] ;
my %majorAllele ;
my $line = <FP> ;
while (<FP>)
{
	chomp ;
	$line = $_ ;
	my @cols = split ;
	#my $s = $cols[0] ; #GetMajorAllele($cols[0]) ;
	my $s = GetMajorAllele($cols[0]) ;
	$s = "HLA-".$s if (substr($s, 0, 4) ne "HLA-");
	$majorAllele{$s} = 1 ;
}
close FP ;

open FP, $ARGV[0] ;
while (<FP>)
{
	my $header = $_ ;
	chomp $header ;
	my $seq = <FP> ;
	chomp $seq ;
	my @cols = split /\s/, substr($header, 1) ;
	my $tmp = $cols[0] ;
	my $s = GetMajorAllele($tmp) ;
	if (defined $majorAllele{$s})
	{
		print "$header\n$seq\n" ;
	}
}
close FP ;
