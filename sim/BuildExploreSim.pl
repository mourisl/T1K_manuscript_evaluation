#!/usr/bin/env perl

use strict ;
use warnings ;

die "usage: a.pl kir.fa haplotype.list\n" if (@ARGV == 0) ;

sub system_call
{
	print STDERR "[".localtime()."] SYSTEM CALL: ".join(" ",@_)."\n";
	system(@_) == 0
		or die "system @_ failed: $?";
	#print STDERR " finished\n";
} 

my %kirSeq ;
open FP, $ARGV[0] ;
while (<FP>)
{
	chomp ;
	my $gene = (split /\s+/, substr($_, 1))[0] ;
	my $seq = <FP> ;
	chomp $seq ;
	$kirSeq{$gene} = $seq ;
}
close FP ;

open FP, $ARGV[1] ;
my $i = 0 ;
my $mason = "mason-0.1.2-Linux-x86_64/bin/mason" ;
while (<FP>)
{
	chomp ;
	my @cols = split ;
	open FPout, ">ref.fa" ;
	foreach my $gene (@cols)
	{
		print FPout ">$gene\n".$kirSeq{$gene}."\n" ;
	}
	close FPout ;
	
	my $readCnt = 50 * scalar(@cols) ;
	my $otherOption = "-hs 0 -hi 0 -pi 0 -pd 0 -pmmb 0.0005 -pmm 0.001 -pmme 0.003 -nN" ;
	system_call("$mason illumina -i -s 17 -sq -mp -n 100 -ll 500 $otherOption -o sample_${i}.fq -N $readCnt ref.fa 2> mason.log") ;
	unlink("sample_${i}.fq.sam") ;
	++$i ;
}
close FP ;
