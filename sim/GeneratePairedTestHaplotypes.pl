#!/usr/bin/env perl

use strict ;
use warnings ;

die "usage: a.pl kir_seq.fa num_of_tests > testcases.tsv \n" if (@ARGV == 0) ;

my %kirSeq ;
my %kirAlleles ;

srand(17) ;

open FP, $ARGV[0] ;
while (<FP>)
{
	chomp ;
	my $gene = (split /\s+/, substr($_, 1))[0] ;
	my $mainGene = (split /\*/, $gene)[0] ;
	my $seq = <FP> ;
	chomp $seq ;
	$kirSeq{$gene} = $seq ;
	push @{$kirAlleles{$mainGene}}, $gene ;
}
close FP ;

# Select between 6 to 10 genes for each haplotype
my $t = 0 ;
for ($t = 0 ; $t < $ARGV[1] ; ++$t)
{
	my $n = int(rand(5)) + 6 ;
	my $i = 0 ;
	my @list = keys %kirAlleles ;
	my $total = scalar(@list) ;
	for ($i = 0 ; $i < $n ; ++$i) 
	{
		my $j = $i + int(rand($total - $i)) ;
		($list[$i], $list[$j]) = ($list[$j], $list[$i]) ;
	}

	# Select two allele from each gene
	my @selectedGene ;
	for ($i = 0 ; $i < $n ; ++$i)
	{
		my $j = int(rand(scalar(@{$kirAlleles{$list[$i]}}))) ;
		push @selectedGene, ${$kirAlleles{$list[$i]}}[$j]	;
		$j = int(rand(scalar(@{$kirAlleles{$list[$i]}}))) ;
		push @selectedGene, ${$kirAlleles{$list[$i]}}[$j]	;
	}
	print join("\t", @selectedGene), "\n" ;
}
