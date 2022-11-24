#!/bin/perl

use strict ;
use warnings ;

die "usage: a.pl haplotype_list [sample_id]\n" if (@ARGV == 0) ;

open FP, $ARGV[0] ;
my $i = 0 ;
my $j ;
while (<FP>)
{
	chomp ;
	my %genotype ;
	my @cols = split ;
	foreach my $c (@cols)
	{
		my $gene = (split /\*/, $c)[0] ;
		my $majorAllele = substr($c, 0, length($gene) + 4) ;
		push @{$genotype{$gene}}, $majorAllele ;
	}

	if (defined $ARGV[1] && $i != $ARGV[1])
	{ 
		++$i ;
		next ;
	}

	open FPtype, "sample_${i}_genotype.tsv" ;

	my $correctCnt = 0 ;
	my $wrongCnt = 0 ;
	my %prediction ;
	my $FN = 0 ;
	my $FP = 0 ;
	my $T = scalar(@cols) ;
	my $P = 0 ;
	while (<FPtype>)
	{
		chomp ;
		my @cols = split ;
		#print join("_", @cols), "\n" ;
		#if (scalar(@cols) >= 5)
		@{$prediction{$cols[0]}} = () ;
		if ($cols[1] >= 1)
		{
			push @{$prediction{$cols[0]}}, ($cols[2], $cols[4]) ;
		}
		if ($cols[1] >= 2)
		{
			push @{$prediction{$cols[0]}}, ($cols[5], $cols[7]) ;
		}
		$P += $cols[1] ;
	}

	foreach my $gene (keys %prediction)
	{
		my @pred = @{$prediction{$gene}} ;
		my @truth ;
		if (defined $genotype{$gene})
		{
			@truth = @{$genotype{$gene}} ;
		}
		my $foundBad = 0 ;
		my @predCovered ;
		for ($j = 1 ; $j <= 3 ; $j += 2)
		{
			if (defined $pred[$j])
			{	
				if ($pred[$j] > 1)
				{
					push @predCovered, 0 ;
				}
				else
				{
					push @predCovered, 1 ;
				}
			}
		}

		foreach my $truthAllele (@truth)
		{
			my $tmp = $truthAllele ;
			$tmp =~ s/\*/\\\*/ ;
			my $match = 0 ;
			for ($j = 0 ; $j <= 2 ; $j += 2)
			{
				if (defined $pred[$j] && $pred[$j] =~ /$tmp/) 
				{
					++$match ;
					$predCovered[$j / 2] = 1 ;
				}
			}
			if ($match == 0) 
			{
				print("FN: ", $truthAllele, "\n") ;
				++$FN ;
				$foundBad = 1 ;
			}
		}

		foreach my $p (@predCovered)
		{
			if ($p == 0) 
			{
				$foundBad = 1 ;
				++$FP ;
				print("FP: ", $gene, "\n");
			}
		}

		if ($foundBad == 1)
		{
			++$wrongCnt ;
		}
		else
		{
			++$correctCnt ;
		}
	}

	close FPtype ;
	#print "$i $correctCnt $wrongCnt ", $correctCnt / ($correctCnt + $wrongCnt), "\n";
	print "$i $FP $P $FN $T ", 1 - $FN/$T, " ", 1 - $FP/$P, "\n";
	++$i ;
}
close FP ;
