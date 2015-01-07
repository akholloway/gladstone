#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Long;

######################### MAIN  ######################

#print OUT "FILE\tAVGLINES\tSITES\tPI\n";

my ($SEQ);

open IN, $ARGV[0];
$/='>';
($SEQ) = readDATA($ARGV[0]);
$/ = '\n';
calcPI();
close IN;
exit;


#########################################################################
				
							# SUBROUTINES #

#########################################################################

sub readDATA{
	
	my $data;
	my $SEQ;
	
	<IN>;
	while($data = <IN>){
		
		$data =~ s/.*?\n//;	# remove defline
		$data =~ s/\s+//g;	# remove whitespace
		$data =~ s/>//g;	# remove trailing >
		$data =~ tr/a-z/A-Z/;
		push(@{$SEQ}, $data);
	}
	return($SEQ);
}
		


#########################################################################

sub calcPI{		##pi = (sum((1-(p^2 + q^2 + r^2 + s^2)))(n/(n-1)) / number of sites
	
	my $POS=0;
	my $numNT = length $SEQ->[0];
	my $pi=0;
	my $numSEQS = scalar @{$SEQ};
	my $numsites=0;
	my $pipersite=0;
	my $lines=0;
	my $het=0;
	
	while($POS < $numNT){
		my ($A, $C, $G, $T, $N) = (0,0,0,0,0);
		for(my $x=0;$x<$numSEQS;$x++){
			my $NT = substr($SEQ->[$x], $POS, 1);
			if($NT eq "A"){
				$A++;
			}
			elsif($NT eq "C"){
				$C++;
			}
			elsif($NT eq "G"){
				$G++;
			}
			elsif($NT eq "T"){
				$T++;
			}
			else{	#$NT is an N or -
				$N++;
			}
			#print $NT, "\t";
		}
		my $numLINES = $numSEQS-$N; ## gets the number of lines with data
		if($numLINES>1){
			$numsites++;
			
			my $freqA = $A/$numLINES;
			my $freqC = $C/$numLINES;
			my $freqG = $G/$numLINES;
			my $freqT = $T/$numLINES;
			$het = 1 - (($freqA*$freqA) + ($freqC*$freqC) + ($freqG*$freqG) + ($freqT*$freqT));
			$lines += $numLINES;
			$pi += ($het)*($numLINES/($numLINES-1));
			
		}
		$POS++;	
	}
	
	if($numsites){
		my $avglines = $lines/$numsites;
		my $cov = (($avglines*($avglines-1))/2)*$numsites;
		if($cov >= 100){
			$pipersite = $pi/$numsites;
			print $ARGV[0], "\t", $avglines, "\t", $numsites, "\t", $pipersite, "\n";
		}
	}
}



################ EOF ###################
