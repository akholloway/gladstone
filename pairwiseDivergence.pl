#!/usr/bin/perl

use strict;
use warnings;

######################### MAIN  ######################

		
open(IN, $ARGV[0]);
$/='>';
my ($ID,$SEQ) = readDATA();
$/ = '\n';
calcDiv();
close IN;
exit;


#########################################################################
				
							# SUBROUTINES #

#########################################################################

sub readDATA{
	
	my $data;
	my ($ID,$SEQ);
	
	<IN>;
	while($data = <IN>){
		
		($data =~ /([+A-Za-z0-9]+)/);
		if($1){
			push(@{$ID}, $1);
			$data =~ s/.*?\n//;	# remove defline
			$data =~ s/\s+//g;	# remove whitespace
			$data =~ s/>//g;	# remove trailing >
			$data =~ tr/[a-z]/[A-Z]/;
			push(@{$SEQ}, $data);
		}	
	}
	return($ID,$SEQ);
}
		


#########################################################################

sub calcDiv{		##

	my $numNT = length $SEQ->[0];
		
	for(my $i=0; $i < scalar @{$SEQ}; $i++){
		for(my $j=$i+1; $j < scalar @{$SEQ}; $j++){
			my $POS=0;
			my $numsites=0;
			my $diff=0;
			while($POS < $numNT){
				my $nt1 = substr($SEQ->[$i], $POS, 1);
				my $nt2 = substr($SEQ->[$j], $POS, 1);
				unless($nt1 eq "N" || $nt2 eq "N"){
					if($nt1 ne $nt2){
						$diff++;
					}
					$numsites++;
				}
				$POS++;	
			}
			if($numsites){
				my $div = $diff/$numsites;
				print $ID->[$i], "\t", $ID->[$j], "\t", $numsites, "\t", $div, "\n";
			}
		}
	}
}

###############################################################
# Sub: Usage
############################################################### 

sub usage {
	my $u = <<END;
	
calculates pairwise divergence for all seqs in a fasta file
includes indels in estimate of divergence
no Jukes-Cantor correction

perl pairwiseDiv.pl	file

END

 print $u;

 exit(1);
 
}


################ EOF ###################
