#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Long;

## could sort the seqs by total number of nt covered
## then would use the "best" sequences most often

######################### MAIN  ######################

my ($opt_d, $opt_o, $opt_s);

GetOptions('d:s' => \$opt_d,
			'o:s' => \$opt_o,
			's:s' => \$opt_s);

unless (-e $opt_d) {
        usage();
}

open (OUT, ">$opt_o") || die;
print OUT "FILE\tAVGLINES\tSITES\tPI\tTHETA\tTajimasD\n";

my ($SEQ);

##OPEN DIRECTORY WITH FASTA FILES
opendir(DIR, $opt_d);
my @allfiles = readdir DIR;
my $numfiles = scalar @allfiles;
my $a=0;
while ($a<$numfiles){			
	if($allfiles[$a] =~ /.fa$/){
		my $datafile = $opt_d.$allfiles[$a];
		open IN, $datafile;
		$/='>';
		($SEQ) = readDATA($datafile);
		$/ = '\n';
		calcPI();
		close IN;
	}
	$a++;
}

close OUT;
exit;


#########################################################################
				
							# SUBROUTINES #

#########################################################################

sub readDATA{
	
	my $data;
	my $SEQ;
	
	<IN>;
	while($data = <IN>){
		
		($data =~ /^(\S+);/);
		if($1 && $1 eq $opt_s){
			
			$data =~ s/.*?\n//;	# remove defline
			$data =~ s/\s+//g;	# remove whitespace
			$data =~ s/>//g;	# remove trailing >
			
			push(@{$SEQ}, $data);
		}	
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
	my $theta=0;
	my $thetaPerSite=0;
	my $lines=0;
	my $het=0;
	my $S=0;
	
	my @sampleSize;
	my $maxData=0;
	while($POS < $numNT){
		my $num=0;
		for(my $x=0;$x<$numSEQS;$x++){
			my $NT = substr($SEQ->[$x], $POS, 1);
			if($NT eq "A" || $NT eq "C" || $NT eq "G" || $NT eq "T"){
				$num++;
			}
		}
		$sampleSize[$num]++;
		$POS++;
	}
	for( my $ss = 0; $ss < @sampleSize; $ss++){
		if($sampleSize[$ss] && $sampleSize[$ss] > $maxData){ $maxData = $ss; }
	}
	
	$POS=0;
	$numNT = length $SEQ->[0];
	
	while($POS < $numNT){
		my ($A, $C, $G, $T, $N) = (0,0,0,0,0);
		my ($x, $numLINES)= (0,0);
		while($x < $numSEQS && $numLINES < $maxData){
			my $NT = substr($SEQ->[$x], $POS, 1);
			if($NT eq "A"){
				$A++;
				$numLINES++;
			}
			elsif($NT eq "C"){
				$C++;
				$numLINES++;
			}
			elsif($NT eq "G"){
				$G++;
				$numLINES++;
			}
			elsif($NT eq "T"){
				$T++;
				$numLINES++;
			}
			else{	#$NT is an N or -
				$N++;
			}
			$x++;
			#print $NT, "\t";
		}
		if($numLINES == $maxData){
			$numsites++;
			
			my $freqA = $A/$numLINES;
			my $freqC = $C/$numLINES;
			my $freqG = $G/$numLINES;
			my $freqT = $T/$numLINES;
			$het = 1 - (($freqA*$freqA) + ($freqC*$freqC) + ($freqG*$freqG) + ($freqT*$freqT));
			$lines += $numLINES;
			$pi += ($het)*($numLINES/($numLINES-1));
			if($het){
				$S++;
			}
		}
		$POS++;	
	}
	
	if($numsites){
		
		my $cov = (($maxData*($maxData-1))/2)*$numsites;
		if($cov >= 50 && $maxData > 2){
			my $D=0;
			if($S){
				my $a1=0;
				my $a2=0;
				my $n = $maxData;
				for(my $i=1; $i < $n; $i++){
					$a1 += 1/$i;
					$a2 += 1/($i*$i);
				}
				$theta = $S/$a1;
				my $b1 = ($n+1)/(3*($n-1));
				my $e1 = ($b1 - (1/$a1))/$a1;
				my $b2 = (2*(($n*$n) + $n + 3))/(9*$n*($n-1));
				my $c2 = $b2 - (($n+2)/($a1*$n)) + ($a2/($a1*$a1));
				my $e2 = $c2/(($a1*$a1) + $a2);
				my $d = $pi-$theta;
				my $sd = (sqrt(($e1*$S) + ($e2*$S*($S-1))));
#				if($sd == 0){ 
#					print $allfiles[$a]; 
#					print "\tpi-$pi\ttheta\t$theta\t$e1\t$e2\t$S\n";
#				}
				unless($sd == 0){ $D = $d / $sd; }
				$pipersite = $pi/$numsites;
				$thetaPerSite = $theta/$numsites;
			}
			print OUT $allfiles[$a], "\t", $maxData, "\t", $numsites, "\t", $pipersite, "\t", $thetaPerSite, "\t", $D, "\n";
		}
	}
}

###############################################################
# Sub: Usage
############################################################### 

sub usage {
	my $u = <<END;
	
processes multiple fasta files in a directory (-d)
calculates pi, theta, and Tajima's D on each file

pi_theta_TajimasD_maxData.pl	
		-d directory of fasta files with ending *.fa
		-s species with polymorphism data (fasta file format >species [whatever else on id line]
    -o output file

END

 print $u;

 exit(1);
 
}


################ EOF ###################
