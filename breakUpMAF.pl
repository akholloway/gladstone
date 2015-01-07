#!/usr/bin/perl

use strict;
use warnings;

my ($chr) = ($ARGV[0] =~ /\/maf\/(\S+)\.maf/);
my $end = 10000000;
my $start = 1;
my $len = 0;
system `mkdir $ARGV[1]$chr`;

open(IN, $ARGV[0]);
my $info = <IN>;

$/ = "a score";
<IN>;
while(my $data = <IN>){
#	print $data;
	(my $score) = ($data =~ /(\=\-?\d+.+)/);
	$data =~ s/.*?\n//;	# remove defline
	$data =~ s/a score//;
	my (@d) = split(/\s+/, $data);
	my $beg = $d[2];
	unless ($len){ 
		$len=$d[5]; 
		open (OUT, ">$ARGV[1]$chr/$chr.$start-$end.maf");
		print OUT $info;
	}
	if($beg < $end){
		print OUT "\na score";
		print OUT $score, "\n", $data;
	}
	else{
		close OUT;
		$start = $end+1;
		$end += 10000000;
		if($end > $len){ $end = $len; }
		open(OUT, ">$ARGV[1]$chr/$chr.$start-$end.maf");
		print OUT $info;
		print OUT "\na score";
		print OUT $score, "\n", $data;
	}
}
close OUT;
exit;



################	EOF
