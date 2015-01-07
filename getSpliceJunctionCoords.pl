#!/usr/bin/perl

use strict;
use warnings;

my ($c, $s, $e, $t);

# input is bed format with transcript ID

open(IN, $ARGV[0]);	#sorted by transcript, then by coord
while(my $data = <IN>){
	chomp $data;
	my @info = split(/\t/, $data);
	my ($chr, $start, $end, $transc) = ($info[0], $info[1], $info[2], $info[3]);
	if($c){
		if($t eq $transc){
			print "$chr\t$e\t$start\t$transc\n";
		}
	}
	$c=$chr; $s=$start; $e=$end; $t=$transc;
}

close IN;
exit;
