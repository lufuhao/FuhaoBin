#!/usr/bin/perl
use warnings;
use strict;

our $input=$ARGV[0];
our $output=$ARGV[1];

open (INPUT, "$input") || die "Can not open input file\n";
open (OUTPUT, ">>$output") || die "Can not write data to output file\n";
our $temp = 0;
our $count = 0;
our $pre_blast;
our $pre_target;

while (our $blasts= <INPUT>) {
	chomp $blasts;
#	print $blasts."\n";                             ###testline###
	our @blast=split (/\t/, $blasts);
#	print "@blast\n";                                ###testline###
	if ($temp == 0){
		$pre_blast=$blast[0];
		$pre_target=$blast[1];
#		print $pre_blast."\n".$pre_target."\n";   ###testline###
		$temp++;
		$count++;
	}
	elsif($pre_blast eq $blast[0]) {
		if ($blast[1] ne $pre_target) {
		$count++;}
	}
	else {
		print OUTPUT $pre_blast."\t".$count."\n";
		$pre_blast=$blast[0];
		$pre_target=$blast[1];
		$temp=1;
		$count=1;
	}
}

close INPUT;
close OUTPUT;
