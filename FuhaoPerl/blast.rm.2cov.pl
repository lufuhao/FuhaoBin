#!/usr/bin/perl
use strict;
use warnings;

our $input=$ARGV[0];
chomp $input;
open (INPUT, "$input") || die "Can not input\n";
our $output=$ARGV[1];
chomp $output;
open (OUTPUT, ">>$output") || die "Can not output\n";

our ($input_line, @input_comp, $input_line_rm);
while ($input_line=<INPUT>){
	chomp $input_line;
	@input_comp=split(/\t/, $input_line);
	pop @input_comp;
	pop @input_comp;
	$input_line_rm=join ("\t", @input_comp);
	print OUTPUT $input_line_rm."\n";
} 

close INPUT;
close OUTPUT;
