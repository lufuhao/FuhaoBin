#!/usr/bin/env perl
use strict;
use warnings;
use FuhaoPerl5Lib::MiscKit qw/ListMergerMultiCols/;
use constant USAGE =><<EOH;

usage: $0 out.sum "undefined_value" in1 in2 ..

Desc: 
    Merge multiple (>=2) list

Example
    $0 out.sum "[undef]" In1 In2 ..

### In1
Name1	Value1
Name2	value2

### In2
Name1	Value3
Name4	value4

### out.sum
     	[In1] 	[IN2]
Name1	Value1	Value3
Name2	value2	[undef]
Name4	[undef]	value4



v20180924

EOH
die USAGE if (scalar(@ARGV) <4 or $ARGV[0] eq '-h' or $ARGV[0] eq '--help');

my $output=shift @ARGV;
my $undef_mark=shift @ARGV;
unless (ListMergerMultiCols($output, $undef_mark, \@ARGV)) {
	die "Error: ListMerger running";
}
