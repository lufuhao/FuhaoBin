#!/usr/bin/env perl
use warnings;
use strict;
use FuhaoPerl5Lib::FastqKit qw/FastqCompareIds/;
use constant USAGE =><<EOH;

perl $0 MateR1 MateR2 PATTERN prefix_output

Version: 20150508

Description:
    Compare Fastq1 and Fastq2, output ids for share unique1 and unique2

Pattern:
    Pattern1: '\\@(\\S+)\\/[12]\\s*\\S*'
    Pattern2: '\\@(\\S+)\\s*\\S*'

Output
    Shared:     prefix.shared
    R1 unique:  prefix.uniqueR1
    R2 unique:  prefix.uniqueR2

EOH
die USAGE if (scalar(@ARGV) != 4 or $ARGV[0] eq '-h' or $ARGV[0] eq '--help');


##input
my $in1=$ARGV[0];
my $in2=$ARGV[1];
my $pattern=$ARGV[2];
my $output_prefix=$ARGV[3];
#print $pattern, "\n"; #for test
###Output
my $output_share=$output_prefix.'.shared';
my $output_unique1=$output_prefix.'.uniqueR1';
my $output_unique2=$output_prefix.'.uniqueR2';
###

unless (FastqCompareIds($in1, $in2, $pattern, $output_share, $output_unique1, $output_unique2)) {
	die "Error: FastqCompareIds running error";
}
