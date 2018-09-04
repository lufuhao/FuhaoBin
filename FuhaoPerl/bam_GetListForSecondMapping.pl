#!/usr/bin/env perl
use strict;
use warnings;
use FuhaoPerl5Lib::BamKit qw/GetListForSecondMapping/;
use constant USAGE =><<EOH;

usage: $0 BAM_list_file CHROM_keep_list_file CHROM_exclude_list_file read.out.list ref.out.list

v20180904

EOH
die USAGE if (scalar(@ARGV) !=5 or $ARGV[0] eq '-h' or $ARGV[0] eq '--help');

unless (GetListForSecondMapping($ARGV[0],$ARGV[1],$ARGV[2],$ARGV[3],$ARGV[4])) {
	die "Error: GetListForSecondMapping running error\n";
}
