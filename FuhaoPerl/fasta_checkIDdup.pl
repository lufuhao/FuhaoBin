#!/usr/bin/env perl
use strict;
use warnings;
use FuhaoPerl5Lib::FastaKit qw/CheckFastaIdDup/;
use constant USAGE =><<EOH;

perl $0 fasta

v20171205

Description:
    Check if seq IDs are duplicated in a multi-fasta file

EOH
die USAGE if (scalar(@ARGV) !=1 or $ARGV[0] eq '-h' or $ARGV[0] eq '--help');

my $fasta=$ARGV[0];
unless (CheckFastaIdDup($fasta)) {
	die "Error: errors detected\n";
}
else {
	print "Info: All seq IDs are unique\n";
	exit 0;
}
