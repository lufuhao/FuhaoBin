#!/usr/bin/env perl
use strict;
use warnings;
use FuhaoPerl5Lib::FileKit qw/SubsetExtraction/;
use constant USAGE =><<EOH;

usage: $0 infile ID.file outfile

v20160805

EOH
die USAGE if (scalar(@ARGV) !=3 or $ARGV[0] eq '-h' or $ARGV[0] eq '--help');



my $input=$ARGV[0];
my $idfile=$ARGV[1];
my $output=$ARGV[2];
unless (SubsetExtraction($input, $idfile, $output)) {
	print STDERR "Error: extraction filed\n";
	exit 1;
}
