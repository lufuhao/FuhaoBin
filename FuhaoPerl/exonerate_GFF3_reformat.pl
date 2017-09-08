#!/usr/bin/env perl
use strict;
use warnings;
use Data::Dumper qw/Dumper/;
use FuhaoPerl5Lib::GffKit qw /ExonerateGff3Reformat/;
use constant USAGE =><<EOH;

usage: $0 cdna.gff3 cds.gff3 final.gff3

v20170503

EOH
die USAGE if (scalar(@ARGV) !=3 or $ARGV[0] eq '-h' or $ARGV[0] eq '--help');
my $GffKit_success=1; my $GffKit_failure=0; my $GffKit_debug=0;



my $cdnain=shift @ARGV;
my $cdsin=shift @ARGV;
my $gff3out=shift @ARGV;
unless (ExonerateGff3Reformat($cdnain, $cdsin, $gff3out)) {
	die "Error:\n";
}
