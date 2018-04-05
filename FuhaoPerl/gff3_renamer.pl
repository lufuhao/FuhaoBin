#!/usr/bin/env perl
use strict;
use warnings;
use FuhaoPerl5Lib::GffKit qw/Gff3Renamer/;
use constant USAGE =><<EOH;


usage: $0 in.gff3 in.idlist out.gff3

v20180406

EOH
die USAGE if (scalar(@ARGV) !=3 or $ARGV[0] eq '-h' or $ARGV[0] eq '--help');


unless (Gff3Renamer($ARGV[0], $ARGV[1], $ARGV[2])) {
	die "Error: Gff3Renamer running failed\n";
}
