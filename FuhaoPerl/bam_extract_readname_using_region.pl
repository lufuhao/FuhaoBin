#!/usr/bin/env perl
use strict;
use warnings;
use FuhaoPerl5Lib::BamKit qw/BamExtractReadsUsingBed/;
use constant USAGE =><<EOH;

usage: $0 input.bam input.bed readname.out

v20161103

EOH
die USAGE if (scalar(@ARGV) !=3 or $ARGV[0] eq '-h' or $ARGV[0] eq '--help');



die "Error: invalid input.bam\n" unless (defined $ARGV[0] and -s $ARGV[0]);
my $bamin=$ARGV[0];
die "Error: invalid input.bed\n" unless (defined $ARGV[1] and -s $ARGV[1]);
my $bedin=$ARGV[1];
die "Error: invalid readname.out\n" unless (defined $ARGV[2] and $ARGV[2]=~/^\S+$/);
my $readout=$ARGV[2];
unlink $readout if (-e $readout);


unless (BamExtractReadsUsingBed($bamin, $bedin, $readout)) {
	die "Error: failed to extract read names\n";
}
