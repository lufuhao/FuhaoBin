#!/usr/bin/env perl
use strict;
use warnings;
use FuhaoPerl5Lib::BamKit qw/BamFilterReadsByNames/;
use constant USAGE =><<EOH;

usage: $0 input.bam readname.list code output.bam

code
  1=keep the reads in list
  0=Remove the reads in list

v20170123

EOH
die USAGE if (scalar(@ARGV) !=4 or $ARGV[0] eq '-h' or $ARGV[0] eq '--help');



die "Error: invalid input.bam\n" unless (defined $ARGV[0] and -s $ARGV[0]);
my $bamin=$ARGV[0];
die "Error: invalid input.bed\n" unless (defined $ARGV[1] and -s $ARGV[1]);
my $readin=$ARGV[1];
die "Error: invalid code\n" unless (defined $ARGV[2] and $ARGV[2]=~/^[01]{1}$/);
my $code=$ARGV[2];
die "Error: invalid output.bam\n" unless (defined $ARGV[3] and $ARGV[3]=~/^\S+$/);
my $bamout=$ARGV[3];
unlink $bamout if (-e $bamout);


unless (BamFilterReadsByNames($bamin, $readin, $code, $bamout)) {
	die "Error: failed to filter BAM\n";
}
