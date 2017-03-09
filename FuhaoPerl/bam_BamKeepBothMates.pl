#!/usr/bin/env perl
use strict;
use warnings;
use FuhaoPerl5Lib::BamKit qw/BamKeepBothMates/;
use constant USAGE =><<EOH;

usage: $0 in.R1.bam in.R2.bam paied.R1.bam paired.R2.bam [unpaied.R1.bam paired.R2.bam]

v20170217

EOH
die USAGE if ($ARGV[0] eq '-h' or $ARGV[0] eq '--help');
die USAGE unless (scalar(@ARGV) ==4 or scalar(@ARGV) ==6);


die "Error: invalid in.R1.bam\n" unless (defined $ARGV[0] and -s $ARGV[0]);
my $bamin1=$ARGV[0];
die "Error: invalid in.R2.bam\n" unless (defined $ARGV[1] and -s $ARGV[1]);
my $bamin2=$ARGV[1];
die "Error: invalid out.R1.bam\n" unless (defined $ARGV[2] and $ARGV[2]=~/^\S+$/);
my $bamout1=$ARGV[2];
unlink ($bamout1) if (-e $bamout1);
die "Error: invalid out.R2.bam\n" unless (defined $ARGV[3] and $ARGV[3]=~/^\S+$/);
my $bamout2=$ARGV[3];
unlink ($bamout2) if (-e $bamout2);


unless (BamKeepBothMates($bamin1, $bamin2, $bamout1, $bamout2, $ARGV[4], $ARGV[5])) {
	print STDERR "Error: Failed\n";
	exit 1;
}
else {
	print "Info: Success\n";
	exit 0;
}
