#!/usr/bin/env perl
use strict;
use warnings;
use FuhaoPerl5Lib::BamKit qw/BamMarkPairs/;
use constant USAGE =><<EOH;

usage: $0 in.R1.bam in.R2.bam out.R1.bam out.R2.bam

Mark flags in Paired BAM
	*Need to make sure both mates are paired
	R1 will add 67
	R2 will add 131

v20161114

EOH
die USAGE if (scalar(@ARGV) !=4 or $ARGV[0] eq '-h' or $ARGV[0] eq '--help');

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


unless (BamMarkPairs($bamin1, $bamin2, $bamout1, $bamout2)) {
	print STDERR "Failed\n";
	exit 1;
}
else {
	print "Success\n";
	exit 0;
}
