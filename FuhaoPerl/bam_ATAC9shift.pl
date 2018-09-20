#!/usr/bin/env perl
use strict;
use warnings;
use FuhaoPerl5Lib::BamKit qw/BamAtacShift/;
use constant USAGE =><<EOH;

usage: $0 in.bam out.bam

Description: 
    Modify BAM file to adjust 9 bp for ATAC-seq:
        forward strand start +4 and
        reverse strans start -5

Note: got problems when R1 and R2 overlaps
            <---------------------
                 --------------------->
      So used BEDPE instead

Requirements:
    samtools

v20180910

EOH
die USAGE if (scalar(@ARGV) !=2 or $ARGV[0] eq '-h' or $ARGV[0] eq '--help');



unless(BamAtacShift($ARGV[0], $ARGV[1])) {
	die "Error: BamAtacShift running failed\n";
}
