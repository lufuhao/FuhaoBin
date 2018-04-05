#!/usr/bin/env perl
use strict;
use warnings;
use FuhaoPerl5Lib::GffKit qw/AnnotationTransfer/;
use constant USAGE =><<EOH;

usage: $0 GFF3input    configure_file    GFF3output

DESC: transfer GFF from scaffolds to pseudomolecule

* Configure_file [tab-delimited]:

New_scaf	PSstart	Prev	scaffold_ID	strand+/-	start	end
New_scaf1	1	scaffold1	+	1	19900
New_scaf1	20001	scaffold2	-	101	50000
     ### means scaffold2:1-100 not in pseudomolecule

OUTPUT: 
    GFF3output
    GFF3output[.excluded]: not transfered

v20180318

EOH
die USAGE if (scalar(@ARGV) !=3 or $ARGV[0] eq '-h' or $ARGV[0] eq '--help');


unless (AnnotationTransfer($ARGV[0], $ARGV[1], $ARGV[2])) {
	die "Error: AnnotationTransfer running error\n";
}
