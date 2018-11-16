#!/usr/bin/env perl
use warnings;
use strict;
use Bio::Tools::GFF;
use constant USAGE =><<EOH;

Require: Bio::Tools::GFF in BioPerl

usage: $0 file.gff3 > file.gff2

EOH
die USAGE if (scalar(@ARGV)<1 or $ARGV[0] eq '-h' or $ARGV[0] eq '--help' or ! -s $ARGV[0]);
my( $gff3File ) = @ARGV;
my $gffio = Bio::Tools::GFF->new(-file=>"$gff3File", -gff_version=>3);
while( my $feature = $gffio->next_feature() ) {
    my $gff2string = $gffio->_gff2_string( $feature );
    print "$gff2string\n";
}
$gffio->close();
