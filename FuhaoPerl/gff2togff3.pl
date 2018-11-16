#!/usr/bin/env perl
use warnings;
use strict;
use Bio::Tools::GFF;
use constant USAGE =><<EOH;

Require: Bio::Tools::GFF in BioPerl

usage: $0 file.gff3 > file.gff2

EOH
my( $gff2File ) = @ARGV;
my $gffio = Bio::Tools::GFF->new(-file=>"$gff2File", -gff_version=>2);
while( my $feature = $gffio->next_feature() ) {
     my $gff3string = $gffio->_gff3_string( $feature );
     print "$gff3string\n";
}
$gffio->close();
