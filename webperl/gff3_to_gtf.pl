#!/usr/bin/env perl
use strict;
use warnings;
use Data::Dumper;
use Bio::FeatureIO;
use constant USAGE =><<EOH;

usage: $0 infile.gff output.gtf

Desc:
	convert GFF3 to GTF/GFF2.5 format using Bio::FeatureIO;

v20160715

EOH
die USAGE if ($ARGV[0] eq '-h' or $ARGV[0] eq '--help'); 


my $inFile = shift;
my $outFile = shift;
 
my $inGFF = Bio::FeatureIO->new( '-file' => "$inFile",
 '-format' => 'GFF',
 '-version' => 3 );
my $outGTF = Bio::FeatureIO->new( '-file' => ">$outFile",
 '-format' => 'GFF',
 '-version' => 2.5);
 
while (my $feature = $inGFF->next_feature() ) {
 	$outGTF->write_feature($feature); 
}
