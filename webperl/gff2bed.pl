#!/usr/bin/env perl
use warnings;
use strict;
use Bio::Tools::GFF;
use Getopt::Long;
my ($bedfile);
GetOptions(
	'i|input=s' => \$bedfile,
);

my $gffio;
if ($bedfile eq 'stdin') {
	$gffio = Bio::Tools::GFF->new(-fh => \*STDIN, -gff_version => 3);
}
else {
	$gffio = Bio::Tools::GFF->new(-file => $bedfile, -gff_version => 3);
}

while (my $feature = $gffio->next_feature()) {
	my $seq_id = $feature->seq_id();
	my $start = $feature->start() - 1;	
	my $end = $feature->end();
	my $strand = $feature->strand();
	my $id = $feature->primary_id;
	
	next unless ($id =~ /Hr\.Aug\-/);
	next if ($id =~ /TU01/);
	
	if ( $strand == 1 ) {
		$strand="+";
	}
	elsif ( $strand == -1 ) {
		$strand="-";
	}
	print "$seq_id\t$start\t$end\t$id\t0\t$strand\n";
}
