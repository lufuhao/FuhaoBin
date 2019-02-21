#!/usr/bin/env perl
use strict;
use warnings;
use constant USAGE =><<EOH;
use Bio::EnsEMBL::Compara::DBSQL::DBAdaptor;

Requirements:
    Bio::EnsEMBL::Compara::DBSQL::DBAdaptor

usage: $0 <species>

    where <species> is the species name (eg. schistosoma_mansoni)

v20190129

Examples

$0 schistosoma_mansoni

EOH
die USAGE if (scalar(@ARGV) !=1 or $ARGV[0] eq '-h' or $ARGV[0] eq '--help');


# FIND THE NAME OF THE SPECIES OF INTEREST:                     
my $species = $ARGV[0];

#------------------------------------------------------------------#
my $comparadb = new Bio::EnsEMBL::Compara::DBSQL::DBAdaptor (-host => 'xxx', -port => 3378, -user => 'xxx', -pass => 'xxx', -dbname => 'xxx');

# first get the taxon_id for species $species (eg. schistosoma_mansoni)
my $species_taxon_id = 'none';
my $gda = $comparadb->get_adaptor("GenomeDB");
my @genomes = @{$gda->fetch_all()};
foreach my $genome (@genomes){
	my $genome_name = $genome->name(); # eg. homo_sapiens 
	my $taxon = $genome->taxon()->name(); # eg. Homo sapiens
	my $taxon_id = $genome->taxon_id();
	if ($genome_name eq $species){ 
		if ($species_taxon_id ne 'none') {
			print STDERR "ERROR: already have species_taxon_id $species_taxon_id\n"; exit;
		}
		$species_taxon_id = $taxon_id;
	}
}
if ($species_taxon_id eq 'none') {
	print STDERR "ERROR: did not find species_taxon_id for $species\n";
	exit;
}

# get the Bio::EnsEMBL::Compara::SeqMember objects for the species $species 
my $pma = $comparadb->get_adaptor("SeqMember");
my @proteins = @{$pma->fetch_all_by_source_taxon("Ensemblpep", $species_taxon_id)};

foreach my $protein (@proteins) {
	my $stable_id = $protein->stable_id();
	if (defined($protein->sequence())) { 
		my $sequence = $protein->sequence();
		print ">$stable_id\n";
		print "$sequence\n";
	}
	else {
		print STDERR "ERROR: do not know sequence for $stable_id\n";
		exit;
	}
}
