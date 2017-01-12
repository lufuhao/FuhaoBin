#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Std;
use constant USAGE =><<EOH;

BLASTparser.pl -r <query_genome> -t <hit_genome> -i <blout_file> -o <output_file>

	query_genome	name or ID of the query genome.
	hit_genome		name or ID of the hit genome.
	blout_file		name or path to BLAST output file, where query_genome is blasted against hit_genome
	output_file		name of the output file, for e.g., synteny_file.txt
	
Source http://cas-bioinfo.cas.unt.edu/mgsv/tutorial.php
	
EOH
## Usage: BLASTparser.pl -r <query_genome> -t <hit_genome> -i <blout_file> -o <output_file> 

my %options=();
getopts("r:t:i:o:", \%options);
if (!defined $options{r} || !defined $options{t} || !defined $options{i} || !defined $options{o}) {
	error("Usage: BLASTparser.pl -n <name> -i <gff3_file> -c <configuration_file> -o <output_file>");
}


## Open the blout file
open(FILE, $options{i}) or error("Error: blout_file does not exist");
## Output file
open(OUT, ">", $options{o});
print OUT "#Ref\tref_start\tref_end\tTile\ttile_start\ttile_end\te_value\n";
## For each line in blout file
foreach my $line (<FILE>){
	## terminate at the FASTA information
	next if ($line =~ /^#/);
	my @tabs = split(/\t/, $line);
	error("ERROR: BLOUT file is not 12 columns separated by tabs") if (scalar @tabs != 12);
	print OUT $options{r}."\t$tabs[6]\t$tabs[7]\t".$options{t}."\t$tabs[8]\t$tabs[9]\t$tabs[10]\n";
}
close FILE;
close OUT;
exit;

########################################################################

## Print out error message and terminate the command
sub error{
	my $msg = shift;
	print $msg,"\n";
	exit;	
}

