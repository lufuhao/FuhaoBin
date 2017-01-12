#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Std;
use constant USAGE =><<EOH;

Usage: BLASTparser.pl -r <query_genome> -t <hit_genome> -i <blastz__output_file> -o <output_file>

	query_genome: name or ID of the query genome.
	hit_genome: name or ID of the hit genome.
	blastz__output_file: name or path to BLASTZ output file, where query_genome is blasted against hit_genome
	output_file: name of the output file, for e.g., synteny_file.txt

EOH
## 

my %options=();
getopts("r:t:i:o:", \%options);
die USAGE if (!defined $options{r} || !defined $options{t} || !defined $options{i} || !defined $options{o});

$/="}
a {";

## Open the BLASTZ output file
open(FILE, $options{i}) or die "Error: blastz__output does not exist";
## Output file
open(OUT, ">", $options{o});
print OUT "#Ref\tref_start\tred_end\tTile\ttile_start\ttile_end\tscore\n";
## For each line in blastz output file
<FILE>;
foreach my $line (<FILE>){
	my @set = split(/\n/, $line);
	## Extract relevent information
	my ($score) = $set[1] =~ /(\d+)/;
	my ($q_start, $h_start) = $set[2] =~ /(\d+)\ (\d+)/;
	my ($q_end, $h_end) = $set[3] =~ /(\d+)\ (\d+)/;
	print OUT $options{r}."\t$q_start\t$q_end\t".$options{t}."\t$h_start\t$h_end\t$score\n";
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

