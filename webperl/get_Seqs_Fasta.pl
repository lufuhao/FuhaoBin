#!/usr/bin/perl

#######################################################
# Author: Hamid Ashrafi                              
# email: ashrafi@ucdavis.edu
# Pupose: It reads the fasta file and another file with the IDs that you want to extract from the FASTA file
# then print the IDs and seqencing that are matched in the FASTA file.
#
# Requirement. Bio Perl
#
# Modified by Aurelie K
#   -> put the names of files in arguments (and also use the STDOUT for the output file)
#   -> suppr index file generated during the script
######################################################

#usage = path_to/fasta_getseq.pl path_to/fasta_library path_to/IDFILE > log

if ($#ARGV != 1) {
 print "usage:  fasta_getseq.pl IDFILE fasta_library > outputfile\n";
 exit;
}

use strict;
use Bio::DB::Fasta;

my $database;
my $fasta_library = $ARGV[1];   #path to fastalibrary in the second argument
my %records;

open IDFILE, "<$ARGV[0]" or die $!; #first argument is the path of the file containing all the IDs you need to extract
open OUTPUT, <STDOUT>;

# creates the database of the library, based on the file
$database = Bio::DB::Fasta->new("$fasta_library") or die "Failed to creat Fasta DP object on fasta library\n";

# now, it parses the file with the fasta headers you want to get
while (<IDFILE>) {

  my ($id) = (/^>*(\S+)/);  # capture the id string (without the initial ">")
  my $header = $database->header($id);
  #print "$header\n";
  print ">$header\n", $database->seq( $id ), "\n";
  print OUTPUT ">$header\n", $database->seq( $id ), "\n";

  print STDERR ">$header\n",$database->seq( $id), "\n";
}

#remove the index file that is useless for user
unlink "$fasta_library.index";

#close the filehandles
close IDFILE;
close OUTPUT;
exit;
