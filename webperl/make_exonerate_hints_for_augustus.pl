#!/usr/bin/env perl
 
=head1 NAME

    make_exonerate_hints_for_augustus.pl

=head1 SYNOPSIS
 
    make_exonerate_hints_for_augustus.pl input_gff hints_file type
        where input_gff is the input exonerate gff file  
              hints_file is the output hints file,
              type is the exonerate alignment type (prot/est)

=head1 DESCRIPTION

    This script takes an exonerate gff file (<input_gff>), and generates a hints
    file for augustus. The exonerate alignments can be for protein (type="prot") or
    EST (type="est"). 

=head1 VERSION
  
    Perl script last edited 13-Aug-2013.

=head1 CONTACT

    alc@sanger.ac.uk (Avril Coghlan)

=cut
 
# 
# Perl script make_exonerate_hints_for_augustus.pl
# Written by Avril Coghlan (alc@sanger.ac.uk)
# 13-Aug-13.
# Last edited 13-Aug-2013.
# SCRIPT SYNOPSIS: make_exonerate_hints_for_augustus.pl: given an exonerate gff file, generates a hints file for augustus.
#
#------------------------------------------------------------------#
 
# CHECK IF THERE ARE THE CORRECT NUMBER OF COMMAND-LINE ARGUMENTS:
 
use strict;
use warnings;
 
# xxx
# BEGIN {
#     unshift (@INC, '/nfs/users/nfs_a/alc/Documents/git/helminth_scripts/modules');
# }
 
use HelminthGenomeAnalysis::AvrilGenefindingUtils; 
 
my $num_args               = $#ARGV + 1;
if ($num_args != 3) {
    print "Usage of make_exonerate_hints_for_augustus.pl\n\n";
    print "perl make_exonerate_hints_for_augustus.pl <input_gff> <hints_file> <type>\n";
    print "where <input_gff> is the input exonerate gff file,\n";
    print "      <hints_file> is the output hints file,\n";
    print "      <type> is the exonerate alignment type (prot/est)\n";
    print "For example, >perl make_exonerate_hints_for_augustus.pl PTRK_exonerate1.gff PTRK_exonerate_hints.gff prot\n";
    exit;
}
 
# FIND THE PATH TO THE INPUT EXONERATE GFF FILE:
my $input_gff              = $ARGV[0];
 
# FIND THE NAME OF THE OUTPUT HINTS FILE: 
my $hints_file             = $ARGV[1];
 
# FIND WHETHER THE EXONERATE ALIGNMENTS WERE FOR PROTEIN OR EST:
my $type                   = $ARGV[2];
if ($type ne 'prot' && $type ne 'est') { print STDERR "ERROR: type should be prot for protein, or est for ESTs\n"; exit;}
#------------------------------------------------------------------#

# RUN THE MAIN PART OF THE CODE:
 
&run_main_program($input_gff,$hints_file,$type);
 
print STDERR "FINISHED.\n";
 
#------------------------------------------------------------------#
 
# RUN THE MAIN PART OF THE CODE:
 
sub run_main_program {
   my $input_gff           = $_[0]; # INPUT GFF FILE
   my $hints_file          = $_[1]; # OUTPUT GFF FILE
   my $type                = $_[2]; # TYPE OF EXONERATE ALIGNMENT (est/prot)
   my $returnvalue;                 # RETURN VALUE FROM A FUNCTION  
   
   # ASSUME INTRONS WILL BE 15-350000 bp (15 USED IN bam2hints, see 
   # http://avrilomics.blogspot.co.uk/2013/05/running-augustus-gene-finding-software.html AND
   # 350,000 USED IN /nfs/users/nfs_a/alc/Documents/bin/augustus/scripts/exonerate2hints.pl
   $returnvalue = HelminthGenomeAnalysis::AvrilGenefindingUtils::make_exonerate_hints($input_gff,$hints_file,$type,15,350000,9),
 
}
 
#------------------------------------------------------------------#
