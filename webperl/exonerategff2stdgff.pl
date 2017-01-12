#!/usr/bin/env perl
 
=head1 NAME

    convert_exonerate_gff_to_std_gff.pl

=head1 SYNOPSIS
 
    convert_exonerate_gff_to_std_gff.pl input_gff output_gff 
        where input_gff is the input exonerate gff file,
              output_gff is the name of the output gff file.

=head1 DESCRIPTION

    This script takes an input gff file from exonerate (<input_gff>) and converts it to
    a more standard gff file (<output_gff>).

=head1 VERSION
  
    Perl script last edited 12-Aug-2013.

=head1 CONTACT

    alc@sanger.ac.uk (Avril Coghlan)

=cut
 
# 
# Perl script convert_exonerate_gff_to_std_gff.pl
# Written by Avril Coghlan (alc@sanger.ac.uk)
# 12-Aug-13.
# Last edited 12-Aug-2013.
# SCRIPT SYNOPSIS: convert_exonerate_gff_to_std_gff.pl: converts a gff file from exonerate to more standard gff format.
#
#------------------------------------------------------------------#
 
# CHECK IF THERE ARE THE CORRECT NUMBER OF COMMAND-LINE ARGUMENTS:
 
use strict;
use warnings;
 
# xxx
# BEGIN {
#     unshift (@INC, '/nfs/users/nfs_a/alc/Documents/git/helminth_scripts/modules');
# }
#use HelminthGenomeAnalysis::AvrilGffUtils;
 
my $num_args               = $#ARGV + 1;
if ($num_args != 2) {
    print "Usage of convert_exonerate_gff_to_std_gff.pl\n\n";
    print "perl convert_exonerate_gff_to_std_gff.pl <input_gff> <output_gff>\n";
    print "where <input_gff> is the input exonerate gff file,\n";
    print "      <output_gff> is the output gff file\n";
    print "For example, >perl convert_exonerate_gff_to_std_gff.pl PTRK_exonerate1 PTRK_exonerate1.gff\n";
    exit;
}
 
# FIND THE PATH TO THE INPUT EXONERATE GFF FILE:
my $input_gff              = $ARGV[0];
 
# FIND THE NAME OF THE OUTPUT GFF FILE:
my $output_gff             = $ARGV[1];
 
#------------------------------------------------------------------#
 
# RUN THE MAIN PART OF THE CODE:
 
&run_main_program($input_gff,$output_gff);
print STDERR "FINISHED.\n";
 
#------------------------------------------------------------------#
 
# RUN THE MAIN PART OF THE CODE:
 
sub run_main_program {
	my $input_gff           = $_[0]; # INPUT GFF FILE
	my $output_gff          = $_[1]; # OUTPUT GFF FILE
	my $returnvalue;                 # RETURN VALUE FROM A FUNCTION 
	# CONVERT THE $input_gff FILE TO MORE STANDARD GFF FORMAT:
	#   $returnvalue = HelminthGenomeAnalysis::AvrilGffUtils::convert_exonerate_gff_to_standard_gff($input_gff,$output_gff);
	$returnvalue = &convert_exonerate_gff_to_standard_gff($input_gff,$output_gff);
}
 
#------------------------------------------------------------------#



# SUBROUTINE SYNOPSIS: convert_exonerate_gff_to_standard_gff: SUBROUTINE TO CONVERT EXONERATE GFF OUTPUT TO A MORE STANDARD GFF FORMAT.
sub convert_exonerate_gff_to_standard_gff {
   my $input_gff           = $_[0]; # INPUT GFF FILE
   my $output_gff          = $_[1]; # OUTPUT GFF FILE
   my $line;                        # 
   my @temp;                        # 
   my $feature_type;                # FEATURE TYPE
   my $scaffold;                    # SCAFFOLD START
   my $feature_start;               # FEATURE START
   my $feature_end;                 # FEATURE END
   my $strand;                      # FEATURE STRAND 
   my $gene_num            = 0;     # GENE NUMBER IN THE EXONERATE FILE
   my $gene;                        # GENE NAME 
   my $exon_num            = 0;     # EXON NUMBER IN A GENE
   # OPEN THE OUTPUT GFF:
   open(OUTPUT_GFF,">$output_gff") || die "ERROR: convert_exonerate_gff_to_standard_gff: cannot open $output_gff\n";
   # READ IN THE EXONERATE GFF:
   open(INPUT_GFF,"< $input_gff") || die "ERROR: convert_exonerate_gff_to_standard_gff: cannot open $input_gff\n";
   while($line=<INPUT_GFF>) {
      chomp $line;
      next if ($line=~/^#/);
      next if ($line=~/^Command line: \[exonerate/);
      next if ($line=~/^Hostname: \[/);
      next if ($line=~/^vulgar: /);
      next if ($line=~/^>/);
      next if ($line=~/^--/);
      @temp                = split(/\t+/,$line);
      # THROW AN EXCEPTION IF THERE ARE NOT 9 COLUMNS IN THE GFF LINE:
      print STDERR "Warnings: ERRORCODE=1: convert_exonerate_gff_to_standard_gff: do not have 9 columns in $input_gff: line $line\n" if ($#temp != 8); # TESTED FOR
      $scaffold            = $temp[0];
      $feature_type        = $temp[2];
      $feature_start       = $temp[3];
      $feature_end         = $temp[4];
      $strand              = $temp[6];
      if    ($feature_type eq 'gene') {
         # eg. scaffold_000039 exonerate:protein2genome:local  gene    252337  253077  685     +       .       gene_id 0 ; sequence SRAE_2000311000.t1:mRNA ; gene_orientation .
         $gene_num++;
         $exon_num         = 0;
         $gene             = "gene".$gene_num;
         print OUTPUT_GFF "$temp[0]\t$temp[1]\t$temp[2]\t$temp[3]\t$temp[4]\t$temp[5]\t$temp[6]\t$temp[7]\tID=$gene;Parent=$gene\n";
      } 
      elsif ($feature_type eq 'cds') {
         # scaffold_000039 exonerate:protein2genome:local  cds     252337  253077  .       +       .       cds
         print OUTPUT_GFF "$temp[0]\t$temp[1]\tCDS\t$temp[3]\t$temp[4]\t$temp[5]\t$temp[6]\t$temp[7]\tID=$gene:cds;Parent=$gene\n";     
      }
      elsif ($feature_type eq 'exon') {
         # scaffold_000039 exonerate:protein2genome:local  exon    252337  253077  .       +       .       insertions 0 ; deletions 2
         $exon_num++;
         print OUTPUT_GFF "$temp[0]\t$temp[1]\t$temp[2]\t$temp[3]\t$temp[4]\t$temp[5]\t$temp[6]\t$temp[7]\tID=$gene:exon:$exon_num;Parent=$gene\n";
      } 
      elsif ($feature_type eq 'similarity') {
         # scaffold_000039 exonerate:protein2genome:local  similarity      252337  253077  685     +       .       alignment_id 0 ; Query SRAE_2000311000.t1:mRNA ; Align 25175 81 564 ; Align 25739 271 177
         print OUTPUT_GFF "$line\n";
      }
      elsif ($feature_type eq 'splice5') {
         # scaffold_000043 exonerate:protein2genome:local  splice5 45991   45992   .       +       .       intron_id 1 ; splice_site "GT"
         print OUTPUT_GFF "$temp[0]\t$temp[1]\t$temp[2]\t$temp[3]\t$temp[4]\t$temp[5]\t$temp[6]\t$temp[7]\tID=$gene:splice5;Parent=$gene\n";
      }
      elsif ($feature_type eq 'splice3') {
         # scaffold_000043 exonerate:protein2genome:local  splice3 46138   46139   .       +       .       intron_id 0 ; splice_site "AT"
         print OUTPUT_GFF "$temp[0]\t$temp[1]\t$temp[2]\t$temp[3]\t$temp[4]\t$temp[5]\t$temp[6]\t$temp[7]\tID=$gene:splice3;Parent=$gene\n";
      }
      elsif ($feature_type eq 'intron') {
         # scaffold_000043 exonerate:protein2genome:local  intron  45991   46139   .       +       .       intron_id 1
         print OUTPUT_GFF "$temp[0]\t$temp[1]\t$temp[2]\t$temp[3]\t$temp[4]\t$temp[5]\t$temp[6]\t$temp[7]\tID=$gene:intron;Parent=$gene\n";
      }
      else {
         print STDERR "Warnings: ERRORCODE=1: convert_exonerate_gff_to_standard_gff: unknown feature type $feature_type"; # TESTED FOR 
      }
   } 
   close(INPUT_GFF);
   close(OUTPUT_GFF);
 
   return(1);
}
