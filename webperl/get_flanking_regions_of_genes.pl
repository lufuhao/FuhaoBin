#!/usr/bin/env perl
#https://gist.github.com/avrilcoghlan/5917580

=head1 NAME

    get_flanking_regions_of_genes.pl

=head1 SYNOPSIS
 
    get_flanking_regions_of_genes.pl input_gff input_fasta output_gff output_fasta outputdir flank_size
        where input_gff is the input gff file of gene predictions,
              input_fasta is the input assembly fasta file,
              output_gff is the output gff file of gene predictions,
              output_fasta is the output fasta file of genes with flanking regions,
              outputdir is the output directory for writing output files,
              flank_size is the length of flanking region to take.

=head1 DESCRIPTION

    This script takes an input gff file of gene predictions (<input_gff>) for an
    assembly (<input_fasta>), and then makes an output fasta file with flank_size bp
    on either side of each gene (<output_fasta>) and an output gff file with the
    gene prediction with respect to the output fasta (<output_gff>).

=head1 VERSION
  
    Perl script last edited 27-Jun-2013.

=head1 CONTACT

    alc@sanger.ac.uk (Avril Coghlan)

=cut
 
# 
# Perl script get_flanking_regions_of_genes.pl
# Written by Avril Coghlan (alc@sanger.ac.uk)
# 27-Jun-13.
# Last edited 27-Jun-2013.
# SCRIPT SYNOPSIS: get_flanking_regions_of_genes.pl: take an input gff, and take flank_size bp on either side of each gene, and write an output fasta and gff wrt the output fasta 
#
#------------------------------------------------------------------#
 
# CHECK IF THERE ARE THE CORRECT NUMBER OF COMMAND-LINE ARGUMENTS:
 
use strict;
use warnings;
 
# xxx
# BEGIN { 
#     unshift (@INC, '/nfs/users/nfs_a/alc/Documents/git/helminth_scripts/modules'); 
# }
 
use HelminthGenomeAnalysis::AvrilGffUtils;
use HelminthGenomeAnalysis::AvrilFileUtils;
use HelminthGenomeAnalysis::AvrilFastaUtils;
 
my $num_args               = $#ARGV + 1;
if ($num_args != 6)
{
    print "Usage of get_flanking_regions_of_genes.pl\n\n";
    print "perl get_flanking_regions_of_genes.pl <input_gff> <input_fasta> <output_gff> <output_fasta> <outputdir> <flank_size>\n";
    print "where <input_gff> is the input gff file of gene predictions,\n";
    print "      <input_fasta> is the input assembly fasta file,\n";
    print "      <output_gff> is the output gff file of gene predictions,\n";
    print "      <output_fasta> is hte output fasta file of genes with flanking regions,\n";
    print "      <outputdir> is the output directory for writing output files,\n";
    print "      <flank_size> is the length of flanking region to take\n";
    print "For example, >perl get_flanking_regions_of_genes.pl SSTP_training.gff SSTP.v1.fa SSTP_training2.gff SSTP_training2.fa . 1000\n";
    exit;
}
 
# FIND THE PATH TO THE INPUT GFF FILE:                     
 
my $input_gff              = $ARGV[0];
 
# FIND THE PATH TO THE INPUT FASTA FILE:
 
my $input_fasta            = $ARGV[1];
 
# FIND THE PATH TO THE OUTPUT GFF FILE:
 
my $output_gff             = $ARGV[2];
 
# FIND THE PATH TO THE OUTPUT FASTA FILE:
 
my $output_fasta           = $ARGV[3];
 
# FIND THE DIRECTORY TO USE FOR OUTPUT FILES:      
 
my $outputdir              = $ARGV[4];
 
# FIND THE SIZE OF FLANKING REGION TO USE:
 
my $flank_size             = $ARGV[5];
 
#------------------------------------------------------------------#
 
# TEST SUBROUTINES: 
 
my $PRINT_TEST_DATA        = 0;   # SAYS WHETHER TO PRINT DATA USED DURING TESTING.
&test_run_main_program($outputdir);
print STDERR "Finished tests, running main code now...\n";
 
#------------------------------------------------------------------#
 
# RUN THE MAIN PART OF THE CODE:
 
&run_main_program($outputdir,$input_gff,$input_fasta,$output_gff,$output_fasta,$flank_size);
 
print STDERR "FINISHED.\n";
 
#------------------------------------------------------------------#
 
# TEST &run_main_program
 
sub test_run_main_program
{
   my $outputdir           = $_[0]; # DIRECTORY TO WRITE OUTPUT FILES IN
   my $input_gff;                   # INPUT GFF FILE
   my $input_fasta;                 # INPUT FASTA FILE
   my $output_gff;                  # OUTPUT GFF FILE
   my $output_fasta;                # OUTPUT FASTA FILE
   my $expected_output_gff;         # FILE WITH EXPECTED CONTENTS OF $output_gff
   my $expected_output_fasta;       # FILE WITH EXPECTED CONTENTS OF $output_fasta 
   my $differences;                 # DIFFERENCES BETWEEN $output_gff AND $expected_output_gff
   my $length_differences;          # LENGTH OF $differences 
   my $line;                        #  
    
   $input_gff              = HelminthGenomeAnalysis::AvrilFileUtils::make_filename($outputdir); 
   open(INPUT_GFF,">$input_gff") || die "ERROR: test_run_main_program: cannot open $input_gff\n";
   print INPUT_GFF "scaff1  source	gene	20	50	.	+	.	feature1\n";
   print INPUT_GFF "scaff1	source	gene	1	1	.	+	.	feature2\n";
   print INPUT_GFF "scaff1	source	gene	1	2	.	+	.	feature3\n";
   print INPUT_GFF "scaff2	source	gene	5	10	.	+	.	feature4\n";
   close(INPUT_GFF);
   $input_fasta            = HelminthGenomeAnalysis::AvrilFileUtils::make_filename($outputdir); 
   open(INPUT_FASTA,">$input_fasta") || die "ERROR: test_run_main_program: cannot open $input_fasta\n";
   print INPUT_FASTA ">scaff1\n";
   print INPUT_FASTA "ABCDEFGHIJKLMNOPQRSTUVWXYZ\n";
   print INPUT_FASTA "abcdefghijklmnopqrstuvwxyz\n";
   print INPUT_FASTA ">scaff2\n";
   print INPUT_FASTA "ABCDEFGHIJKLMNOPQRSTUVWXYZ\n";
   print INPUT_FASTA "abcdefghijklmnopqrstuvwxyz\n";
   close(INPUT_FASTA); 
   $output_gff             = HelminthGenomeAnalysis::AvrilFileUtils::make_filename($outputdir); 
   $output_fasta           = HelminthGenomeAnalysis::AvrilFileUtils::make_filename($outputdir); 
   &run_main_program($outputdir,$input_gff,$input_fasta,$output_gff,$output_fasta,5);
   # CHECK THAT $output_gff CONTAINS THE CONTENTS WE EXPECT:
   # AFTER TAKING FLANKING REGION OF 5 bp, SHOULD GET:
   #  20--50 ---> 15--52 # 50+5 IS 55, BUT LENGTH OF SCAFFOLD IS 52
   #  1--1   ---> 1--6
   #  1--2   ---> 1--7
   #  5--10  ---> 1--15
   # POSITIONS OF THE FEATURES WITH RESPECT TO THESE SUBSEQUENCES:
   #  20--50 --> 6--36 # 20-15+1 = 6, 6+(50-20) = 36                                     
   #  1--1   --> 1--1  # 1-1+1 = 1, 1+(1-1) = 1 
   #  1--2   --> 1--2  # 1-1+1 = 1, 1+(2-1) = 1+1=2 
   #  5--10  --> 5--10 # 5-1+1 = 5, 5+(10-5) = 5+5=10
   $expected_output_gff    = HelminthGenomeAnalysis::AvrilFileUtils::make_filename($outputdir); 
   open(EXPECTED_OUTPUT_GFF,">$expected_output_gff") || die "ERROR: test_run_main_program: cannot open $expected_output_gff\n";
   print EXPECTED_OUTPUT_GFF "scaff1:15-52	source	gene	6	36	.	+	.	feature1\n";
   print EXPECTED_OUTPUT_GFF "scaff1:1-6	source	gene	1	1	.	+	.	feature2\n";
   print EXPECTED_OUTPUT_GFF "scaff1:1-7	source	gene	1	2	.	+	.	feature3\n";
   print EXPECTED_OUTPUT_GFF "scaff2:1-15	source	gene	5	10	.	+	.	feature4\n";
   close(EXPECTED_OUTPUT_GFF);
   $differences            = "";
   open(TEMP,"diff $output_gff $expected_output_gff |");
   while(<TEMP>)
   {
      $line                = $_;
      $differences         = $differences.$line;
   }
   close(TEMP);
   $length_differences     = length($differences);
   if ($length_differences != 0) { print STDERR "ERROR: test_run_main_program: failed test1 (output_gff $output_gff expected_output_gff $expected_output_gff)\n"; exit;}
   # CHECK THAT $output_fasta CONTAINS THE CONTENTS WE EXPECT:
   $expected_output_fasta  = HelminthGenomeAnalysis::AvrilFileUtils::make_filename($outputdir); 
   open(EXPECTED_OUTPUT_FASTA,">$expected_output_fasta");
   print EXPECTED_OUTPUT_FASTA ">scaff1:15-52\n";
   print EXPECTED_OUTPUT_FASTA "OPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz\n";
   print EXPECTED_OUTPUT_FASTA ">scaff1:1-6\n";
   print EXPECTED_OUTPUT_FASTA "ABCDEF\n";
   print EXPECTED_OUTPUT_FASTA ">scaff1:1-7\n";
   print EXPECTED_OUTPUT_FASTA "ABCDEFG\n";
   print EXPECTED_OUTPUT_FASTA ">scaff2:1-15\n";
   print EXPECTED_OUTPUT_FASTA "ABCDEFGHIJKLMNO\n";
   close(EXPECTED_OUTPUT_FASTA); 
   $differences            = "";
   open(TEMP,"diff $output_fasta $expected_output_fasta |");
   while(<TEMP>)
   {
      $line                = $_;
      $differences         = $differences.$line;
   }
   close(TEMP);
   $length_differences     = length($differences);
   if ($length_differences != 0) { print STDERR "ERROR: test_run_main_program: failed test1 (output_fasta $output_fasta expected_output_fasta $expected_output_fasta)\n"; exit;}
    
}
 
#------------------------------------------------------------------#
 
# RUN THE MAIN PART OF THE CODE:
 
sub run_main_program
{
   my $outputdir           = $_[0]; # DIRECTORY TO PUT OUTPUT FILES IN.
   my $input_gff           = $_[1]; # THE INPUT GFF FILE  
   my $input_fasta         = $_[2]; # THE INPUT FASTA FILE
   my $output_gff          = $_[3]; # THE OUTPUT GFF FILE
   my $output_fasta        = $_[4]; # THE OUTPUT FASTA FILE
   my $flank_size          = $_[5]; # THE SIZE OF FLANKING REGION TO USE
   my $gene_gff;                    # GFF FILE OF GENES
   my $returnvalue;                 # 
   my $feature_types;               # FEATURE TYPES THAT WE ARE INTERESTED IN 
   my $gene_gff2;                   # GFF FILE WITH $flank_size bp ON EITHER SIDE OF GENES 
   my $input_fasta_obj;             # OBJECT FOR THE INPUT FASTA FILE 
   my $contigs2len;                 # HASH TABLE OF THE LENGTHS OF SCAFFOLD/CONTIG SEQUENCES
   my $contigslen_file;             # FILE WITH THE LENGTHS OF SCAFFOLD/CONTIG SEQUENCES 
   my $contig;                      # A CONTIG/SCAFFOLD NAME 
   my $flanks_gff;                  # GFF WITH GENES, PLUS FLANKING REGIONS 
   my $DIFF_IN_START;               # HASH TABLE TO STORE DIFFERENCES IN THE START OF A FEATURE BETWEEN $gene_gff AND $gene_gff3
   my $gene_gff3;                   # GFF WITH POSITIONS OF FEATURES IN $gene_gff, WITH RESPECT TO REGIONS IN $gene_gff2
   my $TRANSCRIPT2GENE;             # HASH TABLE WITH GENE NAMES FOR TRANSCRIPTS
 
   # GET A GFF FILE OF GENES:
   $gene_gff               = HelminthGenomeAnalysis::AvrilFileUtils::make_filename($outputdir); # MAKE A TEMPORARY FILE IN THE CURRENT DIRECTORY
   $feature_types          = ['gene'];
   $returnvalue            = HelminthGenomeAnalysis::AvrilGffUtils::get_gff_lines_for_feature_types($input_gff,$feature_types,$gene_gff); 
 
   # GET THE LENGTHS OF THE SEQUENCES IN THE INPUT FASTA FILE, AND WRITE THEM OUT IN A TEMPORARY FILE:
   $input_fasta_obj        = HelminthGenomeAnalysis::AvrilFastaUtils->new(fasta_file => $input_fasta); 
   $contigs2len            = HelminthGenomeAnalysis::AvrilFastaUtils::build_contigs2len($input_fasta_obj);
   $contigslen_file        = HelminthGenomeAnalysis::AvrilFileUtils::make_filename($outputdir); # MAKE A TEMPORARY FILE IN THE CURRENT DIRECTORY
   open(CONTIGSLENFILE,">$contigslen_file") || die "ERROR: run_main_program: cannot open $contigslen_file\n";
   foreach my $contig (keys %{$contigs2len})
   {
      print CONTIGSLENFILE "$contig\t$contigs2len->{$contig}\n";
   } 
   close(CONTIGSLENFILE);
 
   # MAKE A NEW GFF FILE WITH $flank_size bp ON EITHER SIDE OF EACH GENE: 
   $flanks_gff             = HelminthGenomeAnalysis::AvrilFileUtils::make_filename($outputdir); 
   $returnvalue            = HelminthGenomeAnalysis::AvrilGffUtils::add_flanking_region_to_gff_features($gene_gff,$flank_size,$contigslen_file,$flanks_gff);
   # MAKE A FASTA FILE OF THE SEQUENCES IN $flanks_gff, WITH $flank_size_bp ON EITHER SIDE OF EACH GENE:
   $returnvalue            = HelminthGenomeAnalysis::AvrilGffUtils::make_fasta_for_gff($flanks_gff,$input_fasta,$output_fasta,$outputdir);
 
   # MAKE A GFF FILE WTIH THE COORDINATES OF THE 'gene' FEATURES OF $gene_gff, WITH RESPECT TO THE SUBSEQUENCES IN $flanks_gff:
   $gene_gff3              = HelminthGenomeAnalysis::AvrilFileUtils::make_filename($outputdir); # MAKE A TEMPORARY FILE IN THE CURRENT DIRECTORY
   $DIFF_IN_START          = HelminthGenomeAnalysis::AvrilGffUtils::make_gff_for_features_in_subsequences($gene_gff,$flanks_gff,$gene_gff3),
 
   # ADD THE OTHER FEATURES IN $input_gff TO $gene_gff3, TO MAKE $output_gff:
   $TRANSCRIPT2GENE        = HelminthGenomeAnalysis::AvrilGffUtils::read_genenames_for_transcripts($input_gff),
   $returnvalue            = HelminthGenomeAnalysis::AvrilGffUtils::add_gene_features_to_gff_for_subsequences($input_gff,$gene_gff3,$DIFF_IN_START,$output_gff,$TRANSCRIPT2GENE);
}
 
#------------------------------------------------------------------#
