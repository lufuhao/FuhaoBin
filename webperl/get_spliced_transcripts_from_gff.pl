#!/usr/bin/env perl
#https://gist.github.com/avrilcoghlan/4605556
 
=head1 NAME

    get_spliced_transcripts_from_gff.pl

=head1 SYNOPSIS
 
    get_spliced_transcripts_from_gff.pl input_gff input_fasta output_fasta outputdir ignore_phase ignore_semicolons from_augustus
        where input_gff is the input gff file,
              input_fasta is the input fasta file for scaffolds,
              output_fasta is the output fasta file for transcripts,
              outputdir is the output directory for writing output files,
              ignore_phase tells the program to assume each transcript starts in phase 0, and use the expected phase for each exon (yes/no),
              ignore_semicolons says whether to ignore semicolons at the end of sequences names in input_fasta,
              from_augustus says whether the gff is from augustus.

=head1 DESCRIPTION

    This script takes an input gff file (<input_gff>), and infers the spliced
    sequences of transcripts, and writes them out into an output fasta file
    (<output_fasta>). That is, it takes the coding exons' (CDS) DNA sequences for a gene,
    and splices them together into the spliced DNA sequences (the sequence that's translated
    by the cell into protein). For genes on the minus strand, it finds the reverse
    complement.
       This script takes into account the phase of the first exon - if the first exon does
    not have phase 0, it gets rid of some of the DNA at the start of the sequence.
       Note that the transcript sequences do not necessarily start with ATG and end in
    TAA/TGA/TAG.

=head1 VERSION
  
    Perl script last edited 10-Jan-2013.

=head1 CONTACT

    alc@sanger.ac.uk (Avril Coghlan)

=cut
 
# 
# Perl script get_spliced_transcripts_from_gff.pl
# Written by Avril Coghlan (alc@sanger.ac.uk)
# 10-Jan-13.
# Last edited 10-Jan-2013.
# SCRIPT SYNOPSIS: get_spliced_transcripts_from_gff.pl: infer spliced DNA sequences of transcripts, based on a gff
# 
#------------------------------------------------------------------#
 
# CHECK IF THERE ARE THE CORRECT NUMBER OF COMMAND-LINE ARGUMENTS:
 
use strict;
use warnings;
 
my $num_args               = $#ARGV + 1;
if ($num_args != 7)
{
    print "Usage of get_spliced_transcripts_from_gff.pl\n\n";
    print "perl get_spliced_transcripts_from_gff.pl <input_gff> <input_fasta> <output_fasta> <outputdir> <ignore_phase> <ignore_semicolons> <from_augustus>\n";
    print "where <input_gff> is the input gff file,\n";
    print "      <input_fasta> is the input fasta file for scaffolds,\n";
    print "      <output_fasta> is the output fasta file,\n";
    print "      <outputdir> is the output directory for writing output files,\n";
    print "      <ignore_phase> means tells the program to assume each transcript starts in phase 0, and use the expected phase for each exon (yes/no),\n";
    print "      <ignore_semicolons> says whether to ignore semicolons at the end of sequence names in <input_fasta> (yes/no),\n";
    print "      <from_augustus> says whether the gff is from augustus\n";
    print "For example, >perl get_spliced_transcripts_from_gff.pl Pk_strainH_scaffolds.gff\n";
    print "Pk_strainH_scaffolds.fa Pk_strainH_transcripts.fa\n";
    print "/lustre/scratch108/parasites/alc/RNA-SeQC/Pknowlesi no yes no\n";
    exit;
}
 
# FIND THE PATH TO THE INPUT GFF FILE:                     
 
my $input_gff              = $ARGV[0];
 
# FIND THE PATH TO THE INPUT FASTA FILE:
 
my $input_fasta            = $ARGV[1];
 
# FIND THE PATH TO THE OUTPUT FASTA FILE:
 
my $output_fasta           = $ARGV[2];
 
# FIND THE DIRECTORY TO USE FOR OUTPUT FILES:      
 
my $outputdir              = $ARGV[3];
 
# FIND OUT WHETHER TO IGNORE THE PHASE INFORMATION:
 
my $ignore_phase           = $ARGV[4]; 
 
# FIND OUT WHETHER TO IGNORE SEMICOLONS AT THE END OF SEQUENCES NAMES IN $input_fasta:
 
my $ignore_semicolons      = $ARGV[5];
 
# FIND OUT WHETHER THE GFF FILE IS FROM AUGUSTUS:
 
my $from_augustus          = $ARGV[6];
 
#------------------------------------------------------------------#
 
# TEST SUBROUTINES: 
 
my $PRINT_TEST_DATA        = 0;   # SAYS WHETHER TO PRINT DATA USED DURING TESTING.
&test_get_exon_sequences($outputdir);
&test_infer_spliced_dna($outputdir);
&test_calculate_expected_phase;
&test_correct_first_exon_seq;
&test_correct_end_prev_exon_seq;
&test_reverse_complement;
&test_read_assembly($outputdir);
&test_print_error;
&test_correct_start_exon_seq;
&test_read_exon_gff($outputdir);
&test_make_exon_gff($outputdir);
&test_run_main_program($outputdir);
print STDERR "Finished tests: now running main code...\n";
 
#------------------------------------------------------------------#
 
# RUN THE MAIN PART OF THE CODE:
 
&run_main_program($outputdir,$input_gff,$input_fasta,$output_fasta,$ignore_phase,$ignore_semicolons,0,$from_augustus);
 
print STDERR "FINISHED.\n";
 
#------------------------------------------------------------------#
 
# TEST &run_main_program
 
sub test_run_main_program
{
   my $outputdir           = $_[0]; # DIRECTORY TO WRITE OUTPUT FILES INTO
   my $input_gff;                   # INPUT GFF FILE
   my $input_fasta;                 # INPUT FASTA FILE
   my $output_fasta;                # OUTPUT FASTA FILE
   my $expected_output_fasta;       # FILE WITH EXPECTED CONTENTS OF $output_fasta
   my $errorcode;                   # RETURNED AS 0 FROM A FUNCTION IF THERE IS NO ERROR
   my $errormsg;                    # RETURNED AS 'none' FROM A FUNCTION IF THERE IS NO ERROR
   my $differences;                 # DIFFERENCES BETWEEN $expected_output_fasta AND $output_fasta
   my $length_differences;          # LENGTH OF $differences
   my $line;                        #  
   my @temp;                        #  
 
   ($input_gff,$errorcode,$errormsg) = &make_filename($outputdir);
   if ($errorcode != 0) { ($errorcode,$errormsg) = &print_error($errormsg,$errorcode,0); } 
   open(INPUT_GFF,">$input_gff") || die "ERROR: test_run_main_program: cannot open $input_gff\n";
   print INPUT_GFF "##gff-version 3\n";
   print INPUT_GFF "##sequence-region Pk_strainH_chr01 1 838594\n";
   print INPUT_GFF "Pk_strainH_chr01    chado	gene	108	1421	.	+	.	ID=PKH_010010;isObsolete=false;isFminPartial;feature_id=222341;timelastmodified=10.12.2012+02:43:44+GMT\n";
   print INPUT_GFF "Pk_strainH_chr01	chado	CDS	108	758	.	+	0	ID=PKH_010010.1:exon:1;Parent=PKH_010010.1;isFminPartial;isObsolete=false;timelastmodified=10.12.2012+02:43:44+GMT;colour=7\n";
   print INPUT_GFF "Pk_strainH_chr01	chado	CDS	978	1421	.	+	0	ID=PKH_010010.1:exon:2;Parent=PKH_010010.1;isFminPartial;isObsolete=false;timelastmodified=10.12.2012+02:43:44+GMT;colour=7\n";
   print INPUT_GFF "Pk_strainH_chr01	chado	mRNA	108	1421	.	+	.	ID=PKH_010010.1;Parent=PKH_010010;isObsolete=false;isFminPartial;feature_id=222342;timelastmodified=10.12.2012+02:43:44+GMT;previous_systematic_id=PK00_0600c\n";
   print INPUT_GFF "Pk_strainH_chr01	chado	gene	5697	8325	.	+	.	ID=PKH_010020;isObsolete=false;feature_id=222257;timelastmodified=25.11.2012+01:47:32+GMT\n";
   print INPUT_GFF "Pk_strainH_chr01	chado	CDS	5697	6209	.	+	0	ID=PKH_010020.1:exon:1;Parent=PKH_010020.1;isObsolete=false;timelastmodified=25.11.2012+01:47:32+GMT;colour=7\n";
   print INPUT_GFF "Pk_strainH_chr01	chado	CDS	7891	8325	.	+	0	ID=PKH_010020.1:exon:2;Parent=PKH_010020.1;isObsolete=false;timelastmodified=25.11.2012+01:47:32+GMT;colour=7\n";
   print INPUT_GFF "Pk_strainH_chr01	chado	mRNA	5697	8325	.	+	.	ID=PKH_010020.1;Parent=PKH_010020;isObsolete=false;feature_id=222258;timelastmodified=25.11.2012+01:47:32+GMT;previous_systematic_id=PK9_3420w\n";
   print INPUT_GFF "Pk_strainH_chr01	chado	gene	8988	9593	.	-	.	ID=PKH_010030;isObsolete=false;feature_id=222118;isFmaxPartial;timelastmodified=12.01.2011+05:18:24+GMT\n";
   print INPUT_GFF "Pk_strainH_chr01	chado	CDS	8988	9593	.	-	0	ID=PKH_010030.1:exon:1;Parent=PKH_010030.1;isObsolete=false;timelastmodified=12.01.2011+05:18:24+GMT;colour=12\n";
   print INPUT_GFF "Pk_strainH_chr01	chado	mRNA	8988	9593	.	-	.	ID=PKH_010030.1;Parent=PKH_010030;isObsolete=false;feature_id=222119;timelastmodified=21.10.2007+01:51:33+BST;previous_systematic_id=PK4_2020w\n";
   print INPUT_GFF "Pk_strainH_chr01	chado	gap	10499	12498	.	+	.	ID=gap10499-12498:corrected;isObsolete=false;feature_id=222641;timelastmodified=21.10.2007+01:51:37+BST\n";
   print INPUT_GFF "Pk_strainH_chr01	chado	gene	13925	14846	.	-	.	ID=PKH_010040;isObsolete=false;feature_id=222402;timelastmodified=21.10.2007+01:51:35+BST\n";
   print INPUT_GFF "Pk_strainH_chr01	chado	CDS	14778	14846	.	-	0	ID=PKH_010040.1:exon:2;Parent=PKH_010040.1;isObsolete=false;timelastmodified=21.10.2007+01:51:35+BST;colour=8\n";
   print INPUT_GFF "Pk_strainH_chr01	chado	CDS	13925	14581	.	-	0	ID=PKH_010040.1:exon:1;Parent=PKH_010040.1;isObsolete=false;timelastmodified=21.10.2007+01:51:35+BST;colour=8\n";
   print INPUT_GFF "Pk_strainH_chr01	chado	mRNA	13925	14846	.	-	.	ID=PKH_010040.1;Parent=PKH_010040;isObsolete=false;feature_id=222403;timelastmodified=21.10.2007+01:51:35+BST\n";
   print INPUT_GFF "Pk_strainH_chr01	chado	gene	18394	18738	.	+	.	ID=PKH_010050;isObsolete=false;feature_id=222197;timelastmodified=21.10.2007+01:51:33+BST\n";
   print INPUT_GFF "Pk_strainH_chr01	chado	CDS	18394	18738	.	+	0	ID=PKH_010050.1:exon:1;Parent=PKH_010050.1;isObsolete=false;timelastmodified=21.10.2007+01:51:33+BST;colour=10\n";
   print INPUT_GFF "Pk_strainH_chr01	chado	mRNA	18394	18738	.	+	.	ID=PKH_010050.1;Parent=PKH_010050;isObsolete=false;feature_id=222198;timelastmodified=21.10.2007+01:51:33+BST;previous_systematic_id=PK7_0005w\n";
   print INPUT_GFF "Pk_strainH_chr01	chado	gene	22127	22613	.	-	.	ID=PKH_010080;isObsolete=false;feature_id=222375;timelastmodified=25.11.2012+01:16:50+GMT\n";
   print INPUT_GFF "Pk_strainH_chr01	chado	CDS	22404	22613	.	-	0	ID=PKH_010080.1:exon:2;Parent=PKH_010080.1;isObsolete=false;timelastmodified=25.11.2012+01:16:50+GMT;colour=10\n";
   print INPUT_GFF "Pk_strainH_chr01	chado	CDS	22127	22228	.	-	0	ID=PKH_010080.1:exon:1;Parent=PKH_010080.1;isObsolete=false;timelastmodified=25.11.2012+01:16:50+GMT;colour=10\n";
   print INPUT_GFF "Pk_strainH_chr01	chado	mRNA	22127	22613	.	-	.	ID=PKH_010080.1;Parent=PKH_010080;isObsolete=false;feature_id=222376;timelastmodified=25.11.2012+01:16:50+GMT\n";
   print INPUT_GFF "Pk_strainH_chr01	chado	gene	32144	34243	.	+	.	ID=PKH_010100;isObsolete=false;feature_id=222114;timelastmodified=21.10.2007+01:51:33+BST\n";
   print INPUT_GFF "Pk_strainH_chr01	chado	CDS	32144	34243	.	+	0	ID=PKH_010100.1:exon:1;Parent=PKH_010100.1;isObsolete=false;timelastmodified=21.10.2007+01:51:33+BST;colour=7\n";
   print INPUT_GFF "Pk_strainH_chr01	chado	mRNA	32144	34243	.	+	.	ID=PKH_010100.1;Parent=PKH_010100;isObsolete=false;feature_id=222115;timelastmodified=21.10.2007+01:51:33+BST;previous_systematic_id=PK7_0010w\n";
   print INPUT_GFF "Pk_strainH_chr01	chado	gene	35321	35393	.	-	.	ID=PKH_010102;isObsolete=false;feature_id=222634;timelastmodified=21.10.2007+01:51:37+BST\n";
   print INPUT_GFF "Pk_strainH_chr01	chado	CDS	35321	35393	.	-	0	ID=PKH_010102.1:exon:1;Parent=PKH_010102:tRNA;isObsolete=false;timelastmodified=21.10.2007+01:51:37+BST\n";
   print INPUT_GFF "Pk_strainH_chr01	chado	tRNA	35321	35393	.	-	.	ID=PKH_010102:tRNA;Parent=PKH_010102;isObsolete=false;feature_id=222635;product=term%3DtRNA+Alanine%3B;timelastmodified=07.09.2012+12:23:39+BST;comment=tRNA+Ala+anticodon+AGC%2C+Cove+score+51.85\n";
   close(INPUT_GFF);
   ($input_fasta,$errorcode,$errormsg) = &make_filename($outputdir);
   if ($errorcode != 0) { ($errorcode,$errormsg) = &print_error($errormsg,$errorcode,0); } 
   #  --no-cache ensures it downloads the latest one
   system "wget --quiet --no-cache http://www.maths.tcd.ie/~avrillee/tests/get_spliced_transcripts_from_gff2.fa -O $input_fasta";
   ($output_fasta,$errorcode,$errormsg) = &make_filename($outputdir);
   if ($errorcode != 0) { ($errorcode,$errormsg) = &print_error($errormsg,$errorcode,0); } 
   @temp                  = split(/\//,$output_fasta);
   $output_fasta          = $temp[$#temp];
   &run_main_program($outputdir,$input_gff,$input_fasta,$output_fasta,'no','no',1,'no');
   # MAKE A FILE WITH THE EXPECTED CONTENTS OF $output_fasta:
   ($expected_output_fasta,$errorcode,$errormsg) = &make_filename($outputdir);
   if ($errorcode != 0) { ($errorcode,$errormsg) = &print_error($errormsg,$errorcode,0); } 
   open(EXPECTED,">$expected_output_fasta") || die "ERROR: test_run_main_program: cannot open expected_output_fasta $expected_output_fasta\n"; 
   print EXPECTED ">PKH_010010.1\n";
   print EXPECTED "GACGCCAAAGACGAATTAACGGAACTTCTGGATTATATGATAAAACCGGATAATCAGAAG\n";
   print EXPECTED "GACGTTGCCGAATACTGCAAAGACAAAGAAGATAAATGGAATGCAATGGGCCATAAAGAG\n";
   print EXPECTED "GGCAAAACAAATAAAGCAGCTTGTTTGCTTTTTGCTTCAGGATTAAAGCACATTTATACC\n";
   print EXPECTED "CATGGTATTGTCCAAAAGAATGGCCAGGTTAAGGATCATGTTAAAGGCCCATCGTTTGAA\n";
   print EXPECTED "CAAACGATGGGTTGTTTATTTCTTAAGGAATATGCAAAACAATTGAAAAAAATGGCAGAA\n";
   print EXPECTED "GAACAAAAAAAATATAGGGTACATCCTAATTGTAGTGTAGATTCTGGCATAGATCATGCT\n";
   print EXPECTED "TTTGAAAAAAGTAAAGTCATTATGCAAGCATCATCTCAATGCAAAAAGAATGTTAATAAC\n";
   print EXPECTED "GACTGCTTTGAATGCGACCTAAATAAAGGTTATGACAAATGCTCCATTGGCGATGACGAT\n";
   print EXPECTED "GTAGGAAATAAAGCTAAAAATCTGTTCGAAGGGGAATCGGAACAAAACCATATGCAACAA\n";
   print EXPECTED "ACTTTAGAGAATACAGTCTGTCCCATCCTTCTTACGGATATCCTTACCCCTTTTCTTCCT\n";
   print EXPECTED "TTGGCTCCTGTCTCTATTGGTCTTTCTGCTATGGCTTATTACCTTTGGAAGTATTTTGGT\n";
   print EXPECTED "CCTCTTGGTAAAGGAGGACCACGTTTCAGAAGATCTCCTGCTGAAATTCCTGGTTCATCG\n";
   print EXPECTED "ATACAAGAACACCTACTCGATCATGTGGAAGAAGCTGGTCCACATGAATATCAATTGGTG\n";
   print EXPECTED "AAGGAACGAAAACCTCGTTCTGCTCCAACGAGAACAAAACGTTCTGGTCCCGTGAATCGT\n";
   print EXPECTED "CGCACGATTATTGAAATTCATTTTGAAGTGTTGGATGAATGTCAAAAAGGGGACACACAA\n";
   print EXPECTED "GTGAACCAGAAGGATTTTCTGGAACTTTTGGTTCAAGAGTTCATGGGATCGGAATTAATG\n";
   print EXPECTED "GAAGAAGAACAGGTTCCTAAGGAAGAGCTTTTTATGGAAGGGGTTCCTATGGAAGAGGTT\n";
   print EXPECTED "CCTATGGAAAGTATTCCTTTAGAGCAGGTTCCAATGGAACGTGTTCCAAATTTAGGTTCC\n";
   print EXPECTED "GGGTTTATGGTTTAG\n";
   print EXPECTED ">PKH_010020.1\n";
   print EXPECTED "ATGTCAAGTTCGAGTTCAAGTTCAAGTTCAAGTTCCGCTTCTTCAGGTTCAGGTGGTGGT\n"; 
   print EXPECTED "TGGTGGGACGTAGGTCTTTCAGGTGGTGTTGCTCGTGCTCGTACGCCAGCACCAGCCGCA\n";
   print EXPECTED "CTTTCGGCTAATTCCAGTATAACGAGTAATGTATTTTCCGCTATTCGAGGATTTTGGAAT\n"; 
   print EXPECTED "TCTGGGTCACAAATTATTGAGAACGTGCGAACAAGAATTGCGACTGCGCTGTCTCCAGTT\n";
   print EXPECTED "ACATTGATCACATCAGCTTCATCCGCCATTAATAGTGTTTCTAATGTAGCTCCCTCAGTT\n";
   print EXPECTED "ATAAAGACAGTTCATACAGTTGTAGGAGGTAAACTAGAAACTTCTGCTCCAGCACCTACT\n"; 
   print EXPECTED "CCTCCTCAACCAGCCCAACCGGCTGCACGGCCTGAACCAGTTACTCCTGCTCCCCCCACG\n";
   print EXPECTED "AAGCCTACAAAACCATGGGATCGGCTTATCCCTTTTGTTGCGTTGGCTCCTGCTACAGTT\n";
   print EXPECTED "GCTATTTCTATTTTTTCCTACTTAATCTGGAAGTATTTTGCTCAACTGCGTAAGATAAGA\n";
   print EXPECTED "CTCTACAGAAGAGCTCCTTTAAGAATTCCCGGTCCATCCGTACAGGAACAACTCCTTGAT\n"; 
   print EXPECTED "CATGTGGAAGAAGCTGGTCCACATGAATATCAATTGGTGAAGGAACGAAAACCTCCATCT\n";
   print EXPECTED "GTTCCAGCGAGAACAAAGCGTTCTGGTCGCGATCCTGCTGATGGTGGTCGCGTAAATCGT\n";
   print EXPECTED "CGAACGATTATTAAAATTCATTTTGAACTGGTGGACGAATGTCAAAAAGGGAACACAAAA\n";
   print EXPECTED "TTGACTCAGAATGATTTTCTGGAACTTTTGGTTCAAGAGTTTATGGGATCCGAATTAATG\n";
   print EXPECTED "GAAGAAGAAGAACAGGTTTCTAAGGAAGAGGTTCTTATGCAAGGGGTTCCTATGGAAAGT\n";
   print EXPECTED "GTTCCAATGGAACGTGTTCCAATTTTAGGTTCCGTGTTTATGGTTTAG\n";
   print EXPECTED ">PKH_010030.1\n";
   print EXPECTED "CACGAAGATGATTTTTCAAGACAGGGTGTGAGTACTTTATTACAAGTACATAATAAAGTA\n";
   print EXPECTED "GGTAGAAATTTAGCAAGCAAAGAAGAAAAAAAAGAAGAAAAAAAAGAAGAAAAGAAAGAA\n";
   print EXPECTED "GAAAAGAAAGAAGAAAAAAAAGAAGAAAAGAAAAAAACAAATAAAACATTATTAGGGAAA\n";
   print EXPECTED "AAAGGGAAAAAGGAATATAAGGATGGTAGATCCACATTTGTAGGAGCTAGAGATATTAAC\n"; 
   print EXPECTED "AAAGAGAAGAGTGAAATGGATAGTTATCGACAGAGAAAATTAGATTTCTGGGAACATTTC\n";
   print EXPECTED "GAACCAACAATGACAGCAAATTTTGAGAAGGTATTAAAAAGGTGTTTTGCACGTAAAGCA\n";
   print EXPECTED "GGAGAAGAAGACATGGATGAAGATTATTCAGTTGTACTTCCAAAGCTGGGATGGAATGCG\n";
   print EXPECTED "GATCCATATGGTGTATTAAAGAAAACAAAGAAAGGGACGGTGCAGGAAGGATATAGGAAA\n"; 
   print EXPECTED "ATGCTAGAAAACAATTTCCGCAATGTTCCATACGTACATGATTCACAGAATGCTCATAAA\n";
   print EXPECTED "TCAAACGAACCATACCTAGAAAGAGAATACGGACGCGTAGAGTTGGGCGCTGATGCAAAA\n";
   print EXPECTED "ACATGA\n"; 
   print EXPECTED ">PKH_010040.1\n";
   print EXPECTED "ATGATATCTCTTAGGAAACTTCCTCTATTCATCTTTGTTGTATACTCCTGGCAATGTTTC\n";
   print EXPECTED "GGCAACTATGTGAGTGCATTAAGCACGTTAGATGATAGCACTTCCATGGAAAGTCTATTT\n"; 
   print EXPECTED "ACGAAAGTGCAAAAGACCAAATTTCCCCACCATAATGGTCATCATCATCATAATGATGGT\n";
   print EXPECTED "TATGGTGGTCATAATGACAGCATGTCTACCCTTTCGAGCAGTGTAGACACCATTTTTAAA\n"; 
   print EXPECTED "AAAGAAGGAATTCGTACAAATTATAAAGCCTCATACAAGGCACCACATAAAACACTCTAT\n";
   print EXPECTED "GAAGCACTACATGAAATCCCATATAGAGCCCCACATGAAGTGTCATATAGAGGCCCACAT\n";
   print EXPECTED "AAAGTCCCATATAAAGCCCTACATGAATTCTCATATAGAGCCCCACATGAAGTGTCATAT\n";
   print EXPECTED "AGAGGCCCACATAAAACTCTCTTCGATGCACCGCATATGGCACGATATGGACGAATGAAA\n";
   print EXPECTED "ATTCGACCTAGAGGATTTTTTAATAAATTGCTATACAATATAGAGAAAATTAAAAGAACA\n";
   print EXPECTED "GGAGCATTGTTATCCCCTAAGGTAATTCGCATATCGCTCTTATATATATTAGCTTTGTTT\n"; 
   print EXPECTED "TGTATAGGTTTGATTGTGGTCCTTCCTATAGGTCCCCAAACACCTTTTGCATGGATGACT\n";
   print EXPECTED "TGTGGTATGGCAATTTCTGTATTAATAGCATGTATATATTCTAGTGTTCGGATTGCAATA\n"; 
   print EXPECTED "AAATAA\n";
   print EXPECTED ">PKH_010050.1\n";
   print EXPECTED "ATGTATCAGCAAAGGAATTTCGACAAAGGCGACCCCACTTCGACGCACCAACGCCAATTT\n"; 
   print EXPECTED "GAAGAGTCAGATGACGGTATACCAGGTTATGGAATCCCTCCAAACCCAACCATGATAAAC\n"; 
   print EXPECTED "CTCACTGGTAACCAAGACTCACGATCAAACCTAATGCAACAATTTGGAATAAACAACAAA\n";
   print EXPECTED "ACTGTATCGCAGTTTTTAGTAAACATGTTCGTCTACGTTGCTGCTATATTAATTAGTTTA\n";
   print EXPECTED "AAAATATGGGACTACATGTCTTATAGCAAATGTGATTATTATAAAGATTTATTATTAAGA\n";
   print EXPECTED "ATTGTAAGATACCAGTCACATATGAATGATGTTAAGATGGCCTAA\n"; 
   print EXPECTED ">PKH_010080.1\n";
   print EXPECTED "ATGTTAGTGAGGATGTTCCGAATCACAGGTAGCACAAACCCCTATACAACAACTGCAGAT\n";
   print EXPECTED "GATATTAGTGCGATTAAGCATGCTCTTGGAGGATGCCATTTAAACAGGGAAGGTTATTTT\n";
   print EXPECTED "TTCACGGACCCTCGATTCTTGCAAGGAGTTTTTGCTGGTATCCTTGTGTTTTATGTAATG\n"; 
   print EXPECTED "CTATACTATACTCGCCGATACAAAAATAGGGTTAATGACGTGCCAAGACATGAATCATTT\n";
   print EXPECTED "TCCTATGTTAACCTGCCAACTCACGGTCCATATAGAGAGTTCCATGGAAGCGAAATGAAT\n";
   print EXPECTED "CAAATACATTAG\n";
   print EXPECTED ">PKH_010100.1\n";
   print EXPECTED "ATGAATACCTTCAAGTCTCTCTTCTTAATTTTTGCTTCTGTGCTTTATTGCACACATTCC\n";
   print EXPECTED "GTAAGTCTTAAGGGCAGAAATACCAGAGGAGATAATCATGAACTAAATATGATCATGGAT\n";
   print EXPECTED "GGCGTAAGTAAACTACAACTGACCAAAGCGCACAACCAATCATTCATTTCTGGCTTATAC\n";
   print EXPECTED "ACCGACGATTCGAAAAAAGTCATTTGCGAATGTACATGCACAGGAAGTAATAACATTGAG\n";
   print EXPECTED "GGTGGAAGCGGAAGCGGAAGCGGTTCCGGAAGTGGCAGCGGTTCCGGAAGTGGCAGCGGT\n";
   print EXPECTED "TCCGGAAGTGGCAGCGGTTCAGGAAGTGGAAGCGGAAGCGGAAGCGGTTCCGGAAGTGGC\n";
   print EXPECTED "AGCGGTTCCGGAAGTGGCAGCGGTTCAGGAAGTGGCAGCGGTTCAGGAAGTGGAAGCGGA\n";
   print EXPECTED "ACAGGAAGTGGAAGCGGAACAGGAAGTGGCAGCGGAACAGGAAGTGGCAGCGGAACAGGA\n";
   print EXPECTED "AGTGGCAGCGGAACAGGAAGTGGCAGCGGTTCAGGAAGTGGCAGCGGTTCAGGAAGTGGC\n";
   print EXPECTED "AGCGGTTCAGGAAGCGGCAGCGGTTCAGGAAGCGGCAGCGGTTCAGGAAGTGGTAGCCAT\n";
   print EXPECTED "GATCGAAAACCTCCAAGAGAAATTTTAGAGGAATATAAAAGAAGAAAACAAGGAATAGTT\n";
   print EXPECTED "GCAGGATATTATGGTTCATGGAATAGTCAAGGAGATAGGGCCAAAAATATGACAGATTCT\n";
   print EXPECTED "AGCTCAATGGTATCAATACTGTACATCGCGTTTGCTCGTATAAATATGTTTTATGATGTA\n";
   print EXPECTED "TCCAGACCATTTAACGGAAAGCAGAAATTTTTAGTAAGGAAACATGGACTGGAATATGAA\n";
   print EXPECTED "ACTTACGGAGTTATGCTAAACGAAATCAGGCGTATAAGGAAAGCTCGCCCAGACATCATA\n";
   print EXPECTED "TTAATTCTATCATTAGGAGGAGAGACGTACATGATAGATATTACAAAAGATATTGACTAT\n";
   print EXPECTED "ATGGAACAAATAGTTAAGCTTGTGAAAGATTTTGACTTAGATGGTGTAGATATTGACTGG\n";
   print EXPECTED "GAACCACATGGCAGCTTTAACAATCTAAACGAACTGAACTATTCAGAATATTATATTAAA\n";
   print EXPECTED "TTAATCAACTTAGTAAGAAGTAGCATTCCAGAGGAGAAATTAATCTCCATATCGGGATCA\n";
   print EXPECTED "TCCAACGCAGCATTATCATGCGTGTCCGTGAATGAAACATTTTGCAAAGACGATGACTCT\n";
   print EXPECTED "CCATATAACACGAACTACTTGTCTGAACAAATGGGTAAGAACGATGAGTTATACAAAGCC\n";
   print EXPECTED "TCAACTATGTTATCTACAGGAACATTTGTCAATATTTTTAACAGAGCAAAGGACAAAATT\n";
   print EXPECTED "GACCTTGTATTTATTCAGACGTACAACTTAGAGACGACGAACCCAAGTATTATGGTAGAC\n";
   print EXPECTED "ATGTATCTATCCCATCTCTATTTCGGATTGAAATATAACATCACAGTTTTGCTAGGATTT\n";
   print EXPECTED "TCACTGGAACATAACCGAGGTGGCTTCAGCCCAGATGACAGAGCATTAGTCGAATTAGTA\n";
   print EXPECTED "TCAAAAACTATTCATGACGAAAATCATAAGCACAATCGGGCAGACGGAGTAGGAATATGG\n";
   print EXPECTED "CACTTGTTCATGAAAGAACAGATGCCAAGTGGATCATATGATATAGATGCGTTCCTCACG\n";
   print EXPECTED "AATGTATGGAAACATTTAAATCCCCAGGTACAAGTACCAAAAGATGTAGTTACTACTCAG\n";
   print EXPECTED "AACCCTGACGACTGTAACAGTATAGATGAATATATTTCTGGACTTGTTGTGTCTAAATCA\n";
   print EXPECTED "GGCGTGTATTACAAACACAATGGTGCAATATGGAAAACTAGGTCATACTCAACCCGTGCT\n";
   print EXPECTED "CCTGGTGTTGACAGATATGAATGGGACTTGGTCAAAATATGTTATGAAAAAGCGTGCAAC\n";
   print EXPECTED "GGTATGGCAGCCCATTACTTTAACACTGATTATCAAAATGGATCTATAGTCATATGGAAA\n";
   print EXPECTED "GGAGAAGCTTTCACTATTAAATGGTGGCAATCTGGACCTCCTGAAGGAGCGGCGTTGGAA\n";
   print EXPECTED "GCATATGAAAAATTGGAAGCATCCGCATGTCCAGGACTCTCGGAATGGAACAAGGAGCAC\n";
   print EXPECTED "CCACACAAGCCACTTGAGGAAGACATTCCCTATGAGCAAGAGCAAGACGCCCCATTATAA\n";
   print EXPECTED ">PKH_010102:tRNA\n";
   print EXPECTED "GGGCGACTAGCTCAAGTGGTAGAGCGCTCGCTTAGCATGCGAGAGGTACGGGGATCGATA\n";
   print EXPECTED "CCCCGGTCGTCCA\n"; 
   close(EXPECTED);
   $differences            = "";
   $output_fasta           = $outputdir."/".$output_fasta;
   open(TEMP,"diff $output_fasta $expected_output_fasta |");
   while(<TEMP>)
   {
      $line                = $_;
      $differences         = $differences.$line;
   }
   close(TEMP);  
   $length_differences     = length($differences);
   if ($length_differences != 0) { print STDERR "ERROR: test_run_main_program: failed test1 (output_fasta $output_fasta expected_output_fasta $expected_output_fasta)\n"; exit;}
   system "rm -f $expected_output_fasta"; 
   system "rm -f $output_fasta";
   system "rm -f $input_gff"; 
   system "rm -f $input_fasta";
 
   ($input_gff,$errorcode,$errormsg) = &make_filename($outputdir);
   if ($errorcode != 0) { ($errorcode,$errormsg) = &print_error($errormsg,$errorcode,0); }   
   open(INPUT_GFF,">$input_gff") || die "ERROR: test_run_main_program: cannot open $input_gff\n";
   print INPUT_GFF "CHROMOSOME_X\tCoding_transcript\tCDS\t2446067\t2446593\t.\t-\t2\tID=CDS:T07D1.4;Note=RNA-binding protein;Parent=Transcript:T07D1.4.1,Transcript:T07D1.4.2;locus=fox-1;status=Partially_confirmed;wormpep=CE:CE25105\n";
   print INPUT_GFF "CHROMOSOME_X\tCoding_transcript\tCDS\t2446639\t2446783\t.\t-\t0\tID=CDS:T07D1.4;Note=RNA-binding protein;Parent=Transcript:T07D1.4.1,Transcript:T07D1.4.2;locus=fox-1;status=Partially_confirmed;wormpep=CE:CE25105\n";
   print INPUT_GFF "CHROMOSOME_X\tCoding_transcript\tCDS\t2446839\t2446975\t.\t-\t2\tID=CDS:T07D1.4;Note=RNA-binding protein;Parent=Transcript:T07D1.4.1,Transcript:T07D1.4.2;locus=fox-1;status=Partially_confirmed;wormpep=CE:CE25105\n";
   print INPUT_GFF "CHROMOSOME_X\tCoding_transcript\tCDS\t2447335\t2447530\t.\t-\t0\tID=CDS:T07D1.4;Note=RNA-binding protein;Parent=Transcript:T07D1.4.1,Transcript:T07D1.4.2;locus=fox-1;status=Partially_confirmed;wormpep=CE:CE25105\n";
   print INPUT_GFF "CHROMOSOME_X\tCoding_transcript\tCDS\t2448754\t2448842\t.\t-\t2\tID=CDS:T07D1.4;Note=RNA-binding protein;Parent=Transcript:T07D1.4.1,Transcript:T07D1.4.2;locus=fox-1;status=Partially_confirmed;wormpep=CE:CE25105\n";
   print INPUT_GFF "CHROMOSOME_X\tCoding_transcript\tCDS\t2450334\t2450487\t.\t-\t0\tID=CDS:T07D1.4;Note=RNA-binding protein;Parent=Transcript:T07D1.4.1,Transcript:T07D1.4.2;locus=fox-1;status=Partially_confirmed;wormpep=CE:CE25105\n";
   print INPUT_GFF "CHROMOSOME_X\tCoding_transcript\tgene\t2445794\t2450500\t.\t-\t.\tID=Gene:T07D1.4\n";
   print INPUT_GFF "CHROMOSOME_X\tCoding_transcript\tmRNA\t2445794\t2450496\t.\t-\t.\tID=Transcript:T07D1.4.1;Parent=Gene:T07D1.4;cds=T07D1.4;wormpep=CE:CE25105\n";
   print INPUT_GFF "CHROMOSOME_X\tCoding_transcript\tmRNA\t2445794\t2450500\t.\t-\t.\tID=Transcript:T07D1.4.2;Parent=Gene:T07D1.4;cds=T07D1.4;wormpep=CE:CE25105\n";
   close(INPUT_GFF);
   ($input_fasta,$errorcode,$errormsg) = &make_filename($outputdir);
   if ($errorcode != 0) { ($errorcode,$errormsg) = &print_error($errormsg,$errorcode,0); }  
   system "wget --quiet http://www.maths.tcd.ie/~avrillee/tests/Ce_WS226_CHRX.fa -O $input_fasta";
   ($output_fasta,$errorcode,$errormsg) = &make_filename($outputdir);
   if ($errorcode != 0) { ($errorcode,$errormsg) = &print_error($errormsg,$errorcode,0); }
   @temp                  = split(/\//,$output_fasta);
   $output_fasta          = $temp[$#temp];
   &run_main_program($outputdir,$input_gff,$input_fasta,$output_fasta,'no','no',1,'no');
   # MAKE A FILE WITH THE EXPECTED CONTENTS OF $output_fasta:
   ($expected_output_fasta,$errorcode,$errormsg) = &make_filename($outputdir);
   if ($errorcode != 0) { ($errorcode,$errormsg) = &print_error($errormsg,$errorcode,0); }
   open(EXPECTED,">$expected_output_fasta") || die "ERROR: test_run_main_program: cannot open expected_output_fasta $expected_output_fasta\n";
   print EXPECTED ">Transcript:T07D1.4.1\n";
   print EXPECTED "ATGCAAGCCCTGTACCAACTGTCTGCCACAGGGGCTCAACAACAGAATCAACAAATTCCG\n";
   print EXPECTED "ATTGGATTGAGCAACTCATTGCTCTATCAGCAATTGGCGGCACATCAGCAAATTGCTGCT\n";
   print EXPECTED "CAACAGCATCAGCAACAGCTCGCTGTCTCTGCAGCTCATCAAACACAAAACAATATAATG\n";
   print EXPECTED "CTTGCTACCAGCGCTCCATCATTGATTAATCATATGGAGAACTCGACGGACGGTAAAGTC\n";
   print EXPECTED "AAGGACGATCCGAACAGTGATTACGATTTGCAACTCTCGATTCAGCAACAATTGGCGGCG\n";
   print EXPECTED "GCAGCGCAAGCTGCTCAGATGGGACAAACGCAGATTGGTCCACAAATTGTGGGACAACAA\n";
   print EXPECTED "GGGCAGCCGGTAGTCGCCACAACGGCCGGTTCGACGAATGGCTCGGCGGCCGTCACACAG\n";
   print EXPECTED "CCCGATCCCAGCACTAGCTCCGGACCCGATGGACCAAAAAGACTTCACGTATCGAATATC\n";
   print EXPECTED "CCATTCAGATTCAGAGACCCAGACCTCAAAACGATGTTTGAGAAATTTGGTGTGGTTTCC\n";
   print EXPECTED "GACGTAGAAATTATTTTCAATGAGCGAGGATCAAAGGGGTTTGGATTTGTGACAATGGAG\n";
   print EXPECTED "AGACCGCAGGATGCTGAGAGAGCTAGACAAGAGCTTCATGGATCAATGATTGAAGGACGG\n";
   print EXPECTED "AAAATTGAAGTCAACTGCGCTACAGCTCGTGTTCACTCGAAGAAAGTTAAACCAACTGGA\n";
   print EXPECTED "GGAATCCTGGACCAAATGAACCCACTGATGGCCCAATCAGCTCTTGCCGCACAGGCTCAG\n";
   print EXPECTED "ATGAACAGAGCCCTATTGCTCCGTAGTCCATTGGTAGCACAATCTCTTCTTGGTCGCGGA\n";
   print EXPECTED "GCCGCTTTAATCCCCGGAATGCAACAGCCCGCGTTCCAATTGCAAGCGGCTTTGGCTGGC\n";
   print EXPECTED "AATCCGCTGGCACAACTTCAAGGTCAACCTTTGCTGTTCAACGCTGCAGCCCTTCAAACG\n";
   print EXPECTED "AACGCACTCCAACAGTCGGCGTTTGGAATGGATCCGGCCGCCGTTCAAGCTGCTCTGCTT\n";
   print EXPECTED "GCCAATGAGCAAGCTCGGTTTCAACTCGCTGCTGCCGCTGCTCAAGGTAACGAGTACATT\n";
   print EXPECTED "ATGTACCATCAAGCCAAACAACAAGAATTACCAGGACGGATCCCATCATCTGGAAATGCC\n";
   print EXPECTED "TCGGCGTTTGGCGAACAATACCTTAGCAACGCGTTGGCCACCGCCTCACTCCCCTCATAT\n";
   print EXPECTED "CAAATGAACCCGGCGCTTAGAACGTTAAATCGATTTACTCCGTATTGA\n";
   print EXPECTED ">Transcript:T07D1.4.2\n";
   print EXPECTED "ATGCAAGCCCTGTACCAACTGTCTGCCACAGGGGCTCAACAACAGAATCAACAAATTCCG\n";
   print EXPECTED "ATTGGATTGAGCAACTCATTGCTCTATCAGCAATTGGCGGCACATCAGCAAATTGCTGCT\n";
   print EXPECTED "CAACAGCATCAGCAACAGCTCGCTGTCTCTGCAGCTCATCAAACACAAAACAATATAATG\n";
   print EXPECTED "CTTGCTACCAGCGCTCCATCATTGATTAATCATATGGAGAACTCGACGGACGGTAAAGTC\n";
   print EXPECTED "AAGGACGATCCGAACAGTGATTACGATTTGCAACTCTCGATTCAGCAACAATTGGCGGCG\n";
   print EXPECTED "GCAGCGCAAGCTGCTCAGATGGGACAAACGCAGATTGGTCCACAAATTGTGGGACAACAA\n";
   print EXPECTED "GGGCAGCCGGTAGTCGCCACAACGGCCGGTTCGACGAATGGCTCGGCGGCCGTCACACAG\n";
   print EXPECTED "CCCGATCCCAGCACTAGCTCCGGACCCGATGGACCAAAAAGACTTCACGTATCGAATATC\n";
   print EXPECTED "CCATTCAGATTCAGAGACCCAGACCTCAAAACGATGTTTGAGAAATTTGGTGTGGTTTCC\n";
   print EXPECTED "GACGTAGAAATTATTTTCAATGAGCGAGGATCAAAGGGGTTTGGATTTGTGACAATGGAG\n";
   print EXPECTED "AGACCGCAGGATGCTGAGAGAGCTAGACAAGAGCTTCATGGATCAATGATTGAAGGACGG\n";
   print EXPECTED "AAAATTGAAGTCAACTGCGCTACAGCTCGTGTTCACTCGAAGAAAGTTAAACCAACTGGA\n";
   print EXPECTED "GGAATCCTGGACCAAATGAACCCACTGATGGCCCAATCAGCTCTTGCCGCACAGGCTCAG\n";
   print EXPECTED "ATGAACAGAGCCCTATTGCTCCGTAGTCCATTGGTAGCACAATCTCTTCTTGGTCGCGGA\n";
   print EXPECTED "GCCGCTTTAATCCCCGGAATGCAACAGCCCGCGTTCCAATTGCAAGCGGCTTTGGCTGGC\n";
   print EXPECTED "AATCCGCTGGCACAACTTCAAGGTCAACCTTTGCTGTTCAACGCTGCAGCCCTTCAAACG\n";
   print EXPECTED "AACGCACTCCAACAGTCGGCGTTTGGAATGGATCCGGCCGCCGTTCAAGCTGCTCTGCTT\n";
   print EXPECTED "GCCAATGAGCAAGCTCGGTTTCAACTCGCTGCTGCCGCTGCTCAAGGTAACGAGTACATT\n";
   print EXPECTED "ATGTACCATCAAGCCAAACAACAAGAATTACCAGGACGGATCCCATCATCTGGAAATGCC\n";
   print EXPECTED "TCGGCGTTTGGCGAACAATACCTTAGCAACGCGTTGGCCACCGCCTCACTCCCCTCATAT\n";
   print EXPECTED "CAAATGAACCCGGCGCTTAGAACGTTAAATCGATTTACTCCGTATTGA\n";
   close(EXPECTED);
   $differences            = "";
   $output_fasta           = $outputdir."/".$output_fasta;
   open(TEMP,"diff $output_fasta $expected_output_fasta |");
   while(<TEMP>)
   {
      $line                = $_;
      $differences         = $differences.$line;
   }
   close(TEMP);
   $length_differences     = length($differences);
   if ($length_differences != 0) { print STDERR "ERROR: test_run_main_program: failed test2 (output_fasta $output_fasta expected_output_fasta $expected_output_fasta)\n"; exit;}
   system "rm -f $expected_output_fasta";
   system "rm -f $output_fasta";
   system "rm -f $input_gff";
   system "rm -f $input_fasta";
 
}
 
#------------------------------------------------------------------#
 
# RUN THE MAIN PART OF THE CODE:
 
sub run_main_program
{
   my $outputdir           = $_[0]; # DIRECTORY TO PUT OUTPUT FILES IN.
   my $input_gff           = $_[1]; # THE INPUT GFF FILE  
   my $input_fasta         = $_[2]; # THE INPUT FASTA FILE
   my $output_fasta        = $_[3]; # THE OUTPUT FASTA FILE 
   my $ignore_phase        = $_[4]; # SAYS WHETHER TO IGNORE PHASE INFORMATION, AND ASSUME THE PHASE OF THE FIRST EXON IS 0 
   my $ignore_semicolons   = $_[5]; #
   my $testing             = $_[6]; # SAYS WHETHER THIS IS BEING CALLED FROM A TESTING SUBROUTINE
   my $from_augustus       = $_[7]; # SAYS WHETHER THE GFF FILE IS FROM AUGUSTUS
   my $errorcode;                   # RETURNED AS 0 IF THERE IS NO ERROR.
   my $errormsg;                    # RETURNED AS 'none' IF THERE IS NO ERROR. 
   my $exon_gff;                    # GFF FILE OF EXONS 
   my $SEQ;                         # HASH TABLE TO STORE THE SEQUENCES OF SCAFFOLDS/CHROMOSOMES
   my $EXONSEQ;                     # HASH TABLE OF EXON SEQUENCES
   my $EXONS;                       # HASH TABLE OF EXONS IN EACH GENE
   my $EXON_START;                  # HASH TABLE WITH THE START POSITIONS OF EXONS
   my $EXON_FRAME;                  # HASH TABLE WITH THE PHASES OF EXONS 
 
   # CHECK THAT $ignore_phase IS 'yes' OR 'no':
   if ($ignore_phase ne 'yes' && $ignore_phase ne 'no')
   {
      print STDERR "ERROR: run_main_program: ignore_phase $ignore_phase (should be yes/no)\n";
      exit;
   }
 
   # CHECK THAT $from_augustus IS 'yes' OR 'no':
   if ($from_augustus ne 'yes' && $from_augustus ne 'no')
   {
      print STDERR "ERROR: run_main_program: from_augustus $from_augustus (should be yes/no)\n";
      exit;
   }
 
   # MAKE A GFF FILE OF THE EXONS:
   if ($testing == 0) { print STDERR "Making exon gff file...\n";}
   ($exon_gff,$errorcode,$errormsg) = &make_exon_gff($input_gff,$outputdir);
   if ($errorcode != 0) { ($errorcode,$errormsg) = &print_error($errormsg,$errorcode,0); }
 
   # READ IN THE SEQUENCES OF THE SCAFFOLDS/CHROMOSOMES:
   if ($testing == 0) { print STDERR "Reading in sequences of scaffolds...\n";}
   ($SEQ,$errorcode,$errormsg) = &read_assembly($input_fasta,$ignore_semicolons);
   if ($errorcode != 0) { ($errorcode,$errormsg) = &print_error($errormsg,$errorcode,0); }
 
   # GET THE DNA SEQUENCES FOR THE EXONS IN THE EXON GFF FILE:
   if ($testing == 0) { print STDERR "Getting DNA sequences for the exons...\n";}
   ($EXONSEQ,$errorcode,$errormsg) = &get_exon_sequences($exon_gff,$SEQ,$outputdir);
   if ($errorcode != 0) { ($errorcode,$errormsg) = &print_error($errormsg,$errorcode,0); }
 
   # READ IN THE EXONS FROM EACH GENE FROM THE GFF FILE:
   if ($testing == 0) { print STDERR "Reading in the exons for each gene from the gff...\n";}
   ($EXONS,$EXON_START,$EXON_FRAME,$errorcode,$errormsg) = &read_exon_gff($exon_gff,$ignore_phase,$from_augustus);
   if ($errorcode != 0) { ($errorcode,$errormsg) = &print_error($errormsg,$errorcode,0); }
  
   # INFER THE SPLICED DNA SEQUENCE FOR EACH GENE:
   if ($testing == 0) { print STDERR "Inferring spliced DNA sequences for transcripts...\n";}
   ($errorcode,$errormsg)  = &infer_spliced_dna($exon_gff,$EXONS,$EXON_START,$EXON_FRAME,$EXONSEQ,$output_fasta,$outputdir,$ignore_phase,$from_augustus);
   if ($errorcode != 0) { ($errorcode,$errormsg) = &print_error($errormsg,$errorcode,0); }
 
   # DELETE TEMPORARY FILES:
   system "rm -f $exon_gff";
}
 
#------------------------------------------------------------------#
 
# TEST &infer_spliced_dna
 
sub test_infer_spliced_dna
{
   my $outputdir           = $_[0]; # DIRECTORY TO WRITE OUTPUT FILES INTO
   my $exon_gff;                    # EXON GFF FILE
   my $errorcode;                   # RETURNED AS 0 FROM A FUNCTION IF THERE IS NO ERROR
   my $errormsg;                    # RETURNED AS 'none' FROM A FUNCTION IF THERE IS NO ERROR 
   my %EXONS               = ();    # HASH TABLE TO RECORD THE EXONS IN A GENE
   my %EXON_START          = ();    # HASH TABLE TO RECORD THE START POSITION OF EXONS
   my %EXON_FRAME          = ();    # HASH TABLE TO RECORD THE SEQUENCE OF EXONS  
   my %EXONSEQ             = ();    # HASH TABLE TO STORE THE SEQUENCE OF EXONS 
   my $output_fasta;                # OUTPUT FASTA FILE WITH THE SPLICED DNA SEQUENCE
   my $expected_output_fasta;       # FILE CONTAINING THE EXPECTED CONTENTS OF $output_fasta 
   my $differences;                 # DIFFERENCES BETWEEN $output_fasta AND $expected_output_fasta 
   my $length_differences;          # LENGTH OF $differences
   my $line;                        #  
   my @temp;                        #  
 
   ($exon_gff,$errorcode,$errormsg) = &make_filename($outputdir);
   if ($errorcode != 0) { ($errorcode,$errormsg) = &print_error($errormsg,$errorcode,0); } 
   open(EXON_GFF,">$exon_gff") || die "ERROR: test_get_exon_sequences: cannot open $exon_gff\n";
   print EXON_GFF "chr1\tchado\tCDS\t11\t20\t.\t+\t0\tID=gene1:exon1;Parent=gene1;\n";
   print EXON_GFF "chr1\tchado\tCDS\t31\t40\t.\t+\t2\tID=gene1:exon2;Parent=gene1;\n";
   close(EXON_GFF);
   $EXONS{"gene1"} = "ID=gene1:exon1*,*ID=gene1:exon2";
   $EXON_START{"ID=gene1:exon1"} = 11;
   $EXON_START{"ID=gene1:exon2"} = 31;
   $EXON_FRAME{"ID=gene1:exon1"} = "zero";
   $EXON_FRAME{"ID=gene1:exon2"} = 2;
   $EXONSEQ{"ID=gene1:exon1"} = "CCCCCGGGGG"; # 10 BASES LONG, SO THE NEXT EXON SHOULD START IN PHASE 2
   $EXONSEQ{"ID=gene1:exon2"} = "CGCGCGCGCG";
   ($output_fasta,$errorcode,$errormsg) = &make_filename($outputdir);
   if ($errorcode != 0) { ($errorcode,$errormsg) = &print_error($errormsg,$errorcode,0); } 
   @temp                  = split(/\//,$output_fasta);
   $output_fasta          = $temp[$#temp];
   # CALL &infer_spliced_dna:
   ($errorcode,$errormsg) = &infer_spliced_dna($exon_gff,\%EXONS,\%EXON_START,\%EXON_FRAME,\%EXONSEQ,$output_fasta,$outputdir,'no','no');
   if ($errorcode != 0) { print STDERR "ERROR: test_infer_spliced_dna: failed test1 (errorcode $errorcode errormsg $errormsg)\n"; exit;}
   ($expected_output_fasta,$errorcode,$errormsg) = &make_filename($outputdir);
   if ($errorcode != 0) { ($errorcode,$errormsg) = &print_error($errormsg,$errorcode,0); } 
   open(EXPECTED,">$expected_output_fasta") || die "ERROR: test_infer_spliced_dna: cannot open $expected_output_fasta\n";
   print EXPECTED ">gene1\n";
   print EXPECTED "CCCCCGGGGGCGCGCGCGCG\n";
   close(EXPECTED);
   $differences            = "";
   $output_fasta           = $outputdir."/".$output_fasta;
   open(TEMP,"diff $output_fasta $expected_output_fasta |");
   while(<TEMP>)
   {
      $line                = $_;
      $differences         = $differences.$line;
   }
   close(TEMP);  
   $length_differences     = length($differences);
   if ($length_differences != 0) { print STDERR "ERROR: test_infer_spliced_dna: failed test1 (output_fasta $output_fasta expected_output_fasta $expected_output_fasta)\n"; exit;}
   system "rm -f $exon_gff";
   system "rm -f $output_fasta";
   system "rm -f $expected_output_fasta";
 
   # TEST ERRORCODE=29 (A GENE IN THE GFF IS NOT IN %EXONS):
   ($exon_gff,$errorcode,$errormsg) = &make_filename($outputdir);
   if ($errorcode != 0) { ($errorcode,$errormsg) = &print_error($errormsg,$errorcode,0); } 
   open(EXON_GFF,">$exon_gff") || die "ERROR: test_get_exon_sequences: cannot open $exon_gff\n";
   print EXON_GFF "chr1\tchado\tCDS\t11\t20\t.\t+\t0\tID=gene1:exon1;Parent=gene2;\n";
   print EXON_GFF "chr1\tchado\tCDS\t31\t40\t.\t+\t2\tID=gene1:exon2;Parent=gene2;\n";
   close(EXON_GFF);
   $EXONS{"gene1"} = "ID=gene1:exon1*,*ID=gene1:exon2";
   $EXON_START{"ID=gene1:exon1"} = 11;
   $EXON_START{"ID=gene1:exon2"} = 31;
   $EXON_FRAME{"ID=gene1:exon1"} = "zero";
   $EXON_FRAME{"ID=gene1:exon2"} = 2;
   $EXONSEQ{"ID=gene1:exon1"} = "CCCCCGGGGG"; # 10 BASES LONG, SO THE NEXT EXON SHOULD START IN PHASE 2
   $EXONSEQ{"ID=gene1:exon2"} = "CGCGCGCGCG";
   ($output_fasta,$errorcode,$errormsg) = &make_filename($outputdir);
   if ($errorcode != 0) { ($errorcode,$errormsg) = &print_error($errormsg,$errorcode,0); } 
   @temp                  = split(/\//,$output_fasta);
   $output_fasta          = $temp[$#temp];
   # CALL &infer_spliced_dna:
   ($errorcode,$errormsg) = &infer_spliced_dna($exon_gff,\%EXONS,\%EXON_START,\%EXON_FRAME,\%EXONSEQ,$output_fasta,$outputdir,'no','no');
   if ($errorcode != 29) { print STDERR "ERROR: test_infer_spliced_dna: failed test2 (errorcode $errorcode errormsg $errormsg)\n"; exit;}
   system "rm -f $exon_gff";
   system "rm -f $output_fasta";
 
   # TEST ERRORCODE=23 (STRAND IS NOT + OR -):
   ($exon_gff,$errorcode,$errormsg) = &make_filename($outputdir);
   if ($errorcode != 0) { ($errorcode,$errormsg) = &print_error($errormsg,$errorcode,0); } 
   open(EXON_GFF,">$exon_gff") || die "ERROR: test_get_exon_sequences: cannot open $exon_gff\n";
   print EXON_GFF "chr1\tchado\tCDS\t11\t20\t.\t*\t0\tID=gene1:exon1;Parent=gene1;\n";
   print EXON_GFF "chr1\tchado\tCDS\t31\t40\t.\t*\t2\tID=gene1:exon2;Parent=gene1;\n";
   close(EXON_GFF);
   $EXONS{"gene1"} = "ID=gene1:exon1*,*ID=gene1:exon2";
   $EXON_START{"ID=gene1:exon1"} = 11;
   $EXON_START{"ID=gene1:exon2"} = 31;
   $EXON_FRAME{"ID=gene1:exon1"} = "zero";
   $EXON_FRAME{"ID=gene1:exon2"} = 2;
   $EXONSEQ{"ID=gene1:exon1"} = "CCCCCGGGGG"; # 10 BASES LONG, SO THE NEXT EXON SHOULD START IN PHASE 2
   $EXONSEQ{"ID=gene1:exon2"} = "CGCGCGCGCG";
   ($output_fasta,$errorcode,$errormsg) = &make_filename($outputdir);
   if ($errorcode != 0) { ($errorcode,$errormsg) = &print_error($errormsg,$errorcode,0); } 
   @temp                  = split(/\//,$output_fasta);
   $output_fasta          = $temp[$#temp];
   # CALL &infer_spliced_dna:
   ($errorcode,$errormsg) = &infer_spliced_dna($exon_gff,\%EXONS,\%EXON_START,\%EXON_FRAME,\%EXONSEQ,$output_fasta,$outputdir,'no','no');
   if ($errorcode != 23) { print STDERR "ERROR: test_infer_spliced_dna: failed test3 (errorcode $errorcode errormsg $errormsg)\n"; exit;}
   system "rm -f $exon_gff";
   system "rm -f $output_fasta";
 
   # TEST ERRORCODE=24/27/28 (DO NOT HAVE SEQUENCE OF ALL EXONS IN THE GFF):
   ($exon_gff,$errorcode,$errormsg) = &make_filename($outputdir);
   if ($errorcode != 0) { ($errorcode,$errormsg) = &print_error($errormsg,$errorcode,0); } 
   open(EXON_GFF,">$exon_gff") || die "ERROR: test_get_exon_sequences: cannot open $exon_gff\n";
   print EXON_GFF "chr1\tchado\tCDS\t11\t20\t.\t+\t0\tID=gene1:exon1;Parent=gene1;\n";
   print EXON_GFF "chr1\tchado\tCDS\t31\t40\t.\t+\t2\tID=gene1:exon2;Parent=gene1;\n";
   close(EXON_GFF);
   $EXONS{"gene1"} = "ID=gene1:exon1*,*ID=gene1:exon2";
   $EXON_START{"ID=gene1:exon1"} = 11;
   $EXON_START{"ID=gene1:exon2"} = 31;
   $EXON_FRAME{"ID=gene1:exon1"} = "zero";
   $EXON_FRAME{"ID=gene1:exon2"} = 2;
   %EXONSEQ                   = ();
   ($output_fasta,$errorcode,$errormsg) = &make_filename($outputdir);
   if ($errorcode != 0) { ($errorcode,$errormsg) = &print_error($errormsg,$errorcode,0); } 
   @temp                  = split(/\//,$output_fasta);
   $output_fasta          = $temp[$#temp];
   # CALL &infer_spliced_dna:
   ($errorcode,$errormsg) = &infer_spliced_dna($exon_gff,\%EXONS,\%EXON_START,\%EXON_FRAME,\%EXONSEQ,$output_fasta,$outputdir,'no','no');
   if ($errorcode != 24 && $errorcode != 27 && $errorcode != 28) { print STDERR "ERROR: test_infer_spliced_dna: failed test4 (errorcode $errorcode errormsg $errormsg)\n"; exit;}
   system "rm -f $exon_gff";
   system "rm -f $output_fasta";
 
   # TEST ERRORCODE=25/26 (DO NOT KNOW PHASE OF AN EXON):
   ($exon_gff,$errorcode,$errormsg) = &make_filename($outputdir);
   if ($errorcode != 0) { ($errorcode,$errormsg) = &print_error($errormsg,$errorcode,0); } 
   open(EXON_GFF,">$exon_gff") || die "ERROR: test_get_exon_sequences: cannot open $exon_gff\n";
   print EXON_GFF "chr1\tchado\tCDS\t11\t20\t.\t+\t0\tID=gene1:exon1;Parent=gene1;\n";
   print EXON_GFF "chr1\tchado\tCDS\t31\t40\t.\t+\t2\tID=gene1:exon2;Parent=gene1;\n";
   close(EXON_GFF);
   $EXONS{"gene1"} = "ID=gene1:exon1*,*ID=gene1:exon2";
   $EXON_START{"ID=gene1:exon1"} = 11;
   $EXON_START{"ID=gene1:exon2"} = 31;
   %EXON_FRAME            = ();
   $EXONSEQ{"ID=gene1:exon1"} = "CCCCCGGGGG"; # 10 BASES LONG, SO THE NEXT EXON SHOULD START IN PHASE 2
   $EXONSEQ{"ID=gene1:exon2"} = "CGCGCGCGCG";
   ($output_fasta,$errorcode,$errormsg) = &make_filename($outputdir);
   if ($errorcode != 0) { ($errorcode,$errormsg) = &print_error($errormsg,$errorcode,0); } 
   @temp                  = split(/\//,$output_fasta);
   $output_fasta          = $temp[$#temp];
   # CALL &infer_spliced_dna:
   ($errorcode,$errormsg) = &infer_spliced_dna($exon_gff,\%EXONS,\%EXON_START,\%EXON_FRAME,\%EXONSEQ,$output_fasta,$outputdir,'no','no');
   if ($errorcode != 25 && $errorcode != 26) { print STDERR "ERROR: test_infer_spliced_dna: failed test5 (errorcode $errorcode errormsg $errormsg)\n"; exit;}
   system "rm -f $exon_gff";
   system "rm -f $output_fasta";
 
   # TEST USING AN AUGUSTUS GFF FILE:
   ($exon_gff,$errorcode,$errormsg) = &make_filename($outputdir);
   if ($errorcode != 0) { ($errorcode,$errormsg) = &print_error($errormsg,$errorcode,0); } 
   open(EXON_GFF,">$exon_gff") || die "ERROR: test_get_exon_sequences: cannot open $exon_gff\n";
   print EXON_GFF "chr1\tAUGUSTUS\tCDS\t11\t20\t0.5\t+\t0\ttranscript_id \"gene1.t1\"; gene_id \"gene1\"\n";
   print EXON_GFF "chr1\tAUGUSTUS\tCDS\t31\t40\t0.53\t+\t2\ttranscript_id \"gene1.t1\"; gene_id \"gene1\"\n";
   close(EXON_GFF);
   $EXONS{"gene1"} = "ID=gene1:exon1*,*ID=gene1:exon2";
   $EXON_START{"ID=gene1:exon1"} = 11;
   $EXON_START{"ID=gene1:exon2"} = 31;
   $EXON_FRAME{"ID=gene1:exon1"} = "zero";
   $EXON_FRAME{"ID=gene1:exon2"} = 2;
   $EXONSEQ{"ID=gene1:exon1"} = "CCCCCGGGGG"; # 10 BASES LONG, SO THE NEXT EXON SHOULD START IN PHASE 2
   $EXONSEQ{"ID=gene1:exon2"} = "CGCGCGCGCG";
   ($output_fasta,$errorcode,$errormsg) = &make_filename($outputdir);
   if ($errorcode != 0) { ($errorcode,$errormsg) = &print_error($errormsg,$errorcode,0); } 
   @temp                  = split(/\//,$output_fasta);
   $output_fasta          = $temp[$#temp];
   # CALL &infer_spliced_dna:
   ($errorcode,$errormsg) = &infer_spliced_dna($exon_gff,\%EXONS,\%EXON_START,\%EXON_FRAME,\%EXONSEQ,$output_fasta,$outputdir,'no','yes');
   if ($errorcode != 0) { print STDERR "ERROR: test_infer_spliced_dna: failed test6 (errorcode $errorcode errormsg $errormsg)\n"; exit;}
   ($expected_output_fasta,$errorcode,$errormsg) = &make_filename($outputdir);
   if ($errorcode != 0) { ($errorcode,$errormsg) = &print_error($errormsg,$errorcode,0); } 
   open(EXPECTED,">$expected_output_fasta") || die "ERROR: test_infer_spliced_dna: cannot open $expected_output_fasta\n";
   print EXPECTED ">gene1\n";
   print EXPECTED "CCCCCGGGGGCGCGCGCGCG\n";
   close(EXPECTED);
   $differences            = "";
   $output_fasta           = $outputdir."/".$output_fasta;
   open(TEMP,"diff $output_fasta $expected_output_fasta |");
   while(<TEMP>)
   {
      $line                = $_;
      $differences         = $differences.$line;
   }
   close(TEMP);  
   $length_differences     = length($differences);
   if ($length_differences != 0) { print STDERR "ERROR: test_infer_spliced_dna: failed test6 (output_fasta $output_fasta expected_output_fasta $expected_output_fasta)\n"; exit;}
   system "rm -f $exon_gff";
   system "rm -f $output_fasta";
   system "rm -f $expected_output_fasta";
 
}
 
#------------------------------------------------------------------#
 
# INFER THE SPLICED DNA SEQUENCE FOR EACH GENE:
 
sub infer_spliced_dna
{
   my $exon_gff            = $_[0]; # EXON GFF
   my $EXONS               = $_[1]; # HASH TABLE WITH EXONS IN EACH GENE
   my $EXON_START          = $_[2]; # HASH TABLE WITH START POINTS OF EXONS
   my $EXON_FRAME          = $_[3]; # HASH TABLE WITH FRAME OF EXONS
   my $EXONSEQ             = $_[4]; # HASH TABLE WITH DNA SEQUENCES OF EXONS
   my $output_fasta        = $_[5]; # OUTPUT FASTA FILE FOR SPLICED SEQUENCES
   my $outputdir           = $_[6]; # DIRECTORY TO WRITE OUTPUT FILES INTO
   my $ignore_phase        = $_[7]; # SAYS WHETHER TO IGNORE PHASE INFORMATION IN THE GFF, AND ASSUME THE FIRST EXON IS PHASE 0 
   my $from_augustus       = $_[8]; # SAYS WHETHER THE GFF FILE IS FROM AUGUSTUS 
   my %FOUND_SPLICED_DNA   = ();    # HASH TABLE TO RECORD WHETHER WE HAVE ALREADY FOUND THE SPLICED 
                                    # DNA SEQUENCE FOR A GENE
   my $errorcode           = 0;     # RETURNED AS 0 IF THERE IS NO ERROR
   my $errormsg            = "none";# RETURNED AS 'none' IF THERE IS NO ERROR
   my $line;                        # 
   my @temp;                        # 
   my $strand;                      # STRAND OF A GENE
   my $exon;                        # EXON NAME
   my $gene;                        # GENE NAME 
   my $exons;                       # EXONS IN A GENE
   my @exons;                       # EXONS IN A GENE 
   my $no_exons;                    # NUMBER OF EXONS IN A GENE 
   my $i;                           # 
   my $expected_phase;              # EXPECTED PHASE OF AN EXON 
   my $exon_phase;                  # PHASE OF AN EXON
   my $exon_seq;                    # SEQUENCE OF AN EXON
   my $exon_seq_length;             # LENGTH OF AN EXON'S SEQUENCE 
   my $prev_exon;                   # THE PREVIOUS EXON
   my $prev_exon_seq;               # SEQUENCE OF THE PREVIOUS EXON 
   my $prev_exon_seq_length;        # LENGTH OF $prev_exon_seq 
   my $spliced_dna;                 # THE SPLICED DNA SEQUENCE FOR A GENE 
   my $length;                      # LENGTH OF $spliced_dna
   my $offset;                      # COUNTER USED FOR PRINTING OUT THE SEQUENCE
   my @genes;                       # TRANSCRIPTS IN A GENE (THERE CAN BE MULTIPLE TRANSCRIPTS PER GENE, eg. FOR C. ELEGANS)
   my $k;                           #  
 
   # OPEN THE OUTPUT FASTA FILE:
   $output_fasta           = $outputdir."/".$output_fasta;
   open(OUTPUT,">$output_fasta") || die "ERROR: infer_spliced_dna: cannot open output_fasta $output_fasta\n";
 
   # READ THROUGH THE EXON GFF FILE:
   open(EXON_GFF,"$exon_gff") || die "ERROR: infer_spliced_dna: cannot open $exon_gff.\n";
   while(<EXON_GFF>)
   {
      $line                = $_;
      chomp $line;
      @temp                = split(/\t+/,$line);
      $strand              = $temp[6];
      $exon                = $temp[8];
      # FIND THE NAME OF THE GENE:
      if ($from_augustus eq 'no')
      {
         @temp             = split(/Parent=/,$exon); # eg. ID=0:SSTP.scaffold.00025.352065-TF352220-GW-3:cds;Parent=0:SSTP.scaffold.00025.352065-TF352220-GW-3
         $gene             = $temp[1];               # eg. 0:SSTP.scaffold.00025.352065-TF352220-GW-3 
         @temp             = split(/\;/,$gene);
         $gene             = $temp[0];
         @genes            = split(/\,/,$gene);
      }
      elsif ($from_augustus eq 'yes') # GFF IS FROM AUGUSTUS
      {
         # eg. transcript_id "g1.t1"; gene_id "g1"; 
         @temp             = split(/gene_id/,$exon);
         $gene             = $temp[1]; # eg.  "g1"; 
         @temp             = split(/\"/,$gene);
         $gene             = $temp[1]; # eg. g1
         @genes            = split(/\,/,$gene); 
      }
      for ($k = 0; $k <= $#genes; $k++) # THERE CAN BE MULTIPLE TRANSCRIPTS PER GENE, eg. FOR C. ELEGANS
      {
         $gene             = $genes[$k]; 
         # IF WE HAVEN'T ALREADY FOUND THE SPLICED DNA FOR THIS GENE:
         if (!($FOUND_SPLICED_DNA{$gene}))
         {
            if (!($EXONS->{$gene}))
            {
               system "rm -f $output_fasta";
               system "rm -f $exon_gff";
               $errormsg   = "ERROR: infer_spliced_dna: do not know exons in gene $gene\n"; 
               $errorcode  = 29; # ERRORCODE=29 (TESTED FOR)
               return($errorcode,$errormsg);
            }
            $exons         = $EXONS->{$gene};
            @exons         = split(/\*\,\*/,$exons);
            # SORT THE EXONS BY START-POINT IN THE DNA:
            @exons         = sort { $EXON_START->{$a} <=> $EXON_START->{$b} } @exons;
            # GO THROUGH THE EXONS, AND CHECK THAT THEY ALL HAVE THE SAME PHASE AS EXPECTED; IF
            # NOT SOME DNA MAY HAVE TO BE REMOVED FROM THE EXON'S ENDS:
            $no_exons      = $#exons + 1;
            for ($i = 0; $i < $no_exons; $i++)
            { 
               # GO THROUGH THE EXONS IN ORDER IN THE PROTEIN:
               if    ($strand eq '+') { $exon = $exons[$i];        }
               elsif ($strand eq '-') { $exon = $exons[$#exons-$i];}
               else
               {
                  system "rm -f $output_fasta";
                  system "rm -f $exon_gff";
                  $errormsg= "ERROR: infer_spliced_dna: strand $strand for gene $gene\n";
                  $errorcode = 23; # ERRORCODE=23 (TESTED FOR)
                  return($errorcode,$errormsg);
               }   
               # GET THE DNA SEQUENCE OF THE EXON:
               if (!($EXONSEQ->{$exon}))
               {
                  system "rm -f $output_fasta";
                  system "rm -f $exon_gff";
                  $errormsg= "ERROR: infer_spliced_dna: do not know DNA sequence for exon $exon\n";
                  $errorcode = 24; # ERRORCODE=24 (TESTED FOR)
                  return($errorcode,$errormsg);
               }
               $exon_seq   = $EXONSEQ->{$exon};
               if ($exon_seq eq "none") { $exon_seq = "";}
               # IF IT IS THE FIRST EXON, CHECK WHETHER THE PHASE IS ZERO:
               if ($i == 0)
               {
                  if (!($EXON_FRAME->{$exon}))
                  {
                     system "rm -f $output_fasta";
                     system "rm -f $exon_gff";
                     $errormsg= "ERROR: infer_spliced_dna: do not know frame for exon $exon\n";
                     $errorcode = 25; # ERRORCODE=25 (TESTED FOR)
                     return($errorcode,$errormsg);
                  }
                  $exon_phase = $EXON_FRAME->{$exon};
                  if ($exon_phase eq "zero") { $exon_phase = 0;}
                  if ($ignore_phase eq 'yes') { $exon_phase = 0;} # ASSUME THE FIRST EXON IS PHASE 0 IF 'ignore_phase' IS 'yes'
                  if ($exon_phase != 0) # IF THE PHASE OF THE FIRST EXON IN THE PROTEIN IS NOT ZERO.
                  {
                     ($exon_seq,$errorcode,$errormsg) = &correct_first_exon_seq($exon_seq,$strand,$exon_phase,$EXONSEQ,$exon); 
                     if ($errorcode != 0) { ($errorcode,$errormsg) = &print_error($errormsg,$errorcode,0); }
                  }
               }   
               else # IF IT IS NOT THE FIRST EXON, CHECK THAT THE START PHASE OF THIS EXON IS THE END PHASE OF THE PREVIOUS EXON:
               {
                  if (!($EXON_FRAME->{$exon})) 
                  {
                     system "rm -f $output_fasta";
                     system "rm -f $exon_gff";
                     $errormsg= "ERROR: infer_spliced_dna: do not know expected phase of exon $exon\n";
                     $errorcode = 26; # ERRORCODE=26 (TESTED FOR)
                     return($errorcode,$errormsg);
                  }
                  $exon_phase = $EXON_FRAME->{$exon};
                  if ($exon_phase eq "zero") { $exon_phase = 0;}
 
                  # CALCULATE THE EXPECTED FRAME OF THIS EXON:
                  ($expected_phase,$errorcode,$errormsg) = &calculate_expected_phase($i,$strand,\@exons,$EXONSEQ); 
                  if ($errorcode != 0) { ($errorcode,$errormsg) = &print_error($errormsg,$errorcode,0); }
                  if ($ignore_phase eq 'yes') { $exon_phase = $expected_phase;} # IF 'ignore_phase' IS 'yes', SET THE PHASE TO BE THE EXPECTED PHASE.
 
                  # CHECK IF THE EXON PHASE IS EQUAL TO THE EXPECTED EXON PHASE:
                  if ($exon_phase != $expected_phase)
                  {
                     # ADD A BIT TO THE START OF THIS EXON TO MAKE IT STARTS IN PHASE 0:
                     ($exon_seq,$errorcode,$errormsg) = &correct_start_exon_seq($exon_seq,$strand,$exon_phase,$EXONSEQ,$exon); 
                     if ($errorcode != 0) { ($errorcode,$errormsg) = &print_error($errormsg,$errorcode,0); }
 
                     # ADD A BIT ONTO THE END OF THE LAST EXON TO MAKE IT END IN PHASE 0:
                     if     ($strand eq '+') { $prev_exon = $exons[$i-1];        }
                     elsif  ($strand eq '-') { $prev_exon = $exons[$#exons-$i+1];}
                     if (!($EXONSEQ->{$prev_exon})) 
                     { 
                        system "rm -f $output_fasta";
                        system "rm -f $exon_gff";
                        $errormsg = "ERROR: infer_spliced_dna: do not know the sequence of prev_exon $prev_exon\n";
                        $errorcode = 27; # ERRORCODE=27 (TESTED FOR)
                        return($errorcode,$errormsg);
                     } 
                     $prev_exon_seq = $EXONSEQ->{$prev_exon};
                     ($prev_exon_seq,$errorcode,$errormsg) = &correct_end_prev_exon_seq($prev_exon_seq,$strand,$expected_phase,$EXONSEQ,$prev_exon);
                     if ($errorcode != 0) { ($errorcode,$errormsg) = &print_error($errormsg,$errorcode,0); }
 
                  } # END OF if ($exon_phase != $expected_phase)
               } # END OF IF IT IS NOT THE FIRST EXON
            } # END OF LOOP THROUGH ALL THE EXONS.
 
            # SPLICE TOGETHER THE DNA FOR ALL THE EXONS:
            $no_exons      = $#exons + 1;
            $spliced_dna   = "";
            for ($i = 0; $i < $no_exons; $i++)
            {
               $exon       = $exons[$i];
               if (!($EXONSEQ->{$exon}))
               {
                  system "rm -f $output_fasta";
                  system "rm -f $exon_gff";
                  $errormsg = "ERROR: infer_spliced_dna: do not know sequence of exon $exon\n";
                  $errorcode = 28; # ERRORCODE=28 (TESTED FOR)
                  return($errorcode,$errormsg);
               }
               $exon_seq   = $EXONSEQ->{$exon};
               if ($exon_seq eq "none") { $exon_seq = "";}
               $spliced_dna= $spliced_dna.$exon_seq;
            }
            $spliced_dna   =~ tr/[a-z]/[A-Z]/;
 
            # GET THE REVERSE COMPLEMENT:
            if ($strand eq '-')
            {
               ($spliced_dna,$errorcode,$errormsg) = &reverse_complement($spliced_dna); 
               if ($errorcode != 0) { ($errorcode,$errormsg) = &print_error($errormsg,$errorcode,0); }
            }
 
            # PRINT OUT TO THE OUTPUT FILE:
            print OUTPUT ">$gene\n";
            $length        = length($spliced_dna);
            $offset        = 0;
            while ($offset < $length)
            {
               $line       = substr($spliced_dna,$offset,60);
               print OUTPUT "$line\n";
               $offset     = $offset + 60;
            }    
            $FOUND_SPLICED_DNA{$gene} = 1;
         } # END OF IF WE DON'T HAVE THE SPLICED DNA FOR THE GENE ALREADY   
      }
   }
   close(GFF);
   close(OUTPUT);
 
   return($errorcode,$errormsg); 
}
 
#------------------------------------------------------------------#
 
# TEST &calculate_expected_phase
 
sub test_calculate_expected_phase
{
   my $exon_phase;                  # PHASE OF THE EXON
   my $errorcode;                   # RETURNED BY 0 AS A FUNCTION IF THERE IS NO ERROR
   my $errormsg;                    # RETURNED BY 'none' AS A FUNCTION IF THERE IS NO ERROR
   my %EXONSEQ             = ();    # HASH TABLE WITH THE SEQUENCES OF EXONS
   my @exons;                       # ARRAY OF EXONS
   
   @exons                  = ("exon1","exon2","exon3");
   $EXONSEQ{"exon1"}       = "ATCGATCG";
   $EXONSEQ{"exon2"}       = "ATATATAT";
   $EXONSEQ{"exon3"}       = "CGACGACG";   
   ($exon_phase,$errorcode,$errormsg) = &calculate_expected_phase(0,"+",\@exons,\%EXONSEQ);
   if ($errorcode != 0 || $exon_phase != 0) { print STDERR "ERROR: test_calculate_expected_exon_phase: failed test1\n"; exit;}
   
   @exons                  = ("exon1","exon2","exon3");
   $EXONSEQ{"exon1"}       = "ATCGATCG"; # 2nd EXON STARTS ON 3rd BASE OF A CODON (PHASE 1)
   $EXONSEQ{"exon2"}       = "ATATATAT";
   $EXONSEQ{"exon3"}       = "CGACGACG";   
   ($exon_phase,$errorcode,$errormsg) = &calculate_expected_phase(1,"+",\@exons,\%EXONSEQ);
   if ($errorcode != 0 || $exon_phase != 1) { print STDERR "ERROR: test_calculate_expected_exon_phase: failed test2 (errorcode $errorcode errormsg $errormsg exon_phase $exon_phase)\n"; exit;}
   
   @exons                  = ("exon1","exon2","exon3");
   $EXONSEQ{"exon1"}       = "ATCGATCG"; # 2nd EXON STARTS ON 3rd BASE OF A CODON (PHASE 1)
   $EXONSEQ{"exon2"}       = "ATATATAT"; # 3rd EXON STARTS ON 2nd BASE OF A CODON (PHASE 2)
   $EXONSEQ{"exon3"}       = "CGACGACG";   
   ($exon_phase,$errorcode,$errormsg) = &calculate_expected_phase(2,"+",\@exons,\%EXONSEQ);
   if ($errorcode != 0 || $exon_phase != 2) { print STDERR "ERROR: test_calculate_expected_exon_phase: failed test3 (errorcode $errorcode errormsg $errormsg exon_phase $exon_phase)\n"; exit;}
 
   # TEST FOR ERRORCODE=30 (DON'T KNOW AN EXON'S SEQUENCE):
   @exons                  = ("exon1","exon2","exon3");
   %EXONSEQ                = ();
   ($exon_phase,$errorcode,$errormsg) = &calculate_expected_phase(1,"+",\@exons,\%EXONSEQ);
   if ($errorcode != 30) { print STDERR "ERROR: test_calculate_expected_exon_phase: failed test4 (errorcode $errorcode errormsg $errormsg exon_phase $exon_phase)\n"; exit;}
 
}
 
#------------------------------------------------------------------#
 
# SUBROUTINE TO CALCULATE THE EXPECTED PHASE OF EXON $i IN A GENE:
 
sub calculate_expected_phase
{
   my $exon_index          = $_[0]; # INDEX OF THE EXON TO CALCULATE THE PHASE OF.
   my $strand              = $_[1]; # STRAND OF THE EXON
   my $exons               = $_[2]; # ARRAY OF EXONS IN THE GENE 
   my $EXONSEQ             = $_[3]; # HASH TABLE WITH THE SEQUENCES OF EXONS 
   my $no_exons;                    # NUMBER OF EXONS IN THE GENE
   my $i;                           # 
   my $exon;                        # EXON NAME
   my $exon_seq;                    # SEQUENCE OF EXON
   my $length_of_previous_exons;    # LENGTH OF PREVIOUS EXONS BEFORE THIS ONE IN THE GENE
   my $exon_length;                 # LENGTH OF EXON
   my $first_exon;                  # FIRST EXON IN THE GENE
   my $phase_first_exon;            # PHASE OF THE FIRST EXON IN THE GENE
   my $exon_phase;                  # PHASE OF THE EXON
   my $errorcode           = 0;     # RETURNED AS 0 IF THERE IS NO ERROR
   my $errormsg            = "none";# RETURNED AS 'none' IF THERE IS NO ERROR 
 
   # FIND THE SUMMED LENGTHS OF ALL THE EXONS PREVIOUS TO THE EXON
   # $exon_index IN THE PROTEIN:
   if ($strand eq '-') { $exon_index = (@$exons - 1) - $exon_index;}
   $length_of_previous_exons= 0;
   $no_exons               = @$exons;      
   if    ($strand eq '+')
   {
      for ($i = 0; $i < $exon_index; $i++)
      {
         $exon             = $exons->[$i];
         if (!($EXONSEQ->{$exon})) 
         {
            $errormsg      = "ERROR: calculate_expected_phase : do not know the sequence of $exon.\n"; 
            $errorcode     = 30; # ERRORCODE=30 (TESTED FOR)
            return($exon_phase,$errorcode,$errormsg);
         }
         $exon_seq         = $EXONSEQ->{$exon};
         if ($exon_seq eq "none") { $exon_seq = "";}
         $exon_length      = length($exon_seq);
         $length_of_previous_exons = $length_of_previous_exons + $exon_length;
      }
   }
   elsif ($strand eq '-')
   {
      for ($i = $no_exons-1; $i > $exon_index; $i--)
      {
         $exon             = $exons->[$i];
         if (!($EXONSEQ->{$exon})) 
         {
            $errormsg      = "ERROR: calculate_expected_phase : do not know the sequence of $exon.\n"; 
            $errorcode     = 30; # ERRORCODE=30 (TESTED FOR)
            return($exon_phase,$errorcode,$errormsg);
         } 
         $exon_seq         = $EXONSEQ->{$exon};
         if ($exon_seq eq "none") { $exon_seq = "";}
         $exon_length      = length($exon_seq);
         $length_of_previous_exons = $length_of_previous_exons + $exon_length;
      }
   }
 
   # CALCULATE THE EXPECTED PHASE OF THE EXON:
   # FOR EXAMPLE, IF THE LENGTH OF THE PREVIOUS EXONS IS 4, $exon_phase = 2:  XXX|X=====XX|XXX...
   #              IF THE LENGTH OF THE PREVIOUS EXONS IS 5, $exon_phase = 1:  XXX|XX=====X|XXX...
   $exon_phase             = $length_of_previous_exons % 3;
   if    ($exon_phase == 2) { $exon_phase = 1;}
   elsif ($exon_phase == 1) { $exon_phase = 2;}
   # ELSE IF $exon_phase == 0, THEN $exon_phase = 0
 
   return ($exon_phase,$errorcode,$errormsg);
}
 
#------------------------------------------------------------------#
 
# TEST &correct_first_exon_seq
 
sub test_correct_first_exon_seq
{
   my $errorcode;                   # RETURNED AS 0 BY A FUNCTION IF THERE IS NO ERROR
   my $errormsg;                    # RETURNED AS 'none' BY A FUNCTION IF THERE IS NO ERROR
   my $exon_seq;                    # SEQUENCE OF AN EXON  
   my %EXONSEQ             = ();    # HASH TABLE WITH SEQUENCES OF EXONS
 
   ($exon_seq,$errorcode,$errormsg) = &correct_first_exon_seq("AATGCAG","+",1,\%EXONSEQ,"exon1");
   if ($exon_seq ne 'ATGCAG' || $errorcode != 0 || $EXONSEQ{"exon1"} ne 'ATGCAG') { print STDERR "ERROR: test_correct_first_exon_seq: failed test1\n"; exit;}
  
   %EXONSEQ                = ();
   ($exon_seq,$errorcode,$errormsg) = &correct_first_exon_seq("ATGCAG","+",0,\%EXONSEQ,"exon1");
   if ($exon_seq ne 'ATGCAG' || $errorcode != 0 || $EXONSEQ{"exon1"} ne 'ATGCAG') { print STDERR "ERROR: test_correct_first_exon_seq: failed test2\n"; exit;}
  
   %EXONSEQ                = ();
   ($exon_seq,$errorcode,$errormsg) = &correct_first_exon_seq("AAATGCAG","+",2,\%EXONSEQ,"exon1");
   if ($exon_seq ne 'ATGCAG' || $errorcode != 0 || $EXONSEQ{"exon1"} ne 'ATGCAG') { print STDERR "ERROR: test_correct_first_exon_seq: failed test3\n"; exit;}
 
}
 
#------------------------------------------------------------------#
 
# IF THE PHASE OF THE FIRST EXON IN THE PROTEIN IS NOT ZERO.
 
sub correct_first_exon_seq
{
   my $exon_seq            = $_[0]; # SEQUENCE OF THE FIRST EXON 
   my $strand              = $_[1]; # STRAND OF THE FIRST EXON
   my $exon_phase          = $_[2]; # PHASE OF THE FIRST EXON
   my $EXONSEQ             = $_[3]; # HASH TABLE WITH THE SEQUENCES OF EXONS
   my $exon                = $_[4]; # EXON NAME
   my $exon_seq_length;             # LENGTH OF AN EXON'S SEQUENCE
   my $errorcode           = 0;     # RETURNED AS 0 IF THERE IS NO ERROR
   my $errormsg            = "none";# RETURNED AS 'none' IF THERE IS NO ERROR
 
   if    ($exon_phase == 1 && $strand eq '+')
   {
      # IGNORE THE FIRST BASE:
      $exon_seq            = substr($exon_seq,1,length($exon_seq)-1);
   }
   elsif ($exon_phase == 1 && $strand eq '-')
   {
      # IGNORE THE LAST BASE:
      $exon_seq            = substr($exon_seq,0,length($exon_seq)-1);
   }
   elsif ($exon_phase == 2 && $strand eq '+')
   {
      # IGNORE THE FIRST TWO BASES:
      $exon_seq            = substr($exon_seq,2,length($exon_seq)-2);
   }
   elsif ($exon_phase == 2 && $strand eq '-')
   {
      # IGNORE THE LAST TWO BASES:
      $exon_seq            = substr($exon_seq,0,length($exon_seq)-2);
   }
 
   # REPLACE THE SEQUENCE OF THE FIRST EXON:
   $exon_seq_length        = length($exon_seq);
   if ($exon_seq_length == 0) { $EXONSEQ->{$exon} = "none";   }
   else                       { $EXONSEQ->{$exon} = $exon_seq;}
 
   return($exon_seq,$errorcode,$errormsg);
}
 
#------------------------------------------------------------------#
 
# TEST &correct_end_prev_exon_seq
 
sub test_correct_end_prev_exon_seq
{
   my $errorcode;                   # RETURNED AS 0 BY A FUNCTION IF THERE IS NO ERROR
   my $errormsg;                    # RETURNED AS 'none' BY A FUNCTION IF THERE IS NO ERROR
   my $prev_exon_seq;               # SEQUENCE OF THE PREVIOUS EXON
   my %EXONSEQ             = ();    # HASH TABLE TO STORE THE SEQUENCES OF EXONS
 
   ($prev_exon_seq,$errorcode,$errormsg) = &correct_end_prev_exon_seq("ATGCAG","+",0,\%EXONSEQ,"exon1");
   if ($errorcode != 0) { print STDERR "ERROR: test_correct_end_prev_exon_seq: failed test1\n"; exit;}
   if ($prev_exon_seq ne 'ATGCAGNNN') { print STDERR "ERROR: test_correct_end_prev_exon_seq: failed test1b (prev_exon_seq $prev_exon_seq)\n"; exit;}
   if ($EXONSEQ{"exon1"} ne 'ATGCAGNNN') { print STDERR "ERROR: test_correct_end_prev_exon_seq: failed test1c\n"; exit;}
 
   ($prev_exon_seq,$errorcode,$errormsg) = &correct_end_prev_exon_seq("ATGCA","+",1,\%EXONSEQ,"exon1"); 
   if ($errorcode != 0) { print STDERR "ERROR: test_correct_end_prev_exon_seq: failed test2\n"; exit;}
   if ($prev_exon_seq ne 'ATGCAN') { print STDERR "ERROR: test_correct_prev_exon_seq: failed test2b (prev_exon_seq $prev_exon_seq)\n"; exit;}
   if ($EXONSEQ{"exon1"} ne 'ATGCAN') { print STDERR "ERROR: test_correct_end_prev_exon_seq: failed test2c\n"; exit;}
 
   ($prev_exon_seq,$errorcode,$errormsg) = &correct_end_prev_exon_seq("ATGC","+",2,\%EXONSEQ,"exon1");
   if ($errorcode != 0) { print STDERR "ERROR: test_correct_end_prev_exon_seq: failed test3\n"; exit;}
   if ($prev_exon_seq ne 'ATGCNN') { print STDERR "ERROR: test_correct_prev_exon_seq: failed test3b (prev_exon_seq $prev_exon_seq)\n"; exit;}
   if ($EXONSEQ{"exon1"} ne 'ATGCNN') { print STDERR "ERROR: test_correct_end_prev_exon_seq: failed test3c\n"; exit;}
   
}
 
#------------------------------------------------------------------#
 
# ADD A BIT ONTO THE END OF THE LAST EXON TO MAKE IT END IN PHASE 0:
 
sub correct_end_prev_exon_seq
{
   my $prev_exon_seq       = $_[0]; # SEQUENCE OF THE PREVIOUS EXON
   my $strand              = $_[1]; # STRAND
   my $expected_phase      = $_[2]; # EXPECTED PHASE OF THE PREVIOUS EXON
   my $EXONSEQ             = $_[3]; # HASH TABLE WITH THE SEQUENCES OF EXONS
   my $prev_exon           = $_[4]; # NAME OF THE PREVIOUS EXON
   my $prev_exon_seq_length;        # LENGTH OF THE PREVIOUS EXON'S SEQUENCE 
   my $errorcode           = 0;     # RETURNED AS 0 IF THERE IS NO ERROR
   my $errormsg            = "none";# RETURNED AS 'none' IF THERE IS NO ERROR 
 
   if ($prev_exon_seq ne "none")
   {
      # ADD SOME BASES ONTO THE END OF THE PREVIOUS EXON SO IT ENDS IN PHASE 0:
      if    ($expected_phase == 1 && $strand eq '+')
      {
         # WE HAVE ...123|123|12---intron---3|123|123...
         # ADD ANOTHER BASE TO THE END OF THE PREVIOUS EXON:
         $prev_exon_seq    = $prev_exon_seq."N";
      }
      elsif ($expected_phase == 2 && $strand eq '+')
      {
         # WE HAVE ...123|123|1---intron--23|123|123...
         # ADD TWO BASES TO THE END OF THE PREVIOUS EXON:
         $prev_exon_seq    = $prev_exon_seq."NN";
      }
      elsif ($expected_phase == 1 && $strand eq '-')
      {
         # WE HAVE ...321|321|321|3---intron---21|321|321...
         # ADD ANOTHER BASE TO THE START OF THE PREVIOUS EXON:
         $prev_exon_seq    = "N".$prev_exon_seq;
      }
      elsif ($expected_phase == 2 && $strand eq '-')
      {
         # WE HAVE ...321|321|321|32---intron---1|321|321...
         # ADD ANOTHER TWO BASE TO THE START OF THE PREVIOUS EXON:
         $prev_exon_seq    = "NN".$prev_exon_seq;
      }  
      elsif ($expected_phase == 0 && $strand eq '+')
      {
         # ADD ANOTHER THREE BASES TO THE END OF THE PREVIOUS EXON:
         $prev_exon_seq    = $prev_exon_seq."NNN";
      }
      elsif ($expected_phase == 0 && $strand eq '-')
      {
         # ADD ANOTHER THREE BASES TO THE END OF THE PREVIOUS EXON:
         $prev_exon_seq    = "NNN".$prev_exon_seq;
      } 
      # REPLACE THE SEQUENCE FOR THE PREVIOUS EXON:
      $prev_exon_seq_length= length($prev_exon_seq);
      if ($prev_exon_seq_length == 0) { $EXONSEQ->{$prev_exon} = "none";        }
      else                            { $EXONSEQ->{$prev_exon} = $prev_exon_seq;}
   }
 
   return($prev_exon_seq,$errorcode,$errormsg);
}
 
#------------------------------------------------------------------#
 
# TEST &correct_start_exon_seq
 
sub test_correct_start_exon_seq
{
   my $errorcode;                   # RETURNED AS 0 FROM A FUNCTION IF THERE IS NO ERROR
   my $errormsg;                    # RETURNED AS 'none' FROM A FUNCTION IF THERE IS NO ERROR
   my $exonseq;                     # SEQUENCE OF THE EXON
   my %EXONSEQ             = ();    # HASH TABLE TO STORE EXON SEQUENCES IN
   
   ($exonseq,$errorcode,$errormsg) = &correct_start_exon_seq("ATGCAG","+",0,\%EXONSEQ,"exon1");
   if ($errorcode != 0) { print STDERR "ERROR: test_correct_start_exon_seq: failed test1\n"; exit;}
   if ($exonseq ne 'NNNATGCAG') { print STDERR "ERROR: test_correct_start_exon_seq: failed test1b\n"; exit;}
   if ($EXONSEQ{"exon1"} ne "NNNATGCAG") { print STDERR "ERROR: test_correct_start_exon_seq: failed test1c\n"; exit;}
 
   ($exonseq,$errorcode,$errormsg) = &correct_start_exon_seq("GCAG","+",1,\%EXONSEQ,"exon1");
   if ($errorcode != 0) { print STDERR "ERROR: test_correct_start_exon_seq: failed test2\n"; exit;}
   if ($exonseq ne 'NNGCAG') { print STDERR "ERROR: test_correct_start_exon_seq: failed test2b\n"; exit;} 
   if ($EXONSEQ{"exon1"} ne "NNGCAG") { print STDERR "ERROR: test_correct_start_exon_seq: failed test2c\n"; exit;}
 
   ($exonseq,$errorcode,$errormsg) = &correct_start_exon_seq("TGCAG","+",2,\%EXONSEQ,"exon1");
   if ($errorcode != 0) { print STDERR "ERROR: test_correct_start_exon_seq: failed test3\n"; exit;}
   if ($exonseq ne "NTGCAG") { print STDERR "ERROR: test_correct_start_exon_seq: failed test3b\n"; exit;}
   if ($EXONSEQ{"exon1"} ne "NTGCAG") { print STDERR "ERROR: test_correct_start_exon_seq: failed test3c\n"; exit;}
}
 
#------------------------------------------------------------------#
 
# ADD A BIT TO THE START OF THIS EXON TO MAKE IT STARTS IN PHASE 0:
 
sub correct_start_exon_seq
{
   my $exon_seq            = $_[0]; # SEQUENCE OF THE EXON
   my $strand              = $_[1]; # STRAND OF THE EXON
   my $exon_phase          = $_[2]; # PHASE OF THE EXON
   my $EXONSEQ             = $_[3]; # HASH TABLE CONTAINING EXON SEQUENCES
   my $exon                = $_[4]; # NAME OF THE EXON
   my $exon_seq_length;             # LENGTH OF THE EXON'S SEQUENCE 
   my $errorcode           = 0;     # RETURNED AS 0 IF THERE IS NO ERROR
   my $errormsg            = "none";# RETURNED AS 'none' IF THERE IS NO ERROR 
 
   # ADD A BIT TO THE START OF THIS EXON TO MAKE IT STARTS IN PHASE 0:
   if    ($exon_phase == 1 && $strand eq '+')
   {
      # WE HAVE : 3|123|123|...
      # ADD TWO BASES BEFORE THE FIRST BASE:
      $exon_seq            = "NN".$exon_seq;
   }
   elsif ($exon_phase == 1 && $strand eq '-')
   {
      # WE HAVE : ...|321|321|3
      # ADD TWO BASES AFTER THE LAST BASE:
      $exon_seq            = $exon_seq."NN";
   }
   elsif ($exon_phase == 2 && $strand eq '+')
   {
      # WE HAVE : 23|123|123|...
      # ADD ONE BASE BEFORE THE FIRST BASE:
      $exon_seq            = "N".$exon_seq;
   }
   elsif ($exon_phase == 2 && $strand eq '-')
   {
      # WE HAVE : ...|321|321|32
      # ADD ONE BASE AFTER THE LAST BASE:
      $exon_seq            = $exon_seq."N";
   }
   elsif ($exon_phase == 0 && $strand eq '+')
   {
      # ADD THREE BASES BEFORE THE FIRST BASE:
      $exon_seq            = "NNN".$exon_seq;
   }   
   elsif ($exon_phase == 0 && $strand eq '-')
   {
      # ADD THREE BASES AFTER THE LAST BASE:
      $exon_seq            = $exon_seq."NNN";
   }
 
   # REPLACE THE SEQUENCE FOR THIS EXON:
   $exon_seq_length = length($exon_seq);
   if ($exon_seq_length == 0) { $EXONSEQ->{$exon} = "none";   }
   else                       { $EXONSEQ->{$exon} = $exon_seq;}
 
   return($exon_seq,$errorcode,$errormsg);
}
 
#------------------------------------------------------------------#
 
# TEST &read_exon_gff
 
sub test_read_exon_gff
{
   my $outputdir           = $_[0]; # DIRECTORY TO WRITE OUTPUT FILES INTO
   my $errorcode;                   # RETURNED AS 0 FROM A FUNCTION IF THERE IS NO ERROR 
   my $errormsg;                    # RETURNED AS 'none' FROM A FUNCTION IF THERE IS NO ERROR
   my $EXONS;                       # HASH TABLE WITH EXONS IN A GENE
   my $EXON_START;                  # HASH TABLE WITH START POSITIONS OF EXONS IN A GENE
   my $EXON_FRAME;                  # HASH TABLE WITH FRAME OF EXONS IN A GENE
   my $exon_gff;                    # EXON GFF FILE
   my $input_gff;                   # THE INPUT GFF FILE 
 
   ($input_gff,$errorcode,$errormsg) = &make_filename($outputdir);
   if ($errorcode != 0) { ($errorcode,$errormsg) = &print_error($errormsg,$errorcode,0); } 
   open(INPUT_GFF,">$input_gff") || die "ERROR: test_read_exon_gff: cannot open $input_gff\n";
   print INPUT_GFF "##gff-version 3\n";
   print INPUT_GFF "##sequence-region Pk_strainH_chr01 1 838594\n";
   print INPUT_GFF "Pk_strainH_chr01	chado	gene	108	1421	.	+	.	ID=PKH_010010;isObsolete=false;isFminPartial;feature_id=222341;timelastmodified=10.12.2012+02:43:44+GMT\n";
   print INPUT_GFF "Pk_strainH_chr01	chado	CDS	108	758	.	+	0	ID=PKH_010010.1:exon:1;Parent=PKH_010010.1;isFminPartial;isObsolete=false;timelastmodified=10.12.2012+02:43:44+GMT;colour=7\n";
   print INPUT_GFF "Pk_strainH_chr01	chado	CDS	978	1421	.	+	0	ID=PKH_010010.1:exon:2;Parent=PKH_010010.1;isFminPartial;isObsolete=false;timelastmodified=10.12.2012+02:43:44+GMT;colour=7\n";
   print INPUT_GFF "Pk_strainH_chr01	chado	mRNA	108	1421	.	+	.	ID=PKH_010010.1;Parent=PKH_010010;isObsolete=false;isFminPartial;feature_id=222342;timelastmodified=10.12.2012+02:43:44+GMT;previous_systematic_id=PK00_0600c\n";
   print INPUT_GFF "Pk_strainH_chr01	chado	gene	5697	8325	.	+	.	ID=PKH_010020;isObsolete=false;feature_id=222257;timelastmodified=25.11.2012+01:47:32+GMT\n";
   print INPUT_GFF "Pk_strainH_chr01	chado	CDS	5697	6209	.	+	0	ID=PKH_010020.1:exon:1;Parent=PKH_010020.1;isObsolete=false;timelastmodified=25.11.2012+01:47:32+GMT;colour=7\n";
   print INPUT_GFF "Pk_strainH_chr01	chado	CDS	7891	8325	.	+	0	ID=PKH_010020.1:exon:2;Parent=PKH_010020.1;isObsolete=false;timelastmodified=25.11.2012+01:47:32+GMT;colour=7\n";
   print INPUT_GFF "Pk_strainH_chr01	chado	mRNA	5697	8325	.	+	.	ID=PKH_010020.1;Parent=PKH_010020;isObsolete=false;feature_id=222258;timelastmodified=25.11.2012+01:47:32+GMT;previous_systematic_id=PK9_3420w\n";
   print INPUT_GFF "Pk_strainH_chr01	chado	gene	8988	9593	.	-	.	ID=PKH_010030;isObsolete=false;feature_id=222118;isFmaxPartial;timelastmodified=12.01.2011+05:18:24+GMT\n";
   print INPUT_GFF "Pk_strainH_chr01	chado	CDS	8988	9593	.	-	0	ID=PKH_010030.1:exon:1;Parent=PKH_010030.1;isObsolete=false;timelastmodified=12.01.2011+05:18:24+GMT;colour=12\n";
   print INPUT_GFF "Pk_strainH_chr01	chado	mRNA	8988	9593	.	-	.	ID=PKH_010030.1;Parent=PKH_010030;isObsolete=false;feature_id=222119;timelastmodified=21.10.2007+01:51:33+BST;previous_systematic_id=PK4_2020w\n";
   print INPUT_GFF "Pk_strainH_chr01	chado	gap	10499	12498	.	+	.	ID=gap10499-12498:corrected;isObsolete=false;feature_id=222641;timelastmodified=21.10.2007+01:51:37+BST\n";
   print INPUT_GFF "Pk_strainH_chr01	chado	gene	13925	14846	.	-	.	ID=PKH_010040;isObsolete=false;feature_id=222402;timelastmodified=21.10.2007+01:51:35+BST\n";
   print INPUT_GFF "Pk_strainH_chr01	chado	CDS	14778	14846	.	-	0	ID=PKH_010040.1:exon:2;Parent=PKH_010040.1;isObsolete=false;timelastmodified=21.10.2007+01:51:35+BST;colour=8\n";
   print INPUT_GFF "Pk_strainH_chr01	chado	CDS	13925	14581	.	-	0	ID=PKH_010040.1:exon:1;Parent=PKH_010040.1;isObsolete=false;timelastmodified=21.10.2007+01:51:35+BST;colour=8\n";
   print INPUT_GFF "Pk_strainH_chr01	chado	mRNA	13925	14846	.	-	.	ID=PKH_010040.1;Parent=PKH_010040;isObsolete=false;feature_id=222403;timelastmodified=21.10.2007+01:51:35+BST\n";
   print INPUT_GFF "Pk_strainH_chr01	chado	gene	18394	18738	.	+	.	ID=PKH_010050;isObsolete=false;feature_id=222197;timelastmodified=21.10.2007+01:51:33+BST\n";
   print INPUT_GFF "Pk_strainH_chr01	chado	CDS	18394	18738	.	+	0	ID=PKH_010050.1:exon:1;Parent=PKH_010050.1;isObsolete=false;timelastmodified=21.10.2007+01:51:33+BST;colour=10\n";
   print INPUT_GFF "Pk_strainH_chr01	chado	mRNA	18394	18738	.	+	.	ID=PKH_010050.1;Parent=PKH_010050;isObsolete=false;feature_id=222198;timelastmodified=21.10.2007+01:51:33+BST;previous_systematic_id=PK7_0005w\n";
   print INPUT_GFF "Pk_strainH_chr01	chado	gene	22127	22613	.	-	.	ID=PKH_010080;isObsolete=false;feature_id=222375;timelastmodified=25.11.2012+01:16:50+GMT\n";
   print INPUT_GFF "Pk_strainH_chr01	chado	CDS	22404	22613	.	-	0	ID=PKH_010080.1:exon:2;Parent=PKH_010080.1;isObsolete=false;timelastmodified=25.11.2012+01:16:50+GMT;colour=10\n";
   print INPUT_GFF "Pk_strainH_chr01	chado	CDS	22127	22228	.	-	0	ID=PKH_010080.1:exon:1;Parent=PKH_010080.1;isObsolete=false;timelastmodified=25.11.2012+01:16:50+GMT;colour=10\n";
   print INPUT_GFF "Pk_strainH_chr01	chado	mRNA	22127	22613	.	-	.	ID=PKH_010080.1;Parent=PKH_010080;isObsolete=false;feature_id=222376;timelastmodified=25.11.2012+01:16:50+GMT\n";
   print INPUT_GFF "Pk_strainH_chr01	chado	gene	32144	34243	.	+	.	ID=PKH_010100;isObsolete=false;feature_id=222114;timelastmodified=21.10.2007+01:51:33+BST\n";
   print INPUT_GFF "Pk_strainH_chr01	chado	CDS	32144	34243	.	+	0	ID=PKH_010100.1:exon:1;Parent=PKH_010100.1;isObsolete=false;timelastmodified=21.10.2007+01:51:33+BST;colour=7\n";
   print INPUT_GFF "Pk_strainH_chr01	chado	mRNA	32144	34243	.	+	.	ID=PKH_010100.1;Parent=PKH_010100;isObsolete=false;feature_id=222115;timelastmodified=21.10.2007+01:51:33+BST;previous_systematic_id=PK7_0010w\n";
   print INPUT_GFF "Pk_strainH_chr01	chado	gene	35321	35393	.	-	.	ID=PKH_010102;isObsolete=false;feature_id=222634;timelastmodified=21.10.2007+01:51:37+BST\n";
   print INPUT_GFF "Pk_strainH_chr01	chado	CDS	35321	35393	.	-	0	ID=PKH_010102.1:exon:1;Parent=PKH_010102:tRNA;isObsolete=false;timelastmodified=21.10.2007+01:51:37+BST\n";
   print INPUT_GFF "Pk_strainH_chr01	chado	tRNA	35321	35393	.	-	.	ID=PKH_010102:tRNA;Parent=PKH_010102;isObsolete=false;feature_id=222635;product=term%3DtRNA+Alanine%3B;timelastmodified=07.09.2012+12:23:39+BST;comment=tRNA+Ala+anticodon+AGC%2C+Cove+score+51.85\n";
   close(INPUT_GFF);
   # MAKE A GFF FILE OF THE EXONS:
   ($exon_gff,$errorcode,$errormsg) = &make_exon_gff($input_gff,$outputdir);
   if ($errorcode != 0) { ($errorcode,$errormsg) = &print_error($errormsg,$errorcode,0); }
   # CALL &read_exon_gff: 
   ($EXONS,$EXON_START,$EXON_FRAME,$errorcode,$errormsg) = &read_exon_gff($exon_gff,'no','no');
   if ($errorcode != 0) { print STDERR "ERROR: test_read_exon_gff: failed test1\n"; exit;}
   if ($EXONS->{"PKH_010010.1"} ne "Pk_strainH_chr01#108#758#+#ID=PKH_010010.1:exon:1;Parent=PKH_010010.1;isFminPartial;isObsolete=false;timelastmodified=10.12.2012+02:43:44+GMT;colour=7*,*Pk_strainH_chr01#978#1421#+#ID=PKH_010010.1:exon:2;Parent=PKH_010010.1;isFminPartial;isObsolete=false;timelastmodified=10.12.2012+02:43:44+GMT;colour=7") { print STDERR "ERROR: test_read_exon_gff: failed test1b\n"; exit;}
   if ($EXONS->{"PKH_010020.1"} ne "Pk_strainH_chr01#5697#6209#+#ID=PKH_010020.1:exon:1;Parent=PKH_010020.1;isObsolete=false;timelastmodified=25.11.2012+01:47:32+GMT;colour=7*,*Pk_strainH_chr01#7891#8325#+#ID=PKH_010020.1:exon:2;Parent=PKH_010020.1;isObsolete=false;timelastmodified=25.11.2012+01:47:32+GMT;colour=7") { print STDERR "ERROR: test_read_exon_gff: failed test1c\n"; exit;}
   if ($EXONS->{"PKH_010030.1"} ne "Pk_strainH_chr01#8988#9593#-#ID=PKH_010030.1:exon:1;Parent=PKH_010030.1;isObsolete=false;timelastmodified=12.01.2011+05:18:24+GMT;colour=12") { print STDERR "ERROR: test_read_exon_gff: failed test1d\n"; exit;}
   if ($EXONS->{"PKH_010040.1"} ne "Pk_strainH_chr01#13925#14581#-#ID=PKH_010040.1:exon:1;Parent=PKH_010040.1;isObsolete=false;timelastmodified=21.10.2007+01:51:35+BST;colour=8*,*Pk_strainH_chr01#14778#14846#-#ID=PKH_010040.1:exon:2;Parent=PKH_010040.1;isObsolete=false;timelastmodified=21.10.2007+01:51:35+BST;colour=8") { print STDERR "ERROR: test_read_exon_gff: failed test1e\n"; exit;}
   if ($EXONS->{"PKH_010050.1"} ne "Pk_strainH_chr01#18394#18738#+#ID=PKH_010050.1:exon:1;Parent=PKH_010050.1;isObsolete=false;timelastmodified=21.10.2007+01:51:33+BST;colour=10") { print STDERR "ERROR: test_read_exon_gff: failed test1f\n"; exit;}
   if ($EXONS->{"PKH_010080.1"} ne "Pk_strainH_chr01#22127#22228#-#ID=PKH_010080.1:exon:1;Parent=PKH_010080.1;isObsolete=false;timelastmodified=25.11.2012+01:16:50+GMT;colour=10*,*Pk_strainH_chr01#22404#22613#-#ID=PKH_010080.1:exon:2;Parent=PKH_010080.1;isObsolete=false;timelastmodified=25.11.2012+01:16:50+GMT;colour=10") { print STDERR "ERROR: test_read_exon_gff: failed test1g\n"; exit;}
   if ($EXONS->{"PKH_010100.1"} ne "Pk_strainH_chr01#32144#34243#+#ID=PKH_010100.1:exon:1;Parent=PKH_010100.1;isObsolete=false;timelastmodified=21.10.2007+01:51:33+BST;colour=7") { print STDERR "ERROR: test_read_exon_gff: failed test1h\n"; exit;}
   if ($EXONS->{"PKH_010102:tRNA"} ne "Pk_strainH_chr01#35321#35393#-#ID=PKH_010102.1:exon:1;Parent=PKH_010102:tRNA;isObsolete=false;timelastmodified=21.10.2007+01:51:37+BST") { print STDERR "ERROR: test_read_exon_gff: failed test1i\n"; exit;}
   if ($EXON_START->{"Pk_strainH_chr01#108#758#+#ID=PKH_010010.1:exon:1;Parent=PKH_010010.1;isFminPartial;isObsolete=false;timelastmodified=10.12.2012+02:43:44+GMT;colour=7"} != 108) { print STDERR "ERROR: test_read_exon_gff: failed test1j\n"; exit;}
   if ($EXON_FRAME->{"Pk_strainH_chr01#108#758#+#ID=PKH_010010.1:exon:1;Parent=PKH_010010.1;isFminPartial;isObsolete=false;timelastmodified=10.12.2012+02:43:44+GMT;colour=7"} ne "zero") { print STDERR "ERROR: test_read_exon_gff: failed test1k\n"; exit;}
   if ($EXON_START->{"Pk_strainH_chr01#978#1421#+#ID=PKH_010010.1:exon:2;Parent=PKH_010010.1;isFminPartial;isObsolete=false;timelastmodified=10.12.2012+02:43:44+GMT;colour=7"} != 978) { print STDERR "ERROR: test_read_exon_gff: failed test1l\n"; exit;}
   if ($EXON_FRAME->{"Pk_strainH_chr01#978#1421#+#ID=PKH_010010.1:exon:2;Parent=PKH_010010.1;isFminPartial;isObsolete=false;timelastmodified=10.12.2012+02:43:44+GMT;colour=7"} ne "zero") { print STDERR "ERROR: test_read_exon_gff: failed test1m\n"; exit;}
   if ($EXON_START->{"Pk_strainH_chr01#5697#6209#+#ID=PKH_010020.1:exon:1;Parent=PKH_010020.1;isObsolete=false;timelastmodified=25.11.2012+01:47:32+GMT;colour=7"} != 5697) { print STDERR "ERROR: test_read_exon_gff: failed test1n\n"; exit;}
   if ($EXON_FRAME->{"Pk_strainH_chr01#5697#6209#+#ID=PKH_010020.1:exon:1;Parent=PKH_010020.1;isObsolete=false;timelastmodified=25.11.2012+01:47:32+GMT;colour=7"} ne "zero") { print STDERR "ERROR: test_read_exon_gff: failed test1o\n"; exit;}
   if ($EXON_START->{"Pk_strainH_chr01#35321#35393#-#ID=PKH_010102.1:exon:1;Parent=PKH_010102:tRNA;isObsolete=false;timelastmodified=21.10.2007+01:51:37+BST"} != 35321) { print STDERR "ERROR: test_read_exon_gff: failed test1p\n"; exit;}
   if ($EXON_FRAME->{"Pk_strainH_chr01#35321#35393#-#ID=PKH_010102.1:exon:1;Parent=PKH_010102:tRNA;isObsolete=false;timelastmodified=21.10.2007+01:51:37+BST"} ne "zero") { print STDERR "ERROR: test_read_exon_gff: failed test1q\n"; exit;}
   # DELETE TEMPORARY FILES:
   system "rm -f $input_gff";
   system "rm -f $exon_gff";
 
   # TEST ERRORCODE=20 (ALREADY KNOW START OF EXON):
   ($input_gff,$errorcode,$errormsg) = &make_filename($outputdir);
   if ($errorcode != 0) { ($errorcode,$errormsg) = &print_error($errormsg,$errorcode,0); } 
   open(INPUT_GFF,">$input_gff") || die "ERROR: test_read_exon_gff: cannot open $input_gff\n";
   print INPUT_GFF "##gff-version 3\n";
   print INPUT_GFF "##sequence-region Pk_strainH_chr01 1 838594\n";
   print INPUT_GFF "Pk_strainH_chr01	chado	gene	108	1421	.	+	.	ID=PKH_010010;isObsolete=false;isFminPartial;feature_id=222341;timelastmodified=10.12.2012+02:43:44+GMT\n";
   print INPUT_GFF "Pk_strainH_chr01	chado	CDS	108	758	.	+	0	ID=PKH_010010.1:exon:1;Parent=PKH_010010.1;isFminPartial;isObsolete=false;timelastmodified=10.12.2012+02:43:44+GMT;colour=7\n";
   print INPUT_GFF "Pk_strainH_chr01	chado	CDS	978	1421	.	+	0	ID=PKH_010010.1:exon:2;Parent=PKH_010010.1;isFminPartial;isObsolete=false;timelastmodified=10.12.2012+02:43:44+GMT;colour=7\n";
   print INPUT_GFF "Pk_strainH_chr01	chado	CDS	978	1421	.	+	0	ID=PKH_010010.1:exon:2;Parent=PKH_010010.1;isFminPartial;isObsolete=false;timelastmodified=10.12.2012+02:43:44+GMT;colour=7\n";
   print INPUT_GFF "Pk_strainH_chr01	chado	mRNA	108	1421	.	+	.	ID=PKH_010010.1;Parent=PKH_010010;isObsolete=false;isFminPartial;feature_id=222342;timelastmodified=10.12.2012+02:43:44+GMT;previous_systematic_id=PK00_0600c\n";
   print INPUT_GFF "Pk_strainH_chr01	chado	gene	5697	8325	.	+	.	ID=PKH_010020;isObsolete=false;feature_id=222257;timelastmodified=25.11.2012+01:47:32+GMT\n";
   print INPUT_GFF "Pk_strainH_chr01	chado	CDS	5697	6209	.	+	0	ID=PKH_010020.1:exon:1;Parent=PKH_010020.1;isObsolete=false;timelastmodified=25.11.2012+01:47:32+GMT;colour=7\n";
   print INPUT_GFF "Pk_strainH_chr01	chado	CDS	7891	8325	.	+	0	ID=PKH_010020.1:exon:2;Parent=PKH_010020.1;isObsolete=false;timelastmodified=25.11.2012+01:47:32+GMT;colour=7\n";
   print INPUT_GFF "Pk_strainH_chr01	chado	mRNA	5697	8325	.	+	.	ID=PKH_010020.1;Parent=PKH_010020;isObsolete=false;feature_id=222258;timelastmodified=25.11.2012+01:47:32+GMT;previous_systematic_id=PK9_3420w\n";
   print INPUT_GFF "Pk_strainH_chr01	chado	gene	8988	9593	.	-	.	ID=PKH_010030;isObsolete=false;feature_id=222118;isFmaxPartial;timelastmodified=12.01.2011+05:18:24+GMT\n";
   print INPUT_GFF "Pk_strainH_chr01	chado	CDS	8988	9593	.	-	0	ID=PKH_010030.1:exon:1;Parent=PKH_010030.1;isObsolete=false;timelastmodified=12.01.2011+05:18:24+GMT;colour=12\n";
   print INPUT_GFF "Pk_strainH_chr01	chado	mRNA	8988	9593	.	-	.	ID=PKH_010030.1;Parent=PKH_010030;isObsolete=false;feature_id=222119;timelastmodified=21.10.2007+01:51:33+BST;previous_systematic_id=PK4_2020w\n";
   print INPUT_GFF "Pk_strainH_chr01	chado	gap	10499	12498	.	+	.	ID=gap10499-12498:corrected;isObsolete=false;feature_id=222641;timelastmodified=21.10.2007+01:51:37+BST\n";
   print INPUT_GFF "Pk_strainH_chr01	chado	gene	13925	14846	.	-	.	ID=PKH_010040;isObsolete=false;feature_id=222402;timelastmodified=21.10.2007+01:51:35+BST\n";
   print INPUT_GFF "Pk_strainH_chr01	chado	CDS	14778	14846	.	-	0	ID=PKH_010040.1:exon:2;Parent=PKH_010040.1;isObsolete=false;timelastmodified=21.10.2007+01:51:35+BST;colour=8\n";
   print INPUT_GFF "Pk_strainH_chr01	chado	CDS	13925	14581	.	-	0	ID=PKH_010040.1:exon:1;Parent=PKH_010040.1;isObsolete=false;timelastmodified=21.10.2007+01:51:35+BST;colour=8\n";
   print INPUT_GFF "Pk_strainH_chr01	chado	mRNA	13925	14846	.	-	.	ID=PKH_010040.1;Parent=PKH_010040;isObsolete=false;feature_id=222403;timelastmodified=21.10.2007+01:51:35+BST\n";
   print INPUT_GFF "Pk_strainH_chr01	chado	gene	18394	18738	.	+	.	ID=PKH_010050;isObsolete=false;feature_id=222197;timelastmodified=21.10.2007+01:51:33+BST\n";
   print INPUT_GFF "Pk_strainH_chr01	chado	CDS	18394	18738	.	+	0	ID=PKH_010050.1:exon:1;Parent=PKH_010050.1;isObsolete=false;timelastmodified=21.10.2007+01:51:33+BST;colour=10\n";
   print INPUT_GFF "Pk_strainH_chr01	chado	mRNA	18394	18738	.	+	.	ID=PKH_010050.1;Parent=PKH_010050;isObsolete=false;feature_id=222198;timelastmodified=21.10.2007+01:51:33+BST;previous_systematic_id=PK7_0005w\n";
   print INPUT_GFF "Pk_strainH_chr01	chado	gene	22127	22613	.	-	.	ID=PKH_010080;isObsolete=false;feature_id=222375;timelastmodified=25.11.2012+01:16:50+GMT\n";
   print INPUT_GFF "Pk_strainH_chr01	chado	CDS	22404	22613	.	-	0	ID=PKH_010080.1:exon:2;Parent=PKH_010080.1;isObsolete=false;timelastmodified=25.11.2012+01:16:50+GMT;colour=10\n";
   print INPUT_GFF "Pk_strainH_chr01	chado	CDS	22127	22228	.	-	0	ID=PKH_010080.1:exon:1;Parent=PKH_010080.1;isObsolete=false;timelastmodified=25.11.2012+01:16:50+GMT;colour=10\n";
   print INPUT_GFF "Pk_strainH_chr01	chado	mRNA	22127	22613	.	-	.	ID=PKH_010080.1;Parent=PKH_010080;isObsolete=false;feature_id=222376;timelastmodified=25.11.2012+01:16:50+GMT\n";
   print INPUT_GFF "Pk_strainH_chr01	chado	gene	32144	34243	.	+	.	ID=PKH_010100;isObsolete=false;feature_id=222114;timelastmodified=21.10.2007+01:51:33+BST\n";
   print INPUT_GFF "Pk_strainH_chr01	chado	CDS	32144	34243	.	+	0	ID=PKH_010100.1:exon:1;Parent=PKH_010100.1;isObsolete=false;timelastmodified=21.10.2007+01:51:33+BST;colour=7\n";
   print INPUT_GFF "Pk_strainH_chr01	chado	mRNA	32144	34243	.	+	.	ID=PKH_010100.1;Parent=PKH_010100;isObsolete=false;feature_id=222115;timelastmodified=21.10.2007+01:51:33+BST;previous_systematic_id=PK7_0010w\n";
   print INPUT_GFF "Pk_strainH_chr01	chado	gene	35321	35393	.	-	.	ID=PKH_010102;isObsolete=false;feature_id=222634;timelastmodified=21.10.2007+01:51:37+BST\n";
   print INPUT_GFF "Pk_strainH_chr01	chado	CDS	35321	35393	.	-	0	ID=PKH_010102.1:exon:1;Parent=PKH_010102:tRNA;isObsolete=false;timelastmodified=21.10.2007+01:51:37+BST\n";
   print INPUT_GFF "Pk_strainH_chr01	chado	tRNA	35321	35393	.	-	.	ID=PKH_010102:tRNA;Parent=PKH_010102;isObsolete=false;feature_id=222635;product=term%3DtRNA+Alanine%3B;timelastmodified=07.09.2012+12:23:39+BST;comment=tRNA+Ala+anticodon+AGC%2C+Cove+score+51.85\n";
   close(INPUT_GFF);
   # MAKE A GFF FILE OF THE EXONS:
   ($exon_gff,$errorcode,$errormsg) = &make_exon_gff($input_gff,$outputdir);
   if ($errorcode != 0) { ($errorcode,$errormsg) = &print_error($errormsg,$errorcode,0); }
   # CALL &read_exon_gff: 
   ($EXONS,$EXON_START,$EXON_FRAME,$errorcode,$errormsg) = &read_exon_gff($exon_gff,'no','no');
   if ($errorcode != 20) { print STDERR "ERROR: test_read_exon_gff: failed test2 (errorcode $errorcode errormsg $errormsg)\n"; exit;}
   # DELETE TEMPORARY FILES:
   system "rm -f $input_gff";
   system "rm -f $exon_gff";
 
   # TEST ERRORCODE=31 (PHASE IS NOT 0, 1 OR 2):
   ($input_gff,$errorcode,$errormsg) = &make_filename($outputdir);
   if ($errorcode != 0) { ($errorcode,$errormsg) = &print_error($errormsg,$errorcode,0); } 
   open(INPUT_GFF,">$input_gff") || die "ERROR: test_read_exon_gff: cannot open $input_gff\n";
   print INPUT_GFF "##gff-version 3\n";
   print INPUT_GFF "##sequence-region Pk_strainH_chr01 1 838594\n";
   print INPUT_GFF "Pk_strainH_chr01	chado	gene	108	1421	.	+	.	ID=PKH_010010;isObsolete=false;isFminPartial;feature_id=222341;timelastmodified=10.12.2012+02:43:44+GMT\n";
   print INPUT_GFF "Pk_strainH_chr01	chado	CDS	108	758	.	+	0	ID=PKH_010010.1:exon:1;Parent=PKH_010010.1;isFminPartial;isObsolete=false;timelastmodified=10.12.2012+02:43:44+GMT;colour=7\n";
   print INPUT_GFF "Pk_strainH_chr01	chado	CDS	978	1421	.	+	3	ID=PKH_010010.1:exon:2;Parent=PKH_010010.1;isFminPartial;isObsolete=false;timelastmodified=10.12.2012+02:43:44+GMT;colour=7\n";
   print INPUT_GFF "Pk_strainH_chr01	chado	mRNA	108	1421	.	+	.	ID=PKH_010010.1;Parent=PKH_010010;isObsolete=false;isFminPartial;feature_id=222342;timelastmodified=10.12.2012+02:43:44+GMT;previous_systematic_id=PK00_0600c\n";
   print INPUT_GFF "Pk_strainH_chr01	chado	gene	5697	8325	.	+	.	ID=PKH_010020;isObsolete=false;feature_id=222257;timelastmodified=25.11.2012+01:47:32+GMT\n";
   print INPUT_GFF "Pk_strainH_chr01	chado	CDS	5697	6209	.	+	0	ID=PKH_010020.1:exon:1;Parent=PKH_010020.1;isObsolete=false;timelastmodified=25.11.2012+01:47:32+GMT;colour=7\n";
   print INPUT_GFF "Pk_strainH_chr01	chado	CDS	7891	8325	.	+	0	ID=PKH_010020.1:exon:2;Parent=PKH_010020.1;isObsolete=false;timelastmodified=25.11.2012+01:47:32+GMT;colour=7\n";
   print INPUT_GFF "Pk_strainH_chr01	chado	mRNA	5697	8325	.	+	.	ID=PKH_010020.1;Parent=PKH_010020;isObsolete=false;feature_id=222258;timelastmodified=25.11.2012+01:47:32+GMT;previous_systematic_id=PK9_3420w\n";
   print INPUT_GFF "Pk_strainH_chr01	chado	gene	8988	9593	.	-	.	ID=PKH_010030;isObsolete=false;feature_id=222118;isFmaxPartial;timelastmodified=12.01.2011+05:18:24+GMT\n";
   print INPUT_GFF "Pk_strainH_chr01	chado	CDS	8988	9593	.	-	0	ID=PKH_010030.1:exon:1;Parent=PKH_010030.1;isObsolete=false;timelastmodified=12.01.2011+05:18:24+GMT;colour=12\n";
   print INPUT_GFF "Pk_strainH_chr01	chado	mRNA	8988	9593	.	-	.	ID=PKH_010030.1;Parent=PKH_010030;isObsolete=false;feature_id=222119;timelastmodified=21.10.2007+01:51:33+BST;previous_systematic_id=PK4_2020w\n";
   print INPUT_GFF "Pk_strainH_chr01	chado	gap	10499	12498	.	+	.	ID=gap10499-12498:corrected;isObsolete=false;feature_id=222641;timelastmodified=21.10.2007+01:51:37+BST\n";
   print INPUT_GFF "Pk_strainH_chr01	chado	gene	13925	14846	.	-	.	ID=PKH_010040;isObsolete=false;feature_id=222402;timelastmodified=21.10.2007+01:51:35+BST\n";
   print INPUT_GFF "Pk_strainH_chr01	chado	CDS	14778	14846	.	-	0	ID=PKH_010040.1:exon:2;Parent=PKH_010040.1;isObsolete=false;timelastmodified=21.10.2007+01:51:35+BST;colour=8\n";
   print INPUT_GFF "Pk_strainH_chr01	chado	CDS	13925	14581	.	-	0	ID=PKH_010040.1:exon:1;Parent=PKH_010040.1;isObsolete=false;timelastmodified=21.10.2007+01:51:35+BST;colour=8\n";
   print INPUT_GFF "Pk_strainH_chr01	chado	mRNA	13925	14846	.	-	.	ID=PKH_010040.1;Parent=PKH_010040;isObsolete=false;feature_id=222403;timelastmodified=21.10.2007+01:51:35+BST\n";
   print INPUT_GFF "Pk_strainH_chr01	chado	gene	18394	18738	.	+	.	ID=PKH_010050;isObsolete=false;feature_id=222197;timelastmodified=21.10.2007+01:51:33+BST\n";
   print INPUT_GFF "Pk_strainH_chr01	chado	CDS	18394	18738	.	+	0	ID=PKH_010050.1:exon:1;Parent=PKH_010050.1;isObsolete=false;timelastmodified=21.10.2007+01:51:33+BST;colour=10\n";
   print INPUT_GFF "Pk_strainH_chr01	chado	mRNA	18394	18738	.	+	.	ID=PKH_010050.1;Parent=PKH_010050;isObsolete=false;feature_id=222198;timelastmodified=21.10.2007+01:51:33+BST;previous_systematic_id=PK7_0005w\n";
   print INPUT_GFF "Pk_strainH_chr01	chado	gene	22127	22613	.	-	.	ID=PKH_010080;isObsolete=false;feature_id=222375;timelastmodified=25.11.2012+01:16:50+GMT\n";
   print INPUT_GFF "Pk_strainH_chr01	chado	CDS	22404	22613	.	-	0	ID=PKH_010080.1:exon:2;Parent=PKH_010080.1;isObsolete=false;timelastmodified=25.11.2012+01:16:50+GMT;colour=10\n";
   print INPUT_GFF "Pk_strainH_chr01	chado	CDS	22127	22228	.	-	0	ID=PKH_010080.1:exon:1;Parent=PKH_010080.1;isObsolete=false;timelastmodified=25.11.2012+01:16:50+GMT;colour=10\n";
   print INPUT_GFF "Pk_strainH_chr01	chado	mRNA	22127	22613	.	-	.	ID=PKH_010080.1;Parent=PKH_010080;isObsolete=false;feature_id=222376;timelastmodified=25.11.2012+01:16:50+GMT\n";
   print INPUT_GFF "Pk_strainH_chr01	chado	gene	32144	34243	.	+	.	ID=PKH_010100;isObsolete=false;feature_id=222114;timelastmodified=21.10.2007+01:51:33+BST\n";
   print INPUT_GFF "Pk_strainH_chr01	chado	CDS	32144	34243	.	+	0	ID=PKH_010100.1:exon:1;Parent=PKH_010100.1;isObsolete=false;timelastmodified=21.10.2007+01:51:33+BST;colour=7\n";
   print INPUT_GFF "Pk_strainH_chr01	chado	mRNA	32144	34243	.	+	.	ID=PKH_010100.1;Parent=PKH_010100;isObsolete=false;feature_id=222115;timelastmodified=21.10.2007+01:51:33+BST;previous_systematic_id=PK7_0010w\n";
   print INPUT_GFF "Pk_strainH_chr01	chado	gene	35321	35393	.	-	.	ID=PKH_010102;isObsolete=false;feature_id=222634;timelastmodified=21.10.2007+01:51:37+BST\n";
   print INPUT_GFF "Pk_strainH_chr01	chado	CDS	35321	35393	.	-	0	ID=PKH_010102.1:exon:1;Parent=PKH_010102:tRNA;isObsolete=false;timelastmodified=21.10.2007+01:51:37+BST\n";
   print INPUT_GFF "Pk_strainH_chr01	chado	tRNA	35321	35393	.	-	.	ID=PKH_010102:tRNA;Parent=PKH_010102;isObsolete=false;feature_id=222635;product=term%3DtRNA+Alanine%3B;timelastmodified=07.09.2012+12:23:39+BST;comment=tRNA+Ala+anticodon+AGC%2C+Cove+score+51.85\n";
   close(INPUT_GFF);
   # MAKE A GFF FILE OF THE EXONS:
   ($exon_gff,$errorcode,$errormsg) = &make_exon_gff($input_gff,$outputdir);
   if ($errorcode != 0) { ($errorcode,$errormsg) = &print_error($errormsg,$errorcode,0); }
   # CALL &read_exon_gff: 
   ($EXONS,$EXON_START,$EXON_FRAME,$errorcode,$errormsg) = &read_exon_gff($exon_gff,'no','no');
   if ($errorcode != 31) { print STDERR "ERROR: test_read_exon_gff: failed test3\n"; exit;}
   # DELETE TEMPORARY FILES:
   system "rm -f $input_gff";
   system "rm -f $exon_gff";
 
   # TEST ERRORCODE=32 ($ignore_phase IS NOT yes/no):
   ($input_gff,$errorcode,$errormsg) = &make_filename($outputdir);
   if ($errorcode != 0) { ($errorcode,$errormsg) = &print_error($errormsg,$errorcode,0); } 
   open(INPUT_GFF,">$input_gff") || die "ERROR: test_read_exon_gff: cannot open $input_gff\n";
   print INPUT_GFF "##gff-version 3\n";
   print INPUT_GFF "##sequence-region Pk_strainH_chr01 1 838594\n";
   print INPUT_GFF "Pk_strainH_chr01	chado	gene	108	1421	.	+	.	ID=PKH_010010;isObsolete=false;isFminPartial;feature_id=222341;timelastmodified=10.12.2012+02:43:44+GMT\n";
   print INPUT_GFF "Pk_strainH_chr01	chado	CDS	108	758	.	+	0	ID=PKH_010010.1:exon:1;Parent=PKH_010010.1;isFminPartial;isObsolete=false;timelastmodified=10.12.2012+02:43:44+GMT;colour=7\n";
   print INPUT_GFF "Pk_strainH_chr01	chado	CDS	978	1421	.	+	1	ID=PKH_010010.1:exon:2;Parent=PKH_010010.1;isFminPartial;isObsolete=false;timelastmodified=10.12.2012+02:43:44+GMT;colour=7\n";
   print INPUT_GFF "Pk_strainH_chr01	chado	mRNA	108	1421	.	+	.	ID=PKH_010010.1;Parent=PKH_010010;isObsolete=false;isFminPartial;feature_id=222342;timelastmodified=10.12.2012+02:43:44+GMT;previous_systematic_id=PK00_0600c\n";
   print INPUT_GFF "Pk_strainH_chr01	chado	gene	5697	8325	.	+	.	ID=PKH_010020;isObsolete=false;feature_id=222257;timelastmodified=25.11.2012+01:47:32+GMT\n";
   print INPUT_GFF "Pk_strainH_chr01	chado	CDS	5697	6209	.	+	0	ID=PKH_010020.1:exon:1;Parent=PKH_010020.1;isObsolete=false;timelastmodified=25.11.2012+01:47:32+GMT;colour=7\n";
   print INPUT_GFF "Pk_strainH_chr01	chado	CDS	7891	8325	.	+	0	ID=PKH_010020.1:exon:2;Parent=PKH_010020.1;isObsolete=false;timelastmodified=25.11.2012+01:47:32+GMT;colour=7\n";
   print INPUT_GFF "Pk_strainH_chr01	chado	mRNA	5697	8325	.	+	.	ID=PKH_010020.1;Parent=PKH_010020;isObsolete=false;feature_id=222258;timelastmodified=25.11.2012+01:47:32+GMT;previous_systematic_id=PK9_3420w\n";
   print INPUT_GFF "Pk_strainH_chr01	chado	gene	8988	9593	.	-	.	ID=PKH_010030;isObsolete=false;feature_id=222118;isFmaxPartial;timelastmodified=12.01.2011+05:18:24+GMT\n";
   print INPUT_GFF "Pk_strainH_chr01	chado	CDS	8988	9593	.	-	0	ID=PKH_010030.1:exon:1;Parent=PKH_010030.1;isObsolete=false;timelastmodified=12.01.2011+05:18:24+GMT;colour=12\n";
   print INPUT_GFF "Pk_strainH_chr01	chado	mRNA	8988	9593	.	-	.	ID=PKH_010030.1;Parent=PKH_010030;isObsolete=false;feature_id=222119;timelastmodified=21.10.2007+01:51:33+BST;previous_systematic_id=PK4_2020w\n";
   print INPUT_GFF "Pk_strainH_chr01	chado	gap	10499	12498	.	+	.	ID=gap10499-12498:corrected;isObsolete=false;feature_id=222641;timelastmodified=21.10.2007+01:51:37+BST\n";
   print INPUT_GFF "Pk_strainH_chr01	chado	gene	13925	14846	.	-	.	ID=PKH_010040;isObsolete=false;feature_id=222402;timelastmodified=21.10.2007+01:51:35+BST\n";
   print INPUT_GFF "Pk_strainH_chr01	chado	CDS	14778	14846	.	-	0	ID=PKH_010040.1:exon:2;Parent=PKH_010040.1;isObsolete=false;timelastmodified=21.10.2007+01:51:35+BST;colour=8\n";
   print INPUT_GFF "Pk_strainH_chr01	chado	CDS	13925	14581	.	-	0	ID=PKH_010040.1:exon:1;Parent=PKH_010040.1;isObsolete=false;timelastmodified=21.10.2007+01:51:35+BST;colour=8\n";
   print INPUT_GFF "Pk_strainH_chr01	chado	mRNA	13925	14846	.	-	.	ID=PKH_010040.1;Parent=PKH_010040;isObsolete=false;feature_id=222403;timelastmodified=21.10.2007+01:51:35+BST\n";
   print INPUT_GFF "Pk_strainH_chr01	chado	gene	18394	18738	.	+	.	ID=PKH_010050;isObsolete=false;feature_id=222197;timelastmodified=21.10.2007+01:51:33+BST\n";
   print INPUT_GFF "Pk_strainH_chr01	chado	CDS	18394	18738	.	+	0	ID=PKH_010050.1:exon:1;Parent=PKH_010050.1;isObsolete=false;timelastmodified=21.10.2007+01:51:33+BST;colour=10\n";
   print INPUT_GFF "Pk_strainH_chr01	chado	mRNA	18394	18738	.	+	.	ID=PKH_010050.1;Parent=PKH_010050;isObsolete=false;feature_id=222198;timelastmodified=21.10.2007+01:51:33+BST;previous_systematic_id=PK7_0005w\n";
   print INPUT_GFF "Pk_strainH_chr01	chado	gene	22127	22613	.	-	.	ID=PKH_010080;isObsolete=false;feature_id=222375;timelastmodified=25.11.2012+01:16:50+GMT\n";
   print INPUT_GFF "Pk_strainH_chr01	chado	CDS	22404	22613	.	-	0	ID=PKH_010080.1:exon:2;Parent=PKH_010080.1;isObsolete=false;timelastmodified=25.11.2012+01:16:50+GMT;colour=10\n";
   print INPUT_GFF "Pk_strainH_chr01	chado	CDS	22127	22228	.	-	0	ID=PKH_010080.1:exon:1;Parent=PKH_010080.1;isObsolete=false;timelastmodified=25.11.2012+01:16:50+GMT;colour=10\n";
   print INPUT_GFF "Pk_strainH_chr01	chado	mRNA	22127	22613	.	-	.	ID=PKH_010080.1;Parent=PKH_010080;isObsolete=false;feature_id=222376;timelastmodified=25.11.2012+01:16:50+GMT\n";
   print INPUT_GFF "Pk_strainH_chr01	chado	gene	32144	34243	.	+	.	ID=PKH_010100;isObsolete=false;feature_id=222114;timelastmodified=21.10.2007+01:51:33+BST\n";
   print INPUT_GFF "Pk_strainH_chr01	chado	CDS	32144	34243	.	+	0	ID=PKH_010100.1:exon:1;Parent=PKH_010100.1;isObsolete=false;timelastmodified=21.10.2007+01:51:33+BST;colour=7\n";
   print INPUT_GFF "Pk_strainH_chr01	chado	mRNA	32144	34243	.	+	.	ID=PKH_010100.1;Parent=PKH_010100;isObsolete=false;feature_id=222115;timelastmodified=21.10.2007+01:51:33+BST;previous_systematic_id=PK7_0010w\n";
   print INPUT_GFF "Pk_strainH_chr01	chado	gene	35321	35393	.	-	.	ID=PKH_010102;isObsolete=false;feature_id=222634;timelastmodified=21.10.2007+01:51:37+BST\n";
   print INPUT_GFF "Pk_strainH_chr01	chado	CDS	35321	35393	.	-	0	ID=PKH_010102.1:exon:1;Parent=PKH_010102:tRNA;isObsolete=false;timelastmodified=21.10.2007+01:51:37+BST\n";
   print INPUT_GFF "Pk_strainH_chr01	chado	tRNA	35321	35393	.	-	.	ID=PKH_010102:tRNA;Parent=PKH_010102;isObsolete=false;feature_id=222635;product=term%3DtRNA+Alanine%3B;timelastmodified=07.09.2012+12:23:39+BST;comment=tRNA+Ala+anticodon+AGC%2C+Cove+score+51.85\n";
   close(INPUT_GFF);
   # MAKE A GFF FILE OF THE EXONS:
   ($exon_gff,$errorcode,$errormsg) = &make_exon_gff($input_gff,$outputdir);
   if ($errorcode != 0) { ($errorcode,$errormsg) = &print_error($errormsg,$errorcode,0); }
   # CALL &read_exon_gff: 
   ($EXONS,$EXON_START,$EXON_FRAME,$errorcode,$errormsg) = &read_exon_gff($exon_gff,'hello','no');
   if ($errorcode != 32) { print STDERR "ERROR: test_read_exon_gff: failed test4\n"; exit;}
   # DELETE TEMPORARY FILES:
   system "rm -f $input_gff";
   system "rm -f $exon_gff";
 
   # TEST READING A C. ELEGANS EXON GFF FORMAT (GFF3 DOWNLOADED FROM WORMBASE, PROCESSED USING make_small_gff_for_rnaseqqc.pl):
   ($input_gff,$errorcode,$errormsg) = &make_filename($outputdir);   
   if ($errorcode != 0) { ($errorcode,$errormsg) = &print_error($errormsg,$errorcode,0); }   
   open(INPUT_GFF,">$input_gff") || die "ERROR: test_read_exon_gff: cannot open $input_gff\n";
   print INPUT_GFF "CHROMOSOME_X\tCoding_transcript\tCDS\t2446067\t2446593\t.\t-\t2\tID=CDS:T07D1.4;Note=RNA-binding protein;Parent=Transcript:T07D1.4.1,Transcript:T07D1.4.2;locus=fox-1;status=Partially_confirmed;wormpep=CE:CE25105\n";
   print INPUT_GFF "CHROMOSOME_X\tCoding_transcript\tCDS\t2446639\t2446783\t.\t-\t0\tID=CDS:T07D1.4;Note=RNA-binding protein;Parent=Transcript:T07D1.4.1,Transcript:T07D1.4.2;locus=fox-1;status=Partially_confirmed;wormpep=CE:CE25105\n";
   print INPUT_GFF "CHROMOSOME_X\tCoding_transcript\tCDS\t2446839\t2446975\t.\t-\t2\tID=CDS:T07D1.4;Note=RNA-binding protein;Parent=Transcript:T07D1.4.1,Transcript:T07D1.4.2;locus=fox-1;status=Partially_confirmed;wormpep=CE:CE25105\n";
   print INPUT_GFF "CHROMOSOME_X\tCoding_transcript\tCDS\t2447335\t2447530\t.\t-\t0\tID=CDS:T07D1.4;Note=RNA-binding protein;Parent=Transcript:T07D1.4.1,Transcript:T07D1.4.2;locus=fox-1;status=Partially_confirmed;wormpep=CE:CE25105\n";
   print INPUT_GFF "CHROMOSOME_X\tCoding_transcript\tCDS\t2448754\t2448842\t.\t-\t2\tID=CDS:T07D1.4;Note=RNA-binding protein;Parent=Transcript:T07D1.4.1,Transcript:T07D1.4.2;locus=fox-1;status=Partially_confirmed;wormpep=CE:CE25105\n"; 
   print INPUT_GFF "CHROMOSOME_X\tCoding_transcript\tCDS\t2450334\t2450487\t.\t-\t0\tID=CDS:T07D1.4;Note=RNA-binding protein;Parent=Transcript:T07D1.4.1,Transcript:T07D1.4.2;locus=fox-1;status=Partially_confirmed;wormpep=CE:CE25105\n";
   print INPUT_GFF "CHROMOSOME_X\tCoding_transcript\tgene\t2445794\t2450500\t.\t-\t.\tID=Gene:T07D1.4\n";
   print INPUT_GFF "CHROMOSOME_X\tCoding_transcript\tmRNA\t2445794\t2450496\t.\t-\t.\tID=Transcript:T07D1.4.1;Parent=Gene:T07D1.4;cds=T07D1.4;wormpep=CE:CE25105\n";
   print INPUT_GFF "CHROMOSOME_X\tCoding_transcript\tmRNA\t2445794\t2450500\t.\t-\t.\tID=Transcript:T07D1.4.2;Parent=Gene:T07D1.4;cds=T07D1.4;wormpep=CE:CE25105\n";
   close(INPUT_GFF);
   # MAKE A GFF FILE OF THE EXONS:
   ($exon_gff,$errorcode,$errormsg) = &make_exon_gff($input_gff,$outputdir);
   if ($errorcode != 0) { ($errorcode,$errormsg) = &print_error($errormsg,$errorcode,0); }
   # CALL &read_exon_gff:
   ($EXONS,$EXON_START,$EXON_FRAME,$errorcode,$errormsg) = &read_exon_gff($exon_gff,'no','no');
   if ($errorcode != 0) { print STDERR "ERROR: test_read_exon_gff: failed test5 (errorcode $errorcode errormsg $errormsg)\n"; exit;}
   if ($EXONS->{"Transcript:T07D1.4.1"} ne "CHROMOSOME_X#2446067#2446593#-#ID=CDS:T07D1.4;Note=RNA-binding protein;Parent=Transcript:T07D1.4.1,Transcript:T07D1.4.2;locus=fox-1;status=Partially_confirmed;wormpep=CE:CE25105*,*CHROMOSOME_X#2446639#2446783#-#ID=CDS:T07D1.4;Note=RNA-binding protein;Parent=Transcript:T07D1.4.1,Transcript:T07D1.4.2;locus=fox-1;status=Partially_confirmed;wormpep=CE:CE25105*,*CHROMOSOME_X#2446839#2446975#-#ID=CDS:T07D1.4;Note=RNA-binding protein;Parent=Transcript:T07D1.4.1,Transcript:T07D1.4.2;locus=fox-1;status=Partially_confirmed;wormpep=CE:CE25105*,*CHROMOSOME_X#2447335#2447530#-#ID=CDS:T07D1.4;Note=RNA-binding protein;Parent=Transcript:T07D1.4.1,Transcript:T07D1.4.2;locus=fox-1;status=Partially_confirmed;wormpep=CE:CE25105*,*CHROMOSOME_X#2448754#2448842#-#ID=CDS:T07D1.4;Note=RNA-binding protein;Parent=Transcript:T07D1.4.1,Transcript:T07D1.4.2;locus=fox-1;status=Partially_confirmed;wormpep=CE:CE25105*,*CHROMOSOME_X#2450334#2450487#-#ID=CDS:T07D1.4;Note=RNA-binding protein;Parent=Transcript:T07D1.4.1,Transcript:T07D1.4.2;locus=fox-1;status=Partially_confirmed;wormpep=CE:CE25105") { print STDERR "ERROR: test_read_exon_gff: failed test5b\n"; exit;} 
   if ($EXONS->{"Transcript:T07D1.4.2"} ne "CHROMOSOME_X#2446067#2446593#-#ID=CDS:T07D1.4;Note=RNA-binding protein;Parent=Transcript:T07D1.4.1,Transcript:T07D1.4.2;locus=fox-1;status=Partially_confirmed;wormpep=CE:CE25105*,*CHROMOSOME_X#2446639#2446783#-#ID=CDS:T07D1.4;Note=RNA-binding protein;Parent=Transcript:T07D1.4.1,Transcript:T07D1.4.2;locus=fox-1;status=Partially_confirmed;wormpep=CE:CE25105*,*CHROMOSOME_X#2446839#2446975#-#ID=CDS:T07D1.4;Note=RNA-binding protein;Parent=Transcript:T07D1.4.1,Transcript:T07D1.4.2;locus=fox-1;status=Partially_confirmed;wormpep=CE:CE25105*,*CHROMOSOME_X#2447335#2447530#-#ID=CDS:T07D1.4;Note=RNA-binding protein;Parent=Transcript:T07D1.4.1,Transcript:T07D1.4.2;locus=fox-1;status=Partially_confirmed;wormpep=CE:CE25105*,*CHROMOSOME_X#2448754#2448842#-#ID=CDS:T07D1.4;Note=RNA-binding protein;Parent=Transcript:T07D1.4.1,Transcript:T07D1.4.2;locus=fox-1;status=Partially_confirmed;wormpep=CE:CE25105*,*CHROMOSOME_X#2450334#2450487#-#ID=CDS:T07D1.4;Note=RNA-binding protein;Parent=Transcript:T07D1.4.1,Transcript:T07D1.4.2;locus=fox-1;status=Partially_confirmed;wormpep=CE:CE25105") { print STDERR "ERROR: test_read_exon_gff: failed test5c\n"; exit;}
   if ($EXON_START->{"CHROMOSOME_X#2446067#2446593#-#ID=CDS:T07D1.4;Note=RNA-binding protein;Parent=Transcript:T07D1.4.1,Transcript:T07D1.4.2;locus=fox-1;status=Partially_confirmed;wormpep=CE:CE25105"} != 2446067) { print STDERR "ERROR: test_read_exon_gff: failed test5c\n"; exit;}
   if ($EXON_FRAME->{"CHROMOSOME_X#2446067#2446593#-#ID=CDS:T07D1.4;Note=RNA-binding protein;Parent=Transcript:T07D1.4.1,Transcript:T07D1.4.2;locus=fox-1;status=Partially_confirmed;wormpep=CE:CE25105"} != 2) { print STDERR "ERROR: test_read_exon_gff: failed test5d\n"; exit;}
   system "rm -f $input_gff";   
   system "rm -f $exon_gff";
}
 
#------------------------------------------------------------------#
 
# READ IN THE EXONS FROM EACH GENE FROM THE GFF FILE:
 
sub read_exon_gff
{
   my $exon_gff            = $_[0]; # EXON GFF
   my $ignore_phase        = $_[1]; # SAYS WHETHER TO IGNORE THE PHASE INFORMATION IN THE GFF FILE 
   my $from_augustus       = $_[2]; # SAYS WHETHER THE GFF FILE IS FROM AUGUSTUS 
   my %EXONS               = ();    # HASH TABLE OF THE EXONS IN EACH GENE
   my %EXON_START          = ();    # HASH TABLE OF THE START POINTS OF EXONS
   my %EXON_FRAME          = ();    # HASH TABLE OF THE FRAME OF EXONS
   my $errorcode           = 0;     # RETURNED AS 0 IF THERE IS NO ERROR
   my $errormsg            = "none";# RETURNED AS 'none' IF THERE IS NO ERROR
   my $line;                        # 
   my @temp;                        #  
   my $scaffold;                    # NAME OF THE SCAFFOLD
   my $start;                       # START POSITION OF THE EXON
   my $end;                         # END OF THE EXON
   my $strand;                      # STRAND OF THE EXON 
   my $frame;                       # FRAME OF THE EXON
   my $exon;                        # NAME OF THE EXON
   my $gene;                        # NAME OF THE GENE 
   my @genes;                       # NAME OF THE GENES (ACTUALLY TRANSCRIPTS, IF THERE ARE SEVERAL TRANSCRIPTS FOR ONE GENE)
   my $i;                           # 
 
   # CHECK THAT $ignore_phase IS yes/no:
   if ($ignore_phase ne 'yes' && $ignore_phase ne 'no')
   {
      $errormsg            = "ERROR: read_exon_gff: ignore_phase $ignore_phase\n";
      $errorcode           = 32; # ERRORCODE=32 (TESTED FOR)
      return(\%EXONS,\%EXON_START,\%EXON_FRAME,$errorcode,$errormsg);
   }
 
   # READ IN THE GFF FILE:
   open(EXON_GFF,"$exon_gff") || die "ERROR: read_exon_gff: cannot open $exon_gff.\n";
   while(<EXON_GFF>)
   {
      $line                = $_;
      chomp $line;
      @temp                = split(/\t+/,$line);
      $scaffold            = $temp[0];
      $start               = $temp[3];
      $end                 = $temp[4]; 
      $strand              = $temp[6];
      $frame               = $temp[7];
      # IF IT IS A CDS FROM a rRNA GENE, TREAT THE PHASE AS 0: (THIS ALLOWS US TO GET THE SEQUENCE USING get_spliced_transcripts_from_gff.pl
      # eg. CHROMOSOME_V    rRNA    CDS     17121920        17122038        .       +       .       Parent=Transcript:ZK218.20
      if ($temp[1] eq 'rRNA') { $frame = '0';}
      $exon                = $temp[8];
      $exon                = $scaffold."#".$start."#".$end."#".$strand."#".$exon; 
      # FIND THE NAME OF THE GENE:
      if ($from_augustus eq 'no')
      {
         @temp             = split(/Parent=/,$exon); # eg. ID=0:SSTP.scaffold.00025.352065-TF352220-GW-3:cds;Parent=0:SSTP.scaffold.00025.352065-TF352220-GW-3
         $gene             = $temp[1];               # eg. 0:SSTP.scaffold.00025.352065-TF352220-GW-3 
         # NOTE IN WORMBASE GFF3 FORMAT, THERE COULD BE AN EXON IN MORE THAN ONE TRANSCRIPT eg. Parent=Transcript:T07D1.4.1,Transcript:T07D1.4.2 
         @temp             = split(/\;/,$gene); 
         $gene             = $temp[0];               # eg. 0:SSTP.scaffold.00025.352065-TF352220-GW-3
         @genes            = split(/\,/,$gene);
      }
      elsif ($from_augustus eq 'yes') # GFF IS FROM AUGUSTUS 
      {
         # eg. transcript_id "g1.t1"; gene_id "g1"; 
         @temp             = split(/gene_id/,$exon);
         $gene             = $temp[1]; # eg.  "g1"; 
         @temp             = split(/\"/,$gene);
         $gene             = $temp[1]; # eg. g1
         @genes            = split(/\,/,$gene); 
      }
      for ($i = 0; $i <= $#genes; $i++)
      {
         $gene             = $genes[$i];
         # RECORD THE EXONS IN THE GENE:
         if (!($EXONS{$gene})) { $EXONS{$gene} = $exon;                    }
         else                  { $EXONS{$gene} = $EXONS{$gene}."*,*".$exon;} # HAD TO DO THIS, AS SOMETIMES THE LAST COLUMN HAS ,S IN IT.
         # RECORD THE START OF THE EXON:
         if ($EXON_START{$exon} && $i == 0) # $i==0 MEANS IT IS THE FIRST TRANSCRIPT WE SEE FOR THIS EXON
         {
            $errormsg      = "ERROR: read_exon_gff: already know start of exon $exon\n";
            $errorcode     = 20; # ERRORCODE=20 (TESTED FOR)
            return(\%EXONS,\%EXON_START,\%EXON_FRAME,$errorcode,$errormsg);
         }
         $EXON_START{$exon}= $start;
         # IF $ignore_phase = 'yes', SET THE EXON PHASE TO 0:
         if ($ignore_phase eq 'yes') { $frame = 0; } 
         # RECORD THE EXON FRAME:
         if ($frame ne '0' && $frame ne '1' && $frame ne '2' && $frame ne 'zero')
         {
            $errormsg      = "ERROR: read_exon_gff: frame is $frame (line $line)\n";
            $errorcode     = 31; # ERRORCODE=31 (TESTED FOR)
            return(\%EXONS,\%EXON_START,\%EXON_FRAME,$errorcode,$errormsg); 
         }
         if ($frame eq '0') { $frame = "zero";} 
         if ($EXON_FRAME{$exon} && $i == 0) # $i==0 MEANS IT IS THE FIRST TRANSCRIPT WE SEE FOR THIS EXON
         {
            $errormsg      = "ERROR: read_exon_gff: already know frame of exon $exon\n";
            $errorcode     = 21; # ERRORCODE=21 (SHOULDN'T HAPPEN, SHOULD BE PICKED UP BY ERRORCODE=20)
            return(\%EXONS,\%EXON_START,\%EXON_FRAME,$errorcode,$errormsg);
         }
         $EXON_FRAME{$exon}= $frame;
      }
   }
   close(GFF);
 
   return(\%EXONS,\%EXON_START,\%EXON_FRAME,$errorcode,$errormsg);   
}
 
#------------------------------------------------------------------#
 
# TEST &get_exon_sequences
 
sub test_get_exon_sequences
{
   my $outputdir           = $_[0]; # DIRECTORY TO WRITE OUTPUT FILES INTO
   my $exon_gff;                    # EXON GFF FILE 
   my %SEQ                 = ();    # HASH TABLE WITH THE SEQUENCES OF SCAFFOLDS
   my $errorcode;                   # RETURNED AS 0 FROM A FUNCTION IF THERE IS NO ERROR
   my $errormsg;                    # RETURNED AS 'none' FROM A FUNCTION IF THERE IS NO ERROR
   my $EXONSEQ;                     # HASH TABLE TO STORE THE SEQUENCES OF EXONS
 
   ($exon_gff,$errorcode,$errormsg) = &make_filename($outputdir);
   if ($errorcode != 0) { ($errorcode,$errormsg) = &print_error($errormsg,$errorcode,0); } 
   open(EXON_GFF,">$exon_gff") || die "ERROR: test_get_exon_sequences: cannot open $exon_gff\n";
   print EXON_GFF "chr1\tchado\tCDS\t11\t20\t.\t+\t0\tID=gene1:exon1;Parent=gene1;\n";
   print EXON_GFF "chr1\tchado\tCDS\t31\t40\t.\t+\t0\tID=gene1:exon2;Parent=gene1;\n";
   close(EXON_GFF);
   $SEQ{"chr1"} = "AAAAATTTTTCCCCCGGGGGAAAAATTTTTCGCGCGCGCGAAAAA\n";
   # CALL &get_exon_sequences:
   ($EXONSEQ,$errorcode,$errormsg) = &get_exon_sequences($exon_gff,\%SEQ,$outputdir);
   if ($errorcode != 0) { print STDERR "ERROR: test_get_exon_sequences: failed test1 (errorcode $errorcode errormsg $errormsg)\n"; exit;} 
   if ($EXONSEQ->{"chr1#11#20#+#ID=gene1:exon1;Parent=gene1;"} ne "CCCCCGGGGG") { print STDERR "ERROR: test_get_exon_sequences: failed test1b\n"; exit;}
   if ($EXONSEQ->{"chr1#31#40#+#ID=gene1:exon2;Parent=gene1;"} ne "CGCGCGCGCG") { print STDERR "ERROR: test_get_exon_sequences: failed test1c\n"; exit;}
   system "rm -f $exon_gff\n";
 
   # TEST FOR ERRORCODE=8 (FEATURE IS NOT 'CDS'):
   ($exon_gff,$errorcode,$errormsg) = &make_filename($outputdir);
   if ($errorcode != 0) { ($errorcode,$errormsg) = &print_error($errormsg,$errorcode,0); } 
   open(EXON_GFF,">$exon_gff") || die "ERROR: test_get_exon_sequences: cannot open $exon_gff\n";
   print EXON_GFF "chr1\tchado\tCDS\t11\t20\t.\t+\t0\tID=gene1:exon1;Parent=gene1;\n";
   print EXON_GFF "chr1\tchado\tgene\t31\t40\t.\t+\t0\tID=gene1:exon2;Parent=gene1;\n";
   close(EXON_GFF);
   $SEQ{"chr1"} = "AAAAATTTTTCCCCCGGGGGAAAAATTTTTCGCGCGCGCGAAAAA\n";
   # CALL &get_exon_sequences:
   ($EXONSEQ,$errorcode,$errormsg) = &get_exon_sequences($exon_gff,\%SEQ,$outputdir);
   if ($errorcode != 8) { print STDERR "ERROR: test_get_exon_sequences: failed test2 (errorcode $errorcode errormsg $errormsg)\n"; exit;} 
   system "rm -f $exon_gff\n";
 
   # TEST FOR ERRORCODE=9 (EXON LENGTH IS NEGATIVE OR ZERO):
   ($exon_gff,$errorcode,$errormsg) = &make_filename($outputdir);
   if ($errorcode != 0) { ($errorcode,$errormsg) = &print_error($errormsg,$errorcode,0); } 
   open(EXON_GFF,">$exon_gff") || die "ERROR: test_get_exon_sequences: cannot open $exon_gff\n";
   print EXON_GFF "chr1\tchado\tCDS\t21\t20\t.\t+\t0\tID=gene1:exon1;Parent=gene1;\n";
   print EXON_GFF "chr1\tchado\tCDS\t31\t40\t.\t+\t0\tID=gene1:exon2;Parent=gene1;\n";
   close(EXON_GFF);
   $SEQ{"chr1"} = "AAAAATTTTTCCCCCGGGGGAAAAATTTTTCGCGCGCGCGAAAAA\n";
   # CALL &get_exon_sequences:
   ($EXONSEQ,$errorcode,$errormsg) = &get_exon_sequences($exon_gff,\%SEQ,$outputdir);
   if ($errorcode != 9) { print STDERR "ERROR: test_get_exon_sequences: failed test3 (errorcode $errorcode errormsg $errormsg)\n"; exit;} 
   system "rm -f $exon_gff\n";
 
   # TEST FOR ERRORCODE=10 (DO NOT KNOW SEQUENCE OF CHROMOSOME):
   ($exon_gff,$errorcode,$errormsg) = &make_filename($outputdir);
   if ($errorcode != 0) { ($errorcode,$errormsg) = &print_error($errormsg,$errorcode,0); } 
   open(EXON_GFF,">$exon_gff") || die "ERROR: test_get_exon_sequences: cannot open $exon_gff\n";
   print EXON_GFF "chr1\tchado\tCDS\t11\t20\t.\t+\t0\tID=gene1:exon1;Parent=gene1;\n";
   print EXON_GFF "chr1\tchado\tCDS\t31\t40\t.\t+\t0\tID=gene1:exon2;Parent=gene1;\n";
   close(EXON_GFF);
   %SEQ                    = ();
   # CALL &get_exon_sequences:
   ($EXONSEQ,$errorcode,$errormsg) = &get_exon_sequences($exon_gff,\%SEQ,$outputdir);
   if ($errorcode != 10) { print STDERR "ERROR: test_get_exon_sequences: failed test4 (errorcode $errorcode errormsg $errormsg)\n"; exit;} 
   system "rm -f $exon_gff\n";
 
   # MINUS STRAND EXONS:
   ($exon_gff,$errorcode,$errormsg) = &make_filename($outputdir);
   if ($errorcode != 0) { ($errorcode,$errormsg) = &print_error($errormsg,$errorcode,0); } 
   open(EXON_GFF,">$exon_gff") || die "ERROR: test_get_exon_sequences: cannot open $exon_gff\n";
   print EXON_GFF "chr1\tchado\tCDS\t11\t20\t.\t-\t0\tID=gene1:exon1;Parent=gene1;\n";
   print EXON_GFF "chr1\tchado\tCDS\t31\t40\t.\t-\t0\tID=gene1:exon2;Parent=gene1;\n";
   close(EXON_GFF);
   $SEQ{"chr1"} = "AAAAATTTTTAATAACCGCCAAAAATTTTTCACACACACAAAAAA\n";
   # CALL &get_exon_sequences:
   ($EXONSEQ,$errorcode,$errormsg) = &get_exon_sequences($exon_gff,\%SEQ,$outputdir);
   if ($errorcode != 0) { print STDERR "ERROR: test_get_exon_sequences: failed test6 (errorcode $errorcode errormsg $errormsg)\n"; exit;} 
   if ($EXONSEQ->{"chr1#11#20#-#ID=gene1:exon1;Parent=gene1;"} ne "AATAACCGCC") { print STDERR "ERROR: test_get_exon_sequences: failed test6b\n"; exit;}
   if ($EXONSEQ->{"chr1#31#40#-#ID=gene1:exon2;Parent=gene1;"} ne "CACACACACA") { print STDERR "ERROR: test_get_exon_sequences: failed test6c\n"; exit;}
   system "rm -f $exon_gff\n";
  
}
 
#------------------------------------------------------------------#
 
# GET THE DNA SEQUENCES FOR THE EXONS IN THE EXON GFF FILE:
# SUBROUTINE SYNOPSIS: get_exon_sequences: get the DNA sequences for exons in a gff file, given the scaffold fasta file
 
sub get_exon_sequences
{
   my $exon_gff            = $_[0]; # EXON GFF FILE
   my $SEQ                 = $_[1]; # HASH TABLE OF THE SEQUENCES OF SCAFFOLDS/CHROMOSOMES
   my $outputdir           = $_[2]; # DIRECTORY TO PUT OUTPUT FILES INTO
   my $errorcode           = 0;     # RETURNED AS 0 IF THERE IS NO ERROR
   my $errormsg            = "none";# RETURNED AS 'none' IF THERE IS NO ERROR
   my $line;                        # 
   my @temp;                        #
   my $name;                        # NAME OF THE CHROMOSOME/SCAFFOLD
   my $feature;                     # GFF FEATURE TYPE 
   my $exon_start;                  # START OF THE EXON
   my $exon_end;                    # END OF THE EXON
   my $strand;                      # STRAND OF THE EXON
   my $exon_name;                   # NAME OF THE EXON
   my $exon_length;                 # LENGTH OF THE EXON SEQUENCE
   my $seq;                         # NAME FOR THE SCAFFOLD/CHROMOSOME
   my $seq_length;                  # LENGTH OF THE SCAFFOLD/CHROMOSOME SEQUENCE
   my $exon_seq;                    # SEQUENCE OF THE EXON 
   my %EXONSEQ             = ();    # HASH TABLE OF EXON SEQUENCES 
 
   # READ IN THE EXON GFF FILE:
   open(GFF,"$exon_gff") || die "ERROR: get_exon_sequences: cannot open exon_gff $exon_gff.\n";
   while(<GFF>)
   {
      $line                = $_;
      chomp $line;
      @temp                = split(/\t+/,$line);
      $name                = $temp[0];
      $feature             = $temp[2];
      if ($feature ne 'CDS')
      {
         $errormsg         = "ERROR: get_exon_sequences: feature $feature line $line\n";
         $errorcode        = 8; # ERRORCODE=8 (TESTED FOR)
         return(\%EXONSEQ,$errorcode,$errormsg);
      }
      $exon_start          = $temp[3];
      $exon_end            = $temp[4];
      $strand              = $temp[6];
      $exon_name           = $temp[8];
      $exon_name           = $name."#".$exon_start."#".$exon_end."#".$strand."#".$exon_name; 
      $exon_length         = $exon_end - $exon_start + 1;
      if ($exon_length <= 0) 
      { 
         $errormsg         = "ERROR: get_exon_sequences: exon_length $exon_length: $line.\n"; 
         $errorcode        = 9; # ERRORCODE=9 (TESTED FOR)
         return(\%EXONSEQ,$errorcode,$errormsg);
      }
      if (!($SEQ->{$name}))    
      { 
         $errormsg         = "ERROR: get_exon_sequences: do not know sequence for name $name (exon_gff $exon_gff).\n"; 
         $errorcode        = 10; # ERRORCODE=10 (TESTED FOR)
         return(\%EXONSEQ,$errorcode,$errormsg);
      }
      $seq                 = $SEQ->{$name};
      $seq_length          = length($seq);
      if ($seq_length == 0)  
      { 
         $errormsg         = "ERROR: get_exon_sequences: seq_length $seq_length for name $name.\n"; 
         $errorcode        = 14; # ERRORCODE=14 (SHOULDN'T OCCUR, CAN'T TEST FOR)
         return(\%EXONSEQ,$errorcode,$errormsg);
      }
 
      # CHECK IF THE $exon_seq WILL BE OUTSIDE OF $name:
      if ($exon_end > $seq_length)
      {
         print STDERR "WARNING: cannot take $name $exon_start $exon_end because $name is $seq_length bp\n";
      }
      else
      {
         $exon_seq         = substr($seq,$exon_start-1,$exon_length);
         if ($exon_seq eq '') 
         { 
            $errormsg      = "ERROR: get_exon_sequences: exon_seq $exon_seq exon_start $exon_start exon_end $exon_end exon_length $exon_length seq_length $seq_length name $name line $line\n"; 
            $errorcode     = 17; # ERRORCODE=17 (SHOULDN'T OCCUR, CAN'T TEST FOR)
            return(\%EXONSEQ,$errorcode,$errormsg);
         }
         if ($exon_seq eq '' || $exon_seq eq 'none')
         {
            $errormsg      = "ERROR: get_exon_sequences: exon_seq $exon_seq\n";
            $errorcode     = 22; # ERRORCODE=22 (SHOULDN'T OCCUR, CAN'T TEST FOR)
            return(\%EXONSEQ,$errorcode,$errormsg);
         }
         # RECORD THE EXON SEQUENCE FOR $exon_name:
         if ($EXONSEQ{$exon_name})
         {
            $errormsg      = "ERROR: get_exon_sequences: already have sequence for $exon_name\n";
            $errorcode     = 19; # ERRORCODE=19 (SHOULDN'T OCCUR, CAN'T TEST FOR)
            return(\%EXONSEQ,$errorcode,$errormsg);
         } 
         $EXONSEQ{$exon_name} = $exon_seq;
      }
   }
   close(GFF);
 
   return(\%EXONSEQ,$errorcode,$errormsg);
}
 
#------------------------------------------------------------------#
 
# TEST &reverse_complement
 
sub test_reverse_complement
{
   my $errorcode;                   # RETURNED AS 0 IF THERE IS NO ERROR
   my $errormsg;                    # RETURNED AS 'none' IF THERE IS NO ERROR
   my $seq;                         # SEQUENCE TO FIND THE REVERSE COMPLEMENT OF
   my $revcomp;                     # THE REVERSE COMPLEMENT SEQUENCE
 
   ($revcomp,$errorcode,$errormsg) = &reverse_complement("AACCTTGG");
   if ($errorcode != 0 && $revcomp ne 'CCAAGGTT') { print STDERR "ERROR: test_reverse_complement: failed test1\n"; exit;}
   
   # TEST FOR ERRORCODE=15:
   ($revcomp,$errorcode,$errormsg) = &reverse_complement('');
   if ($errorcode != 15) { print STDERR "ERROR: test_reverse_complement: failed test2\n"; exit;} 
 
}
 
#------------------------------------------------------------------#
 
# GET THE REVERSE COMPLEMENT OF A SEQUENCE:
 
sub reverse_complement
{
   my $seq                 = $_[0]; ## SEQUENCE THAT WE WANT TO FIND THE COMPLEMENT OF
   my $complement          = "";    ## COMPLEMENT OF SEQUENCE $seq
   my $errorcode           = 0;     ## RETURNED AS 0 IF THERE IS NO ERROR
   my $errormsg            = "none";## RETURNED AS 'none' IF THERE IS NO ERROR
 
   if ($seq eq '') 
   {
      $errormsg            = "ERROR: find_complement: seq is $seq\n";
      $errorcode           = 15; # ERRORCODE=15 (TESTED FOR)
      return($complement,$errorcode,$errormsg);
   }
   $seq                    =~ tr/[a-z]/[A-Z]/; # CHANGE TO UPPERCASE
   $complement             = $seq;
   # SWAP As WITH Ts:
   $complement             =~ s/A/1/g;        # SUBSTITUTE '1' FOR 'A'
   $complement             =~ s/T/A/g;        # SUBSTITUTE 'A' FOR 'T'
   $complement             =~ s/1/T/g;        # SUBSTITUTE 'T' FOR '1'
   # SWAP Gs WITH Cs:
   $complement             =~ s/G/1/g;        # SUBSTITUTE '1' FOR 'G'
   $complement             =~ s/C/G/g;        # SUBSITTUTE 'G' FOR 'C'
   $complement             =~ s/1/C/g;        # SUBSTITUTE 'C' FOR '1'
   # FIND THE REVERSE COMPLEMENT:
   $complement             = reverse $complement;
   if ($complement eq '')
   {
      $errormsg            = "ERROR: reverse_complement: complement $complement (seq $seq)\n";
      $errorcode           = 16; # ERRORCODE=16 (SHOULDN'T HAPPEN, CAN'T TEST FOR)
      return($complement,$errorcode,$errormsg);
   }
 
   return($complement,$errorcode,$errormsg);
}
 
#------------------------------------------------------------------#
 
# TEST &make_exon_gff
 
sub test_make_exon_gff
{
   my $outputdir           = $_[0]; # DIRECTORY TO WRITE OUTPUT FILES IN
   my $errorcode;                   # RETURNED AS 0 FROM A FUNCTION IF THERE IS NO ERROR
   my $errormsg;                    # RETURNED AS 'none' FROM A FUNCTION IF THERE IS NO ERROR
   my $input_gff;                   # INPUT GFF FILE
   my $exon_gff;                    # EXON GFF FILE 
   my $expected_exon_gff;           # FILE WITH THE EXPECTED CONTENTS OF $exon_gff
   my $differences;                 # DIFFERENCES BETWEEN $exon_gff AND $expected_exon_gff
   my $length_differences;          # LENGTH OF $differences
   my $line;                        # 
 
   ($input_gff,$errorcode,$errormsg) = &make_filename($outputdir);
   if ($errorcode != 0) { ($errorcode,$errormsg) = &print_error($errormsg,$errorcode,0); } 
   open(INPUT_GFF,">$input_gff") || die "ERROR: test_make_exon_gff: cannot open $input_gff\n";
   print INPUT_GFF "##gff-version 3\n";
   print INPUT_GFF "##sequence-region Pk_strainH_chr01 1 838594\n";
   print INPUT_GFF "Pk_strainH_chr01	chado	gene	108	1421	.	+	.	ID=PKH_010010;isObsolete=false;isFminPartial;feature_id=222341;timelastmodified=10.12.2012+02:43:44+GMT\n";
   print INPUT_GFF "Pk_strainH_chr01	chado	CDS	108	758	.	+	0	ID=PKH_010010.1:exon:1;Parent=PKH_010010.1;isFminPartial;isObsolete=false;timelastmodified=10.12.2012+02:43:44+GMT;colour=7\n";
   print INPUT_GFF "Pk_strainH_chr01	chado	CDS	978	1421	.	+	0	ID=PKH_010010.1:exon:2;Parent=PKH_010010.1;isFminPartial;isObsolete=false;timelastmodified=10.12.2012+02:43:44+GMT;colour=7\n";
   print INPUT_GFF "Pk_strainH_chr01	chado	mRNA	108	1421	.	+	.	ID=PKH_010010.1;Parent=PKH_010010;isObsolete=false;isFminPartial;feature_id=222342;timelastmodified=10.12.2012+02:43:44+GMT;previous_systematic_id=PK00_0600c\n";
   print INPUT_GFF "Pk_strainH_chr01	chado	gene	5697	8325	.	+	.	ID=PKH_010020;isObsolete=false;feature_id=222257;timelastmodified=25.11.2012+01:47:32+GMT\n";
   print INPUT_GFF "Pk_strainH_chr01	chado	CDS	5697	6209	.	+	0	ID=PKH_010020.1:exon:1;Parent=PKH_010020.1;isObsolete=false;timelastmodified=25.11.2012+01:47:32+GMT;colour=7\n";
   print INPUT_GFF "Pk_strainH_chr01	chado	CDS	7891	8325	.	+	0	ID=PKH_010020.1:exon:2;Parent=PKH_010020.1;isObsolete=false;timelastmodified=25.11.2012+01:47:32+GMT;colour=7\n";
   print INPUT_GFF "Pk_strainH_chr01	chado	mRNA	5697	8325	.	+	.	ID=PKH_010020.1;Parent=PKH_010020;isObsolete=false;feature_id=222258;timelastmodified=25.11.2012+01:47:32+GMT;previous_systematic_id=PK9_3420w\n";
   print INPUT_GFF "Pk_strainH_chr01	chado	gene	8988	9593	.	-	.	ID=PKH_010030;isObsolete=false;feature_id=222118;isFmaxPartial;timelastmodified=12.01.2011+05:18:24+GMT\n";
   print INPUT_GFF "Pk_strainH_chr01	chado	CDS	8988	9593	.	-	0	ID=PKH_010030.1:exon:1;Parent=PKH_010030.1;isObsolete=false;timelastmodified=12.01.2011+05:18:24+GMT;colour=12\n";
   print INPUT_GFF "Pk_strainH_chr01	chado	mRNA	8988	9593	.	-	.	ID=PKH_010030.1;Parent=PKH_010030;isObsolete=false;feature_id=222119;timelastmodified=21.10.2007+01:51:33+BST;previous_systematic_id=PK4_2020w\n";
   print INPUT_GFF "Pk_strainH_chr01	chado	gap	10499	12498	.	+	.	ID=gap10499-12498:corrected;isObsolete=false;feature_id=222641;timelastmodified=21.10.2007+01:51:37+BST\n";
   print INPUT_GFF "Pk_strainH_chr01	chado	gene	13925	14846	.	-	.	ID=PKH_010040;isObsolete=false;feature_id=222402;timelastmodified=21.10.2007+01:51:35+BST\n";
   print INPUT_GFF "Pk_strainH_chr01	chado	CDS	14778	14846	.	-	0	ID=PKH_010040.1:exon:2;Parent=PKH_010040.1;isObsolete=false;timelastmodified=21.10.2007+01:51:35+BST;colour=8\n";
   print INPUT_GFF "Pk_strainH_chr01	chado	CDS	13925	14581	.	-	0	ID=PKH_010040.1:exon:1;Parent=PKH_010040.1;isObsolete=false;timelastmodified=21.10.2007+01:51:35+BST;colour=8\n";
   print INPUT_GFF "Pk_strainH_chr01	chado	mRNA	13925	14846	.	-	.	ID=PKH_010040.1;Parent=PKH_010040;isObsolete=false;feature_id=222403;timelastmodified=21.10.2007+01:51:35+BST\n";
   print INPUT_GFF "Pk_strainH_chr01	chado	gene	18394	18738	.	+	.	ID=PKH_010050;isObsolete=false;feature_id=222197;timelastmodified=21.10.2007+01:51:33+BST\n";
   print INPUT_GFF "Pk_strainH_chr01	chado	CDS	18394	18738	.	+	0	ID=PKH_010050.1:exon:1;Parent=PKH_010050.1;isObsolete=false;timelastmodified=21.10.2007+01:51:33+BST;colour=10\n";
   print INPUT_GFF "Pk_strainH_chr01	chado	mRNA	18394	18738	.	+	.	ID=PKH_010050.1;Parent=PKH_010050;isObsolete=false;feature_id=222198;timelastmodified=21.10.2007+01:51:33+BST;previous_systematic_id=PK7_0005w\n";
   print INPUT_GFF "Pk_strainH_chr01	chado	gene	22127	22613	.	-	.	ID=PKH_010080;isObsolete=false;feature_id=222375;timelastmodified=25.11.2012+01:16:50+GMT\n";
   print INPUT_GFF "Pk_strainH_chr01	chado	CDS	22404	22613	.	-	0	ID=PKH_010080.1:exon:2;Parent=PKH_010080.1;isObsolete=false;timelastmodified=25.11.2012+01:16:50+GMT;colour=10\n";
   print INPUT_GFF "Pk_strainH_chr01	chado	CDS	22127	22228	.	-	0	ID=PKH_010080.1:exon:1;Parent=PKH_010080.1;isObsolete=false;timelastmodified=25.11.2012+01:16:50+GMT;colour=10\n";
   print INPUT_GFF "Pk_strainH_chr01	chado	mRNA	22127	22613	.	-	.	ID=PKH_010080.1;Parent=PKH_010080;isObsolete=false;feature_id=222376;timelastmodified=25.11.2012+01:16:50+GMT\n";
   print INPUT_GFF "Pk_strainH_chr01	chado	gene	32144	34243	.	+	.	ID=PKH_010100;isObsolete=false;feature_id=222114;timelastmodified=21.10.2007+01:51:33+BST\n";
   print INPUT_GFF "Pk_strainH_chr01	chado	CDS	32144	34243	.	+	0	ID=PKH_010100.1:exon:1;Parent=PKH_010100.1;isObsolete=false;timelastmodified=21.10.2007+01:51:33+BST;colour=7\n";
   print INPUT_GFF "Pk_strainH_chr01	chado	mRNA	32144	34243	.	+	.	ID=PKH_010100.1;Parent=PKH_010100;isObsolete=false;feature_id=222115;timelastmodified=21.10.2007+01:51:33+BST;previous_systematic_id=PK7_0010w\n";
   print INPUT_GFF "Pk_strainH_chr01	chado	gene	35321	35393	.	-	.	ID=PKH_010102;isObsolete=false;feature_id=222634;timelastmodified=21.10.2007+01:51:37+BST\n";
   print INPUT_GFF "Pk_strainH_chr01	chado	CDS	35321	35393	.	-	0	ID=PKH_010102.1:exon:1;Parent=PKH_010102:tRNA;isObsolete=false;timelastmodified=21.10.2007+01:51:37+BST\n";
   print INPUT_GFF "Pk_strainH_chr01	chado	tRNA	35321	35393	.	-	.	ID=PKH_010102:tRNA;Parent=PKH_010102;isObsolete=false;feature_id=222635;product=term%3DtRNA+Alanine%3B;timelastmodified=07.09.2012+12:23:39+BST;comment=tRNA+Ala+anticodon+AGC%2C+Cove+score+51.85\n";
   close(INPUT_GFF);
   # MAKE A GFF FILE OF THE EXONS:
   ($exon_gff,$errorcode,$errormsg) = &make_exon_gff($input_gff,$outputdir);
   if ($errorcode != 0) { ($errorcode,$errormsg) = &print_error($errormsg,$errorcode,0); }
   # MAKE A FILE WITH THE EXPECTED CONTENTS OF $exon_gff:
   ($expected_exon_gff,$errorcode,$errormsg) = &make_filename($outputdir);
   if ($errorcode != 0) { ($errorcode,$errormsg) = &print_error($errormsg,$errorcode,0); } 
   open(EXPECTED,">$expected_exon_gff") || die "ERROR: test_make_exon_gff: cannot open $expected_exon_gff\n";
   print EXPECTED "Pk_strainH_chr01	chado	CDS	108	758	.	+	0	ID=PKH_010010.1:exon:1;Parent=PKH_010010.1;isFminPartial;isObsolete=false;timelastmodified=10.12.2012+02:43:44+GMT;colour=7\n";
   print EXPECTED "Pk_strainH_chr01	chado	CDS	978	1421	.	+	0	ID=PKH_010010.1:exon:2;Parent=PKH_010010.1;isFminPartial;isObsolete=false;timelastmodified=10.12.2012+02:43:44+GMT;colour=7\n";
   print EXPECTED "Pk_strainH_chr01	chado	CDS	5697	6209	.	+	0	ID=PKH_010020.1:exon:1;Parent=PKH_010020.1;isObsolete=false;timelastmodified=25.11.2012+01:47:32+GMT;colour=7\n";
   print EXPECTED "Pk_strainH_chr01	chado	CDS	7891	8325	.	+	0	ID=PKH_010020.1:exon:2;Parent=PKH_010020.1;isObsolete=false;timelastmodified=25.11.2012+01:47:32+GMT;colour=7\n";
   print EXPECTED "Pk_strainH_chr01	chado	CDS	8988	9593	.	-	0	ID=PKH_010030.1:exon:1;Parent=PKH_010030.1;isObsolete=false;timelastmodified=12.01.2011+05:18:24+GMT;colour=12\n";
   print EXPECTED "Pk_strainH_chr01	chado	CDS	13925	14581	.	-	0	ID=PKH_010040.1:exon:1;Parent=PKH_010040.1;isObsolete=false;timelastmodified=21.10.2007+01:51:35+BST;colour=8\n";
   print EXPECTED "Pk_strainH_chr01	chado	CDS	14778	14846	.	-	0	ID=PKH_010040.1:exon:2;Parent=PKH_010040.1;isObsolete=false;timelastmodified=21.10.2007+01:51:35+BST;colour=8\n";
   print EXPECTED "Pk_strainH_chr01	chado	CDS	18394	18738	.	+	0	ID=PKH_010050.1:exon:1;Parent=PKH_010050.1;isObsolete=false;timelastmodified=21.10.2007+01:51:33+BST;colour=10\n";
   print EXPECTED "Pk_strainH_chr01	chado	CDS	22127	22228	.	-	0	ID=PKH_010080.1:exon:1;Parent=PKH_010080.1;isObsolete=false;timelastmodified=25.11.2012+01:16:50+GMT;colour=10\n";
   print EXPECTED "Pk_strainH_chr01	chado	CDS	22404	22613	.	-	0	ID=PKH_010080.1:exon:2;Parent=PKH_010080.1;isObsolete=false;timelastmodified=25.11.2012+01:16:50+GMT;colour=10\n";
   print EXPECTED "Pk_strainH_chr01	chado	CDS	32144	34243	.	+	0	ID=PKH_010100.1:exon:1;Parent=PKH_010100.1;isObsolete=false;timelastmodified=21.10.2007+01:51:33+BST;colour=7\n";
   print EXPECTED "Pk_strainH_chr01	chado	CDS	35321	35393	.	-	0	ID=PKH_010102.1:exon:1;Parent=PKH_010102:tRNA;isObsolete=false;timelastmodified=21.10.2007+01:51:37+BST\n";
   close(EXPECTED); 
   $differences            = "";
   open(TEMP,"diff $exon_gff $expected_exon_gff |");
   while(<TEMP>)
   {
      $line                = $_;
      $differences         = $differences.$line;
   }
   close(TEMP);  
   $length_differences     = length($differences);
   if ($length_differences != 0) { print STDERR "ERROR: test_make_exon_gff: failed test1 (exon_gff $exon_gff expected_exon_gff $expected_exon_gff)\n"; exit;}
   system "rm -f $expected_exon_gff"; 
   system "rm -f $exon_gff";
   system "rm -f $input_gff"; 
 
}
 
#------------------------------------------------------------------#
 
# MAKE A GFF FILE OF THE EXONS:
 
sub make_exon_gff
{
   my $input_gff           = $_[0]; # INPUT GFF FILE
   my $outputdir           = $_[1]; # DIRECTORY TO PUT OUTPUT FILES INTO
   my $errorcode           = 0;     # RETURNED AS 0 IF THERE IS NO ERROR
   my $errormsg            = "none";# RETURNED AS 'none' IF THERE IS NO ERROR
   my $exon_gff;                    # GFF FILE OF EXONS
   my $cmd;                         # COMMAND TO RUN 
   my $line;                        # 
   my @temp;                        # 
   my $sorted_exon_gff;             # SORTED EXON GFF FILE 
 
   ($exon_gff,$errorcode,$errormsg) = &make_filename($outputdir);
   if ($errorcode != 0) { ($errorcode,$errormsg) = &print_error($errormsg,$errorcode,0); } 
   open(EXON_GFF,">$exon_gff") || die "ERROR: make_exon_gff: cannot open $exon_gff\n";
   open(INPUT_GFF,"$input_gff") || die "ERROR: make_exon_gff: cannot open $input_gff\n";
   while(<INPUT_GFF>)
   {
      $line                = $_;
      chomp $line;
      @temp                = split(/\t+/,$line);
      if ($#temp == 8)
      {
         if ($temp[2] eq 'CDS')
         {
            print EXON_GFF "$line\n";
         }
      }
   }
   close(INPUT_GFF);  
   close(EXON_GFF);
   # SORT THE EXON GFF FILE:
   ($sorted_exon_gff,$errorcode,$errormsg) = &make_filename($outputdir);
   if ($errorcode != 0) { ($errorcode,$errormsg) = &print_error($errormsg,$errorcode,0); } 
   $cmd                    = "sort  -k1,1 -k9,9 -k4,4n $exon_gff > $sorted_exon_gff"; 
   system "$cmd";
   system "rm -f $exon_gff";
 
   return($sorted_exon_gff,$errorcode,$errormsg);
}
 
#------------------------------------------------------------------#
 
# TEST &read_assembly
 
sub test_read_assembly
{
   my $outputdir           = $_[0]; # DIRECTORY WHERE WE CAN WRITE OUTPUT FILES
   my $random_number;               # RANDOM NUMBER TO USE IN TEMPORARY FILE NAMES
   my $assembly;                    # TEMPORARY ASSEMBLY FILE NAME 
   my $SEQ;                         # HASH TABLE WITH SEQUENCES OF SCAFFOLDS
   my $errorcode;                   # RETURNED AS 0 FROM A FUNCTION IF THERE IS NO ERROR
   my $errormsg;                    # RETURNED AS 'none' FROM A FUNCTION IF THERE IS NO ERROR 
 
   $random_number          = rand();
   $assembly               = $outputdir."/tmp".$random_number;
   open(ASSEMBLY,">$assembly") || die "ERROR: test_read_assembly: cannot open $assembly\n";
   print ASSEMBLY ">seq1\n";
   print ASSEMBLY "AAAAA\n";
   print ASSEMBLY ">seq2\n";
   print ASSEMBLY "AAAAATTTTT\n";
   print ASSEMBLY "\n";
   print ASSEMBLY ">seq3\n";
   print ASSEMBLY "AAAAA\n";
   print ASSEMBLY "TTTTT\n";
   print ASSEMBLY ">seq4\n";
   print ASSEMBLY ">seq5\n";
   print ASSEMBLY " AAA AA \n";
   close(ASSEMBLY);
   ($SEQ,$errorcode,$errormsg) = &read_assembly($assembly,'no');
   if ($SEQ->{'seq1'} ne 'AAAAA' || $SEQ->{'seq2'} ne 'AAAAATTTTT' || $SEQ->{'seq3'} ne 'AAAAATTTTT' || defined($SEQ->{'seq4'}) || 
       $SEQ->{'seq5'} ne 'AAAAA' || $errorcode != 0) 
   { 
      print STDERR "ERROR: test_read_assembly: failed test1\n"; 
      exit;
   }
   system "rm -f $assembly";
 
   # TEST ERRORCODE=4:
   $random_number          = rand();
   $assembly               = $outputdir."/tmp".$random_number;
   open(ASSEMBLY,">$assembly") || die "ERROR: test_assembly: cannot open $assembly\n";
   print ASSEMBLY ">seq1\n";
   print ASSEMBLY "AAAAA\n";
   print ASSEMBLY ">seq2\n";
   print ASSEMBLY "AAAAATTTTT\n";
   print ASSEMBLY ">seq1\n";
   print ASSEMBLY "AAAAA\n";
   close(ASSEMBLY);
   ($SEQ,$errorcode,$errormsg) = &read_assembly($assembly,'no');
   if ($errorcode != 4) { print STDERR "ERROR: test_assembly: failed test2\n"; exit;}
   system "rm -f $assembly";
 
   # TEST FOR ERRORCODE=5:
   $random_number          = rand();
   $assembly               = $outputdir."/tmp".$random_number;
   open(ASSEMBLY,">$assembly") || die "ERROR: test_assembly: cannot open $assembly\n";
   print ASSEMBLY "AAAAA\n";
   print ASSEMBLY ">seq2\n";
   print ASSEMBLY "AAAAATTTTT\n";
   print ASSEMBLY ">seq1\n";
   print ASSEMBLY "AAAAA\n";
   close(ASSEMBLY);
   ($SEQ,$errorcode,$errormsg) = &read_assembly($assembly,'no');
   if ($errorcode != 5) { print STDERR "ERROR: test_assembly: failed test2\n"; exit;}
   system "rm -f $assembly";
 
}
 
#------------------------------------------------------------------#
 
# READ IN THE ASSEMBLY FILE:
# SUBROUTINE SYNOPSIS: read_assembly(): read fasta file of scaffold sequences into a hash
 
sub read_assembly       
{
   my $input_assembly      = $_[0]; # THE INPUT ASSEMBLY FILE
   my $ignore_semicolons   = $_[1]; # SAYS WHETHER TO IGNORE SEMICOLONS AT THE END OF SCAFFOLD NAMES 
   my $errorcode           = 0;     # RETURNED AS 0 IF THERE IS NO ERROR
   my $errormsg            = "none";# RETURNED AS 'none' IF THERE IS NO ERROR
   my $line;                        # 
   my $scaffold            = "none";# NAME OF SCAFFOLD
   my @temp;                        # 
   my $seq;                         # SEQUENCE OF SCAFFOLD 
   my %SEQ                 = ();    # HASH TABLE FOR STORING SCAFFOLD SEQUENCE
 
   $seq                    = "";
   open(ASSEMBLY,"$input_assembly") || die "ERROR: read_assembly: cannot open $input_assembly\n";
   while(<ASSEMBLY>)
   {
      $line                = $_;
      chomp $line;
      if (substr($line,0,1) eq '>')
      {
         if ($seq ne "")
         {
            if ($scaffold eq 'none')
            {
               $errormsg   = "ERROR: read_assembly: do not know name of scaffold\n";
               $errorcode  = 5; # ERRORCODE=5 (TESTED FOR)
               return(\%SEQ,$errorcode,$errormsg);
            }
            if ($SEQ{$scaffold}) 
            {
               $errormsg   = "ERROR: read_assembly: already know sequence for scaffold $scaffold\n";
               $errorcode  = 4; # ERRORCODE=4 (TESTED FOR)
               return(\%SEQ,$errorcode,$errormsg);
            }
            $SEQ{$scaffold}= $seq;
            $seq           = "";
         }
         @temp             = split(/\s+/,$line);
         $scaffold         = $temp[0];
         $scaffold         = substr($scaffold,1,length($scaffold)-1);
         if ($ignore_semicolons eq 'yes') # IGNORE SEMICOLONS AT THE END OF SCAFFOLD NAMES
         {
            if (substr($scaffold,length($scaffold)-1,1) eq ';') { chop($scaffold);}
         }
      }
      else
      {
         $line             =~ s/\s+//g; # REMOVE SPACES
         $seq              = $seq.$line; 
      }
   }
   close(ASSEMBLY); 
   if ($seq ne "")
   {
      if ($scaffold eq 'none')
      {
         $errormsg         = "ERROR: read_assembly: do not know name of scaffold\n";
         $errorcode        = 5; # ERRORCODE=5 (TESTED FOR)
         return(\%SEQ,$errorcode,$errormsg);
      }
 
      if ($SEQ{$scaffold}) 
      {
         $errormsg         = "ERROR: read_assembly: already know sequence for scaffold $scaffold\n";
         $errorcode        = 4; # ERRORCODE=4 (TESTED FOR)
         return(\%SEQ,$errorcode,$errormsg);
      }
      $SEQ{$scaffold}      = $seq;
      $seq                 = "";
   }
 
   return(\%SEQ,$errorcode,$errormsg);
 
}
 
#------------------------------------------------------------------#
 
# SUBROUTINE TO MAKE A FILE NAME FOR A TEMPORARY FILE:
 
sub make_filename
{
   my $outputdir             = $_[0]; # OUTPUT DIRECTORY TO WRITE TEMPORARY FILE NAME TO
   my $found_name            = 0;     # SAYS WHETHER WE HAVE FOUND A FILE NAME YET
   my $filename              = "none";# NEW FILE NAME TO USE 
   my $errorcode             = 0;     # RETURNED AS 0 IF THERE IS NO ERROR
   my $errormsg              = "none";# RETURNED AS 'none' IF THERE IS NO ERROR
   my $poss_filename;                 # POSSIBLE FILE NAME TO USE
   my $random_number;                 # RANDOM NUMBER TO USE IN TEMPORARY FILE NAME
 
   while($found_name == 0)
   {
      $random_number         = rand();
      $poss_filename         = $outputdir."/tmp".$random_number;
      if (!(-e $poss_filename))
      {
         $filename           = $poss_filename;
         $found_name         = 1;
      } 
   } 
   if ($found_name == 0 || $filename eq 'none')
   {
      $errormsg              = "ERROR: make_filename: found_name $found_name filename $filename\n";
      $errorcode             = 6; # ERRORCODE=6 (SHOULDN'T HAPPEN, CAN'T TEST FOR)
      return($filename,$errorcode,$errormsg);
   }
 
   return($filename,$errorcode,$errormsg); 
}
 
#------------------------------------------------------------------#
 
# TEST &print_error
 
sub test_print_error
{
   my $errormsg;                    # RETURNED AS 'none' FROM A FUNCTION IF THERE WAS NO ERROR
   my $errorcode;                   # RETURNED AS 0 FROM A FUNCTION IF THERE WAS NO ERROR
 
   ($errormsg,$errorcode)  = &print_error(45,45,1);
   if ($errorcode != 12) { print STDERR "ERROR: test_print_error: failed test1\n"; exit;}
 
   ($errormsg,$errorcode)  = &print_error('My error message','My error message',1);
   if ($errorcode != 11) { print STDERR "ERROR: test_print_error: failed test2\n"; exit;}
 
   ($errormsg,$errorcode)  = &print_error('none',45,1);
   if ($errorcode != 13) { print STDERR "ERROR: test_print_error: failed test3\n"; exit;} 
 
   ($errormsg,$errorcode)  = &print_error('My error message', 0, 1);
   if ($errorcode != 13) { print STDERR "ERROR: test_print_error: failed test4\n"; exit;}
}
 
#------------------------------------------------------------------#
 
# PRINT OUT AN ERROR MESSAGE AND EXIT.
 
sub print_error
{
   my $errormsg            = $_[0]; # THIS SHOULD BE NOT 'none' IF AN ERROR OCCURRED.
   my $errorcode           = $_[1]; # THIS SHOULD NOT BE 0 IF AN ERROR OCCURRED.
   my $called_from_test    = $_[2]; # SAYS WHETHER THIS WAS CALLED FROM test_print_error OR NOT
 
   if ($errorcode =~ /[A-Z]/ || $errorcode =~ /[a-z]/) 
   { 
      if ($called_from_test == 1)
      {
         $errorcode = 11; $errormsg = "ERROR: print_error: the errorcode is $errorcode, should be a number.\n"; # ERRORCODE=11
         return($errormsg,$errorcode);
      }
      else 
      { 
         print STDERR "ERROR: print_error: the errorcode is $errorcode, should be a number.\n"; 
         exit;
      }
   }
 
   if (!($errormsg =~ /[A-Z]/ || $errormsg =~ /[a-z]/)) 
   { 
      if ($called_from_test == 1)
      {
         $errorcode = 12; $errormsg = "ERROR: print_error: the errormessage $errormsg does not seem to contain text.\n"; # ERRORCODE=12
         return($errormsg,$errorcode);
      }
      else
      {
         print STDERR "ERROR: print_error: the errormessage $errormsg does not seem to contain text.\n"; 
         exit;
      }
   }
 
   if    ($errormsg eq 'none' || $errorcode == 0) 
   { 
      if ($called_from_test == 1)
      {
         $errorcode = 13; $errormsg = "ERROR: print_error: errormsg $errormsg, errorcode $errorcode.\n"; # ERRORCODE=13
         return($errormsg,$errorcode);
      }
      else 
      {
         print STDERR "ERROR: print_error: errormsg $errormsg, errorcode $errorcode.\n"; 
         exit;
      }
   }
   else                                           
   { 
      print STDERR "$errormsg"; 
      exit;                                                      
   } 
 
   return($errormsg,$errorcode);
}
 
#------------------------------------------------------------------#
