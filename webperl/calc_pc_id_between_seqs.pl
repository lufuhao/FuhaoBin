#!/usr/bin/env perl
#https://gist.github.com/avrilcoghlan/5311008

=head1 NAME

    calc_pc_id_between_seqs.pl

=head1 SYNOPSIS
 
    calc_pc_id_between_seqs.pl input_fasta output outputdir ggsearch
        where input_fasta is the input fasta file of sequences,
              output is the output file of percent identities,
              outputdir is the output directory for writing output files,
              ggsearch is the path to the ggsearch executable.

=head1 DESCRIPTION

    This script takes an input fasta file of protein sequences (<input_fasta>), and uses ggsearch
    from the fasta package to find the percent identity between each pair of sequences
    over a global pairwise alignment, and writes this to the output file (<output>).

=head1 VERSION
  
    Perl script last edited 4-Apr-2013.

=head1 CONTACT

    alc@sanger.ac.uk (Avril Coghlan)

=cut
 
# 
# Perl script calc_pc_id_between_seqs.pl
# Written by Avril Coghlan (alc@sanger.ac.uk)
# 4-Apr-13.
# Last edited 4-Apr-2013.
# SCRIPT SYNOPSIS: calc_pc_id_between_seqs.pl: calculate (global) percent identity between each pair of protein sequences in a fasta file.
#
#------------------------------------------------------------------#
 
# CHECK IF THERE ARE THE CORRECT NUMBER OF COMMAND-LINE ARGUMENTS:
 
use strict;
use warnings;
 
my $num_args               = $#ARGV + 1;
if ($num_args != 4)
{
    print "Usage of calc_pc_id_between_seqs.pl\n\n";
    print "perl calc_pc_id_between_seqs.pl <input_fasta> <output> <outputdir> <ggsearch>\n";
    print "where <input_fasta> is the input fasta file,\n";
    print "      <output> is the output file,\n";
    print "      <outputdir> is the output directory for writing output files,\n";
    print "      <ggsearch> is the path to the ggsearch executable\n";
    print "For example, >perl calc_pc_id_between_seqs.pl PTRK_training.pep PTRK_training.pep.id\n";
    exit;
}
 
# FIND THE PATH TO THE INPUT FASTA FILE:                     
 
my $input_fasta            = $ARGV[0];
 
# FIND THE PATH TO THE OUTPUT FILE:
 
my $output                 = $ARGV[1];
 
# FIND THE DIRECTORY TO USE FOR OUTPUT FILES:      
 
my $outputdir              = $ARGV[2];
 
# FIND THE PATH TO THE GGSEARCH EXECUTABLE:       
 
my $ggsearch               = $ARGV[3];
 
#------------------------------------------------------------------#
 
# TEST SUBROUTINES: 
 
my $PRINT_TEST_DATA        = 0;   # SAYS WHETHER TO PRINT DATA USED DURING TESTING.
&test_read_assembly($outputdir);
&test_calc_pc_id_between_two_seqs($outputdir,$ggsearch); 
&test_calc_pc_id_between_seqs($outputdir,$ggsearch);
&test_print_error;
print STDERR "Finished tests, now running main code...\n";
 
#------------------------------------------------------------------#
 
# RUN THE MAIN PART OF THE CODE:
 
&run_main_program($outputdir,$input_fasta,$output,$ggsearch);
 
print STDERR "FINISHED.\n";
 
#------------------------------------------------------------------#
 
# RUN THE MAIN PART OF THE CODE:
 
sub run_main_program
{
   my $outputdir           = $_[0]; # DIRECTORY TO PUT OUTPUT FILES IN.
   my $input_fasta         = $_[1]; # THE INPUT FASTA FILE  
   my $output              = $_[2]; # THE OUTPUT FILE 
   my $ggsearch            = $_[3]; # PATH TO THE GGSEARCH EXECUTABLE            
   my $errorcode;                   # RETURNED AS 0 IF THERE IS NO ERROR.
   my $errormsg;                    # RETURNED AS 'none' IF THERE IS NO ERROR. 
   my $SEQ;                         # HASH TABLE OF SEQUENCES 
 
   # READ IN THE INPUT FILE OF SEQUENCES:
   ($SEQ,$errorcode,$errormsg) = &read_assembly($input_fasta);
   if ($errorcode != 0) { ($errorcode,$errormsg) = &print_error($errormsg,$errorcode,0); }
   
   # CALCULATE THE PERCENT IDENTITY BETWEEN EACH PAIR OF SEQUENCES:
   ($errorcode,$errormsg)  = &calc_pc_id_between_seqs($SEQ,$output,$outputdir,$ggsearch,0);
   if ($errorcode != 0) { ($errorcode,$errormsg) = &print_error($errormsg,$errorcode,0); }
 
}
 
#------------------------------------------------------------------#
 
# TEST &calc_pc_id_between_seqs
 
sub test_calc_pc_id_between_seqs
{
   my $outputdir           = $_[0]; # DIRECTORY TO WRITE OUTPUT FILES IN
   my $ggsearch            = $_[1]; # PATH TO THE ggsearch EXECUTABLE 
   my $errorcode;                   # RETURNED AS 0 FROM A FUNCTION IF THERE IS NO ERROR
   my $errormsg;                    # RETURNED AS 'none' FROM A FUNCTION IF THERE IS NO ERROR
   my $output;                      # THE OUTPUT FILE
   my %SEQ                 = ();    # HASH TABLE OF SEQUENCES
   my $expected_output;             # FILE CONTAINING EXPECTED RESULTS OF $output
   my @temp;                        #  
   my $differences;                 # DIFFERENCES BETWEEN $output AND $expected_output
   my $length_differences;          # LENGTH OF $differences
   my $line;                        #  
 
   ($output,$errorcode,$errormsg) = &make_filename($outputdir);
   if ($errorcode != 0) { ($errorcode,$errormsg) = &print_error($errormsg,$errorcode,0); }
   @temp                   = split(/\//,$output);
   $output                 = $temp[$#temp];
   $SEQ{'seq1'}            = "AAAAGAAAAATTTTTCTTTTCACACACACACGGGGTGGGGACCACCACCACCCCCACCCCGAGAGAGAGACCCCGGCCCCAAAACCAAAATTTTCTCTTT";
   $SEQ{'seq2'}            = "AAAAGAAAAATTTTTCTTTTCACACACACACGGGGTGGGGACCACCACCACCCCCACCCCGAGAGAGAGACCCCGGCCCCAAAACCAAAATTTTCTCTTT";
   $SEQ{'seq3'}            = "AAAAGAAAAATTTTTCTTTTCACACACACACGGGGTTGGGACCACCACCACCCCCACCCCGAGAGAGAGACCCCGGCCCCAAAACCAAAATTTTCTCTTT";
   # CALCULATE THE PERCENT IDENTITY BETWEEN EACH PAIR OF SEQUENCES:
   ($errorcode,$errormsg)  = &calc_pc_id_between_seqs(\%SEQ,$output,$outputdir,$ggsearch,1);
   if ($errorcode != 0) { print STDERR "ERROR: test_calc_pc_id_between_seqs: failed test1\n"; exit;}
   ($expected_output,$errorcode,$errormsg) = &make_filename($outputdir);
   if ($errorcode != 0) { ($errorcode,$errormsg) = &print_error($errormsg,$errorcode,0); }
   open(EXPECTED,">$expected_output") || die "ERROR: test_calc_pc_id_between_seqs: failed test1\n";
   print EXPECTED "seq1 seq2 100.0\n";
   print EXPECTED "seq1 seq3 99.0\n";
   print EXPECTED "seq2 seq3 99.0\n";
   close(EXPECTED);
   $differences            = "";
   open(TEMP,"diff $output $expected_output |");
   while(<TEMP>)
   {
      $line                = $_;
      $differences         = $differences.$line;
   }
   close(TEMP);  
   $length_differences     = length($differences);
   if ($length_differences != 0) { print STDERR "ERROR: test_calc_pc_id_between_seqs: failed test1 (output $output expected_output $expected_output)\n"; exit;}
   system "rm -f $output";
   system "rm -f $expected_output"; 
}
 
#------------------------------------------------------------------#
 
# CALCULATE THE PERCENT IDENTITY BETWEEN EACH PAIR OF SEQUENCES:
# 
sub calc_pc_id_between_seqs
{
   my $SEQ                 = $_[0]; # HASH TABLE OF SEQUENCES
   my $output              = $_[1]; # OUTPUT FILE
   my $outputdir           = $_[2]; # DIRECTORY TO WRITE OUTPUT FILES IN
   my $ggsearch            = $_[3]; # PATH TO THE ggsearch EXECUTABLE 
   my $testing             = $_[4]; # SAYS WHETHER THIS IS BEING CALLED FROM A TESTING FUNCTION 
   my $errorcode           = 0;     # RETURNED AS 0 IF THERE IS NO ERROR
   my $errormsg            = "none";# RETURNED AS 'none' IF THERE IS NO ERROR
   my $name1;                       # NAME OF THE FIRST SEQUENCE
   my $name2;                       # NAME OF THE SECOND SEQUENCE
   my $seq1;                        # THE FIRST SEQUENCE
   my $seq2;                        # THE SECOND SEQUENCE
   my $pc_id;                       # PERCENT IDENTITY BETWEEN TWO SEQUENCES
   my $pair1;                       # PAIR OF SEQUENCES
   my $pair2;                       # PAIR OF SEQUENCES
   my %SEEN                = ();    # HASH TABLE OF PAIRS OF SEQUENCES THAT WE HAVE SEEN ALREADY
   my $cnt                 = 0;     # NUMBER OF SEQUENCES SEEN 
   my $num_seqs;                    # TOTAL NUMBER OF SEQUENCES
 
   # OPEN THE OUTPUT FILE:
   $output                 = $outputdir."/".$output;
   open(OUTPUT,">$output") || die "ERROR: calc_pc_id_between_seqs: cannot open $output\n";
   
   # CALCULATE THE PERCENT IDENTITY BETWEEN EACH PAIR OF SEQUENCES:
   $num_seqs               = keys %{$SEQ};
   foreach $name1 (sort keys %{$SEQ}) # SORTING THE KEYS MAKES IT EASIER TO TEST
   {
      $seq1                = $SEQ->{$name1};
      $cnt++;
      if ($testing == 0) { print STDERR "Calculating percent identity for $name1 (sequence $cnt of $num_seqs)...\n";}
      foreach $name2 (sort keys %{$SEQ}) # SORTING THE KEYS MAKES IT EASIER TO TEST
      {
         if ($name1 ne $name2)
         {
            $pair1         = $name1."_".$name2;
            $pair2         = $name2."_".$name1;
            if (!($SEEN{$pair1}) && !($SEEN{$pair2}))
            {
               $SEEN{$pair1} = 1;
               $SEEN{$pair2} = 1;
               $seq2       = $SEQ->{$name2};
               # FIND THE PERCENT IDENTITY BETWEEN THESE TWO SEQUENCES:
               ($pc_id,$errorcode,$errormsg) = &calc_pc_id_between_two_seqs($name1,$name2,$seq1,$seq2,$outputdir,$ggsearch);
               if ($errorcode != 0) { ($errorcode,$errormsg) = &print_error($errormsg,$errorcode,0); }
               # WRITE OUT TO THE OUTPUT FILE:
               print OUTPUT "$name1 $name2 $pc_id\n";
            }
         }
      }
   }
   close(OUTPUT);
 
   return($errorcode,$errormsg);
}
 
#------------------------------------------------------------------#
 
# TEST &calc_pc_id_between_two_seqs
 
sub test_calc_pc_id_between_two_seqs
{
   my $outputdir           = $_[0]; # DIRECTORY TO WRITE OUTPUT FILES IN
   my $ggsearch            = $_[1]; # THE SSERACH EXECUTABLE
   my $seq1;                        # THE FIRST SEQUENCE
   my $seq2;                        # THE SECOND SEQUENCE 
   my $pc_id;                       # PERCENT IDENTITY BETWEEN THE TWO SEQUENCES
   my $errorcode;                   # RETURNED AS 0 FROM A FUNCTION IF THERE IS NO ERROR
   my $errormsg;                    # RETURNED AS 'none' FROM A FUNCTION IF THERE IS NO ERROR
   
   $seq1                   = "AAAAGAAAAATTTTTCTTTTCACACACACACGGGGTGGGGACCACCACCACCCCCACCCCGAGAGAGAGACCCCGGCCCCAAAACCAAAATTTTCTCTTT";
   $seq2                   = "AAAAGAAAAATTTTTCTTTTCACACACACACGGGGTGGGGACCACCACCACCCCCACCCCGAGAGAGAGACCCCGGCCCCAAAACCAAAATTTTCTCTTT";
   ($pc_id,$errorcode,$errormsg) = &calc_pc_id_between_two_seqs("seq1","seq2",$seq1,$seq2,$outputdir,$ggsearch);
   if ($errorcode != 0 || $pc_id != 100) { print STDERR "ERROR: test_calc_pc_id_between_two_seqs: failed test1 (errorcode $errorcode errormsg $errormsg pc_id $pc_id)\n"; exit;}
 
   $seq1                   = "AAAAGAAAAATTTTTCTTTTCACACACACACGGGGTGGGGACCACCACCACCCCCACCCCGAGAGAGAGACCCCGGCCCCAAAACCAAAATTTTCTCTTT";
   $seq2                   = "AAAAGAAAAATTTTTCTTTTCACACACACACGGGGTTGGGACCACCACCACCCCCACCCCGAGAGAGAGACCCCGGCCCCAAAACCAAAATTTTCTCTTT";
   ($pc_id,$errorcode,$errormsg) = &calc_pc_id_between_two_seqs("seq1","seq2",$seq1,$seq2,$outputdir,$ggsearch);
   if ($errorcode != 0 || $pc_id != 99) { print STDERR "ERROR: test_calc_pc_id_between_two_seqs: failed test2 (errorcode $errorcode errormsg $errormsg pc_id $pc_id)\n"; exit;}
   
   $seq1                   = "AAAAGAAAAATTTTTCTTTTCACACACACACGGGGTGGGGACCACCACCACCCCCACCCCGAGAGAGAGACCCCGGCCCCAAAACCAAAATTTTCTCTTT";
   $seq2                   = "AAAAGAACAATTTTTCTTTTCACACACACACGGGGTTGGGACCACCACCACCCCCACCCCGAGAGAGAGACCCCGGCCCCAAAACCAAAATTTTCTCTTT";
   ($pc_id,$errorcode,$errormsg) = &calc_pc_id_between_two_seqs("seq1","seq2",$seq1,$seq2,$outputdir,$ggsearch);
   if ($errorcode != 0 || $pc_id != 98) { print STDERR "ERROR: test_calc_pc_id_between_two_seqs: failed test3 (errorcode $errorcode errormsg $errormsg pc_id $pc_id)\n"; exit;}
  
   $seq1                   = "AAAAGAAAAATTTTTCTTTTCACACACACACGGGGTGGGGACCACCACCACCCCCACCCCGAGAGAGAGACCCCGGCCCCAAAACCAAAATTTTCTCTTT";
   $seq2                   = "AAAAGAACAATTTTTCTTTTCACACACACACGGGGTTGGGACCACCACCACCGCCACCCCGAGAGAGAGACCCCGGCCCCAAAACCAAAATTTTCTCTTT";
   ($pc_id,$errorcode,$errormsg) = &calc_pc_id_between_two_seqs("seq1","seq2",$seq1,$seq2,$outputdir,$ggsearch);
   if ($errorcode != 0 || $pc_id != 97) { print STDERR "ERROR: test_calc_pc_id_between_two_seqs: failed test4 (errorcode $errorcode errormsg $errormsg pc_id $pc_id)\n"; exit;}
 
   $seq1                   = "AAAAGAAAAATTTTTCTTTTCACACACACACGGGGTGGGGACCACCACCACCCCCACCCCGAGAGAGAGACCCCGGCCCCAAAACCAAAATTTTCTCTTT";
   $seq2                   = "AAAAGAACAATTTTTCTTTTCACACACACACGGGGTTGGGACCACCACCACCGCCACCCCGAGAGAGAGACCCCGGCCCCAAAACCAAAATTATCTCTTT";
   ($pc_id,$errorcode,$errormsg) = &calc_pc_id_between_two_seqs("seq1","seq2",$seq1,$seq2,$outputdir,$ggsearch);
   if ($errorcode != 0 || $pc_id != 96) { print STDERR "ERROR: test_calc_pc_id_between_two_seqs: failed test5 (errorcode $errorcode errormsg $errormsg pc_id $pc_id)\n"; exit;}
  
   $seq1                   = "AAAAGAAAAATTTTTCTTTTCACACACACACGGGGTGGGGACCACCACCACCCCCACCCCGAGAGAGAGACCCCGGCCCCAAAACCAAAATTTTCTCTTT";
   $seq2                   = "AAAAGAACAATTTTTCTTTTCACACACACACGGGGTTGGGACCACCACCACCGCCACCCCGAGAGAGAGACCCCGGCCCCAAAACCAAAATTATCTCTCT";
   ($pc_id,$errorcode,$errormsg) = &calc_pc_id_between_two_seqs("seq1","seq2",$seq1,$seq2,$outputdir,$ggsearch);
   if ($errorcode != 0 || $pc_id != 95) { print STDERR "ERROR: test_calc_pc_id_between_two_seqs: failed test6 (errorcode $errorcode errormsg $errormsg pc_id $pc_id)\n"; exit;}
}
 
#------------------------------------------------------------------#
 
# FIND THE PERCENT IDENTITY BETWEEN TWO SEQUENCES:
 
sub calc_pc_id_between_two_seqs
{
   my $name1               = $_[0]; # NAME OF THE FIRST SEQUENCE
   my $name2               = $_[1]; # NAME OF THE SECOND SEQUENCE  
   my $seq1                = $_[2]; # THE FIRST SEQUENCE
   my $seq2                = $_[3]; # THE SECOND SEQUENCE
   my $outputdir           = $_[4]; # DIRECTORY TO WRITE OUTPUT FILES IN
   my $ggsearch            = $_[5]; # PATH TO THE ggsearch EXECUTABLE
   my $pc_id               = -100;  # PERCENT IDENTITY BETWEEN THE TWO SEQUENCES
   my $errorcode           = 0;     # RETURNED AS 0 IF THERE IS NO ERROR
   my $errormsg            = "none";# RETURNED AS 'none' IF THERE IS NO ERROR
   my $seqfile1;                    # FILE CONTAINING $seq1
   my $seqfile2;                    # FILE CONTAINING $seq2 
   my $ggsearch_output;             # FILE CONTAINING THE ggsearch OUTPUT 
   my $cmd;                         # COMMAND TO RUN, TO RUN ggsearch
   my $line;                        # 
   my @temp;                        # 
 
   # PUT THE SEQUENCES INTO TEMPORARY FILES:
   ($seqfile1,$errorcode,$errormsg) = &make_filename($outputdir);
   if ($errorcode != 0) { ($errorcode,$errormsg) = &print_error($errormsg,$errorcode,0); }
   open(SEQFILE1,">$seqfile1") || die "ERROR: calc_pc_id_between_two_seqs: cannot open $seqfile1\n";
   print SEQFILE1 ">$name1\n";
   print SEQFILE1 "$seq1\n";
   close(SEQFILE1);
   ($seqfile2,$errorcode,$errormsg) = &make_filename($outputdir);
   if ($errorcode != 0) { ($errorcode,$errormsg) = &print_error($errormsg,$errorcode,0); }
   open(SEQFILE2,">$seqfile2") || die "ERROR: calc_pc_id_between_two_seqs: cannot open $seqfile2\n";
   print SEQFILE2 ">$name2\n";
   print SEQFILE2 "$seq2\n"; 
   close(SEQFILE2);
 
   # CALCULATE THE PERCENT IDENTITY BETWEEN THE TWO SEQUENCES USING ggsearch:
   ($ggsearch_output,$errorcode,$errormsg) = &make_filename($outputdir);
   if ($errorcode != 0) { ($errorcode,$errormsg) = &print_error($errormsg,$errorcode,0); }
   $cmd                    = "$ggsearch $seqfile1 $seqfile2 > $ggsearch_output";
   system "$cmd"; 
   open(GGSEARCH,"$ggsearch_output") || die "ERROR: calc_pc_id_between_two_seqs: cannot open $ggsearch_output\n";
   while(<GGSEARCH>)
   {
      $line                = $_;
      chomp $line;
      if ($line =~ /\% identity/)
      {
         # eg. global/global (N-W) score: -158; 18.3% identity (45.9% similar) in 987 aa overlap (1-912:1-930)
         @temp             = split(/\s+/,$line);
         if ($pc_id == -100) # THERE MAY BE MULTIPLE ALIGNMENTS IN THE FILE, BUT WE WANT TO JUST TAKE THE TOP (BEST) ALIGNMENT
         {
            $pc_id         = $temp[4]; # 18.3%
            if (substr($pc_id,length($pc_id)-1,1) eq '%') { chop($pc_id);}
            else
            {
               $errormsg   = "ERROR: calc_pc_id_between_two_seqs: pc_id $pc_id\n";
               $errorcode  = 2; # ERRORCODE=2 (THIS SHOULDN'T HAPPEN, SO CAN'T TEST FOR)
               return($pc_id,$errorcode,$errormsg);
            }
         }
      } 
      elsif ($line =~ /!! No sequences with E\(\) < 10/) # THE SEQUENCES CAN'T BE ALIGNED USING ggsearch
      {
         $pc_id            = 0.0;
      }
   }
   close(GGSEARCH); 
   if ($pc_id == -100)
   {
      $errormsg            = "ERROR: calc_pc_id_between_two_seqs: did not find percent identity in ggsearch output $ggsearch_output (seqfile1 $seqfile1 seqfile2 $seqfile2)\n";
      $errorcode           = 1; # ERRORCODE=1 (THIS SHOULD'T HAPPEN, SO CAN'T TEST FOR)
      return($pc_id,$errorcode,$errormsg);
   }
 
   # DELETE TEMPORARY FILES:
   system "rm -f $seqfile1";
   system "rm -f $seqfile2";
   system "rm -f $ggsearch_output";
 
   return($pc_id,$errorcode,$errormsg);
   
}
 
#------------------------------------------------------------------#
 
# TEST &read_assembly
 
sub test_read_assembly
{
   my $outputdir           = $_[0]; # DIRECTORY WHERE WE CAN WRITE OUTPUT FILES
   my $assembly;                    # TEMPORARY ASSEMBLY FILE NAME 
   my $SEQ;                         # HASH TABLE WITH SEQUENCES OF SCAFFOLDS
   my $errorcode;                   # RETURNED AS 0 FROM A FUNCTION IF THERE IS NO ERROR
   my $errormsg;                    # RETURNED AS 'none' FROM A FUNCTION IF THERE IS NO ERROR 
 
   ($assembly,$errorcode,$errormsg) = &make_filename($outputdir);
   if ($errorcode != 0) { ($errorcode,$errormsg) = &print_error($errormsg,$errorcode,0); }
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
   ($SEQ,$errorcode,$errormsg) = &read_assembly($assembly);
   if ($SEQ->{'seq1'} ne 'AAAAA' || $SEQ->{'seq2'} ne 'AAAAATTTTT' || $SEQ->{'seq3'} ne 'AAAAATTTTT' || defined($SEQ->{'seq4'}) || 
       $SEQ->{'seq5'} ne 'AAAAA' || $errorcode != 0) 
   { 
      print STDERR "ERROR: test_read_assembly: failed test1\n"; 
      exit;
   }
   system "rm -f $assembly";
 
   # TEST ERRORCODE=4:
   ($assembly,$errorcode,$errormsg) = &make_filename($outputdir);
   if ($errorcode != 0) { ($errorcode,$errormsg) = &print_error($errormsg,$errorcode,0); }
   open(ASSEMBLY,">$assembly") || die "ERROR: test_assembly: cannot open $assembly\n";
   print ASSEMBLY ">seq1\n";
   print ASSEMBLY "AAAAA\n";
   print ASSEMBLY ">seq2\n";
   print ASSEMBLY "AAAAATTTTT\n";
   print ASSEMBLY ">seq1\n";
   print ASSEMBLY "AAAAA\n";
   close(ASSEMBLY);
   ($SEQ,$errorcode,$errormsg) = &read_assembly($assembly);
   if ($errorcode != 4) { print STDERR "ERROR: test_assembly: failed test2\n"; exit;}
   system "rm -f $assembly";
 
   # TEST FOR ERRORCODE=5:
   ($assembly,$errorcode,$errormsg) = &make_filename($outputdir);
   if ($errorcode != 0) { ($errorcode,$errormsg) = &print_error($errormsg,$errorcode,0); }
   open(ASSEMBLY,">$assembly") || die "ERROR: test_assembly: cannot open $assembly\n";
   print ASSEMBLY "AAAAA\n";
   print ASSEMBLY ">seq2\n";
   print ASSEMBLY "AAAAATTTTT\n";
   print ASSEMBLY ">seq1\n";
   print ASSEMBLY "AAAAA\n";
   close(ASSEMBLY);
   ($SEQ,$errorcode,$errormsg) = &read_assembly($assembly);
   if ($errorcode != 5) { print STDERR "ERROR: test_assembly: failed test2\n"; exit;}
   system "rm -f $assembly";
 
}
 
#------------------------------------------------------------------#
 
# READ IN THE ASSEMBLY FILE:
# SUBROUTINE SYNOPSIS: read_assembly(): read fasta file of scaffold sequences into a hash
 
sub read_assembly       
{
   my $input_assembly      = $_[0]; # THE INPUT ASSEMBLY FILE
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
      $errorcode             = 6; # ERRORCODE=6 
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
