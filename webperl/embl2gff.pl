#!/usr/bin/env perl
#https://gist.github.com/avrilcoghlan/5049468
 
=head1 NAME

    embl_to_gff.pl

=head1 SYNOPSIS
 
    embl_to_gff.pl input_embl output_gff outputdir ratt
        where input_embl is the input embl file,
              output_gff is the output gff file,
              outputdir is the output directory for writing output files,
              ratt says whether the embl files are from ratt (yes/no).

=head1 DESCRIPTION

    This script takes an input embl file (<input_embl>), and converts it to 
    gff format, and writes the output file (<output_gff>) in directory
    <outputdir>.

=head1 VERSION
  
    Perl script last edited 27-Feb-2013.

=head1 CONTACT

    alc@sanger.ac.uk (Avril Coghlan)

=cut
 
# 
# Perl script embl_to_gff.pl
# Written by Avril Coghlan (alc@sanger.ac.uk)
# 27-Jan-13.
# Last edited 27-Jan-2013.
# SCRIPT SYNOPSIS: embl_to_gff.pl: convert an embl file to a gff file.
#
#------------------------------------------------------------------#
 
# CHECK IF THERE ARE THE CORRECT NUMBER OF COMMAND-LINE ARGUMENTS:
 
use strict;
use warnings;
 
my $num_args               = $#ARGV + 1;
if ($num_args != 4)
{
    print "Usage of embl_to_gff.pl\n\n";
    print "perl embl_to_gff.pl <input_embl> <output_gff> <outputdir> <ratt>\n";
    print "where <input_embl> is the input embl file,\n";
    print "      <output_gff> is the output gff file,\n";
    print "      <outputdir> is the output directory for writing output files,\n";
    print "      <ratt> says whether the embl files are from ratt (yes/no)\n";
    print "For example, >perl embl_to_gff.pl my.embl my.gff\n";
    print "/lustre/scratch108/parasites/alc/StrongyloidesTophat/PTRK/SrattiCurated yes\n";
    exit;
}
 
# FIND THE PATH TO THE INPUT EMBL FILE:                     
 
my $input_embl             = $ARGV[0];
 
# FIND THE PATH TO THE OUTPUT GFF FILE:
 
my $output_gff             = $ARGV[1];
 
# FIND THE DIRECTORY TO USE FOR OUTPUT FILES:      
 
my $outputdir              = $ARGV[2];
 
# FIND OUT WHETHER THE EMBL FILES ARE FROM RATT:
 
my $ratt                   = $ARGV[3];
 
#------------------------------------------------------------------#
 
# TEST SUBROUTINES: 
 
my $PRINT_TEST_DATA        = 0;   # SAYS WHETHER TO PRINT DATA USED DURING TESTING.
&test_convert_embl_to_gff($outputdir);
&test_print_error;
print STDERR "Finished tests, now running main code...\n";
 
#------------------------------------------------------------------#
 
# RUN THE MAIN PART OF THE CODE:
 
&run_main_program($outputdir,$input_embl,$output_gff,$ratt);
 
print STDERR "FINISHED.\n";
 
#------------------------------------------------------------------#
 
# RUN THE MAIN PART OF THE CODE:
 
sub run_main_program
{
   my $outputdir           = $_[0]; # DIRECTORY TO PUT OUTPUT FILES IN.
   my $input_embl          = $_[1]; # THE INPUT EMBL FILE  
   my $output_gff          = $_[2]; # THE OUTPUT GFF FILE 
   my $ratt                = $_[3]; # SAYS WHETHER THE EMBL FILES ARE FROM RATT
   my $errorcode;                   # RETURNED AS 0 IF THERE IS NO ERROR.
   my $errormsg;                    # RETURNED AS 'none' IF THERE IS NO ERROR. 
 
   # READ IN THE INPUT EMBL FILE, AND MAKE THE OUTPUT GFF FILE:
   ($errorcode,$errormsg)  = &convert_embl_to_gff($outputdir,$input_embl,$output_gff,$ratt);
   if ($errorcode != 0) { ($errorcode,$errormsg) = &print_error($errormsg,$errorcode,0); }
}
 
#------------------------------------------------------------------#
 
# READ IN THE INPUT EMBL FILE, AND MAKE THE OUTPUT GFF FILE:
 
sub convert_embl_to_gff    
{
   my $outputdir           = $_[0]; # DIRECTORY TO PUT OUTPUT FILES IN
   my $input_embl          = $_[1]; # INPUT EMBL FILE
   my $output_gff          = $_[2]; # OUTPUT GFF FILE
   my $ratt                = $_[3]; # SAYS WHETHER THE EMBL FILES ARE FROM RATT
   my $line;                        # 
   my @temp;                        # 
   my $errorcode           = 0;     # RETURNED AS 0 IF THERE IS NO ERROR
   my $errormsg            = "none";# RETURNED AS 'none' IF THERE IS NO ERROR
   my $name                = "none";# NAME OF THE SEQUENCE
   my $exons;                       # EXONS IN A GENE
   my $strand;                      # STRAND OF THE GENE 
   my @exons;                       # EXONS IN A GENE
   my $i;                           #  
   my $exon;                        # AN EXON IN A GENE
   my $start;                       # START OF AN EXON
   my $end;                         # END OF AN EXON
   my $gene;                        # NAME OF THE GENE
   my %CNT                 = ();    # NUMBER OF TIMES THAT WE HAVE SEEN A PARTICULAR GENE NAME
   my $cnt;                         # NUMBER OF TIMES THAT WE HAVE SEEN A PARTICULAR GENE NAME 
   my $first_exon;                  # FIRST EXON IN A GENE
   my $last_exon;                   # LAST EXON IN A GENE
   my $first_exon_start;            # START POINT OF THE FIRST EXON IN A GENE
   my $last_exon_end;               # END OF THE LAST EXON IN A ENE 
   my $exon_name;                   # NAME OF AN EXON 
   my $transcript_name;             # NAME OF A TRANSCRIPT 
   my $found_seq_start     = 0;     # SAYS WE HAVE FOUND THE START OF THE SEQUENCE IN THE EMBL FILE
   my $prev_line           = "";    # PREVIOUS LINE IN THE EMBL FILE 
   my $prev_prev_line      = "";    # LINE BEFORE $prev_line IN THE EMBL FILE 
   my @temp2;                       #  
 
   # CHECK $ratt IS 'yes' OR 'no':
   if ($ratt ne 'yes' && $ratt ne 'no')
   {
      $errormsg            = "ERROR: convert_embl_to_gff: ratt should be yes/no (is $ratt)\n";
      $errorcode           = 4; # ERRORCODE=4 
      return($errorcode,$errormsg);
   }
 
   # OPEN THE OUTPUT GFF FILE:
   $output_gff             = $outputdir."/".$output_gff;
   open(OUTPUT,">$output_gff") || die "ERROR: convert_embl_to_gff: cannot open $output_gff\n";
 
   # READ IN THE INPUT EMBL FILE:
   open(EMBL,"$input_embl") || die "ERROR: convert_embl_to_gff: cannot open input_embl $input_embl\n";
   while(<EMBL>)
   {
      $line                = $_;
      if ($found_seq_start == 0)
      {
         chomp $line;
         @temp             = split(/\s+/,$line);
         if (substr($line,0,2) eq 'ID')
         {
            $name          = $temp[1]; # eg. Transfer1.PTRK.contig.00177.164151.final
            if ($ratt eq 'yes')
            {
               @temp       = split(/\./,$name);
               $name       = "";
               for ($i = 1; $i <= $#temp-1; $i++) { $name = $name.".".$temp[$i];}
               $name       = substr($name,1,length($name)-1);
            }
            # IF THERE IS A SEMI-COLON AT THE END OF THE NAME, REMOVE IT:
            if (substr($name,length($name)-1,1) eq ';') { chop($name);} 
         } 
         elsif (substr($line,0,2) eq 'FT' && $temp[1] eq 'CDS')
         {
            # eg. FT   CDS             complement(join(309654..309719,309724..310158))
            $exons         = $temp[2]; # eg. complement(join(309654..309719,309724..310158))
            if ($exons =~ /complement/) { $strand = "-";}
            else                        { $strand = "+";}
         }
         elsif (substr($line,0,2) eq 'FT' && $line =~ /gene_orientation/)
         {
            # eg. FT                   /systematic_id="gene_id 0 ; sequence ratti_train136 ; gene_orientation +" 
            # or 
            #     FT                   /systematic_id="gene_id 0 ; sequence
            #     FT                   T1.PTRK.scaffold.00008.1453643 ; gene_orientation +" 
            if (!($line =~ /systematic_id/)) # THE INFORMATION ON THE GENE IS SPREAD OVER TWO LINES
            {
               if ($prev_line =~ /systematic_id/) # THE INFORMATION ON THE GENE IS SPREAD OVER TWO LINES
               { 
                  # eg. prev_line IS: FT                   /systematic_id="gene_id 0 ; sequence
                  # line is:          FT                   T1.PTRK.scaffold.00008.1453643 ; gene_orientation +" 
                  @temp2   = split(/FT\s+/,$line); 
                  $line    = $temp2[1];                   # eg. T1.PTRK.scaffold.00008.1453643 ; gene_orientation +" 
                  $line    = $prev_line." ".$line;        
                  # eg. FT                   /systematic_id="gene_id 0 ; sequence T1.PTRK.scaffold.00008.1453643 ; gene_orientation +"  
               }
               elsif ($prev_prev_line =~ /systematic_id/) # THE INFORMATION ON THE GENE IS SPREAD OVER THREE LINES
               {
                  # IF THE INFORMATION ON THE GENE IS SPREAD OVER THREE LINES
                  # eg. prev_prev_line IS: FT                   /systematic_id="gene_id 0 ;
                  #     prev_line IS:      FT                   T1.PTRK.scaffold.00010.1389350.embl sequence
                  #     line IS:           FT                   T1.PTRK.scaffold.00010.1389350 ; gene_orientation +" 
                  @temp2   = split(/FT\s+/,$line);
                  $line    = $temp2[1];
                  @temp2   = split(/FT\s+/,$prev_line);
                  $prev_line = $temp2[1];
                  $line    = $prev_prev_line." ".$prev_line." ".$line;
               }
               else
               {
                  $errormsg= "ERROR: convert_embl_to_gff: line $line prev_line $prev_line prev_prev_line $prev_prev_line\n";
                  $errorcode= 1; # ERRORCODE=1 (TESTED FOR)
                  return($errorcode,$errormsg);
               } 
            }
            @temp          = split(/sequence\s+/,$line);
            $gene          = $temp[1]; # eg. T1.PTRK.scaffold.00010.1389350 ; gene_orientation +" 
            @temp          = split(/\s+/,$gene);
            $gene          = $temp[0]; # eg. T1.PTRK.scaffold.00010.1389350
            # IF THERE IS A SEMI-COLON AT THE END OF THE GENE NAME, IGNORE IT:
            if (substr($gene,length($gene)-1,1) eq ';') { chop($gene);} 
            if (!($CNT{$gene})) { $CNT{$gene} = 1; }
            else                { $CNT{$gene}++;   }
            $cnt           = $CNT{$gene};
            $gene          = $gene."_".$cnt;
            # PRINT OUT THE EXONS:
            if ($exons =~ /\(/)
            {
               @temp       = split(/\(/,$exons);
               if ($temp[1] =~ /join/) { $exons = $temp[2];} # eg. 309654..309719,309724..310158)
               else                    { $exons = $temp[1];} # eg. 308719..308826
               if (substr($exons,length($exons)-1,1) eq ')') { chop($exons);}
               if (substr($exons,length($exons)-1,1) eq ')') { chop($exons);}
            }
            @exons         = split(/\,/,$exons);
            for ($i = 0; $i <= $#exons; $i++)
            {
               $exon       = $exons[$i];
               @temp       = split(/\.\./,$exon);
               $start      = $temp[0];
               $end        = $temp[1];
               $exon_name  = $gene."_".$i;
               $transcript_name = $gene."_mRNA";
               # PRINT OUT GFF LINE:
               print OUTPUT "$name\tsource\tCDS\t$start\t$end\t.\t$strand\t.\t$exon_name\;Parent=$transcript_name\n";
            }
            # PRINT OUT GFF LINE FOR 'mRNA' AND 'gene' FEATURES:
            $first_exon    = $exons[0];
            @temp          = split(/\.\./,$first_exon);
            $first_exon_start = $temp[0];
            $last_exon     = $exons[$#exons]; 
            @temp          = split(/\.\./,$last_exon);
            $last_exon_end = $temp[1];
            print OUTPUT "$name\tsource\tgene\t$first_exon_start\t$last_exon_end\t.\t$strand\t.\t$gene\n";
            print OUTPUT "$name\tsource\tmRNA\t$first_exon_start\t$last_exon_end\t.\t$strand\t.\t$transcript_name\;Parent=$gene\n";
         }
         elsif (substr($line,0,2) eq 'SQ')
         {
            $found_seq_start  = 1;
         }
      }
      $prev_prev_line      = $prev_line;
      $prev_line           = $line;
      chomp ($prev_line);
   }
   close(EMBL);
   close(OUTPUT);
 
   return($errorcode,$errormsg);
}
 
#------------------------------------------------------------------#
 
# TEST &convert_embl_to_gff    
 
sub test_convert_embl_to_gff      
{
   my $outputdir           = $_[0]; # DIRECTORY TO WRITE OUTPUT FILES IN
   my $input_embl;                  # NAME OF EMBL FILE
   my $output_gff;                  # NAME OF OUTPUT GFF FILE
   my $errorcode           = 0;     # RETURNED AS 0 BY A FUNCTION IF THERE IS NO ERROR
   my $errormsg            = "none";# RETURNED AS 'none' BY A FUNCTION IF THERE IS NO ERROR 
   my $expected_output_gff;         # FILE CONTAINING THE EXPECTED CONTENTS OF $output_gff      
   my $differences;                 # DIFFERENES BETWEEN $output_gff AND $expected_output_gff
   my $length_differences;          # LENGTH OF $differences
   my $line;                        # 
 
   ($output_gff,$errorcode,$errormsg) = &make_filename($outputdir);
   if ($errorcode != 0) { ($errorcode,$errormsg) = &print_error($errormsg,$errorcode,0); } 
   ($input_embl,$errorcode,$errormsg) = &make_filename($outputdir);
   if ($errorcode != 0) { ($errorcode,$errormsg) = &print_error($errormsg,$errorcode,0); } 
   open(EMBL,">$input_embl") || die "ERROR: test_convert_embl_to_gff: cannot open $input_embl\n";
   print EMBL "ID                   Transfer1.PTRK.scaffold.00024.740435.final ; ; ; ; ; 740435 BP.\n";
   print EMBL "FH   Key             Location/Qualifiers\n";
   print EMBL "FH\n";                 
   print EMBL "FT   CDS             complement(join(309654..309719,309724..310158))\n";
   print EMBL "FT                   /systematic_id=\"gene_id 0 ; sequence ratti_train136 ; gene_orientation +\"\n";
   print EMBL "FT                   /colour=7\n";
   print EMBL "FT   CDS             complement(308719..308826)\n";
   print EMBL "FT                   /systematic_id=\"gene_id 0 ; sequence ratti_train136 ; gene_orientation +\"\n";
   print EMBL "FT                   /colour=7\n";
   print EMBL "SQ   Sequence 740435 BP; 276856 A; 85528 C; 86034 G; 273093 T; 18924 other;\n";
   print EMBL "ttagtttttc cagaagcttt agcttttcct ccttttctac gtcctgacat attaataata          60\n";
   print EMBL "aagttataaa atgagtaagt tatttaaaaa ataaaatgtt taaattattt tgttatatag         120\n";
   close(EMBL); 
   ($errorcode,$errormsg)  = &convert_embl_to_gff($outputdir,$input_embl,$output_gff,'no');
   if ($errorcode != 0) { print STDERR "ERROR: test_convert_embl_to_gff: failed test1 (errorcode $errorcode errormsg $errormsg)\n"; exit;}
   ($expected_output_gff,$errorcode,$errormsg) = &make_filename($outputdir);
   if ($errorcode != 0) { ($errorcode,$errormsg) = &print_error($errormsg,$errorcode,0); } 
   open(EXPECTED,">$expected_output_gff") || die "ERROR: test_convert_embl_to_gff: cannot open $expected_output_gff\n";
   print EXPECTED "Transfer1.PTRK.scaffold.00024.740435.final\tsource\tCDS\t309654\t309719\t.\t-\t.\tratti_train136_1_0\;Parent=ratti_train136_1_mRNA\n";
   print EXPECTED "Transfer1.PTRK.scaffold.00024.740435.final\tsource\tCDS\t309724\t310158\t.\t-\t.\tratti_train136_1_1\;Parent=ratti_train136_1_mRNA\n";
   print EXPECTED "Transfer1.PTRK.scaffold.00024.740435.final\tsource\tgene\t309654\t310158\t.\t-\t.\tratti_train136_1\n";
   print EXPECTED "Transfer1.PTRK.scaffold.00024.740435.final\tsource\tmRNA\t309654\t310158\t.\t-\t.\tratti_train136_1_mRNA\;Parent=ratti_train136_1\n";
   print EXPECTED "Transfer1.PTRK.scaffold.00024.740435.final\tsource\tCDS\t308719\t308826\t.\t-\t.\tratti_train136_2_0\;Parent=ratti_train136_2_mRNA\n";
   print EXPECTED "Transfer1.PTRK.scaffold.00024.740435.final\tsource\tgene\t308719\t308826\t.\t-\t.\tratti_train136_2\n";
   print EXPECTED "Transfer1.PTRK.scaffold.00024.740435.final\tsource\tmRNA\t308719\t308826\t.\t-\t.\tratti_train136_2_mRNA\;Parent=ratti_train136_2\n";
   close(EXPECTED);               
   $differences            = "";
   open(TEMP,"diff $output_gff $expected_output_gff |");
   while(<TEMP>)
   {
      $line                = $_;
      $differences         = $differences.$line;
   }
   close(TEMP);  
   $length_differences     = length($differences);
   if ($length_differences != 0) { print STDERR "ERROR: test_convert_embl_to_gff: failed test1 (output_gff $output_gff expected_output_gff $expected_output_gff)\n"; exit;}
   system "rm -f $input_embl";
   system "rm -f $output_gff";
   system "rm -f $expected_output_gff"; 
 
   # TEST ERRORCODE=4 ($ratt IS NOT yes/no):
   ($output_gff,$errorcode,$errormsg) = &make_filename($outputdir);
   if ($errorcode != 0) { ($errorcode,$errormsg) = &print_error($errormsg,$errorcode,0); } 
   ($input_embl,$errorcode,$errormsg) = &make_filename($outputdir);
   if ($errorcode != 0) { ($errorcode,$errormsg) = &print_error($errormsg,$errorcode,0); } 
   open(EMBL,">$input_embl") || die "ERROR: test_convert_embl_to_gff: cannot open $input_embl\n";
   print EMBL "ID                   Transfer1.PTRK.scaffold.00024.740435.final ; ; ; ; ; 740435 BP.\n";
   print EMBL "FH   Key             Location/Qualifiers\n";
   print EMBL "FH\n";                 
   print EMBL "FT   CDS             complement(join(309654..309719,309724..310158))\n";
   print EMBL "FT                   /systematic_id=\"gene_id 0 ; sequence ratti_train136 ; gene_orientation +\"\n";
   print EMBL "FT                   /colour=7\n";
   print EMBL "FT   CDS             complement(308719..308826)\n";
   print EMBL "FT                   /systematic_id=\"gene_id 0 ; sequence ratti_train136 ; gene_orientation +\"\n";
   print EMBL "FT                   /colour=7\n";
   print EMBL "SQ   Sequence 740435 BP; 276856 A; 85528 C; 86034 G; 273093 T; 18924 other;\n";
   print EMBL "ttagtttttc cagaagcttt agcttttcct ccttttctac gtcctgacat attaataata          60\n";
   print EMBL "aagttataaa atgagtaagt tatttaaaaa ataaaatgtt taaattattt tgttatatag         120\n";
   close(EMBL); 
   ($errorcode,$errormsg)  = &convert_embl_to_gff($outputdir,$input_embl,$output_gff,'hello');
   if ($errorcode != 4) { print STDERR "ERROR: test_convert_embl_to_gff: failed test2\n"; exit;}
   system "rm -f $input_embl";
   system "rm -f $output_gff";
 
   # TEST $ratt='yes':
   ($output_gff,$errorcode,$errormsg) = &make_filename($outputdir);
   if ($errorcode != 0) { ($errorcode,$errormsg) = &print_error($errormsg,$errorcode,0); } 
   ($input_embl,$errorcode,$errormsg) = &make_filename($outputdir);
   if ($errorcode != 0) { ($errorcode,$errormsg) = &print_error($errormsg,$errorcode,0); } 
   open(EMBL,">$input_embl") || die "ERROR: test_convert_embl_to_gff: cannot open $input_embl\n";
   print EMBL "ID                   Transfer1.PTRK.scaffold.00024.740435.final ; ; ; ; ; 740435 BP.\n";
   print EMBL "FH   Key             Location/Qualifiers\n";
   print EMBL "FH\n";                 
   print EMBL "FT   CDS             complement(join(309654..309719,309724..310158))\n";
   print EMBL "FT                   /systematic_id=\"gene_id 0 ; sequence ratti_train136 ; gene_orientation +\"\n";
   print EMBL "FT                   /colour=7\n";
   print EMBL "FT   CDS             complement(308719..308826)\n";
   print EMBL "FT                   /systematic_id=\"gene_id 0 ; sequence ratti_train136 ; gene_orientation +\"\n";
   print EMBL "FT                   /colour=7\n";
   print EMBL "SQ   Sequence 740435 BP; 276856 A; 85528 C; 86034 G; 273093 T; 18924 other;\n";
   print EMBL "ttagtttttc cagaagcttt agcttttcct ccttttctac gtcctgacat attaataata          60\n";
   print EMBL "aagttataaa atgagtaagt tatttaaaaa ataaaatgtt taaattattt tgttatatag         120\n";
   close(EMBL); 
   ($errorcode,$errormsg)  = &convert_embl_to_gff($outputdir,$input_embl,$output_gff,'yes');
   if ($errorcode != 0) { print STDERR "ERROR: test_convert_embl_to_gff: failed test3\n"; exit;}
   ($expected_output_gff,$errorcode,$errormsg) = &make_filename($outputdir);
   if ($errorcode != 0) { ($errorcode,$errormsg) = &print_error($errormsg,$errorcode,0); } 
   open(EXPECTED,">$expected_output_gff") || die "ERROR: test_convert_embl_to_gff: cannot open $expected_output_gff\n";
   print EXPECTED "PTRK.scaffold.00024.740435\tsource\tCDS\t309654\t309719\t.\t-\t.\tratti_train136_1_0\;Parent=ratti_train136_1_mRNA\n";
   print EXPECTED "PTRK.scaffold.00024.740435\tsource\tCDS\t309724\t310158\t.\t-\t.\tratti_train136_1_1\;Parent=ratti_train136_1_mRNA\n";
   print EXPECTED "PTRK.scaffold.00024.740435\tsource\tgene\t309654\t310158\t.\t-\t.\tratti_train136_1\n";
   print EXPECTED "PTRK.scaffold.00024.740435\tsource\tmRNA\t309654\t310158\t.\t-\t.\tratti_train136_1_mRNA\;Parent=ratti_train136_1\n";
   print EXPECTED "PTRK.scaffold.00024.740435\tsource\tCDS\t308719\t308826\t.\t-\t.\tratti_train136_2_0\;Parent=ratti_train136_2_mRNA\n";
   print EXPECTED "PTRK.scaffold.00024.740435\tsource\tgene\t308719\t308826\t.\t-\t.\tratti_train136_2\n";
   print EXPECTED "PTRK.scaffold.00024.740435\tsource\tmRNA\t308719\t308826\t.\t-\t.\tratti_train136_2_mRNA\;Parent=ratti_train136_2\n";
   close(EXPECTED);               
   $differences            = "";
   open(TEMP,"diff $output_gff $expected_output_gff |");
   while(<TEMP>)
   {
      $line                = $_;
      $differences         = $differences.$line;
   }
   close(TEMP);  
   $length_differences     = length($differences);
   if ($length_differences != 0) { print STDERR "ERROR: test_convert_embl_to_gff: failed test3 (output_gff $output_gff expected_output_gff $expected_output_gff)\n"; exit;}
   system "rm -f $input_embl";
   system "rm -f $output_gff";
   system "rm -f $expected_output_gff"; 
 
   # EXAMPLE LIKE THE EMBL FILES MADE BY THE SGAs:
   ($output_gff,$errorcode,$errormsg) = &make_filename($outputdir);
   if ($errorcode != 0) { ($errorcode,$errormsg) = &print_error($errormsg,$errorcode,0); } 
   ($input_embl,$errorcode,$errormsg) = &make_filename($outputdir);
   if ($errorcode != 0) { ($errorcode,$errormsg) = &print_error($errormsg,$errorcode,0); } 
   open(EMBL,">$input_embl") || die "ERROR: test_convert_embl_to_gff: cannot open $input_embl\n";
   print EMBL "ID   PTRK.scaffold.00008.1453643; SV 1; linear; genomic DNA; HTG; INV; 1453643 BP.\n";
   print EMBL "FH   Key             Location/Qualifiers\n";
   print EMBL "FH\n";
   print EMBL "FT   CDS             join(115993..116307,116360..116413)\n";
   print EMBL "FT                   /systematic_id=\"gene_id 0 ; sequence\n";
   print EMBL "FT                   T1.PTRK.scaffold.00008.1453643 ; gene_orientation +\"\n";
   print EMBL "FT                   /colour=1\n";
   print EMBL "FT   CDS             join(535391..535554,535599..537591,537642..537971)\n";
   print EMBL "FT                   /systematic_id=\"gene_id 0 ; sequence\n";
   print EMBL "FT                   T2.PTRK.scaffold.00008.1453643.embl ; gene_orientation +\"\n";
   print EMBL "FT                   /colour=1\n";
   print EMBL "SQ   Sequence 1453643 BP; 541339 A; 172091 C; 169611 G; 546324 T; 24278 other;\n";
   print EMBL "atctaaaagc ttcaaaaccg gctttaaagt cattgtgcag caaatagaaa ctaaaaaata        60\n";
   close(EMBL); 
   ($errorcode,$errormsg)  = &convert_embl_to_gff($outputdir,$input_embl,$output_gff,'no');
   if ($errorcode != 0) { print STDERR "ERROR: test_convert_embl_to_gff: failed test4\n"; exit;}
   ($expected_output_gff,$errorcode,$errormsg) = &make_filename($outputdir);
   if ($errorcode != 0) { ($errorcode,$errormsg) = &print_error($errormsg,$errorcode,0); } 
   open(EXPECTED,">$expected_output_gff") || die "ERROR: test_convert_embl_to_gff: cannot open $expected_output_gff\n";
   print EXPECTED "PTRK.scaffold.00008.1453643\tsource\tCDS\t115993\t116307\t.\t+\t.\tT1.PTRK.scaffold.00008.1453643_1_0;Parent=T1.PTRK.scaffold.00008.1453643_1_mRNA\n";
   print EXPECTED "PTRK.scaffold.00008.1453643\tsource\tCDS\t116360\t116413\t.\t+\t.\tT1.PTRK.scaffold.00008.1453643_1_1;Parent=T1.PTRK.scaffold.00008.1453643_1_mRNA\n";
   print EXPECTED "PTRK.scaffold.00008.1453643\tsource\tgene\t115993\t116413\t.\t+\t.\tT1.PTRK.scaffold.00008.1453643_1\n";
   print EXPECTED "PTRK.scaffold.00008.1453643\tsource\tmRNA\t115993\t116413\t.\t+\t.\tT1.PTRK.scaffold.00008.1453643_1_mRNA;Parent=T1.PTRK.scaffold.00008.1453643_1\n";
   print EXPECTED "PTRK.scaffold.00008.1453643\tsource\tCDS\t535391\t535554\t.\t+\t.\tT2.PTRK.scaffold.00008.1453643.embl_1_0;Parent=T2.PTRK.scaffold.00008.1453643.embl_1_mRNA\n";
   print EXPECTED "PTRK.scaffold.00008.1453643\tsource\tCDS\t535599\t537591\t.\t+\t.\tT2.PTRK.scaffold.00008.1453643.embl_1_1;Parent=T2.PTRK.scaffold.00008.1453643.embl_1_mRNA\n";
   print EXPECTED "PTRK.scaffold.00008.1453643\tsource\tCDS\t537642\t537971\t.\t+\t.\tT2.PTRK.scaffold.00008.1453643.embl_1_2;Parent=T2.PTRK.scaffold.00008.1453643.embl_1_mRNA\n";
   print EXPECTED "PTRK.scaffold.00008.1453643\tsource\tgene\t535391\t537971\t.\t+\t.\tT2.PTRK.scaffold.00008.1453643.embl_1\n";
   print EXPECTED "PTRK.scaffold.00008.1453643\tsource\tmRNA\t535391\t537971\t.\t+\t.\tT2.PTRK.scaffold.00008.1453643.embl_1_mRNA;Parent=T2.PTRK.scaffold.00008.1453643.embl_1\n";
   close(EXPECTED);
   $differences            = "";
   open(TEMP,"diff $output_gff $expected_output_gff |");
   while(<TEMP>)
   {
      $line                = $_;
      $differences         = $differences.$line;
   }
   close(TEMP);  
   $length_differences     = length($differences);
   if ($length_differences != 0) { print STDERR "ERROR: test_convert_embl_to_gff: failed test4 (output_gff $output_gff expected_output_gff $expected_output_gff)\n"; exit;}
   system "rm -f $input_embl";
   system "rm -f $output_gff";
   system "rm -f $expected_output_gff"; 
 
   # ANOTHER EXAMPLE LIKE THE EMBL FILES MADE BY THE SGAs:
   ($output_gff,$errorcode,$errormsg) = &make_filename($outputdir);
   if ($errorcode != 0) { ($errorcode,$errormsg) = &print_error($errormsg,$errorcode,0); } 
   ($input_embl,$errorcode,$errormsg) = &make_filename($outputdir);
   if ($errorcode != 0) { ($errorcode,$errormsg) = &print_error($errormsg,$errorcode,0); } 
   open(EMBL,">$input_embl") || die "ERROR: test_convert_embl_to_gff: cannot open $input_embl\n";
   print EMBL "ID   PTRK.scaffold.00010.1389350; SV 1; linear; genomic DNA; HTG; INV; 1389350 BP.\n";
   print EMBL "FH   Key             Location/Qualifiers\n";
   print EMBL "FH\n";
   print EMBL "FT   CDS             join(33177..33365,33583..34248,34731..35867)\n";
   print EMBL "FT                   /systematic_id=\"gene_id 0 ;\n";
   print EMBL "FT                   T1.PTRK.scaffold.00010.1389350.embl sequence\n";
   print EMBL "FT                   T1.PTRK.scaffold.00010.1389350 ; gene_orientation +\"\n";
   print EMBL "FT                   /colour=1\n";
   print EMBL "FT   CDS             join(33177..33365,33583..34248,34731..35867)\n";
   print EMBL "FT                   /systematic_id=\"gene_id 0 ;\n";
   print EMBL "FT                   T1.PTRK.scaffold.00010.1389350.embl sequence\n";
   print EMBL "FT                   T1.PTRK.scaffold.00010.1389350 ; gene_orientation +\"\n";
   print EMBL "FT                   /colour=1\n";
   print EMBL "SQ   Sequence 1389350 BP; 500536 A; 170252 C; 170653 G; 506499 T; 41410 other;\n";
   print EMBL "aaaatatgga cttgaaagat tcatttatat ttggcagtgc aattaccgtt gctttgtttt        60\n";
   close(EMBL);
   ($errorcode,$errormsg)  = &convert_embl_to_gff($outputdir,$input_embl,$output_gff,'no');
   if ($errorcode != 0) { print STDERR "ERROR: test_convert_embl_to_gff: failed test4\n"; exit;}
   ($expected_output_gff,$errorcode,$errormsg) = &make_filename($outputdir);
   if ($errorcode != 0) { ($errorcode,$errormsg) = &print_error($errormsg,$errorcode,0); } 
   open(EXPECTED,">$expected_output_gff") || die "ERROR: test_convert_embl_to_gff: cannot open $expected_output_gff\n";
   print EXPECTED "PTRK.scaffold.00010.1389350\tsource\tCDS\t33177\t33365\t.\t+\t.\tT1.PTRK.scaffold.00010.1389350_1_0;Parent=T1.PTRK.scaffold.00010.1389350_1_mRNA\n";
   print EXPECTED "PTRK.scaffold.00010.1389350\tsource\tCDS\t33583\t34248\t.\t+\t.\tT1.PTRK.scaffold.00010.1389350_1_1;Parent=T1.PTRK.scaffold.00010.1389350_1_mRNA\n";
   print EXPECTED "PTRK.scaffold.00010.1389350\tsource\tCDS\t34731\t35867\t.\t+\t.\tT1.PTRK.scaffold.00010.1389350_1_2;Parent=T1.PTRK.scaffold.00010.1389350_1_mRNA\n";
   print EXPECTED "PTRK.scaffold.00010.1389350\tsource\tgene\t33177\t35867\t.\t+\t.\tT1.PTRK.scaffold.00010.1389350_1\n";
   print EXPECTED "PTRK.scaffold.00010.1389350\tsource\tmRNA\t33177\t35867\t.\t+\t.\tT1.PTRK.scaffold.00010.1389350_1_mRNA;Parent=T1.PTRK.scaffold.00010.1389350_1\n";
   print EXPECTED "PTRK.scaffold.00010.1389350\tsource\tCDS\t33177\t33365\t.\t+\t.\tT1.PTRK.scaffold.00010.1389350_2_0;Parent=T1.PTRK.scaffold.00010.1389350_2_mRNA\n";
   print EXPECTED "PTRK.scaffold.00010.1389350\tsource\tCDS\t33583\t34248\t.\t+\t.\tT1.PTRK.scaffold.00010.1389350_2_1;Parent=T1.PTRK.scaffold.00010.1389350_2_mRNA\n";
   print EXPECTED "PTRK.scaffold.00010.1389350\tsource\tCDS\t34731\t35867\t.\t+\t.\tT1.PTRK.scaffold.00010.1389350_2_2;Parent=T1.PTRK.scaffold.00010.1389350_2_mRNA\n";
   print EXPECTED "PTRK.scaffold.00010.1389350\tsource\tgene\t33177\t35867\t.\t+\t.\tT1.PTRK.scaffold.00010.1389350_2\n";  
   print EXPECTED "PTRK.scaffold.00010.1389350\tsource\tmRNA\t33177\t35867\t.\t+\t.\tT1.PTRK.scaffold.00010.1389350_2_mRNA;Parent=T1.PTRK.scaffold.00010.1389350_2\n";
   close(EXPECTED);
   $differences            = "";
   open(TEMP,"diff $output_gff $expected_output_gff |");
   while(<TEMP>)
   {
      $line                = $_;
      $differences         = $differences.$line;
   }
   close(TEMP);  
   $length_differences     = length($differences);
   if ($length_differences != 0) { print STDERR "ERROR: test_convert_embl_to_gff: failed test4 (output_gff $output_gff expected_output_gff $expected_output_gff)\n"; exit;}
   system "rm -f $input_embl";
   system "rm -f $output_gff";
   system "rm -f $expected_output_gff"; 
 
   # CHECK FOR ERRORCODE=1 (GENE INFORMATION OVER MORE THAN 3 LINES):
   ($output_gff,$errorcode,$errormsg) = &make_filename($outputdir);
   if ($errorcode != 0) { ($errorcode,$errormsg) = &print_error($errormsg,$errorcode,0); } 
   ($input_embl,$errorcode,$errormsg) = &make_filename($outputdir);
   if ($errorcode != 0) { ($errorcode,$errormsg) = &print_error($errormsg,$errorcode,0); } 
   open(EMBL,">$input_embl") || die "ERROR: test_convert_embl_to_gff: cannot open $input_embl\n";
   print EMBL "ID   PTRK.scaffold.00010.1389350; SV 1; linear; genomic DNA; HTG; INV; 1389350 BP.\n";
   print EMBL "FH   Key             Location/Qualifiers\n";
   print EMBL "FH\n";
   print EMBL "FT   CDS             join(33177..33365,33583..34248,34731..35867)\n";
   print EMBL "FT                   /systematic_id=\"gene_id 0 ;\n";
   print EMBL "FT                   T1.PTRK.scaffold.00010.1389350.embl\n";
   print EMBL "FT                   sequence\n";
   print EMBL "FT                   T1.PTRK.scaffold.00010.1389350 ;\n";
   print EMBL "FT                   gene_orientation +\"\n";
   print EMBL "FT                   /colour=1\n";
   print EMBL "FT   CDS             join(33177..33365,33583..34248,34731..35867)\n";
   print EMBL "FT                   /systematic_id=\"gene_id 0 ;\n";
   print EMBL "FT                   T1.PTRK.scaffold.00010.1389350.embl sequence\n";
   print EMBL "FT                   T1.PTRK.scaffold.00010.1389350 ; gene_orientation +\"\n";
   print EMBL "FT                   /colour=1\n";
   print EMBL "SQ   Sequence 1389350 BP; 500536 A; 170252 C; 170653 G; 506499 T; 41410 other;\n";
   print EMBL "aaaatatgga cttgaaagat tcatttatat ttggcagtgc aattaccgtt gctttgtttt        60\n";
   close(EMBL);
   ($errorcode,$errormsg)  = &convert_embl_to_gff($outputdir,$input_embl,$output_gff,'no');
   if ($errorcode != 1) { print STDERR "ERROR: test_convert_embl_to_gff: failed test4 (errorcode $errorcode errormsg $errormsg)\n"; exit;}
   system "rm -f $input_embl";
   system "rm -f $output_gff";
 
   # ANOTHER EXAMPLE LIKE THE EMBL FILES MADE BY THE SGAs:
   ($output_gff,$errorcode,$errormsg) = &make_filename($outputdir);
   if ($errorcode != 0) { ($errorcode,$errormsg) = &print_error($errormsg,$errorcode,0); } 
   ($input_embl,$errorcode,$errormsg) = &make_filename($outputdir);
   if ($errorcode != 0) { ($errorcode,$errormsg) = &print_error($errormsg,$errorcode,0); } 
   open(EMBL,">$input_embl") || die "ERROR: test_convert_embl_to_gff: cannot open $input_embl\n";
   print EMBL "ID   PTRK.scaffold.00010.1389350; SV 1; linear; genomic DNA; HTG; INV; 1389350 BP.\n";
   print EMBL "FH   Key             Location/Qualifiers\n";
   print EMBL "FH\n";
   print EMBL "FT   CDS             join(1186266..1186427,1186482..1186763)\n";
   print EMBL "FT                   /systematic_id=\"gene_id 0 ; sequence\n";
   print EMBL "FT                   T30.PTRK.scaffold.00010.1389350 ; gene_orientation +\"\n";
   print EMBL "FT                   /colour=1\n";
   print EMBL "FT   CDS             join(1193863..1195911,1195975..1196085)\n";
   print EMBL "FT                   /systematic_id=\"gene_id 0 ; sequence\n";
   print EMBL "FT                   T31.PTRK.scaffold.00010.1389350  ; gene_orientation +\"\n";
   print EMBL "FT                   /colour=1\n";
   print EMBL "SQ   Sequence 1389350 BP; 500536 A; 170252 C; 170653 G; 506499 T; 41410 other;\n";
   print EMBL "aaaatatgga cttgaaagat tcatttatat ttggcagtgc aattaccgtt gctttgtttt        60\n";
   close(EMBL);
   ($errorcode,$errormsg)  = &convert_embl_to_gff($outputdir,$input_embl,$output_gff,'no');
   if ($errorcode != 0) { print STDERR "ERROR: test_convert_embl_to_gff: failed test5\n"; exit;}
   ($expected_output_gff,$errorcode,$errormsg) = &make_filename($outputdir);
   if ($errorcode != 0) { ($errorcode,$errormsg) = &print_error($errormsg,$errorcode,0); } 
   open(EXPECTED,">$expected_output_gff") || die "ERROR: test_convert_embl_to_gff: cannot open $expected_output_gff\n";
   print EXPECTED "PTRK.scaffold.00010.1389350\tsource\tCDS\t1186266\t1186427\t.\t+\t.\tT30.PTRK.scaffold.00010.1389350_1_0;Parent=T30.PTRK.scaffold.00010.1389350_1_mRNA\n";
   print EXPECTED "PTRK.scaffold.00010.1389350\tsource\tCDS\t1186482\t1186763\t.\t+\t.\tT30.PTRK.scaffold.00010.1389350_1_1;Parent=T30.PTRK.scaffold.00010.1389350_1_mRNA\n";
   print EXPECTED "PTRK.scaffold.00010.1389350\tsource\tgene\t1186266\t1186763\t.\t+\t.\tT30.PTRK.scaffold.00010.1389350_1\n";
   print EXPECTED "PTRK.scaffold.00010.1389350\tsource\tmRNA\t1186266\t1186763\t.\t+\t.\tT30.PTRK.scaffold.00010.1389350_1_mRNA;Parent=T30.PTRK.scaffold.00010.1389350_1\n";
   print EXPECTED "PTRK.scaffold.00010.1389350\tsource\tCDS\t1193863\t1195911\t.\t+\t.\tT31.PTRK.scaffold.00010.1389350_1_0;Parent=T31.PTRK.scaffold.00010.1389350_1_mRNA\n";
   print EXPECTED "PTRK.scaffold.00010.1389350\tsource\tCDS\t1195975\t1196085\t.\t+\t.\tT31.PTRK.scaffold.00010.1389350_1_1;Parent=T31.PTRK.scaffold.00010.1389350_1_mRNA\n";
   print EXPECTED "PTRK.scaffold.00010.1389350\tsource\tgene\t1193863\t1196085\t.\t+\t.\tT31.PTRK.scaffold.00010.1389350_1\n";
   print EXPECTED "PTRK.scaffold.00010.1389350\tsource\tmRNA\t1193863\t1196085\t.\t+\t.\tT31.PTRK.scaffold.00010.1389350_1_mRNA;Parent=T31.PTRK.scaffold.00010.1389350_1\n";
   close(EXPECTED);
   $differences            = "";
   open(TEMP,"diff $output_gff $expected_output_gff |");
   while(<TEMP>)
   {
      $line                = $_;
      $differences         = $differences.$line;
   }
   close(TEMP);  
   $length_differences     = length($differences);
   if ($length_differences != 0) { print STDERR "ERROR: test_convert_embl_to_gff: failed test5 (output_gff $output_gff expected_output_gff $expected_output_gff)\n"; exit;}
   system "rm -f $input_embl";
   system "rm -f $output_gff";
   system "rm -f $expected_output_gff"; 
 
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
