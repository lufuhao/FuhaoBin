#!/usr/bin/env perl
#https://gist.github.com/avrilcoghlan/5310980

=head1 NAME

    translate_spliced_dna.pl

=head1 SYNOPSIS
 
    translate_spliced_dna.pl spliced_dna translated_dna outputdir 
        where spliced_dna is the input fasta file of spliced DNA sequences for transcripts,
              translated_dna is the output fasta file of amino acid translations of transcripts,
              outputdir is the output directory for writing output files.

=head1 DESCRIPTION

    This script reads in a fasta file of transcript DNA sequences (<spliced_dna>) and translates
    each transcript. It prints out a fasta file of amino acid sequences for transcripts (<translated_dna>).
    Note: this expects different genes to have unique names: two genes on different scaffolds cannot have the
    same name.

=head1 VERSION
  
    Perl script last edited 2-Apr-2013.

=head1 CONTACT

    alc@sanger.ac.uk (Avril Coghlan)

=cut
 
# 
# Perl script translate_spliced_dna.pl
# Written by Avril Coghlan (alc@sanger.ac.uk)
# 2-Apr-13.
# Last edited 2-Apr-2013.
# SCRIPT SYNOPSIS: translate_spliced_dna.pl: given an input file of DNA sequences of transcripts, infers their translations. 
# 
#------------------------------------------------------------------#
 
# CHECK IF THERE ARE THE CORRECT NUMBER OF COMMAND-LINE ARGUMENTS:
 
$num_args                  = $#ARGV + 1;
if ($num_args != 3)
{
   print "Usage of translate_spliced_dna.pl\n\n";
   print "perl translate_spliced_dna_.pl <spliced_dna> <translated_dna> <outputdir>\n";
   print "where <spliced_dna> is the input fasta file of spliced DNA sequences for transcripts,\n";
   print "      <translated_dna> is the output fasta file of amino acid translations of transcripts,\n";
   print "      <outputdir> is the output directory for writing output files.\n";
   print "For example, >perl -w translate_spliced_dna.pl spliced_dna translated_dna /lustre/scratch108/parasites/alc/50HGI_maker_translations\n";
   exit;
}
 
# FIND THE NAME OF THE INPUT FASTA FILE:
 
$spliced_dna               = $ARGV[0];
 
# FIND THE NAME OF THE OUTPUT FASTA FILE:
 
$translated_dna            = $ARGV[1];
 
# FIND THE NAME OF THE OUTPUT DIRECTORY:
 
$outputdir                 = $ARGV[2];
 
# TEST SUBROUTINES: 
 
my $PRINT_TEST_DATA        = 0;   # SAYS WHETHER TO PRINT DATA USED DURING TESTING.
&test_print_to_output($outputdir);
&test_print_error;
&test_read_assembly($outputdir);
&test_write_out_translations($outputdir);
&test_translate;
&test_run_main_program($outputdir);
print STDERR "Finished tests: now running main code...\n";
 
#------------------------------------------------------------------#
 
# RUN THE MAIN PART OF THE CODE:
 
($errorcode,$errormsg)     = &run_main_program($outputdir,$spliced_dna,$translated_dna);
if ($errorcode != 0) { ($errorcode,$errormsg) = &print_error($errormsg,$errorcode,0); } 
 
print STDERR "FINISHED.\n";
 
#------------------------------------------------------------------#
 
# TEST &run_main_program
 
sub test_run_main_program
{
   my $outputdir           = $_[0]; # DIRECTORY TO PUT OUTPUT FILES INTO
   my $input;                       # INPUT FILE OF SPLICED DNA SEQUENCES
   my $output;                      # OUTPUT FILE OF TRANSLATIONS
   my $expected_output;             # FILE CONTAINING EXPECTED CONTENT OF $output
   my $differences;                 # DIFFERENCES BETWEEN $output AND $expected_output
   my $length_differences;          # LENGTH OF $differences
   my $line;                        #  
   my @temp;                        #  
 
   ($input,$errorcode,$errormsg) = &make_filename($outputdir);
   if ($errorcode != 0) { ($errorcode,$errormsg) = &print_error($errormsg,$errorcode,0); } 
   open(INPUT,">$input") || die "ERROR: test_run_main_program: cannot open $input\n";
   print INPUT ">seq1\n";
   print INPUT "ATGCAGTAA\n";
   print INPUT ">seq2\n";
   print INPUT "ATGGGTGCAACATTGTCTTCAACTTCTTCTTCACCATCTTCTATTGATGGAAGAGAAAATTTGATAAGAAGAAATAACACATCAACATACAATCATTTATTTGTCACACCTGGTAATACAAAAAATAAGAATAATAATTTTCTATTTGATTCAAAAAATCGGCACTCATCTCTTTCTTCATCTTCTTCTATGCTTGGTATTGCATGGAACTTTGCCAAAAAAACAACTACGGGAATTCCTAATAAAGAAAAAAGAAATATAAATAGATTATCATCATCTAGTTCAGTATCATCTGCTTTATCAAAATCTCATGATTCTGCAACTTCATTAACTTTTAATTCAAGTTTTTCTATTAATTCATCAAGTGACCAATTATCAGCATACTCTTCTAGTGCATCAAATTCTTCATCAGGTGATAGTGGTATTGGAATATCTAGAAAAATAACAATATCTAATCAATGTTATCATGAAATTGATGGCAAAAATAACAATAATAATATTGTAAATAATGTACCTACTTTAAATACTAATGAAAAATTTATTTGTATTAGAAAAAAAGAATCAATTTTATCAAGTAAAATAATGTCCGTACCATTAAATAATATAAAAAAAGGAAGTTCTAAAGTTCGTGATCATATATTAGCAATGACAAGATCAAAATCAACAAATATGCGTACAGAAAAAAATGGAAACATAATTATAAATGGTAATAATTATAATATTGAAAATAAATTTAATGATAATTTAAAATTATCATCAACAACTTATGATTATGAAAAAAAATTAAATCAAAATATAGAAAATAATATATCTAAATCATTGAAAAATGTTACTTTAACAAATCGAAATAATAATTCAAATTTATATAATACTGCAATATTTAATTTCGGCATTAAAAATGATAGTAATGTTATAGATAAAACAAGAAAAAAAACAATTATACAAGCATCAACAACTGAATTATTAAAAGGTGTTTCATTATTAATTACAACTAAATGTTCCCATAAAGTACCAGATTTTTTTGCTAATCAATTAACAATGTGGTTAAGAAGTGTTGATAGATCTTTAATTGTTCAGGGTTGGCAAGATATCGCATTTTTAAATCCAGCTAACATGGTATTTTTTTTTATGCTTTTAAGATCAATGCTTAATGATGAAGAAACATTTCCTGTTAATAATTTGGAAGACCTTCAAATGATTGTTTTTACTTGTTTATTTATTAGTTACAGTTACATGGGTAATGAAATTAGTTATCCGTTAAAACCATTTATTTCACAACAGGATCGCTCAAAATTTTGGGATATGTGTGTTCAAATAGTTAATAAATATTCTGGGGATATGTTACAATTAAATACATCAGCTTCATTTTTTACACAAGTTTTTAGTGATCTCAAAAATTTTACAACTATTGCTAACAATTAA\n";
   close(INPUT);
   # WRITE OUT THE OUTPUT FILE OF TRANSLATED SEQUENCES:
   ($output,$errorcode,$errormsg) = &make_filename($outputdir);
   if ($errorcode != 0) { ($errorcode,$errormsg) = &print_error($errormsg,$errorcode,0); } 
   @temp                   = split(/\//,$output);
   $output                 = $temp[$#temp];
   ($errorcode,$errormsg)  = &run_main_program($outputdir,$input,$output);
   if ($errorcode != 0) { print STDERR "ERROR: test_run_main_program: failed test1\n"; exit;}
   ($expected_output,$errorcode,$errormsg) = &make_filename($outputdir);
   if ($errorcode != 0) { ($errorcode,$errormsg) = &print_error($errormsg,$errorcode,0); } 
   open(EXPECTED,">$expected_output") || die "ERROR: test_run_main_program: cannot open $expected_output\n"; 
   print EXPECTED ">seq1\n";
   print EXPECTED "MQ*\n";
   print EXPECTED ">seq2\n";
   print EXPECTED "MGATLSSTSSSPSSIDGRENLIRRNNTSTYNHLFVTPGNTKNKNNNFLFDSKNRHSSLSS\n";
   print EXPECTED "SSSMLGIAWNFAKKTTTGIPNKEKRNINRLSSSSSVSSALSKSHDSATSLTFNSSFSINS\n";
   print EXPECTED "SSDQLSAYSSSASNSSSGDSGIGISRKITISNQCYHEIDGKNNNNNIVNNVPTLNTNEKF\n";
   print EXPECTED "ICIRKKESILSSKIMSVPLNNIKKGSSKVRDHILAMTRSKSTNMRTEKNGNIIINGNNYN\n";
   print EXPECTED "IENKFNDNLKLSSTTYDYEKKLNQNIENNISKSLKNVTLTNRNNNSNLYNTAIFNFGIKN\n";
   print EXPECTED "DSNVIDKTRKKTIIQASTTELLKGVSLLITTKCSHKVPDFFANQLTMWLRSVDRSLIVQG\n";
   print EXPECTED "WQDIAFLNPANMVFFFMLLRSMLNDEETFPVNNLEDLQMIVFTCLFISYSYMGNEISYPL\n";
   print EXPECTED "KPFISQQDRSKFWDMCVQIVNKYSGDMLQLNTSASFFTQVFSDLKNFTTIANN*\n";
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
   if ($length_differences != 0) { print STDERR "ERROR: test_run_main_program: failed test1 (output $output expected_output $expected_output)\n"; exit;}
   system "rm -f $output";
   system "rm -f $expected_output";
   system "rm -f $input";
 
}
 
#------------------------------------------------------------------#
 
# RUN THE MAIN PART OF THE CODE:
 
sub run_main_program
{
   my $outputdir           = $_[0]; # DIRECTORY TO PUT OUTPUT FILES IN.
   my $spliced_dna         = $_[1]; # FASTA FILE OF SPLICED DNA SEQUENCES
   my $translated_dna      = $_[2]; # OUTPUT FASTA FILE OF TRANSLATIONS
   my $errorcode           = 0;     # RETURNED AS 0 IF THERE IS NO ERROR.
   my $errormsg            = 'none';# RETURNED AS 'none' IF THERE IS NO ERROR. 
   my $SEQ;                         # HASH TABLE OF SPLICED DNA SEQUENCES 
 
   # READ IN THE FASTA FILE OF SPLICED DNA SEQUENCES:
   ($SEQ,$errorcode,$errormsg) = &read_assembly($spliced_dna);
   if ($errorcode != 0) { ($errorcode,$errormsg) = &print_error($errormsg,$errorcode,0); }
 
   # WRITE OUT THE OUTPUT FILE OF TRANSLATED SEQUENCES:
   ($errorcode,$errormsg)  = &write_out_translations($SEQ,$translated_dna,$outputdir);
   if ($errorcode != 0) { ($errorcode,$errormsg) = &print_error($errormsg,$errorcode,0); }
 
   return($errorcode,$errormsg);
} 
 
#------------------------------------------------------------------#
 
# TEST &write_out_translations
 
sub test_write_out_translations
{
   my $outputdir           = $_[0]; # DIRECTORY TO PUT OUTPUT FILES INTO
   my %SEQ                 = ();    # HASH TABLE OF SPLICED DNA SEQUENCES
   my $output;                      # OUTPUT FILE OF TRANSLATIONS
   my $expected_output;             # FILE CONTAINING EXPECTED CONTENT OF $output
   my $differences;                 # DIFFERENCES BETWEEN $output AND $expected_output
   my $length_differences;          # LENGTH OF $differences
   my $line;                        #  
   my @temp;                        #  
 
   $SEQ{"seq1"}            = "ATGCAGTAA";
   $SEQ{"seq2"}            = "ATGGGTGCAACATTGTCTTCAACTTCTTCTTCACCATCTTCTATTGATGGAAGAGAAAATTTGATAAGAAGAAATAACACATCAACATACAATCATTTATTTGTCACACCTGGTAATACAAAAAATAAGAATAATAATTTTCTATTTGATTCAAAAAATCGGCACTCATCTCTTTCTTCATCTTCTTCTATGCTTGGTATTGCATGGAACTTTGCCAAAAAAACAACTACGGGAATTCCTAATAAAGAAAAAAGAAATATAAATAGATTATCATCATCTAGTTCAGTATCATCTGCTTTATCAAAATCTCATGATTCTGCAACTTCATTAACTTTTAATTCAAGTTTTTCTATTAATTCATCAAGTGACCAATTATCAGCATACTCTTCTAGTGCATCAAATTCTTCATCAGGTGATAGTGGTATTGGAATATCTAGAAAAATAACAATATCTAATCAATGTTATCATGAAATTGATGGCAAAAATAACAATAATAATATTGTAAATAATGTACCTACTTTAAATACTAATGAAAAATTTATTTGTATTAGAAAAAAAGAATCAATTTTATCAAGTAAAATAATGTCCGTACCATTAAATAATATAAAAAAAGGAAGTTCTAAAGTTCGTGATCATATATTAGCAATGACAAGATCAAAATCAACAAATATGCGTACAGAAAAAAATGGAAACATAATTATAAATGGTAATAATTATAATATTGAAAATAAATTTAATGATAATTTAAAATTATCATCAACAACTTATGATTATGAAAAAAAATTAAATCAAAATATAGAAAATAATATATCTAAATCATTGAAAAATGTTACTTTAACAAATCGAAATAATAATTCAAATTTATATAATACTGCAATATTTAATTTCGGCATTAAAAATGATAGTAATGTTATAGATAAAACAAGAAAAAAAACAATTATACAAGCATCAACAACTGAATTATTAAAAGGTGTTTCATTATTAATTACAACTAAATGTTCCCATAAAGTACCAGATTTTTTTGCTAATCAATTAACAATGTGGTTAAGAAGTGTTGATAGATCTTTAATTGTTCAGGGTTGGCAAGATATCGCATTTTTAAATCCAGCTAACATGGTATTTTTTTTTATGCTTTTAAGATCAATGCTTAATGATGAAGAAACATTTCCTGTTAATAATTTGGAAGACCTTCAAATGATTGTTTTTACTTGTTTATTTATTAGTTACAGTTACATGGGTAATGAAATTAGTTATCCGTTAAAACCATTTATTTCACAACAGGATCGCTCAAAATTTTGGGATATGTGTGTTCAAATAGTTAATAAATATTCTGGGGATATGTTACAATTAAATACATCAGCTTCATTTTTTACACAAGTTTTTAGTGATCTCAAAAATTTTACAACTATTGCTAACAATTAA";
   # WRITE OUT THE OUTPUT FILE OF TRANSLATED SEQUENCES:
   ($output,$errorcode,$errormsg) = &make_filename($outputdir);
   if ($errorcode != 0) { ($errorcode,$errormsg) = &print_error($errormsg,$errorcode,0); } 
   @temp                   = split(/\//,$output);
   $output                 = $temp[$#temp];
   ($errorcode,$errormsg)  = &write_out_translations(\%SEQ,$output,$outputdir);
   if ($errorcode != 0) { print STDERR "ERROR: test_write_out_translation: failed test1\n"; exit;}
   ($expected_output,$errorcode,$errormsg) = &make_filename($outputdir);
   if ($errorcode != 0) { ($errorcode,$errormsg) = &print_error($errormsg,$errorcode,0); } 
   open(EXPECTED,">$expected_output") || die "ERROR: test_write_out_translations: cannot open $expected_output\n"; 
   print EXPECTED ">seq1\n";
   print EXPECTED "MQ*\n";
   print EXPECTED ">seq2\n";
   print EXPECTED "MGATLSSTSSSPSSIDGRENLIRRNNTSTYNHLFVTPGNTKNKNNNFLFDSKNRHSSLSS\n";
   print EXPECTED "SSSMLGIAWNFAKKTTTGIPNKEKRNINRLSSSSSVSSALSKSHDSATSLTFNSSFSINS\n";
   print EXPECTED "SSDQLSAYSSSASNSSSGDSGIGISRKITISNQCYHEIDGKNNNNNIVNNVPTLNTNEKF\n";
   print EXPECTED "ICIRKKESILSSKIMSVPLNNIKKGSSKVRDHILAMTRSKSTNMRTEKNGNIIINGNNYN\n";
   print EXPECTED "IENKFNDNLKLSSTTYDYEKKLNQNIENNISKSLKNVTLTNRNNNSNLYNTAIFNFGIKN\n";
   print EXPECTED "DSNVIDKTRKKTIIQASTTELLKGVSLLITTKCSHKVPDFFANQLTMWLRSVDRSLIVQG\n";
   print EXPECTED "WQDIAFLNPANMVFFFMLLRSMLNDEETFPVNNLEDLQMIVFTCLFISYSYMGNEISYPL\n";
   print EXPECTED "KPFISQQDRSKFWDMCVQIVNKYSGDMLQLNTSASFFTQVFSDLKNFTTIANN*\n";
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
   if ($length_differences != 0) { print STDERR "ERROR: test_write_out_translation: failed test1 (output $output expected_output $expected_output)\n"; exit;}
   system "rm -f $output";
   system "rm -f $expected_output";
 
}
 
#------------------------------------------------------------------#
 
# WRITE OUT THE OUTPUT FILE OF TRANSLATED SEQUENCES:
 
sub write_out_translations
{
   my $SEQ                 = $_[0]; # HASH TABLE WITH SPLICED DNA SEQUENCES FOR TRANSCRIPTS
   my $translated_dna      = $_[1]; # OUTPUT FILE
   my $outputdir           = $_[2]; # DIRECTORY TO PUT OUTPUT FILES INTO.
   my $errorcode           = 0;     # RETURNED AS 0 IF THERE IS NO ERROR.
   my $errormsg            = "none";# RETURNED AS 'none' IF THERE IS NO ERROR.
   my $seqname;                     # NAME OF A TRANSCRIPT 
   my $seq;                         # DNA SEQUENCE OF A TRANSCRIPT
 
   # OPEN THE OUTPUT FILE:
   $translated_dna         = $outputdir."/".$translated_dna;
   open(OUTPUT,">$translated_dna") || die "ERROR: write_out_translations: cannot open $translated_dna\n";
   close(OUTPUT);
 
   # LOOP THROUGH ALL THE DNA SEQUENCES:
   foreach $seqname (keys %{$SEQ})
   {
      $seq                 = $SEQ->{$seqname};
      # TRANSLATE THIS SEQUENCE:
      ($seq,$errorcode,$errormsg) = &translate($seq);
      if ($errorcode != 0) { ($errorcode,$errormsg) = &print_error($errormsg,$errorcode,0); }
      # WRITE OUT THE SEQUENCE IN THE OUTPUT FASTA FILE:
      ($errorcode,$errormsg)  = &print_to_output($translated_dna,$seq,$seqname);
      if ($errorcode != 0) { ($errorcode,$errormsg) = &print_error($errormsg,$errorcode,0); }     
   }
 
   return($errorcode,$errormsg);
}
 
#------------------------------------------------------------------#
 
# TEST &translate
 
sub test_translate
{
   my $errorcode;                   # RETURNED AS 0 FROM A FUNCTION IF THERE IS NO ERROR
   my $errormsg;                    # RETURNED AS 'none' FROM A FUNCTION IF THERE IS NO ERROR
   my $seq;                         # TRANSLATED SEQUENCE        
 
   ($seq,$errorcode,$errormsg) = &translate("ATGCAGTAA");
   if ($errorcode != 0 || $seq ne "MQ*") { print STDERR "ERROR: test_translate: failed test1 (seq $seq errorcode $errorcode errormsg $errormsg)\n"; exit;}
 
   ($seq,$errorcode,$errormsg) = &translate("ATGGGTGCAACATTGTCTTCAACTTCTTCTTCACCATCTTCTATTGATGGAAGAGAAAATTTGATAAGAAGAAATAACACATCAACATACAATCATTTATTTGTCACACCTGGTAATACAAAAAATAAGAATAATAATTTTCTATTTGATTCAAAAAATCGGCACTCATCTCTTTCTTCATCTTCTTCTATGCTTGGTATTGCATGGAACTTTGCCAAAAAAACAACTACGGGAATTCCTAATAAAGAAAAAAGAAATATAAATAGATTATCATCATCTAGTTCAGTATCATCTGCTTTATCAAAATCTCATGATTCTGCAACTTCATTAACTTTTAATTCAAGTTTTTCTATTAATTCATCAAGTGACCAATTATCAGCATACTCTTCTAGTGCATCAAATTCTTCATCAGGTGATAGTGGTATTGGAATATCTAGAAAAATAACAATATCTAATCAATGTTATCATGAAATTGATGGCAAAAATAACAATAATAATATTGTAAATAATGTACCTACTTTAAATACTAATGAAAAATTTATTTGTATTAGAAAAAAAGAATCAATTTTATCAAGTAAAATAATGTCCGTACCATTAAATAATATAAAAAAAGGAAGTTCTAAAGTTCGTGATCATATATTAGCAATGACAAGATCAAAATCAACAAATATGCGTACAGAAAAAAATGGAAACATAATTATAAATGGTAATAATTATAATATTGAAAATAAATTTAATGATAATTTAAAATTATCATCAACAACTTATGATTATGAAAAAAAATTAAATCAAAATATAGAAAATAATATATCTAAATCATTGAAAAATGTTACTTTAACAAATCGAAATAATAATTCAAATTTATATAATACTGCAATATTTAATTTCGGCATTAAAAATGATAGTAATGTTATAGATAAAACAAGAAAAAAAACAATTATACAAGCATCAACAACTGAATTATTAAAAGGTGTTTCATTATTAATTACAACTAAATGTTCCCATAAAGTACCAGATTTTTTTGCTAATCAATTAACAATGTGGTTAAGAAGTGTTGATAGATCTTTAATTGTTCAGGGTTGGCAAGATATCGCATTTTTAAATCCAGCTAACATGGTATTTTTTTTTATGCTTTTAAGATCAATGCTTAATGATGAAGAAACATTTCCTGTTAATAATTTGGAAGACCTTCAAATGATTGTTTTTACTTGTTTATTTATTAGTTACAGTTACATGGGTAATGAAATTAGTTATCCGTTAAAACCATTTATTTCACAACAGGATCGCTCAAAATTTTGGGATATGTGTGTTCAAATAGTTAATAAATATTCTGGGGATATGTTACAATTAAATACATCAGCTTCATTTTTTACACAAGTTTTTAGTGATCTCAAAAATTTTACAACTATTGCTAACAATTAA");
   if ($errorcode != 0 || $seq ne "MGATLSSTSSSPSSIDGRENLIRRNNTSTYNHLFVTPGNTKNKNNNFLFDSKNRHSSLSSSSSMLGIAWNFAKKTTTGIPNKEKRNINRLSSSSSVSSALSKSHDSATSLTFNSSFSINSSSDQLSAYSSSASNSSSGDSGIGISRKITISNQCYHEIDGKNNNNNIVNNVPTLNTNEKFICIRKKESILSSKIMSVPLNNIKKGSSKVRDHILAMTRSKSTNMRTEKNGNIIINGNNYNIENKFNDNLKLSSTTYDYEKKLNQNIENNISKSLKNVTLTNRNNNSNLYNTAIFNFGIKNDSNVIDKTRKKTIIQASTTELLKGVSLLITTKCSHKVPDFFANQLTMWLRSVDRSLIVQGWQDIAFLNPANMVFFFMLLRSMLNDEETFPVNNLEDLQMIVFTCLFISYSYMGNEISYPLKPFISQQDRSKFWDMCVQIVNKYSGDMLQLNTSASFFTQVFSDLKNFTTIANN*") { print STDERR "ERROR: test_translate: failed test2\n"; exit;}
 
   # TEST ERRORCODE=1 (NON-ACGTN BASES):
   ($seq,$errorcode,$errormsg) = &translate("ATGCAG-TAA");
   if ($errorcode !=1) { print STDERR "ERROR: test_translate: failed test3 (seq $seq errorcode $errorcode errormsg $errormsg)\n"; exit;}
 
}
  
 
#------------------------------------------------------------------#
 
# SUBROUTINE TO TRANSLATE A DNA SEQUENCE TO A PROTEIN SEQUENCE:
 
sub translate
{
   my $dna                 = $_[0]; # INPUT DNA SEQUENCE
   my $protein             = "";    # PROTEIN SEQUENCE
   my $triplet;                     # A TRIPLET OF THE SEQUENCE
   my $aa;                          # AN AMINO ACID
   my $length;                      # LENGTH OF THE SEQUENCE
   my $i;                           #
   my $errorcode           = 0;     # RETURNED AS 0 IF THERE IS NO ERROR
   my $errormsg            = "none";# RETURNED AS 'none' IF THERE IS NO ERROR
   my $first_base;                  # FIRST BASE OF TRIPLET
   my $second_base;                 # SECOND BASE OF TRIPLET
   my $third_base;                  # THIRD BASE OF TRIPLET 
   my $protein_length_bp;           # LENGTH OF THE PROTEIN SEQUENCE IN BASE-PAIRS
 
   $dna                    =~ tr/a-z/A-Z/;
   $length                 = length($dna);
   $protein                = "";
   for ($i = 0; $i <= $length-3; $i = $i+3)
   {
      $triplet             = substr($dna,$i,3);
      $first_base          = substr($triplet,0,1);
      $second_base         = substr($triplet,1,1);
      $third_base          = substr($triplet,2,1);
      if (($first_base ne 'A'  && $first_base ne 'C'  && $first_base ne 'G'  && $first_base ne 'T'  && $first_base ne 'N')  ||
          ($second_base ne 'A' && $second_base ne 'C' && $second_base ne 'G' && $second_base ne 'T' && $second_base ne 'N') || 
          ($third_base ne 'A'  && $third_base ne 'C'  && $third_base ne 'G'  && $third_base ne 'T' && $third_base ne 'N'))     
      {
         $errormsg         = "ERROR: translate: triplet $triplet contains non-ACGTN bases (first_base $first_base second_base $second_base third_base $third_base)\n";
         $errorcode        = 1; # ERRORCODE=1 (TESTED FOR)
         return($protein,$errorcode,$errormsg);   
      }
 
      $aa                  = "none";
      if ($triplet eq 'TTT') { $aa = 'F';}
      if ($triplet eq 'TTC') { $aa = 'F';}
      if ($triplet eq 'TTA') { $aa = 'L';}
      if ($triplet eq 'TTG') { $aa = 'L';} # 4
 
      if ($triplet eq 'TCT') { $aa = 'S';}
      if ($triplet eq 'TCC') { $aa = 'S';}
      if ($triplet eq 'TCA') { $aa = 'S';}
      if ($triplet eq 'TCG') { $aa = 'S';} # 8
 
      if ($triplet eq 'TAT') { $aa = 'Y';}
      if ($triplet eq 'TAC') { $aa = 'Y';}
      if ($triplet eq 'TAA') { $aa = '*';}
      if ($triplet eq 'TAG') { $aa = '*';} # 12
 
      if ($triplet eq 'TGT') { $aa = 'C';}
      if ($triplet eq 'TGC') { $aa = 'C';}
      if ($triplet eq 'TGA') { $aa = '*';}
      if ($triplet eq 'TGG') { $aa = 'W';} # 16
 
      if ($triplet eq 'CTT') { $aa = 'L';}
      if ($triplet eq 'CTC') { $aa = 'L';}
      if ($triplet eq 'CTA') { $aa = 'L';}
      if ($triplet eq 'CTG') { $aa = 'L';} # 20
 
      if ($triplet eq 'CCT') { $aa = 'P';}
      if ($triplet eq 'CCC') { $aa = 'P';}
      if ($triplet eq 'CCA') { $aa = 'P';}
      if ($triplet eq 'CCG') { $aa = 'P';} # 24
 
      if ($triplet eq 'CAT') { $aa = 'H';}
      if ($triplet eq 'CAC') { $aa = 'H';}
      if ($triplet eq 'CAA') { $aa = 'Q';}
      if ($triplet eq 'CAG') { $aa = 'Q';} # 28
 
      if ($triplet eq 'CGT') { $aa = 'R';}
      if ($triplet eq 'CGC') { $aa = 'R';}
      if ($triplet eq 'CGA') { $aa = 'R';}
      if ($triplet eq 'CGG') { $aa = 'R';} # 32
 
      if ($triplet eq 'ATT') { $aa = 'I';}
      if ($triplet eq 'ATC') { $aa = 'I';}
      if ($triplet eq 'ATA') { $aa = 'I';}
      if ($triplet eq 'ATG') { $aa = 'M';} # 36
 
      if ($triplet eq 'ACT') { $aa = 'T';}
      if ($triplet eq 'ACC') { $aa = 'T';}
      if ($triplet eq 'ACA') { $aa = 'T';}
      if ($triplet eq 'ACG') { $aa = 'T';} # 40
 
      if ($triplet eq 'AAT') { $aa = 'N';}
      if ($triplet eq 'AAC') { $aa = 'N';}
      if ($triplet eq 'AAA') { $aa = 'K';}
      if ($triplet eq 'AAG') { $aa = 'K';} # 44
 
      if ($triplet eq 'AGT') { $aa = 'S';}
      if ($triplet eq 'AGC') { $aa = 'S';}
      if ($triplet eq 'AGA') { $aa = 'R';}
      if ($triplet eq 'AGG') { $aa = 'R';} # 48
 
      if ($triplet eq 'GTT') { $aa = 'V';}
      if ($triplet eq 'GTC') { $aa = 'V';}
      if ($triplet eq 'GTA') { $aa = 'V';}
      if ($triplet eq 'GTG') { $aa = 'V';} # 52
 
      if ($triplet eq 'GCT') { $aa = 'A';}
      if ($triplet eq 'GCC') { $aa = 'A';}
      if ($triplet eq 'GCA') { $aa = 'A';}
      if ($triplet eq 'GCG') { $aa = 'A';} # 56
 
      if ($triplet eq 'GAT') { $aa = 'D';}
      if ($triplet eq 'GAC') { $aa = 'D';}
      if ($triplet eq 'GAA') { $aa = 'E';}
      if ($triplet eq 'GAG') { $aa = 'E';} # 60
 
      if ($triplet eq 'GGT') { $aa = 'G';}
      if ($triplet eq 'GGC') { $aa = 'G';}
      if ($triplet eq 'GGA') { $aa = 'G';}
      if ($triplet eq 'GGG') { $aa = 'G';} # 64
 
      if ($aa eq "none") { $aa = 'X';} # FOR AN UNKNOWN AMINO ACIDS JUST PUT X.
      $protein             = $protein.$aa;
   }
   $protein_length_bp      = length($protein) * 3;
 
   # IF THERE WERE 1 OR 2 OVERLHANGING BASES AT THE END OF THE SPLICED GENE DNA:
   if ($protein_length_bp < $length)
   {
      $protein             = $protein."X";
   }
 
   return($protein,$errorcode,$errormsg);
}
 
#------------------------------------------------------------------#
 
# WRITE OUT THE SEQUENCE IN THE OUTPUT FASTA FILE (THIS APPENDS TO AN EXISTING OUTPUT FILE):
 
sub print_to_output
{
   my $output_fasta        = $_[0]; # OUTPUT FASTA FILE
   my $seq                 = $_[1]; # SEQUENCE
   my $name                = $_[2]; # NAME OF THE SEQUENCE
   my $length;                      # LENGTH OF THE SEQUENCE
   my $offset;                      # 
   my $a_line;                      # A LINE OF THE SEQUENCE
   my $errorcode           = 0;     # RETURNED AS 0 IF THERE IS NO ERROR
   my $errormsg            = "none";# RETURNED AS 'none' IF THERE IS NO ERROR
 
   # MAKE AN OUTPUT FASTA FILE:
   open(OUTPUT,">>$output_fasta") || die "ERROR: print_to_output: cannot open output_fasta $output_fasta\n";
 
   print OUTPUT ">$name\n";
   $length                 = length($seq);
   $offset                 = 0;
   while ($offset < $length)
   {
      $a_line              = substr($seq,$offset,60);
      print OUTPUT "$a_line\n";
      $offset              = $offset + 60;
   }
 
   # CLOSE THE OUTPUT FILE:
   close(OUTPUT);  
 
   return($errorcode,$errormsg);
}
 
#------------------------------------------------------------------#
 
# TEST &print_to_output
 
sub test_print_to_output
{
   my $outputdir           = $_[0]; # DIRECTORY TO WRITE OUTPUT FILES IN
   my $output_fasta;                # OUTPUT FILE
   my $seq;                         # SEQUENCE
   my $errorcode;                   # RETURNED AS 0 BY A FUNCTION IF THERE IS NO ERROR
   my $errormsg;                    # RETURNED AS 'none' BY A FUNCTION IF THERE IS NO ERROR 
   my $expected_output_fasta;       # FILE CONTAINING EXPECTED CONTENTS OF $output_fasta
   my $differences;                 # DIFFERENCES BETWEEN $output_fasta AND $expected_output_fasta
   my $line;                        # 
   my $length_differences;          # LENGTH OF $differences
 
   ($output_fasta,$errorcode,$errormsg) = &make_filename($outputdir);
   open(OUTPUT,">$output_fasta") || die "ERROR: test_print_to_output: cannot open $output_fasta\n";
   close(OUTPUT);
   if ($errorcode != 0) { ($errorcode,$errormsg) = &print_error($errormsg,$errorcode,0); } 
   $seq                    = "AAAAATTTTTCCCCCGGGGGAAAAATTTTTCCCCCGGGGGAAAAATTTTTCCCCCGGGGGAAAAATTTTTCCCCCGGGGGAAAAATTTTTCCCCCGGGGGAAAAATTTTTCCCCCGGGGGAAAAATTTTTCCCCCGGGGG";
   ($errorcode,$errormsg)  = &print_to_output($output_fasta,$seq,"seq1");
   if ($errorcode != 0) { print STDERR "ERROR: test_print_to_output: failed test1\n"; exit;}
   ($expected_output_fasta,$errorcode,$errormsg) = &make_filename($outputdir);
   if ($errorcode != 0) { ($errorcode,$errormsg) = &print_error($errormsg,$errorcode,0); } 
   open(EXPECTED_OUTPUT,">$expected_output_fasta") || die "ERROR: test_print_to_output: cannot open $expected_output_fasta\n";
   print EXPECTED_OUTPUT ">seq1\n";
   print EXPECTED_OUTPUT "AAAAATTTTTCCCCCGGGGGAAAAATTTTTCCCCCGGGGGAAAAATTTTTCCCCCGGGGG\n";
   print EXPECTED_OUTPUT "AAAAATTTTTCCCCCGGGGGAAAAATTTTTCCCCCGGGGGAAAAATTTTTCCCCCGGGGG\n";
   print EXPECTED_OUTPUT "AAAAATTTTTCCCCCGGGGG\n";
   close(EXPECTED_OUTPUT); 
   $differences            = "";
   open(TEMP,"diff $output_fasta $expected_output_fasta |");
   while(<TEMP>)
   {
      $line                = $_;
      $differences         = $differences.$line;
   }
   close(TEMP);  
   $length_differences     = length($differences);
   if ($length_differences != 0) { print STDERR "ERROR: test_print_to_output: failed test1 (output_fasta $output_fasta expected_output_fasta $expected_output_fasta)\n"; exit;}
   system "rm -f $output_fasta";
   system "rm -f $expected_output_fasta"; 
 
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
   ($SEQ,$errorcode,$errormsg) = &read_assembly($assembly);
   if ($SEQ->{'seq1'} ne 'AAAAA' || $SEQ->{'seq2'} ne 'AAAAATTTTT' || $SEQ->{'seq3'} ne 'AAAAATTTTT' || defined($SEQ->{'seq4'}) || 
       $SEQ->{'seq5'} ne 'AAAAA' || $errorcode != 0) 
   { 
      print STDERR "ERROR: test_read_assembly: failed test1\n"; 
      exit;
   }
   system "rm -f $assembly";
 
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
   ($SEQ,$errorcode,$errormsg) = &read_assembly($assembly);
   if ($errorcode != 4) { print STDERR "ERROR: test_assembly: failed test2\n"; exit;}
   system "rm -f $assembly";
 
   $random_number          = rand();
   $assembly               = $outputdir."/tmp".$random_number;
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
# SUBROUTINE SYNOPSIS: read_assembly(): read sequences in a fasta file into a hash
 
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
