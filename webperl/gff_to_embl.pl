#!/usr/bin/env perl
#https://gist.github.com/avrilcoghlan/5031058
=head1 NAME

    gff_to_embl.pl

=head1 SYNOPSIS
 
    gff_to_embl.pl gff fasta outputdir exonerate cegma ratt augustus
        where gff is the gff file of gene predictions,   
              fasta is the fasta file for the genome assembly,
              outputdir is the output directory for writing output files,
              exonerate says whether the gff file is from exonerate (yes/no),
              cegma says whether the gff file is from cegma (yes/no),
              ratt says whether the gff file is from ratt (yes/no),
              augustus says whether the gff file is from augustus (yes/no).

=head1 DESCRIPTION

    
    This script takes an input gff file, and converts it to embl format, and writes the output 
    embl files for each scaffold that has gene predictions, in directory <outputdir>. The output embl 
    files are named after the scaffolds, ie. scaffold1.embl, scaffold2.embl, etc.
    This makes an embl file for each scaffold that has genes on it.

=head1 VERSION
  
    Perl script last edited 25-Feb-2013.

=head1 CONTACT

    alc@sanger.ac.uk (Avril Coghlan)

=cut

# 
# Perl script gff_to_embl.pl
# Written by Avril Coghlan (alc@sanger.ac.uk)
# 25-Feb-13.
# Last edited 25-Feb-2013.
# SCRIPT SYNOPSIS: gff_to_embl.pl: convert a gff file to embl format.
#
#------------------------------------------------------------------#

# CHECK IF THERE ARE THE CORRECT NUMBER OF COMMAND-LINE ARGUMENTS:

use strict;
use warnings;

my $num_args               = $#ARGV + 1;
if ($num_args != 7)
{
    print "Usage of gff_to_embl.pl\n\n";
    print "perl gff_to_embl.pl <gff> <fasta> <outputdir> <exonerate> <cegma> <ratt> <augustus>\n";
    print "where <gff> is the gff file of gene predictions,\n";
    print "      <fasta> is the fasta file for the genome assembly,\n";
    print "      <outputdir> is the output directory for writing output files,\n";
    print "      <exonerate> says whether the gff file is from exonerate (yes/no),\n";
    print "      <cegma> says whether the gff file is from cegma (yes/no),\n";
    print "      <ratt> says whether the gff file is from ratt,\n";
    print "      <augustus> says whether the gff file is from augustus\n"; 
    print "For example, >perl gff_to_embl.pl job1.gff job1.fasta\n";
    print "/lustre/scratch108/parasites/alc/ yes no no no\n";
    exit;
}

# FIND THE PATH TO THE INPUT GFF FILE:                     

my $gff                    = $ARGV[0];

# FIND THE PATH TO THE INPUT FASTA FILE FOR THE GENOME:

my $fasta                  = $ARGV[1];

# FIND THE DIRECTORY TO USE FOR OUTPUT FILES:      

my $outputdir              = $ARGV[2];

# FIND OUT WHETHER THE GFF FILE WAS FROM EXONERATE:

my $exonerate              = $ARGV[3];

# FIND OUT WHETHER THE GFF FILE WAS FROM CEGMA:

my $cegma                  = $ARGV[4];

# FIND OUT WHETHER THE GFF FILE WAS FROM RATT:

my $ratt                   = $ARGV[5];

# FIND OUT WHETHER THE GFF FILE WAS FROM AUGUSTUS:

my $augustus               = $ARGV[6];

#------------------------------------------------------------------#

# TEST SUBROUTINES: 

my $PRINT_TEST_DATA        = 0;   # SAYS WHETHER TO PRINT DATA USED DURING TESTING.
&test_print_error;
&test_read_scaffold_lengths($outputdir);  
&test_read_assembly($outputdir);
&test_read_gene_positions($outputdir); 
&test_count_bases;
print STDERR "Finished tests, proceeding to main program...\n";

#------------------------------------------------------------------#

# RUN THE MAIN PART OF THE CODE:

&run_main_program($outputdir,$gff,$fasta,$exonerate,$cegma,$ratt,$augustus);

print STDERR "FINISHED.\n";

#------------------------------------------------------------------#

# RUN THE MAIN PART OF THE CODE:

sub run_main_program {
   my $outputdir           = $_[0]; # DIRECTORY TO PUT OUTPUT FILES IN.
   my $gff                 = $_[1]; # THE INPUT GFF FILE
   my $fasta               = $_[2]; # THE INPUT GENOME FASTA FILE   
   my $exonerate           = $_[3]; # SAYS WHETHER THE GFF FILE WAS FROM EXONERATE.
   my $cegma               = $_[4]; # SAYS WHETHER THE GFF FILE WAS FROM CEGMA.
   my $ratt                = $_[5]; # SAYS WHETHER THE GFF FILE WAS FROM RATT.
   my $augustus            = $_[6]; # SAYS WHETHER THE GFF FILE WAS FROM AUGUSTUS.
   my $errorcode;                   # RETURNED AS 0 IF THERE IS NO ERROR.
   my $errormsg;                    # RETURNED AS 'none' IF THERE IS NO ERROR. 
   my $LEN;                         # HASH TABLE TO STORE THE LENGTH OF SCAFFOLDS
   my $GENES;                       # HASH TABLE TO STORE THE GENES ON SCAFFOLDS
   my $EXONS;                       # HASH TABLE TO STORE THE EXONS IN GENES 
   my $SEQ;                         # HASH TABLE TO STORE SCAFFOLD SEQUENCES 

   # READ IN THE LENGTHS OF SCAFFOLDS:
   ($LEN,$errorcode,$errormsg) = &read_scaffold_lengths($fasta);
   if ($errorcode != 0) { ($errorcode,$errormsg) = &print_error($errormsg,$errorcode,0); }
 
   # READ IN THE SEQUENCES FROM THE INPUT FASTA FILE:
   ($SEQ,$errorcode,$errormsg) = &read_assembly($fasta);
   if ($errorcode != 0) { ($errorcode,$errormsg) = &print_error($errormsg,$errorcode,0); }

   # READ IN THE POSITIONS OF GENES ON SCAFFOLDS:
   ($GENES,$EXONS,$errorcode,$errormsg) = &read_gene_positions($gff,$exonerate,$cegma,$ratt,$augustus); 
   if ($errorcode != 0) { ($errorcode,$errormsg) = &print_error($errormsg,$errorcode,0); }

   # READ IN THE FILE, AND CONVERT TO EMBL FORMAT:
   ($errorcode,$errormsg)  = &convert_gff_to_embl($outputdir,$GENES,$EXONS,$LEN,$SEQ);
   if ($errorcode != 0) { ($errorcode,$errormsg) = &print_error($errormsg,$errorcode,0); }
}

#------------------------------------------------------------------#

# TEST &read_gene_positions

sub test_read_gene_positions {
   my $outputdir           = $_[0]; # DIRECTORY TO WRITE OUTPUT FILES IN
   my $errorcode;                   # RETURNED AS 0 FROM A FUNCTION IF THERE IS NO ERROR
   my $errormsg;                    # RETURNED AS 'none' FROM A FUNCTION IF THERE IS NO ERROR
   my $GENES;                       # HASH TABLE TO STORE THE GENES ON SCAFFOLDS
   my $EXONS;                       # HASH TABLE TO STORE THE EXONS IN GENES 
   my $gff;                         # GFF FILE 
 
   ($gff,$errorcode,$errormsg) = &make_filename($outputdir);
   if ($errorcode != 0) { ($errorcode,$errormsg) = &print_error($errormsg,$errorcode,0); } 
   open(GFF,">$gff") || die "ERROR: test_read_gene_positions: cannot open $gff\n";
   print GFF "Ratt_Curated.Sratti_scf00001_Chr00	RATT	gene	676104	677702	.	+	.	ID=ratti_train3\n";
   print GFF "Ratt_Curated.Sratti_scf00001_Chr00	RATT	mRNA	676104	677702	.	+	.	ID=ratti_train3;Parent=ratti_train3\n";
   print GFF "Ratt_Curated.Sratti_scf00001_Chr00	RATT	CDS	676104	676526	.	+	.	ID=ratti_train3:exon:1;Parent=ratti_train3\n";
   print GFF "Ratt_Curated.Sratti_scf00001_Chr00	RATT	CDS	676766	677515	.	+	.	ID=ratti_train3:exon:2;Parent=ratti_train3\n";
   print GFF "Ratt_Curated.Sratti_scf00001_Chr00	RATT	CDS	677568	677702	.	+	.	ID=ratti_train3:exon:3;Parent=ratti_train3\n";
   print GFF "Ratt_Curated.Sratti_scf00001_Chr00	RATT	gene	699685	700943	.	-	.	ID=ratti_train4\n";
   print GFF "Ratt_Curated.Sratti_scf00001_Chr00	RATT	mRNA	699685	700943	.	-	.	ID=ratti_train4;Parent=ratti_train4\n";
   print GFF "Ratt_Curated.Sratti_scf00001_Chr00	RATT	CDS	699685	699951	.	-	.	ID=ratti_train4:exon:1;Parent=ratti_train4\n";
   print GFF "Ratt_Curated.Sratti_scf00001_Chr00	RATT	CDS	700002	700943	.	-	.	ID=ratti_train4:exon:2;Parent=ratti_train4\n";
   print GFF "Ratt_Curated.Sratti_scf00001_Chr01	RATT	gene	1650617	1651867	.	+	.	ID=ratti_train414\n";
   print GFF "Ratt_Curated.Sratti_scf00001_Chr01	RATT	mRNA	1650617	1651867	.	+	.	ID=ratti_train414;Parent=ratti_train414\n";
   print GFF "Ratt_Curated.Sratti_scf00001_Chr01	RATT	CDS	1650617	1650755	.	+	.	ID=ratti_train414:exon:1;Parent=ratti_train414\n";
   print GFF "Ratt_Curated.Sratti_scf00001_Chr01	RATT	CDS	1650801	1651867	.	+	.	ID=ratti_train414:exon:2;Parent=ratti_train414\n";
   close(GFF);
   # READ IN THE POSITIONS OF GENES ON SCAFFOLDS:
   ($GENES,$EXONS,$errorcode,$errormsg)  = &read_gene_positions($gff,'no','no','no','no');
   if ($errorcode != 0) { print STDERR "ERROR: test_read_gene_positions: failed test1\n"; exit;} 
   if ($GENES->{"Ratt_Curated.Sratti_scf00001_Chr00"} ne "ratti_train3,ratti_train4") { print STDERR "ERROR: test_read_gene_positions: failed test1b\n"; exit;}
   if ($GENES->{"Ratt_Curated.Sratti_scf00001_Chr01"} ne "ratti_train414") { print STDERR "ERROR: test_read_gene_positions: failed test1c\n"; exit;}
   if ($EXONS->{"ratti_train3"} ne "676104=676526=+,676766=677515=+,677568=677702=+") { print STDERR "ERROR: test_read_gene_positions: failed test1d\n"; exit;}
   if ($EXONS->{"ratti_train4"} ne "699685=699951=-,700002=700943=-") { print STDERR "ERROR: test_read_gene_positions: failed test1e\n"; exit;}
   if ($EXONS->{"ratti_train414"} ne "1650617=1650755=+,1650801=1651867=+") { print STDERR "ERROR: test_read_gene_positions: failed test1f\n"; exit;}
   system "rm -f $gff";

   # TEST ERRORCODE=14 (READING IN A 'cds' LINE BEFORE A 'gene' LINE):
   ($gff,$errorcode,$errormsg) = &make_filename($outputdir);
   if ($errorcode != 0) { ($errorcode,$errormsg) = &print_error($errormsg,$errorcode,0); } 
   open(GFF,">$gff") || die "ERROR: test_read_gene_positions: cannot open $gff\n";
   print GFF "Ratt_Curated.Sratti_scf00001_Chr00	RATT	CDS	676766	677515	.	+	.	ID=ratti_train3:exon:2;Parent=ratti_train3\n";
   print GFF "Ratt_Curated.Sratti_scf00001_Chr00	RATT	gene	676104	677702	.	+	.	ID=ratti_train3\n";
   print GFF "Ratt_Curated.Sratti_scf00001_Chr00	RATT	mRNA	676104	677702	.	+	.	ID=ratti_train3;Parent=ratti_train3\n";
   print GFF "Ratt_Curated.Sratti_scf00001_Chr00	RATT	CDS	676104	676526	.	+	.	ID=ratti_train3:exon:1;Parent=ratti_train3\n";
   print GFF "Ratt_Curated.Sratti_scf00001_Chr00	RATT	CDS	677568	677702	.	+	.	ID=ratti_train3:exon:3;Parent=ratti_train3\n";
   close(GFF);
   # READ IN THE POSITIONS OF GENES ON SCAFFOLDS:
   ($GENES,$EXONS,$errorcode,$errormsg)  = &read_gene_positions($gff,'no','no','no','no');
   if ($errorcode != 14) { print STDERR "ERROR: test_read_gene_positions: failed test2 (errorcode $errorcode errormsg $errormsg)\n"; exit;} 
   system "rm -f $gff";
 
   # TEST ERRORCODE=16 ($exonerate IS NOT yes/no):
   ($gff,$errorcode,$errormsg) = &make_filename($outputdir);
   if ($errorcode != 0) { ($errorcode,$errormsg) = &print_error($errormsg,$errorcode,0); } 
   open(GFF,">$gff") || die "ERROR: test_read_gene_positions: cannot open $gff\n";
   print GFF "Ratt_Curated.Sratti_scf00001_Chr00	RATT	gene	676104	677702	.	+	.	ID=ratti_train3\n";
   print GFF "Ratt_Curated.Sratti_scf00001_Chr00	RATT	mRNA	676104	677702	.	+	.	ID=ratti_train3;Parent=ratti_train3\n";
   print GFF "Ratt_Curated.Sratti_scf00001_Chr00	RATT	CDS	676104	676526	.	+	.	ID=ratti_train3:exon:1;Parent=ratti_train3\n";
   print GFF "Ratt_Curated.Sratti_scf00001_Chr00	RATT	CDS	676766	677515	.	+	.	ID=ratti_train3:exon:2;Parent=ratti_train3\n";
   print GFF "Ratt_Curated.Sratti_scf00001_Chr00	RATT	CDS	677568	677702	.	+	.	ID=ratti_train3:exon:3;Parent=ratti_train3\n";
   print GFF "Ratt_Curated.Sratti_scf00001_Chr00	RATT	gene	699685	700943	.	-	.	ID=ratti_train4\n";
   print GFF "Ratt_Curated.Sratti_scf00001_Chr00	RATT	mRNA	699685	700943	.	-	.	ID=ratti_train4;Parent=ratti_train4\n";
   print GFF "Ratt_Curated.Sratti_scf00001_Chr00	RATT	CDS	699685	699951	.	-	.	ID=ratti_train4:exon:1;Parent=ratti_train4\n";
   print GFF "Ratt_Curated.Sratti_scf00001_Chr00	RATT	CDS	700002	700943	.	-	.	ID=ratti_train4:exon:2;Parent=ratti_train4\n";
   print GFF "Ratt_Curated.Sratti_scf00001_Chr01	RATT	gene	1650617	1651867	.	+	.	ID=ratti_train414\n";
   print GFF "Ratt_Curated.Sratti_scf00001_Chr01	RATT	mRNA	1650617	1651867	.	+	.	ID=ratti_train414;Parent=ratti_train414\n";
   print GFF "Ratt_Curated.Sratti_scf00001_Chr01	RATT	CDS	1650617	1650755	.	+	.	ID=ratti_train414:exon:1;Parent=ratti_train414\n";
   print GFF "Ratt_Curated.Sratti_scf00001_Chr01	RATT	CDS	1650801	1651867	.	+	.	ID=ratti_train414:exon:2;Parent=ratti_train414\n";
   close(GFF);
   # READ IN THE POSITIONS OF GENES ON SCAFFOLDS:
   ($GENES,$EXONS,$errorcode,$errormsg)  = &read_gene_positions($gff,'hello','no','no','no');
   if ($errorcode != 16) { print STDERR "ERROR: test_read_gene_positions: failed test3\n"; exit;} 
   system "rm -f $gff";

   # TEST AN EXONERATE FORMAT GFF:
   ($gff,$errorcode,$errormsg) = &make_filename($outputdir);
   if ($errorcode != 0) { ($errorcode,$errormsg) = &print_error($errormsg,$errorcode,0); } 
   open(GFF,">$gff") || die "ERROR: test_read_gene_positions: cannot open $gff\n";
   print GFF "PTRK.scaffold.00383.64833	exonerate:coding2genome	gene	48175	48813	3992	-	.	gene_id 0 ; sequence PTC00001_1 ; gene_orientation +\n";
   print GFF "PTRK.scaffold.00383.64833	exonerate:coding2genome	cds	48675	48813	.	-	.	\n";
   print GFF "PTRK.scaffold.00383.64833	exonerate:coding2genome	exon	48675	48813	.	-	.	insertions 0 ; deletions 0\n";
   print GFF "PTRK.scaffold.00383.64833	exonerate:coding2genome	splice5	48673	48674	.	-	.	intron_id 1 ; splice_site \"GT\"\n";
   print GFF "PTRK.scaffold.00383.64833	exonerate:coding2genome	intron	48625	48674	.	-	.	intron_id 1\n";
   print GFF "PTRK.scaffold.00383.64833	exonerate:coding2genome	splice3	48625	48626	.	-	.	intron_id 0 ; splice_site \"AG\"\n";
   print GFF "PTRK.scaffold.00383.64833	exonerate:coding2genome	cds	48530	48624	.	-	.	\n";
   print GFF "PTRK.scaffold.00383.64833	exonerate:coding2genome	exon	48530	48624	.	-	.	insertions 0 ; deletions 0\n";
   print GFF "PTRK.scaffold.00383.64833	exonerate:coding2genome	splice5	48528	48529	.	-	.	intron_id 2 ; splice_site \"GT\"\n";
   print GFF "PTRK.scaffold.00383.64833	exonerate:coding2genome	intron	48482	48529	.	-	.	intron_id 2\n";
   print GFF "PTRK.scaffold.00383.64833	exonerate:coding2genome	splice3	48482	48483	.	-	.	intron_id 1 ; splice_site \"AG\"\n";
   print GFF "PTRK.scaffold.00383.64833	exonerate:coding2genome	cds	48334	48481	.	-	.	\n";
   print GFF "PTRK.scaffold.00383.64833	exonerate:coding2genome	exon	48334	48481	.	-	.	insertions 0 ; deletions 0\n";
   print GFF "PTRK.scaffold.00383.64833	exonerate:coding2genome	splice5	48332	48333	.	-	.	intron_id 3 ; splice_site \"GT\"\n";
   print GFF "PTRK.scaffold.00383.64833	exonerate:coding2genome	intron	48282	48333	.	-	.	intron_id 3\n";
   print GFF "PTRK.scaffold.00383.64833	exonerate:coding2genome	splice3	48282	48283	.	-	.	intron_id 2 ; splice_site \"AG\"\n";
   print GFF "PTRK.scaffold.00383.64833	exonerate:coding2genome	cds	48175	48281	.	-	.	\n";
   print GFF "PTRK.scaffold.00383.64833	exonerate:coding2genome	exon	48175	48281	.	-	.	insertions 0 ; deletions 0\n";
   print GFF "PTRK.scaffold.00010.1389350	exonerate:coding2genome	gene	647860	696655	3969	-	.	gene_id 0 ; sequence PTC00423_1 ; gene_orientation +\n";
   print GFF "PTRK.scaffold.00010.1389350	exonerate:coding2genome	cds	696389	696655	.	-	.	\n";
   print GFF "PTRK.scaffold.00010.1389350	exonerate:coding2genome	exon	696389	696655	.	-	.	insertions 0 ; deletions 0\n";
   print GFF "PTRK.scaffold.00010.1389350	exonerate:coding2genome	splice5	696387	696388	.	-	.	intron_id 1 ; splice_site \"GT\"\n";
   print GFF "PTRK.scaffold.00010.1389350	exonerate:coding2genome	intron	696342	696388	.	-	.	intron_id 1\n";
   print GFF "PTRK.scaffold.00010.1389350	exonerate:coding2genome	splice3	696342	696343	.	-	.	intron_id 0 ; splice_site \"AG\"\n";
   print GFF "PTRK.scaffold.00010.1389350	exonerate:coding2genome	cds	696135	696341	.	-	.	\n";	
   print GFF "PTRK.scaffold.00010.1389350	exonerate:coding2genome	exon	696135	696341	.	-	.	insertions 0 ; deletions 0\n";
   print GFF "PTRK.scaffold.00010.1389350	exonerate:coding2genome	splice5	696133	696134	.	-	.	intron_id 2 ; splice_site \"CT\"\n";
   print GFF "PTRK.scaffold.00010.1389350	exonerate:coding2genome	intron	647869	696134	.	-	.	intron_id 2\n";
   print GFF "PTRK.scaffold.00010.1389350	exonerate:coding2genome	splice3	647869	647870	.	-	.	intron_id 1 ; splice_site \"CA\"\n";
   print GFF "PTRK.scaffold.00010.1389350	exonerate:coding2genome	cds	647860	647868	.	-	.	\n";
   print GFF "PTRK.scaffold.00010.1389350	exonerate:coding2genome	exon	647860	647868	.	-	.	insertions 0 ; deletions 0\n";
   close(GFF);
   # READ IN THE POSITIONS OF GENES ON SCAFFOLDS:
   ($GENES,$EXONS,$errorcode,$errormsg)  = &read_gene_positions($gff,'yes','no','no','no');
   if ($errorcode != 0) { print STDERR "ERROR: test_read_gene_positions: failed test4\n"; exit;} 
   if ($GENES->{"PTRK.scaffold.00383.64833"} ne "PTC00001_1") { print STDERR "ERROR: test_read_gene_positions: failed test4b\n"; exit;}
   if ($GENES->{"PTRK.scaffold.00010.1389350"} ne "PTC00423_1") { print STDERR "ERROR: test_read_gene_positions: failed test4c\n"; exit;}
   if ($EXONS->{"PTC00001_1"} ne "48675=48813=-,48530=48624=-,48334=48481=-,48175=48281=-") { print STDERR "ERROR: test_read_gene_positions: failed test4d\n"; exit;}
   if ($EXONS->{"PTC00423_1"} ne "696389=696655=-,696135=696341=-,647860=647868=-") { print STDERR "ERROR: test_read_gene_positions: failed test4e\n"; exit;}
   system "rm -f $gff";

   # TEST ERRORCODE=17 ($cegma IS NOT yes/no):
   ($gff,$errorcode,$errormsg) = &make_filename($outputdir);
   if ($errorcode != 0) { ($errorcode,$errormsg) = &print_error($errormsg,$errorcode,0); } 
   open(GFF,">$gff") || die "ERROR: test_read_gene_positions: cannot open $gff\n";
   print GFF "Ratt_Curated.Sratti_scf00001_Chr00	RATT	gene	676104	677702	.	+	.	ID=ratti_train3\n";
   print GFF "Ratt_Curated.Sratti_scf00001_Chr00	RATT	mRNA	676104	677702	.	+	.	ID=ratti_train3;Parent=ratti_train3\n";
   print GFF "Ratt_Curated.Sratti_scf00001_Chr00	RATT	CDS	676104	676526	.	+	.	ID=ratti_train3:exon:1;Parent=ratti_train3\n";
   print GFF "Ratt_Curated.Sratti_scf00001_Chr00	RATT	CDS	676766	677515	.	+	.	ID=ratti_train3:exon:2;Parent=ratti_train3\n";
   print GFF "Ratt_Curated.Sratti_scf00001_Chr00	RATT	CDS	677568	677702	.	+	.	ID=ratti_train3:exon:3;Parent=ratti_train3\n";
   print GFF "Ratt_Curated.Sratti_scf00001_Chr00	RATT	gene	699685	700943	.	-	.	ID=ratti_train4\n";
   print GFF "Ratt_Curated.Sratti_scf00001_Chr00	RATT	mRNA	699685	700943	.	-	.	ID=ratti_train4;Parent=ratti_train4\n";
   print GFF "Ratt_Curated.Sratti_scf00001_Chr00	RATT	CDS	699685	699951	.	-	.	ID=ratti_train4:exon:1;Parent=ratti_train4\n";
   print GFF "Ratt_Curated.Sratti_scf00001_Chr00	RATT	CDS	700002	700943	.	-	.	ID=ratti_train4:exon:2;Parent=ratti_train4\n";
   print GFF "Ratt_Curated.Sratti_scf00001_Chr01	RATT	gene	1650617	1651867	.	+	.	ID=ratti_train414\n";
   print GFF "Ratt_Curated.Sratti_scf00001_Chr01	RATT	mRNA	1650617	1651867	.	+	.	ID=ratti_train414;Parent=ratti_train414\n";
   print GFF "Ratt_Curated.Sratti_scf00001_Chr01	RATT	CDS	1650617	1650755	.	+	.	ID=ratti_train414:exon:1;Parent=ratti_train414\n";
   print GFF "Ratt_Curated.Sratti_scf00001_Chr01	RATT	CDS	1650801	1651867	.	+	.	ID=ratti_train414:exon:2;Parent=ratti_train414\n";
   close(GFF);
   # READ IN THE POSITIONS OF GENES ON SCAFFOLDS:
   ($GENES,$EXONS,$errorcode,$errormsg)  = &read_gene_positions($gff,'no','hello','no','no');
   if ($errorcode != 17) { print STDERR "ERROR: test_read_gene_positions: failed test5\n"; exit;} 
   system "rm -f $gff";

   # TEST ERRORCODE=51 ($cegma AND $exonerate ARE 'yes'):
   ($gff,$errorcode,$errormsg) = &make_filename($outputdir);
   if ($errorcode != 0) { ($errorcode,$errormsg) = &print_error($errormsg,$errorcode,0); } 
   open(GFF,">$gff") || die "ERROR: test_read_gene_positions: cannot open $gff\n";
   print GFF "Ratt_Curated.Sratti_scf00001_Chr00	RATT	gene	676104	677702	.	+	.	ID=ratti_train3\n";
   print GFF "Ratt_Curated.Sratti_scf00001_Chr00	RATT	mRNA	676104	677702	.	+	.	ID=ratti_train3;Parent=ratti_train3\n";
   print GFF "Ratt_Curated.Sratti_scf00001_Chr00	RATT	CDS	676104	676526	.	+	.	ID=ratti_train3:exon:1;Parent=ratti_train3\n";
   print GFF "Ratt_Curated.Sratti_scf00001_Chr00	RATT	CDS	676766	677515	.	+	.	ID=ratti_train3:exon:2;Parent=ratti_train3\n";
   print GFF "Ratt_Curated.Sratti_scf00001_Chr00	RATT	CDS	677568	677702	.	+	.	ID=ratti_train3:exon:3;Parent=ratti_train3\n";
   print GFF "Ratt_Curated.Sratti_scf00001_Chr00	RATT	gene	699685	700943	.	-	.	ID=ratti_train4\n";
   print GFF "Ratt_Curated.Sratti_scf00001_Chr00	RATT	mRNA	699685	700943	.	-	.	ID=ratti_train4;Parent=ratti_train4\n";
   print GFF "Ratt_Curated.Sratti_scf00001_Chr00	RATT	CDS	699685	699951	.	-	.	ID=ratti_train4:exon:1;Parent=ratti_train4\n";
   print GFF "Ratt_Curated.Sratti_scf00001_Chr00	RATT	CDS	700002	700943	.	-	.	ID=ratti_train4:exon:2;Parent=ratti_train4\n";
   print GFF "Ratt_Curated.Sratti_scf00001_Chr01	RATT	gene	1650617	1651867	.	+	.	ID=ratti_train414\n";
   print GFF "Ratt_Curated.Sratti_scf00001_Chr01	RATT	mRNA	1650617	1651867	.	+	.	ID=ratti_train414;Parent=ratti_train414\n";
   print GFF "Ratt_Curated.Sratti_scf00001_Chr01	RATT	CDS	1650617	1650755	.	+	.	ID=ratti_train414:exon:1;Parent=ratti_train414\n";
   print GFF "Ratt_Curated.Sratti_scf00001_Chr01	RATT	CDS	1650801	1651867	.	+	.	ID=ratti_train414:exon:2;Parent=ratti_train414\n";
   close(GFF);
   # READ IN THE POSITIONS OF GENES ON SCAFFOLDS:
   ($GENES,$EXONS,$errorcode,$errormsg)  = &read_gene_positions($gff,'yes','yes','no','no');
   if ($errorcode != 51) { print STDERR "ERROR: test_read_gene_positions: failed test6\n"; exit;} 
   system "rm -f $gff";

   # TEST A CEGMA FORMAT GFF:
   ($gff,$errorcode,$errormsg) = &make_filename($outputdir);
   if ($errorcode != 0) { ($errorcode,$errormsg) = &print_error($errormsg,$errorcode,0); } 
   open(GFF,">$gff") || die "ERROR: test_read_gene_positions: cannot open $gff\n";
   print GFF "PTRK.scaffold.00789.21789	cegma	Single	18301	18456	72.81	+	0	KOG0002.14\n";
   print GFF "PTRK.scaffold.00114.272491	cegma	First	53275	53413	82.29	-	0	KOG0003.3\n";
   print GFF "PTRK.scaffold.00114.272491	cegma	Terminal	52977	53224	151.50	-	2	KOG0003.3\n";
   print GFF "PTRK.contig.00177.164151	cegma	First	58256	58364	62.85	+	0	KOG0018.7\n";
   print GFF "PTRK.contig.00177.164151	cegma	Internal	58419	61837	1229.98	+	2	KOG0018.7\n";
   print GFF "PTRK.contig.00177.164151	cegma	Terminal	61879	62046	60.54	+	0	KOG0018.7\n";
   close(GFF);
   # READ IN THE POSITIONS OF GENES ON SCAFFOLDS:
   ($GENES,$EXONS,$errorcode,$errormsg)  = &read_gene_positions($gff,'no','yes','no','no');
   if ($errorcode != 0) { print STDERR "ERROR: test_read_gene_positions: failed test7 (errorcode $errorcode errormsg $errormsg)\n"; exit;} 
   if ($GENES->{"PTRK.scaffold.00789.21789"} ne "KOG0002.14") { print STDERR "ERROR: test_read_gene_positions: failed test7b\n"; exit;}
   if ($GENES->{"PTRK.scaffold.00114.272491"} ne "KOG0003.3") { print STDERR "ERROR: test_read_gene_positions: failed test7c\n"; exit;}
   if ($GENES->{"PTRK.contig.00177.164151"} ne "KOG0018.7") { print STDERR "ERROR: test_read_gene_positions: failed test7d\n"; exit;}
   if ($EXONS->{"KOG0002.14"} ne "18301=18456=+") { print STDERR "ERROR: test_read_gene_positions: failed test7e\n"; exit;}
   if ($EXONS->{"KOG0003.3"} ne "53275=53413=-,52977=53224=-") { print STDERR "ERROR: test_read_gene_positions: failed test7f\n"; exit;}
   if ($EXONS->{"KOG0018.7"} ne "58256=58364=+,58419=61837=+,61879=62046=+") { print STDERR "ERROR: test_read_gene_positions: failed test7g\n"; exit;}
   system "rm -f $gff";

   # TEST ERRORCODE=21 ($ratt IS NOT 'yes'/'no'):
   ($gff,$errorcode,$errormsg) = &make_filename($outputdir);
   if ($errorcode != 0) { ($errorcode,$errormsg) = &print_error($errormsg,$errorcode,0); }
   open(GFF,">$gff") || die "ERROR: test_read_gene_positions: cannot open $gff\n";
   print GFF "PTRK.scaffold.00789.21789 cegma   Single  18301   18456   72.81   +       0       KOG0002.14\n";
   print GFF "PTRK.scaffold.00114.272491        cegma   First   53275   53413   82.29   -       0       KOG0003.3\n";
   print GFF "PTRK.scaffold.00114.272491        cegma   Terminal        52977   53224   151.50  -       2       KOG0003.3\n";
   print GFF "PTRK.contig.00177.164151  cegma   First   58256   58364   62.85   +       0       KOG0018.7\n";
   print GFF "PTRK.contig.00177.164151  cegma   Internal        58419   61837   1229.98 +       2       KOG0018.7\n";
   print GFF "PTRK.contig.00177.164151  cegma   Terminal        61879   62046   60.54   +       0       KOG0018.7\n";
   close(GFF);
   # READ IN THE POSITIONS OF GENES ON SCAFFOLDS:
   ($GENES,$EXONS,$errorcode,$errormsg)  = &read_gene_positions($gff,'no','yes','hello','no');
   if ($errorcode != 21) { print STDERR "ERROR: test_read_gene_positions: failed test8 (errorcode $errorcode errormsg $errormsg)\n"; exit;}
   system "rm -f $gff";

   # TEST ERRORCODE=51 ($ratt IS 'yes', $cemga IS 'yes'):
   ($gff,$errorcode,$errormsg) = &make_filename($outputdir);
   if ($errorcode != 0) { ($errorcode,$errormsg) = &print_error($errormsg,$errorcode,0); }
   open(GFF,">$gff") || die "ERROR: test_read_gene_positions: cannot open $gff\n";
   print GFF "PTRK.scaffold.00789.21789 cegma   Single  18301   18456   72.81   +       0       KOG0002.14\n";
   print GFF "PTRK.scaffold.00114.272491        cegma   First   53275   53413   82.29   -       0       KOG0003.3\n";
   print GFF "PTRK.scaffold.00114.272491        cegma   Terminal        52977   53224   151.50  -       2       KOG0003.3\n";
   print GFF "PTRK.contig.00177.164151  cegma   First   58256   58364   62.85   +       0       KOG0018.7\n";
   print GFF "PTRK.contig.00177.164151  cegma   Internal        58419   61837   1229.98 +       2       KOG0018.7\n";
   print GFF "PTRK.contig.00177.164151  cegma   Terminal        61879   62046   60.54   +       0       KOG0018.7\n";
   close(GFF);
   # READ IN THE POSITIONS OF GENES ON SCAFFOLDS:
   ($GENES,$EXONS,$errorcode,$errormsg)  = &read_gene_positions($gff,'no','yes','yes','no');
   if ($errorcode != 51) { print STDERR "ERROR: test_read_gene_positions: failed test9 (errorcode $errorcode errormsg $errormsg)\n"; exit;}
   system "rm -f $gff";

   # TEST ERRORCODE=51 ($ratt IS 'yes', $exonerate IS 'yes'):
   ($gff,$errorcode,$errormsg) = &make_filename($outputdir);
   if ($errorcode != 0) { ($errorcode,$errormsg) = &print_error($errormsg,$errorcode,0); }
   open(GFF,">$gff") || die "ERROR: test_read_gene_positions: cannot open $gff\n";
   print GFF "PTRK.scaffold.00789.21789 cegma   Single  18301   18456   72.81   +       0       KOG0002.14\n";
   print GFF "PTRK.scaffold.00114.272491        cegma   First   53275   53413   82.29   -       0       KOG0003.3\n";
   print GFF "PTRK.scaffold.00114.272491        cegma   Terminal        52977   53224   151.50  -       2       KOG0003.3\n";
   print GFF "PTRK.contig.00177.164151  cegma   First   58256   58364   62.85   +       0       KOG0018.7\n";
   print GFF "PTRK.contig.00177.164151  cegma   Internal        58419   61837   1229.98 +       2       KOG0018.7\n";
   print GFF "PTRK.contig.00177.164151  cegma   Terminal        61879   62046   60.54   +       0       KOG0018.7\n";
   close(GFF);
   # READ IN THE POSITIONS OF GENES ON SCAFFOLDS:
   ($GENES,$EXONS,$errorcode,$errormsg)  = &read_gene_positions($gff,'yes','no','yes','no');
   if ($errorcode != 51) { print STDERR "ERROR: test_read_gene_positions: failed test10 (errorcode $errorcode errormsg $errormsg)\n"; exit;}
   system "rm -f $gff";

   # TEST WITH AUGUSTUS EXAMPLE:
   ($gff,$errorcode,$errormsg) = &make_filename($outputdir);
   if ($errorcode != 0) { ($errorcode,$errormsg) = &print_error($errormsg,$errorcode,0); } 
   open(GFF,">$gff") || die "ERROR: test_read_gene_positions: cannot open $gff\n";
   print GFF "Sratt_Chr00_000001\tAUGUSTUS\tgene\t26596\t27840\t3.21\t-\t.\tg10772\n";
   print GFF "Sratt_Chr00_000001\tAUGUSTUS\ttranscript\t26596\t27348\t0.98\t-\t.\tg10772.t1\n";
   print GFF "Sratt_Chr00_000001\tAUGUSTUS\tstop_codon\t26596\t26598\t.\t-\t0\ttranscript_id \"g10772.t1\"; gene_id \"g10772\";\n";
   print GFF "Sratt_Chr00_000001\tAUGUSTUS\tCDS\t26596\t27348\t0.98\t-\t0\ttranscript_id \"g10772.t1\"; gene_id \"g10772\";\n";
   print GFF "Sratt_Chr00_000001\tAUGUSTUS\tstart_codon\t27346\t27348\t.\t-\t0\ttranscript_id \"g10772.t1\"; gene_id \"g10772\";\n";
   print GFF "Sratt_Chr00_000001\tAUGUSTUS\ttranscript\t26596\t27840\t0.7\t-\t.\tg10772.t2\n";
   print GFF "Sratt_Chr00_000001\tAUGUSTUS\tstop_codon\t26596\t26598\t.\t-\t0\ttranscript_id \"g10772.t2\"; gene_id \"g10772\";\n";
   print GFF "Sratt_Chr00_000001\tAUGUSTUS\tintron\t27355\t27391\t1\t-\t.\ttranscript_id \"g10772.t2\"; gene_id \"g10772\";\n";
   print GFF "Sratt_Chr00_000001\tAUGUSTUS\tintron\t27518\t27558\t1\t-\t.\ttranscript_id \"g10772.t2\"; gene_id \"g10772\";\n";
   print GFF "Sratt_Chr00_000001\tAUGUSTUS\tCDS\t26596\t27354\t0.83\t-\t0\ttranscript_id \"g10772.t2\"; gene_id \"g10772\";\n";
   print GFF "Sratt_Chr00_000001\tAUGUSTUS\tCDS\t27392\t27517\t1\t-\t0\ttranscript_id \"g10772.t2\"; gene_id \"g10772\";\n";
   print GFF "Sratt_Chr00_000001\tAUGUSTUS\tCDS\t27559\t27840\t0.83\t-\t0\ttranscript_id \"g10772.t2\"; gene_id \"g10772\";\n";
   print GFF "Sratt_Chr00_000001\tAUGUSTUS\tstart_codon\t27838\t27840\t.\t-\t0\ttranscript_id \"g10772.t2\"; gene_id \"g10772\";\n";
   print GFF "Sratt_Chr00_000001\tAUGUSTUS\ttranscript\t26596\t27315\t0.76\t-\t.\tg10772.t3\n";
   print GFF "Sratt_Chr00_000001\tAUGUSTUS\tstop_codon\t26596\t26598\t.\t-\t0\ttranscript_id \"g10772.t3\"; gene_id \"g10772\";\n";
   print GFF "Sratt_Chr00_000001\tAUGUSTUS\tCDS\t26596\t27315\t0.76\t-\t0\ttranscript_id \"g10772.t3\"; gene_id \"g10772\";\n";
   print GFF "Sratt_Chr00_000001\tAUGUSTUS\tstart_codon\t27313\t27315\t.\t-\t0\ttranscript_id \"g10772.t3\"; gene_id \"g10772\";\n";
   print GFF "Sratt_Chr00_000001\tAUGUSTUS\ttranscript\t26596\t27840\t0.77\t-\t.\tg10772.t4\n";
   print GFF "Sratt_Chr00_000001\tAUGUSTUS\tstop_codon\t26596\t26598\t.\t-\t0\ttranscript_id \"g10772.t4\"; gene_id \"g10772\";\n";
   print GFF "Sratt_Chr00_000001\tAUGUSTUS\tintron\t27252\t27305\t0.96\t-\t.\ttranscript_id \"g10772.t4\"; gene_id \"g10772\";\n";
   print GFF "Sratt_Chr00_000001\tAUGUSTUS\tintron\t27355\t27391\t1\t-\t.\ttranscript_id \"g10772.t4\"; gene_id \"g10772\";\n";
   print GFF "Sratt_Chr00_000001\tAUGUSTUS\tintron\t27518\t27558\t1\t-\t.\ttranscript_id \"g10772.t4\"; gene_id \"g10772\";\n";
   print GFF "Sratt_Chr00_000001\tAUGUSTUS\tCDS\t26596\t27251\t0.96\t-\t2\ttranscript_id \"g10772.t4\"; gene_id \"g10772\";\n";
   print GFF "Sratt_Chr00_000001\tAUGUSTUS\tCDS\t27306\t27354\t0.96\t-\t0\ttranscript_id \"g10772.t4\"; gene_id \"g10772\";\n";
   print GFF "Sratt_Chr00_000001\tAUGUSTUS\tCDS\t27392\t27517\t1\t-\t0\ttranscript_id \"g10772.t4\"; gene_id \"g10772\";\n";
   print GFF "Sratt_Chr00_000001\tAUGUSTUS\tCDS\t27559\t27840\t0.8\t-\t0\ttranscript_id \"g10772.t4\"; gene_id \"g10772\";\n";
   print GFF "Sratt_Chr00_000001\tAUGUSTUS\tstart_codon\t27838\t27840\t.\t-\t0\ttranscript_id \"g10772.t4\"; gene_id \"g10772\";\n";
   close(GFF);
   # READ IN THE POSITIONS OF GENES ON SCAFFOLDS:
   ($GENES,$EXONS,$errorcode,$errormsg)  = &read_gene_positions($gff,'no','no','no','yes');
   if ($errorcode != 0) { print STDERR "ERROR: test_read_gene_positions: failed test11\n"; exit;} 
   if ($GENES->{"Sratt_Chr00_000001"} ne "g10772.t1,g10772.t2,g10772.t3,g10772.t4") { print STDERR "ERROR: test_read_gene_positions: failed test11b\n"; exit;}
   if ($EXONS->{"g10772.t1"} ne "26596=27348=-") { print STDERR "ERROR: test_read_gene_positions: failed test11c\n"; exit;}
   if ($EXONS->{"g10772.t2"} ne "26596=27354=-,27392=27517=-,27559=27840=-") { print STDERR "ERROR: test_read_gene_positions: failed test11d\n"; exit;}
   if ($EXONS->{"g10772.t3"} ne "26596=27315=-") { print STDERR "ERROR: test_read_gene_positions: failed test11e\n"; exit;}
   if ($EXONS->{"g10772.t4"} ne "26596=27251=-,27306=27354=-,27392=27517=-,27559=27840=-") { print STDERR "ERROR: test_read_gene_positions: failed test11f\n"; exit;}
   system "rm -f $gff";

}

#------------------------------------------------------------------#

# READ IN THE POSITIONS OF GENES ON SCAFFOLDS:

sub read_gene_positions
{
   my $gff                 = $_[0]; # GFF FILE
   my $exonerate           = $_[1]; # SAYS WHETHER THE GFF FILE IS FROM EXONERATE 
   my $cegma               = $_[2]; # SAYS WHETHER THE GFF FILE IS FROM CEGMA 
   my $ratt                = $_[3]; # SAYS WHETHER THE GFF FILE IS FROM RATT
   my $augustus            = $_[4]; # SAYS WHETHER THE GFF FILE IS FROM AUGUSTUS
   my $errorcode           = 0;     # RETURNED AS 0 IF THERE IS NO ERROR
   my $errormsg            = "none";# RETURNED AS 'none' IF THERE IS NO ERROR
   my $line;                        # 
   my @temp;                        # 
   my $scaffold;                    # SCAFFOLD NAME
   my $feature;                     # FEATURE TYPE
   my $start;                       # START 
   my $end;                         # END 
   my $strand;                      # STRAND
   my $name                = "none";# FEATURE NAME
   my $score;                       # FEATURE SCORE
   my %EXONS               = ();    # HASH TABLE TO STORE EXONS IN A GENE
   my %GENES               = ();    # HASH TABLE TO STORE GENES ON A SCAFFOLD 
   my $exon;                        # COORDINATES OF AN EXON 
   my %SEEN                = ();    # HASH TABLE TO RECORD WHETHER WE SAW A GENE BEFORE 
   my $yeses               = 0;     # NUMBER OF $cegma, $ratt, $exonerate, $augustus THAT ARE 'yes'

   # CHECK IF THE GFF FILE WAS FROM EXONERATE:
   if ($exonerate ne 'yes' && $exonerate ne 'no')
   {
      $errormsg            = "ERROR: convert_gff_to_embl: exonerate should be yes/no (is $exonerate)\n";
      $errorcode           = 16; # ERRORCODE=16 (TESTED FOR) 
      return(\%GENES,\%EXONS,$errorcode,$errormsg);
   }

   # CHECK IF THE GFF FILE WAS FROM CEGMA:     
   if ($cegma ne 'yes' && $cegma ne 'no')
   {
      $errormsg            = "ERROR: convert_gff_to_embl: cegma should be yes/no (is $cegma)\n";
      $errorcode           = 17; # ERRORCODE=17 (TESTED FOR) 
      return(\%GENES,\%EXONS,$errorcode,$errormsg);
   }
   
   # CHECK IF THE GFF FILE WAS FROM RATT:
   if ($ratt ne 'yes' && $ratt ne 'no')
   {
      $errormsg            = "ERROR: convert_gff_to_embl: ratt should be yes/no (is $ratt)\n";
      $errorcode           = 21; # ERRORCODE=21 (TESTED FOR)
      return(\%GENES,\%EXONS,$errorcode,$errormsg);
   }

   # CHECK IF THE GFF FILE WAS FROM AUGUSTUS:
   if ($augustus ne 'yes' && $augustus ne 'no')
   {
      $errormsg            = "ERROR: convert_gff_to_embl: augustus should be yes/no (is $augustus)\n";
      $errorcode           = 50; # ERRORCODE=50 
      return(\%GENES,\%EXONS,$errorcode,$errormsg);
   }

   # CHECK THAT ONLY ONE OF $cegma, $ratt, $exonerate, $augustus IS 'yes':
   if ($cegma eq 'yes')     { $yeses++;}
   if ($ratt eq 'yes')      { $yeses++;}
   if ($exonerate eq 'yes') { $yeses++;}
   if ($augustus eq 'yes')  { $yeses++;}  
   if ($yeses > 1)
   {
      $errormsg             = "ERROR: convert_gff_to_embl: only one of cegma/exonerate/ratt/augustus can be set to yes\n";
      $errorcode            = 51; # ERRORCODE=51 (TESTED FOR)
      return(\%GENES,\%EXONS,$errorcode,$errormsg);
   }

   # READ IN THE GFF FILE:
   # NOTE: THE GFF FILE HAS A 'gene' FEATURE FIRST FOR THE GENE, AND THIS
   # IS FOLLOWED BY CDS & mRNA FEATURES FOR THE GENE.
   open(GFF,"$gff") || die "ERROR: read_gene_positions: cannot open $gff\n";
   while(<GFF>)
   {
      $line                = $_;   
      chomp $line;
      if (substr($line,0,1) ne '#') # IF IT'S NOT A COMMENT LINE
      {
         @temp             = split(/\t+/,$line);
         $scaffold         = $temp[0];
         $feature          = $temp[2];
         $start            = $temp[3];
         $end              = $temp[4];
         $score            = $temp[5];
         $strand           = $temp[6];
         if    ($feature eq 'gene') {
            # FIND THE GENE NAME:
            $name          = $temp[8];
            if ($exonerate eq 'yes') {
               @temp       = split(/sequence\s/,$name);  # eg. gene_id 0 ; sequence PTC00487_1 ; gene_orientation +
               $name       = $temp[1]; 
               @temp       = split(/\s+/,$name);
               $name       = $temp[0];                   # eg. PTC00487_1
            }
            else {
               if ($name =~ /ID=/) {
                  @temp    = split(/ID=/, $name); # eg. ID=ratti_train414
                  $name    = $temp[1];            # eg. ratti_train414 
               }
            }
            if ($augustus ne 'yes') {
               if (!($GENES{$scaffold})) {
                  $GENES{$scaffold} = $name;
               }
               else {
                  $GENES{$scaffold} = $GENES{$scaffold}.",".$name;
               }
            }
         }
         elsif ($feature eq 'CDS' || $feature eq 'cds' || $feature eq 'First' || $feature eq 'Internal' || $feature eq 'Terminal' || $feature eq 'Single')
         # exonerate USES 'cds', cegma 'First'/'Terminal'/'Internal'/'Single'
         {
            if ($cegma eq 'yes') {
               $name       = $temp[8];
               if (!($SEEN{$name}))
               {
                  if (!($GENES{$scaffold}))
                  {
                     $GENES{$scaffold} = $name;
                  }
                  else
                  {
                    $GENES{$scaffold} = $GENES{$scaffold}.",".$name;
                  }
                  $SEEN{$name} = 1;
               }
            }
            elsif ($augustus eq 'yes') {
               # eg. transcript_id "g1.t1"; gene_id "g1";
               $name       = $temp[8];
               @temp       = split(/transcript_id/,$name);
               $name       = $temp[1]; #   "g1.t1"; gene_id "g1";
               @temp       = split(/\"/,$name);
               $name       = $temp[1]; # g1.t1
               # PUT EACH TRANSCRIPT OF THE GENE AS A SEPARATE GENE IN THE EMBL FILE, SO THAT THEY WILL BE SEEN AS SEPARATE GENES BY ARTEMIS:
               if (!($SEEN{$name})) {
                  if (!($GENES{$scaffold})) {
                     $GENES{$scaffold} = $name;
                  }
                  else {
                    $GENES{$scaffold} = $GENES{$scaffold}.",".$name;
                  }
                  $SEEN{$name} = 1;
               }
            }
            elsif ($ratt eq 'yes') {
               $name       = $temp[8];
               @temp       = split(/\;Parent=/,$name); # eg. ratti_train155_1_0;Parent=ratti_train155_1_mRNA
               $name       = $temp[1];                 # eg. ratti_train155_1_mRNA
               @temp       = split(/\_mRNA/,$name);     
               $name       = $temp[0];                 # eg. ratti_train155_1
            }
            if ($name eq 'none') {
               $errormsg   = "ERROR: read_gene_positions: name $name line $line\n";
               $errorcode  = 14; # ERRORCODE=14 (TESTED FOR)
               return(\%GENES,\%EXONS,$errorcode,$errormsg);
            }
            $exon       = $start."=".$end."=".$strand;
            if (!($EXONS{$name})) {
               $EXONS{$name} = $exon;
            }
            else {
               $EXONS{$name} = $EXONS{$name}.",".$exon;
            }
         }
      }
   }
   close(GFF);

   return(\%GENES,\%EXONS,$errorcode,$errormsg);
}

#------------------------------------------------------------------#

# TEST &read_assembly

sub test_read_assembly {
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
   if ($SEQ->{'seq1'} ne 'AAAAA' || $SEQ->{'seq2'} ne 'AAAAATTTTT' || $SEQ->{'seq3'} ne 'AAAAATTTTT' || defined($SEQ->{'seq4'}) || $SEQ->{'seq5'} ne 'AAAAA' || $errorcode != 0) { 
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

# TEST &convert_gff_to_embl

sub test_convert_gff_to_embl
{
   my $outputdir           = $_[0]; # DIRECTORY TO PUT OUTPUT FILES INTO
   my $errorcode;                   # RETURNED AS 0 FROM A FUNCTION IF THERE IS NO ERROR
   my $errormsg;                    # RETURNED AS 'none' FROM A FUNCTION IF THERE IS NO ERROR  
   my %GENES               = ();    # HASH TABLE OF GENES IN SCAFFOLDS
   my %EXONS               = ();    # HASH TABLE OF EXONS IN GENES
   my %LEN                 = ();    # HASH TABLE WITH THE LENGTHS OF SCAFFOLDS 
   my $expected_output;             # EXPECTED OUTPUT FILE
   my $differences;                 # DIFFERENCES BETWEEN AN OUTPUT FILE AND EXPECTED OUTPUT FILE
   my $length_differences;          # LENGTH OF $differences
   my $line;                        #  
   my $output;                      # OUTPUT FILE
   my %SEQ                 = ();    # HASH TABLE TO STORE SCAFFOLD SEQUENCES 

   $GENES{"scaffold1"}     = "gene1,gene2,gene4";
   $SEQ{"scaffold1"}       = "AAAAAAAAAAGGGGGGGGGGCCCCCCCCCCTTTTTTTTTTAAAAAAAAAAGGGGGGGGGGCCCCCCCCCCTTTTTTTTTT";
   $EXONS{"gene1"}         = "100=200=+,250=350=+,450=490=+";
   $EXONS{"gene2"}         = "5000=5500=+";
   $EXONS{"gene3"}         = "10=20=-,25=30=-";
   $EXONS{"gene4"}         = "100=200=-,250=350=-,450=490=-"; 
   $LEN{"scaffold1"}       = 80;
   ($errorcode,$errormsg)  = &convert_gff_to_embl($outputdir,\%GENES,\%EXONS,\%LEN,\%SEQ);
   if ($errorcode != 0) { print STDERR "ERROR: test_convert_gff_to_embl: failed test1 (errorcode $errorcode errormsg $errormsg)\n"; exit;}
   ($expected_output,$errorcode,$errormsg) = &make_filename($outputdir);
   if ($errorcode != 0) { ($errorcode,$errormsg) = &print_error($errormsg,$errorcode,0); } 
   open(EXPECTED,">$expected_output") || die "ERROR: test_convert_gff_to_embl: cannot open expected_output $expected_output\n";
   print EXPECTED "ID   scaffold1\; SV 1\; linear\; genomic DNA\; HTG\; INV\; 80 BP.\n";
   print EXPECTED "FH   Key             Location/Qualifiers\n";
   print EXPECTED "FH\n";
   print EXPECTED "FT   CDS              join(100..200,250..350,450..490)\n";
   print EXPECTED "FT                   \/systematic_id=\"gene_id 0 \; sequence gene1 \; gene_orientation +\"\n";
   print EXPECTED "FT                   \/colour=7\n";
   print EXPECTED "FT   CDS              join(5000..5500)\n";
   print EXPECTED "FT                   \/systematic_id=\"gene_id 0 \; sequence gene2 \; gene_orientation +\"\n";
   print EXPECTED "FT                   \/colour=7\n";
   print EXPECTED "FT   CDS              complement(join(100..200,250..350,450..490))\n";
   print EXPECTED "FT                   \/systematic_id=\"gene_id 0 \; sequence gene4 \; gene_orientation +\"\n";
   print EXPECTED "FT                   \/colour=7\n";
   print EXPECTED "SQ   Sequence 80 BP; 20 A; 20 C; 20 G; 20 T; 0 other;\n";
   print EXPECTED "     AAAAAAAAAA GGGGGGGGGG CCCCCCCCCC TTTTTTTTTT AAAAAAAAAA GGGGGGGGGG 60\n";
   print EXPECTED "     CCCCCCCCCC TTTTTTTTTT                                             80\n";
   print EXPECTED "//\n";
   close(EXPECTED);
   $differences            = "";
   open(TEMP,"diff scaffold1.embl $expected_output |");
   while(<TEMP>)
   {
      $line                = $_;
      $differences         = $differences.$line;
   }
   close(TEMP);  
   $length_differences     = length($differences);
   if ($length_differences != 0) { print STDERR "ERROR: test_convert_gff_to_embl: failed test1 (output scaffold1.embl expected_output $expected_output)\n"; exit;}
   system "rm -f scaffold1.embl";
   system "rm -f $expected_output";
   ($expected_output,$errorcode,$errormsg) = &make_filename($outputdir);
   if ($errorcode != 0) { ($errorcode,$errormsg) = &print_error($errormsg,$errorcode,0); } 
   open(EXPECTED,">$expected_output") || die "ERROR: test_convert_gff_to_embl: cannot open expected_output $expected_output\n";
   print EXPECTED "ID   scaffold2\; SV 1\; linear\; genomic DNA\; HTG\; INV\; 100 BP.\n";
   print EXPECTED "FH   Key             Location/Qualifiers\n";
   print EXPECTED "FH\n";
   print EXPECTED "FT   CDS              complement(join(10..20))\n";
   print EXPECTED "FT                   \/systematic_id=\"gene_id 0 \; sequence gene3 \; gene_orientation +\"\n";
   print EXPECTED "FT                   \/colour=7\n";
   print EXPECTED "FT   CDS              complement(join(25..30))\n";
   print EXPECTED "FT                   \/systematic_id=\"gene_id 0 \; sequence gene3 \; gene_orientation +\"\n";
   print EXPECTED "FT                   \/colour=7\n";
   print EXPECTED "//\n";
   close(EXPECTED);
   $differences            = "";
   open(TEMP,"diff scaffold2.embl $expected_output |");
   while(<TEMP>)
   {
      $line                = $_;
      $differences         = $differences.$line;
   }
   close(TEMP);  
   $length_differences     = length($differences);
   if ($length_differences != 0) { print STDERR "ERROR: test_convert_gff_to_embl: failed test1 (output scaffold2.embl expected_output $expected_output)\n"; exit;}
   system "rm -f scaffold2.embl";
   system "rm -f $expected_output";

   # TEST ERRORCODE=7 (DO NOT KNOW THE LENGTH OF A SCAFFOLD):
   $GENES{"scaffold1"}     = "gene1,gene2";
   $GENES{"scaffold2"}     = "gene3"; 
   $EXONS{"gene1"}         = "100=200=+,250=350=+,450=490=+";
   $EXONS{"gene2"}         = "5000=5500=+";
   $EXONS{"gene3"}         = "10=20=-,25=30=-";
   %LEN                    = ();
   ($errorcode,$errormsg)  = &convert_gff_to_embl($outputdir,\%GENES,\%EXONS,\%LEN,\%SEQ);
   if ($errorcode != 7) { print STDERR "ERROR: test_convert_gff_to_embl: failed test2 (errorcode $errorcode errormsg $errormsg)\n"; exit;}
   system "rm -f scaffold1.embl";
   system "rm -f scaffold2.embl";

   # TEST ERRORCODE=15 (SCAFFOLD FILE EXISTS ALREADY):
   $GENES{"scaffold1"}     = "gene1,gene2";
   $GENES{"scaffold2"}     = "gene3"; 
   $EXONS{"gene1"}         = "100=200=+,250=350=+,450=490=+";
   $EXONS{"gene2"}         = "5000=5500=+";
   $EXONS{"gene3"}         = "10=20=-,25=30=-";
   $LEN{"scaffold1"}       = 6000;
   $LEN{"scaffold2"}       = 100;
   $output                 = $outputdir."/scaffold1.embl";
   open(OUTPUT,">$output") || die "ERROR: test_convert_gff_to_embl: cannot open output $output\n";
   close(OUTPUT);
   ($errorcode,$errormsg)  = &convert_gff_to_embl($outputdir,\%GENES,\%EXONS,\%LEN,\%SEQ);
   if ($errorcode != 4) { print STDERR "ERROR: test_convert_gff_to_embl: failed test3 (errorcode $errorcode errormsg $errormsg)\n"; exit;}
   system "rm -f $output"; 
   system "rm -f scaffold1.embl";
   system "rm -f scaffold2.embl";

   # TEST ERRORCODE=3 (DO NOT KNOW EXONS IN A GENE):
   $GENES{"scaffold1"}     = "gene1,gene2";
   $GENES{"scaffold2"}     = "gene3"; 
   %EXONS                  = ();
   $LEN{"scaffold1"}       = 6000;
   $LEN{"scaffold2"}       = 100;
   ($errorcode,$errormsg)  = &convert_gff_to_embl($outputdir,\%GENES,\%EXONS,\%LEN,\%SEQ);
   if ($errorcode != 3) { print STDERR "ERROR: test_convert_gff_to_embl: failed test4 (errorcode $errorcode errormsg $errormsg)\n"; exit;}
   system "rm -f scaffold1.embl";
   system "rm -f scaffold2.embl";

   # TEST ERRORCODE=9 (DO NOT KNOW SEQUENCE OF A SCAFFOLD):
   $GENES{"scaffold1"}     = "gene1,gene2";
   %SEQ                    = ();
   $EXONS{"gene1"}         = "100=200=+,250=350=+,450=490=+";
   $EXONS{"gene2"}         = "5000=5500=+";
   $EXONS{"gene3"}         = "10=20=-,25=30=-";
   $LEN{"scaffold1"}       = 80;
   ($errorcode,$errormsg)  = &convert_gff_to_embl($outputdir,\%GENES,\%EXONS,\%LEN,\%SEQ);
   if ($errorcode != 9) { print STDERR "ERROR: test_convert_gff_to_embl: failed test5 (errorcode $errorcode errormsg $errormsg)\n"; exit;}
   system "rm -f scaffold1.embl";
   system "rm -f scaffold2.embl";
}

#------------------------------------------------------------------#

# READ IN THE GFF FILE, AND CONVERT TO EMBL FORMAT:

sub convert_gff_to_embl
{
   my $outputdir           = $_[0]; # DIRECTORY TO WRITE OUTPUT FILES TO
   my $GENES               = $_[1]; # HASH TABLE TO RECORD GENES ON SCAFFOLDS
   my $EXONS               = $_[2]; # HASH TABLE TO RECORD EXONS IN GENES
   my $LEN                 = $_[3]; # HASH TABLE TO RECORD THE LENGTH OF SCAFFOLDS
   my $SEQ                 = $_[4]; # HASH TABLE WITH SCAFFOLD SEQUENCES 
   my $errorcode           = 0;     # RETURNED AS 0 IF THERE IS NO ERROR
   my $errormsg            = "none";# RETURNED AS 'none' IF THERE IS NO ERROR
   my $scaffold;                    # SCAFFOLD NAME
   my $scaffold_len;                # LENGTH OF THE SCAFFOLD
   my $genes;                       # GENES ON A SCAFFOLD
   my @genes;                       # GENES ON A SCAFFOLD
   my $i;                           # 
   my $gene;                        # A GENE
   my $exons;                       # EXONS IN A GENE
   my @exons;                       # EXONS IN A GENE
   my $exon;                        # AN EXON
   my $j;                           # 
   my $exon_start;                  # EXON START
   my $exon_end;                    # EXON END
   my $exon_strand;                 # EXON STRAND
   my @temp;                        #
   my $output;                      # OUTPUT EMBL FILE FOR A SCAFFOLD 
   my $seq;                         # SEQUENCE FOR THE SCAFFOLD
   my $num_as;                      # NUMBER OF As IN THE SCAFFOLD
   my $num_cs;                      # NUMBER OF Cs IN THE SCAFFOLD
   my $num_gs;                      # NUMBER OF Gs IN THE SCAFFOLD
   my $num_ts;                      # NUMBER OF Ts IN THE SCAFFOLD 
   my $num_other;                   # NUMBER OF OTHER BASES IN THE SCAFFOLD
   my $offset;                      # COUNTER USED IN PRINTING OUT THE SEQUENCE
   my $a_line;                      # A LINE OF SEQUENCE TO PRINT OUT
   my $a_chunk;                     # A CHUNK OF SEQUENCE TO PRINT OUT 
   my $new_exons;                   # LINE FOR NEW EXONS TO PRINT TO EMBL FILE 
   my $leftover;                    # LEFTOVER SPACE AFTER PRINTOUT OUT SEQUENCE 
   my $printed;                     # AMOUNT OF SEQUENCE PRINTED ALREADY 

   # MAKE AN EMBL FILE FOR EACH SCAFFOLD THAT HAS GENES ON IT: 
   foreach $scaffold (keys %{$GENES})
   {
      # FIND THE LENGTH OF THE SCAFFOLD:
      if (!($LEN->{$scaffold}))
      {
         $errormsg         = "ERROR: convert_gff_to_embl: do not know length of scaffold $scaffold\n";
         $errorcode        = 7; # ERRORCODE=7 (TESTED FOR)
         return($errorcode,$errormsg);
      }
      $scaffold_len        = $LEN->{$scaffold};

      # MAKE AN OUTPUT FILE FOR THIS SCAFFOLD:
      $output              = $scaffold.".embl";
      $output              = $outputdir."/".$output;
      if (-e $output)
      {
         $errormsg         = "ERROR: convert_gff_to_embl: output file $output exists already\n";
         $errorcode        = 15; # ERRORCODE=15 (TESTED FOR)
         return($errorcode,$errormsg);
      }
      print STDERR "Opening output file $output...\n";
      open(OUTPUT,">$output") || die "ERROR: convert_gff_to_embl: cannot open $output\n";
      print OUTPUT "ID   $scaffold\; SV 1\; linear\; genomic DNA\; HTG\; INV\; $scaffold_len BP.\n";
      print OUTPUT "FH   Key             Location/Qualifiers\n";
      print OUTPUT "FH\n";
      # NOTE: FROM http://www.ebi.ac.uk/help/formats.html#EMBL, THE COLUMNS ARE:
      # ACCESSION NUMBER
      # SV = SEQUENCE VERSION NUMBER
      # DATA CLASS, eg. 'linear'
      # MOLECULE, eg. 'mRNA'
      # DATA CLASS, eg. 'STD'
      # TAXONOMIC DIVISION, eg. 'INV'  
      # SEQUENCE LENGTH. 

      # FIND THE GENES IN THE SCAFFOLD:
      if (!($GENES->{$scaffold}))
      {
         $errormsg         = "ERROR: convert_gff_to_embl: do not know genes for scaffold $scaffold\n";
         $errorcode        = 8; # ERRORCODE=8 (SHOULDN'T GET HERE, SO CAN'T TEST FOR THIS)
         return($errorcode,$errormsg);
      }
      $genes               = $GENES->{$scaffold};
      @genes               = split(/\,/,$genes);
      for ($i = 0; $i <= $#genes; $i++)
      {
         $gene             = $genes[$i];
         # GET THE EXONS IN THIS GENE:
         if (!($EXONS->{$gene}))
         {
            $errormsg      = "ERROR: convert_gff_to_embl: do not know exons in gene $gene\n";
            $errorcode     = 3; # ERRORCODE=3 (TESTED FOR)
            return($errorcode,$errormsg);
         }
         $exons            = $EXONS->{$gene};
         @exons            = split(/\,/,$exons);
         $new_exons        = "";
         for ($j = 0; $j <= $#exons; $j++)
         {
            $exon          = $exons[$j];
            @temp          = split(/=/,$exon);
            $exon_start    = $temp[0];
            $exon_end      = $temp[1];
            $exon_strand   = $temp[2];
            $new_exons     = $new_exons.",".$exon_start."..".$exon_end;
         }
         $new_exons        = substr($new_exons,1,length($new_exons)-1);
         if ($#exons > 0) { $new_exons = "join(".$new_exons.")"; }  
         if ($exon_strand eq '-') { $new_exons = "complement(".$new_exons.")";}
         print OUTPUT "FT   CDS             $new_exons\n";
         print OUTPUT "FT                   \/systematic_id=\"gene_id 0 \; sequence $gene \; gene_orientation +\"\n";
         print OUTPUT "FT                   \/colour=7\n";
      }
      # PRINT OUT THE FASTA SEQUENCE:
      if (!($SEQ->{$scaffold}))
      {
         $errormsg         = "ERROR: convert_gff_to_embl: do not know sequence for $scaffold\n";
         $errorcode        = 9; # ERRORCODE=9 (TESTED FOR)
         return($errorcode,$errormsg);
      }
      $seq                 = $SEQ->{$scaffold};
      ($num_as,$num_cs,$num_gs,$num_ts,$num_other,$errorcode,$errormsg) = &count_bases($seq);
      if ($errorcode != 0) { ($errorcode,$errormsg) = &print_error($errormsg,$errorcode,0); }
      print OUTPUT "SQ   Sequence $scaffold_len BP\; $num_as A\; $num_cs C\; $num_gs G\; $num_ts T\; $num_other other\;\n"; 
      # PRINT OUT THE SEQUENCE, eg.:
      # SQ   Sequence 5957 BP; 1968 A; 916 C; 902 G; 2171 T; 0 other;
      #      agcttttggt atctcgtaat gtagaacaat aatgggttga ccgtcaaatc cttttcatta        60
      #      aaagagtgga ttaatgccca gccaacacca tccaatttga ttgggaataa tctgtgttac       120
      #      aatacttttt gatcccaggc tggtaaaaaa tgtaaacttt tagcccgaaa gaatagaaac       180
      #      agatgccagg ccaaaaaccc aaaaatagag ctatgacgct atcattttaa caagttagat       240
      #      aaattctttc atagaactta acgtttcatc ttccatacat aaataaaacg gtagataggg       300
      $offset              = 0;
      while ($offset < $scaffold_len)
      {
         print OUTPUT "     "; 
         $a_line           = substr($seq,$offset,60); 
         my $length_a_line    = length($a_line); 
         my $offset2          = 0; 
         while ($offset2 < $length_a_line)
         {
            $a_chunk       = substr($a_line,$offset2,10);
            print OUTPUT "$a_chunk ";
            $offset2       = $offset2 + 10; 
         }
         if ($offset + 60 < $scaffold_len)
         {
            $printed       = $offset + 60;
            print OUTPUT "   $printed\n";
         }
         else
         {
            $leftover      = 60 - length($a_line);
            print STDERR "leftover $leftover\n";
            for ($i = 1; $i <= $leftover; $i++)
            {
               print OUTPUT " ";
            }
            print OUTPUT "    $scaffold_len\n";
         }
         $offset           = $offset + 60;
      }
      print OUTPUT "//\n"; 
      close(OUTPUT);
   }   
 
   return($errorcode,$errormsg); 
}

#------------------------------------------------------------------#

# TEST &count_bases

sub test_count_bases
{
   my $seq;                         # SEQUENCE
   my $num_as;                      # NUMBER OF As IN THE SEQUENCE
   my $num_cs;                      # NUMBER OF Cs IN THE SEQUENCE
   my $num_gs;                      # NUMBER OF Gs IN THE SEQUENCE
   my $num_ts;                      # NUMBER OF Ts IN THE SEQUENCE
   my $num_other;                   # NUMBER OF OTHER LETTERS IN THE SEQUENCE
   my $errorcode;                   # RETURNED AS 0 FROM A FUNCTION IF THERE IS NO ERROR
   my $errormsg;                    # RETURNED AS 'none' FROM A FUNCTION IF THERE IS NO ERROR
 
   $seq                    = "AAACCGGTYYAACCGT";
   ($num_as,$num_cs,$num_gs,$num_ts,$num_other,$errorcode,$errormsg) = &count_bases($seq); 
   if ($errorcode != 0 || $num_as != 5 || $num_cs != 4 || $num_gs != 3 || $num_ts != 2 || $num_other != 2) { print STDERR "ERROR: test_count_bases: failed test1\n"; exit;}
   
}

#------------------------------------------------------------------#

# SUBROUTINE SYNOPSIS: count_bases(): COUNTS THE NUMBER OF As, Cs, Gs AND Ts IN A SEQUENCE

sub count_bases
{
   my $seq                 = $_[0]; # SEQUENCE
   my $errorcode           = 0;     # RETURNED AS 0 IF THERE IS NO ERROR
   my $errormsg            = "none";# RETURNED AS 'none' IF THERE IS NO ERROR
   my $num_as              = 0;     # NUMBER OF As IN THE SEQUENCE
   my $num_cs              = 0;     # NUMBER OF Cs IN THE SEQUENCE
   my $num_gs              = 0;     # NUMBER OF Gs IN THE SEQUENCE
   my $num_ts              = 0;     # NUMBER OF Ts IN THE SEQUENCE
   my $num_other           = 0;     # NUMBER OF OTHER LETTERS IN THE SEQUENCE
   my $length;                      # LENGTH OF THE SEQUENC
   my $i;                           # 
   my $base;                        # A BASE IN THE SEQUENCE

   $length                 = length($seq);
   $seq                    =~ tr/[a-z]/[A-Z]/;
   for ($i = 1; $i <= $length; $i++)
   {
      $base                = substr($seq,$i-1,1);
      if    ($base eq 'A') { $num_as++;   }
      elsif ($base eq 'C') { $num_cs++;   }
      elsif ($base eq 'G') { $num_gs++;   }
      elsif ($base eq 'T') { $num_ts++;   }
      else                 { $num_other++;} 
   }

   return($num_as,$num_cs,$num_gs,$num_ts,$num_other,$errorcode,$errormsg); 
}

#------------------------------------------------------------------#

# READ IN THE LENGTHS OF CONTIGS IN THE ASSEMBLY:
# SUBROUTINE SYNOPSIS: read_scaffold_length(): store lengths of sequences from a fasta file in a hash

sub read_scaffold_lengths
{
   my $assembly            = $_[0]; # FILE CONTAINING THE ASSEMBLY
   my %SCAFFOLDLEN         = ();    # HASH TABLE FOR STORING LENGTHS OF SCAFFOLDS
   my $line;                        #  
   my @temp;                        # 
   my $length;                      # LENGTH OF A SCAFFOLD
   my $name                = "none";# NAME OF A SCAFFOLD
   my $seq;                         # SEQUENCE OF A SCAFFOLD
   my $errorcode           = 0;     # THIS IS RETURNED AS 0 IF NO ERROR OCCURRED. 
   my $errormsg            = "none";# THIS IS RETURNED AS 'none' IF NO ERROR OCCURRED.
   my %SEEN                = ();    # HASH TABLE TO RECORD WHICH SCAFFOLD NAMES WE HAVE SEEN.
 
   open(ASSEMBLY,"$assembly") || die "ERROR: read_scaffold_lengths: cannot open $assembly\n"; 
   while(<ASSEMBLY>)
   {
      $line                = $_;
      chomp $line;
      if (substr($line,0,1) eq '>')
      {
         @temp             = split(/\s+/,$line);
         $name             = $temp[0];
         $name             = substr($name,1,length($name)-1);
         if ($SEEN{$name})
         {
            $errormsg      = "ERROR: read_scaffold_lengths: seen scaffold $name already\n"; 
            $errorcode     = 1; # ERRORCODE=1 (TESTED FOR)
            return(\%SCAFFOLDLEN, $errorcode, $errormsg);
         }
         $SEEN{$name}      = 1;
      } 
      else
      {
         $seq              = $line;
         # REMOVAL SPACES FROM THE LINE, EITHER INTERNAL SPACES OR SPACES AT EITHER END:
         $seq              = $line;
         $seq              =~ s/\s+//g;
         # STORE THE LENGTH OF THE SEQUENCE, OR UPDATE THE STORED LENGTH:
         if ($seq eq '') { $length = 0;           }
         else            { $length = length($seq);}
         if ($name eq 'none') 
         { 
            $errormsg      = "ERROR: read_scaffold_lengths: name is $name\n"; 
            $errorcode     = 2; # ERRORCODE=2 (TESTED FOR)
            return(\%SCAFFOLDLEN, $errorcode, $errormsg);
         }
         if (!($SCAFFOLDLEN{$name})){ $SCAFFOLDLEN{$name} = $length;}
         else {$SCAFFOLDLEN{$name} =  $SCAFFOLDLEN{$name} + $length;}   
      }
   }
   close(ASSEMBLY);

   return(\%SCAFFOLDLEN,$errorcode,$errormsg);
}

#------------------------------------------------------------------#

# TEST &read_scaffold_lengths

sub test_read_scaffold_lengths
{
   my $outputdir           = $_[0];  # DIRECTORY TO WRITE OUTPUT IN.
   my $assembly;                     # FILE CONTAINING ASSEMBLY.      
   my $LEN;                          # HASH TABLE CONTAINING LENGTHS OF SEQUENCES.
   my $errorcode;                    # RETURNED AS 0 FROM A FUNCTION IF THERE WAS NO ERROR. 
   my $errormsg;                     # RETURNED AS 'none' FROM A FUNCTION IF THERE WAS NO ERROR.

   ($assembly,$errorcode,$errormsg) = &make_filename($outputdir);
   if ($errorcode != 0) { ($errorcode,$errormsg) = &print_error($errormsg,$errorcode,0); } 
   open(ASSEMBLY,">$assembly") || die "ERROR: test_read_scaffold_lengths: cannot open $assembly\n";
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
   ($LEN,$errorcode,$errormsg) = &read_scaffold_lengths($assembly);
   if ($LEN->{'seq1'} != 5 || $LEN->{'seq2'} != 10 || $LEN->{'seq3'} != 10 || defined($LEN->{'seq4'}) || $LEN->{'seq5'} != 5 || $errorcode != 0) { print STDERR "ERROR: test_read_scaffold_lengths: failed test1\n"; exit;}
   system "rm -f $assembly";

   ($assembly,$errorcode,$errormsg) = &make_filename($outputdir);
   if ($errorcode != 0) { ($errorcode,$errormsg) = &print_error($errormsg,$errorcode,0); } 
   open(ASSEMBLY,">$assembly") || die "ERROR: test_read_scaffold_lengths: cannot open $assembly\n";
   print ASSEMBLY ">seq1\n";
   print ASSEMBLY "AAAAA\n";
   print ASSEMBLY ">seq2\n";
   print ASSEMBLY "AAAAATTTTT\n";
   print ASSEMBLY ">seq1\n";
   print ASSEMBLY "AAAAA\n";
   close(ASSEMBLY);
   ($LEN,$errorcode,$errormsg) = &read_scaffold_lengths($assembly);
   if ($errorcode != 1) { print STDERR "ERROR: test_read_scaffold_lengths: failed test2\n"; exit;}
   system "rm -f $assembly";

   ($assembly,$errorcode,$errormsg) = &make_filename($outputdir);
   if ($errorcode != 0) { ($errorcode,$errormsg) = &print_error($errormsg,$errorcode,0); } 
   open(ASSEMBLY,">$assembly") || die "ERROR: test_read_scaffold_lengths: cannot open $assembly\n";
   print ASSEMBLY "AAAAA\n";
   print ASSEMBLY ">seq2\n";
   print ASSEMBLY "AAAAATTTTT\n";
   print ASSEMBLY ">seq1\n";
   print ASSEMBLY "AAAAA\n";
   close(ASSEMBLY);
   ($LEN,$errorcode,$errormsg) = &read_scaffold_lengths($assembly);
   if ($errorcode != 2) { print STDERR "ERROR: test_read_scaffold_lengths: failed test3\n"; exit;}
   system "rm -f $assembly";
   
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


