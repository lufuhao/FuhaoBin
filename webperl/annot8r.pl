#!/usr/bin/perl -w
# annot8r.pl version 1.1.2 luf
# annotates GO terms, EC numbers and KEGG entries based on blast searches 
# Perl Modules - DBD:Pg and BioPerl
# DBD::Pg is part of the DBI is available from cpan - http://www.cpan.org/
# Bioperl is available from their site - http://www.bioperl.org/
#
# Updated v1.1.1 06/10/2008 Mark Blaxter
# Institute Evolutionary Biology, University of Edinburgh
# Updated to cope with changes in the delivery format of sequences from EBI
# 
# Core program v1.1 25/11/2007 by Ralf Schmid
# Department of Biochemistry, University of Leicester
# Copyright (C) 2006 Ralf Schmid

# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.

use strict;
use File::stat;
use Term::ANSIColor;
use Bio::Tools::Run::StandAloneBlast;
use Bio::SearchIO;
use Bio::SeqIO;
use DBI;
use DBD::Pg;  
use Time::localtime;

# How to get rid of the longish cut-off sermon? Change the value for $quiet from 0 to 1   
my $quiet=0; 
# my $quiet= 1;  


# the following defines the type of evidence for GO-entries which will be used for the uniprot subset
# change if necessary, remove"#" to make active
# The evidence can be thought of in a loose hierarchy of reliability (according to the GOA curators):
#  TAS / IDA >  IMP / IGI / IPI > ISS / IEP > NAS > IEA



my %go_hash;
$go_hash{IC} = 1; 	# IC: Inferred by Curator
$go_hash{TAS} = 1;	# TAS: Traceable Author Statement
$go_hash{IDA} = 1;	# IDA: Inferred from Direct Assay
$go_hash{IMP} = 1;	# IMP: Inferred from Mutant Phenotype
$go_hash{IGI} = 1;	# IGI: Inferred from Genetic Interaction
$go_hash{IPI} = 1;	# IPI: Inferred from Physical Interaction
$go_hash{ISS} = 1;	# ISS: Inferred from Sequence or Structural Similarity
$go_hash{IEP} = 1;	# IEP: Inferred from Expression Pattern
$go_hash{NAS} = 1;	# NAS: Non-traceable Author Statement
#$go_hash{IEA} = 1;	# IEA: Inferred from Electronic Annotation
#$go_hash{ND} = 1;	# ND: No biological Data available
$go_hash{RCA} = 1;	# RCA: inferred from Reviewed Computational Analysis
#$go_hash{NR} = 1;	# NR: Not Recorded 

my $version_number = "1.1  ";
my @PATH=split(":","$ENV{'PATH'}");   


#############################################################################
#############################################################################

# wget & ftp commands, may change if providers move their resources - e.g. below
my $getdata1="wget --passive-ftp --output-document=uniprot_sprot.fasta.gz ftp://ftp.ebi.ac.uk/pub/databases/uniprot/knowledgebase/uniprot_sprot.fasta.gz";
my $getdata2="wget --passive-ftp --output-document=uniprot_trembl.fasta.gz ftp://ftp.ebi.ac.uk/pub/databases/uniprot/knowledgebase/uniprot_trembl.fasta.gz";


#my $getdata3="wget --passive-ftp --output-document=gene_association.goa_uniprot.gz ftp://ftp.ebi.ac.uk/pub/databases/GO/goa/UNIPROT/gene_association.goa_uniprot.gz";
### NEW UPDATE luf 20171019
### goa_uniprot_all.gaf.gz - renamed gene_association.goa_uniprot.gz; identical content
my $getdata3="wget --passive-ftp --output-document=gene_association.goa_uniprot.gz ftp://ftp.ebi.ac.uk/pub/databases/GO/goa/UNIPROT/goa_uniprot_all.gaf.gz";


my $getdata4="wget --passive-ftp --output-document=goaslim.map ftp://ftp.ebi.ac.uk/pub/databases/GO/goa/goslim/goaslim.map";
my $getdata5="wget --passive-ftp --output-document=GO.terms_and_ids http://www.geneontology.org/doc/GO.terms_and_ids";

# my $getdata6="wget --passive-ftp --output-document=enzyme.dat ftp://www.expasy.org/databases/enzyme/release_with_updates/enzyme.dat";
# apparently changed in 06/2006
#my $getdata6="wget --passive-ftp --output-document=enzyme.dat ftp://www.expasy.org/databases/enzyme/enzyme.dat";
### NEW UPDATE luf 20171019
my $getdata6="wget --passive-ftp --output-document=enzyme.dat ftp://ftp.expasy.org/databases/enzyme/enzyme.dat";


# my $getdata7="wget --passive-ftp --output-document=genes_pathway.list ftp://ftp.genome.ad.jp/pub/kegg/linkdb/genes/genes_pathway.list";
# my $getdata7="wget --passive-ftp --output-document=ko ftp://ftp.genome.ad.jp/pub/kegg/tarfiles/ko";
# yet another change in 07/2007
my $getdata7="wget --passive-ftp --output-document=ko ftp://ftp.genome.ad.jp/pub/kegg/genes/ko";
my $getdata8="wget --passive-ftp --output-document=genes_uniprot.list ftp://ftp.genome.ad.jp/pub/kegg/linkdb/genes/genes_uniprot.list";
my $getdata9="wget --passive-ftp --output-document=genes_ko.list ftp://ftp.genome.ad.jp/pub/kegg/linkdb/genes/genes_ko.list";
my $getdata0="wget --passive-ftp --output-document=map_title.tab ftp://ftp.genome.ad.jp/pub/kegg/pathway/map_title.tab";

# comma separated outputfiles
my $goout = "GO.csv";
my $ecout = "EC.csv";
my $keggout = "KEGG.csv";

# databases
my $go_db = "a8r_gobase"; # GO-database created from flatfile
my $ec_db = "a8r_ecbase"; # EC-database created from flatfile
my $kegg_db = "a8r_keggbase"; # KEGG-database created from flatfile
my $pg_database; # PartiGene (or other database) to be decorated

# downloaded files, 
my $goasso_file = "gene_association.goa_uniprot"; # flatfile containing GO terms associated entries
my $gomap_file = "goaslim.map"; # GO_slim file
my $godef_file = "GO.terms_and_ids"; # GO_definitions file
my $ecdata_file = "enzyme.dat"; # EC-stuff
my $sprot_file = "uniprot_sprot.fasta"; # blast database file to be processed to $GO_blastfile
my $trembl_file = "uniprot_trembl.fasta";# blast database file to be processed to $GO_blastfile
my $kegguniprot_file = "genes_uniprot.list";
my $keggko_file = "genes_ko.list";
my $keggpathway_file = "ko";
my $keggmap_file = "map_title.tab";

#other files
my $GO_blastfile = "uniprot_GO.fsa"; # modified blastfile, sequences without GO-terms are removed
my $EC_blastfile = "uniprot_EC.fsa"; # modified blastfile, sequences without EC-terms are removed
my $KEGG_blastfile = "uniprot_KEGG.fsa"; # modified blastfile, sequences without KEGG-terms are removed
my $seqs2blast; # sequences to be blasted and annotated
my $blast_output; # blast output to be parsed

&options();

#######################################################################################################################
#### pre subs
#######################################################################################################################
sub options () {
#### option selection
   my $answer=0;
   while($answer!=6) {
      $answer=&title_page();
      if($answer==1)  { &download_options(); } # downloads
      if($answer==2)  { &annodb_options(); } # get dbs ready
      if($answer==3)  { &prepare_blast_options(); } # BLAST prep
      if($answer==4)  { &blast_options(); } # BLAST      
      if($answer==5)  { &anno_options(); } #do annotations
     
   }
   system("clear");
   exit();   #exit program
}
############################################################################################################################

   
###################################################################################
sub title_page() {
#### intro & sub-menu selection
    print_title();
    print "\n\tPlease select the option you want to run.\n\n";
    print "\t\t1. Download relevant files.\n";
    print "\t\t2. Extract data from files into databases.\n";
    print "\t\t3. Prepare your Blast searches.\n";
    print "\t\t4. Blast your sequences.\n";    
    print "\t\t5. Annotate your sequences.\n";
    print "\t\t6. Quit.\n";

    my $flag=0; my $answer;
    while($flag==0) {
        $answer=<>;
	if($answer=~/^[1|2|3|4|5|6]$/) { $flag=1; next; }
        else {print " You have entered $answer This is not an option. Please try again\n";}
    }
    return $answer;
}
####################################################################################


############################################################################################################################
sub print_title() {
##### displays title
    print colored("\n\n\n", "white bold", "on_black"); 
    print colored("\t###########################################################\n","white bold", "on_black");
    print colored("\t###                                                     ###\n","white bold", "on_black");
    print colored("\t###                    annot8r.pl                       ###\n","white bold", "on_black");
    print colored("\t###           a tool for sequence annotation            ###\n","white bold", "on_black");
    print colored("\t###            based on BLAST results Vs $version_number          ###\n","white bold", "on_black");
    print colored("\t###                                                     ###\n","white bold", "on_black");
    print colored("\t###     BANG 2007                                       ###\n","white bold", "on_black");
    print colored("\t###                                                     ###\n","white bold", "on_black");
    print colored("\t###     For news and upgrades and help:                 ###\n","white bold", "on_black");
    print colored("\t###     nematodes.bioinf\@ed.ac.uk                       ###\n","white bold", "on_black");
    print colored("\t###                                                     ###\n","white bold", "on_black");
    print colored("\t###########################################################\n\n\n","white bold", "on_black");
   
}
###########################################################################################################################

########################################################################################################################
sub download_options() {
   print colored("\n\t\t##### DOWNLOAD OPTIONS #####\n","white bold" , "on_black");
   my $answer=0;
   while($answer!=5) {
      $answer=&download_page();
      if($answer==1)  { &download("1"); } #go stuff
      if($answer==2)  { &download("2"); } #ec stuff
      if($answer==3)  { &download("3"); } #kegg stuff       
      if($answer==4)  { &download("4"); } #all of the above     
   }
   system("clear");
   options();   #exit program
}
######################################################################################################################

#######################################################################################
sub download_page() {
    print "\n\tThis option allows you to download the relevant files for\n";
    print "\tannot8r.pl.\n\n";
    print "\t\t1. Download files required for GO-annotation. \n";
    print "\t\t2. Download files required for EC-annotation. \n";
    print "\t\t3. Download files required for KEGG-annotation. \n";    
    print "\t\t4. Download everything.\n";
    print "\t\t5. Back to main menu.\n";

    my $flag=0; my $answer;
    while($flag==0) {
        $answer=<>;
	if($answer=~/^[1|2|3|4|5]$/) { $flag=1; next; }
        else {print " You have entered $answer This is not an option. Please try again\n";}
    }
    return $answer;
}
####################################################################################


########################################################################################################################
sub annodb_options() {
   print colored("\n\t\t##### DATABASE SETUP OPTIONS #####\n","white bold" , "on_black");
   my $answer=0;
   while($answer!=5) {
      $answer=&annodb_page();
      if($answer==1)  { &annodb("1"); } #go stuff
      if($answer==2)  { &annodb("2"); } #ec stuff
      if($answer==3)  { &annodb("3"); } #kegg stuff       
      if($answer==4)  { &annodb("4"); } #all of the above     
   }
   system("clear");
   options();   #exit program
}
######################################################################################################################


#######################################################################################
sub annodb_page() {
    print "\n\tThis option reads the relevant data from the files downloaded\n";
    print "\tin the previous option into three postgresql databases:\n";
    print "\ta8r_go, a8r_ec, a8r_kegg\n\n";
    print "\t\t1. Create or update a8r_go database.\n";
    print "\t\t2. Create or update a8r_ec database.\n";
    print "\t\t3. Create or update a8r_kegg database.\n";    
    print "\t\t4. Create or update all of them.\n";
    print "\t\t5. Back to main menu.\n";

    my $flag=0; my $answer;
    while($flag==0) {
        $answer=<>;
	if($answer=~/^[1|2|3|4|5]$/) { $flag=1; next; }
        else {print " You have entered $answer This is not an option. Please try again\n";}
    }
    return $answer;
}
####################################################################################


########################################################################################################################
sub prepare_blast_options() {
   print colored("\n\t\t##### PREPARE BLAST OPTIONS #####\n","white bold" , "on_black");
   my $answer=0;
   while($answer!=5) {
      $answer=&prepare_blast_page();
      if($answer==1)  { &prepare_blast("1"); } #go stuff
      if($answer==2)  { &prepare_blast("2"); } #ec stuff
      if($answer==3)  { &prepare_blast("3"); } #kegg stuff       
      if($answer==4)  { &prepare_blast("4"); } #all of the above     
   }
   system("clear");
   options();   #exit program
}
######################################################################################################################


#######################################################################################
sub prepare_blast_page() {
    print "\n\tThis option prepares specific BLAST databases required for\n";
    print "\tthe annotation of GO, EC and KEGG terms.\n\n";
    print "\t\t1. Create BLAST/GO database.\n";
    print "\t\t2. Create BLAST/EC database.\n";
    print "\t\t3. Create BLAST/KEGG database.\n";    
    print "\t\t4. Create all of them.\n";
    print "\t\t5. Back to main menu.\n";

    my $flag=0; my $answer;
    while($flag==0) {
        $answer=<>;
	if($answer=~/^[1|2|3|4|5]$/) { $flag=1; next; }
        else {print " You have entered $answer This is not an option. Please try again.\n";}
    }
    return $answer;
}
####################################################################################

########################################################################################################################
sub blast_options() {
   print colored("\n\t\t##### RUN BLASTS #####\n","white bold" , "on_black");
   my $answer=0;
   while($answer!=5) {
      $answer=&blast_page();
      if($answer==1)  { &blast("1"); } #go stuff
      if($answer==2)  { &blast("2"); } #ec stuff
      if($answer==3)  { &blast("3"); } #kegg stuff       
      if($answer==4)  { &blast("4"); } #all of the above     
   }
   system("clear");
   options();   #exit program
}
######################################################################################################################

#######################################################################################
sub blast_page() {
    print "\n\tThis option prepares specific BLAST databases required for\n";
    print "\tthe annotation of GO, EC and KEGG terms.\n\n";
    print "\t\t1. Run BLAST/GO database.\n";
    print "\t\t2. Run  BLAST/EC database.\n";
    print "\t\t3. Run BLAST/KEGG database.\n";    
    print "\t\t4. Run all of them.\n";
    print "\t\t5. Back to main menu.\n";

    my $flag=0; my $answer;
    while($flag==0) {
        $answer=<>;
	if($answer=~/^[1|2|3|4|5]$/) { $flag=1; next; }
        else {print " You have entered $answer This is not an option. Please try again.\n";}
    }
    return $answer;
}
####################################################################################

########################################################################################################################
sub anno_options() {
   print colored("\n\t\t##### Annotation #####\n","white bold" , "on_black");
   my $answer=0;
   while($answer!=5) {
      $answer=&anno_page();
      if($answer==1)  { &anno("1"); } #go stuff
      if($answer==2)  { &anno("2"); } #ec stuff
      if($answer==3)  { &anno("3"); } #kegg stuff       
      if($answer==4)  { &anno("4"); } #all of the above     
   }
   system("clear");
   options();   #exit program
}
######################################################################################################################

#######################################################################################
sub anno_page() {
    print "\n\tThis option extracts GO, EC and KEGG annotation\n";
    print "\tfrom BLAST results.\n\n";
    print "\t\t1. Annotate GO terms.\n";
    print "\t\t2. Annotate EC numbers.\n";
    print "\t\t3. Annotate KEGG pathway entries.\n";    
    print "\t\t4. Annotate all of above.\n";
    print "\t\t5. Back to main menu.\n";

    my $flag=0; my $answer;
    while($flag==0) {
        $answer=<>;
	if($answer=~/^[1|2|3|4|5]$/) { $flag=1; next; }
        else {print " You have entered $answer This is not an option. Please try again.\n";}
    }
    return $answer;
}
####################################################################################


#### from here on the main subs


##################################################################################
#################  1-download  ###########################
##################################################################################

sub download() {
  my $d_flag = $_[0]; my $flag = 1;
  unless (-e "downloads") { system ("mkdir downloads");}
  chdir ("downloads");

#### all stuff first
  print "Downloading relevant files required for all types of annotation:\n\n";
  $flag = 1;
  if (-e "$sprot_file") {$flag = &delete_file("$sprot_file"); }
  if (-e "uniprot_sprot.fasta.gz") {$flag = &delete_file("uniprot_sprot.fasta.gz")}
  if ($flag ==1) {  
    print "Please wait, downloading $sprot_file now ...\n";
    system("$getdata1");   
    print "Please wait, unpacking $sprot_file now ...\n";
    system ("gunzip uniprot_sprot.fasta.gz");    
    &filecheck($sprot_file);
    print "Done.\n";
  }
  $flag = 1;
  if (-e "$trembl_file") {$flag = &delete_file("$trembl_file"); }
  if (-e "uniprot_trembl.fasta.gz") {$flag = &delete_file("uniprot_trembl.fasta.gz")}
  if ($flag ==1) {  
    print "Please wait, downloading $trembl_file now ...\n";
    system("$getdata2");   
    print "Please wait, unpacking $trembl_file now ...\n";
    system ("gunzip uniprot_trembl.fasta.gz");    
    &filecheck($trembl_file);
    print "Done.\n";
  }
    
#### next GO stuff  
  if (($d_flag==1) || ($d_flag==4)) {
    print "Downloading relevant files for GO-annotation:\n\n";
    $flag = 1;
    if (-e "gene_association.goa_uniprot.gz") {$flag = &delete_file("gene_association.goa_uniprot.gz")}
    if  (-e "$goasso_file") {$flag = &delete_file("$goasso_file")}
    if ($flag ==1) {  
      print "Please wait, downloading $goasso_file now ...\n";
      system("$getdata3");   
      print "Please wait, unpacking $goasso_file now ...\n";
      system ("gunzip gene_association.goa_uniprot.gz");    
      &filecheck($goasso_file);
      print "Done.\n";
    }
    $flag = 1;
    if (-e "$gomap_file") {
      $flag = &delete_file("$gomap_file");
    }
    if ($flag ==1) {  
      print "Please wait, downloading $gomap_file now ...\n";
      system("$getdata4");   
      &filecheck($gomap_file);
      print "Done.\n";
    }
    $flag = 1;
    if (-e "$godef_file") {
      $flag = &delete_file("$godef_file");
    }
    if ($flag ==1) {  
      print "Please wait, downloading $godef_file now ...\n";
      system("$getdata5");   
      &filecheck($godef_file);
      print "Done.\n";  
    }
  }

#### now EC stuff
  if ($d_flag==2 || ($d_flag==4)) {
    print "Downloading relevant files for EC-annotation:\n\n";
    $flag=1;
    if (-e "$ecdata_file") {$flag = &delete_file("$ecdata_file"); }
    if ($flag ==1) {  
      print "Please wait, downloading $ecdata_file now ...\n";
      system("$getdata6");   
      &filecheck($ecdata_file);
      print "Done.\n";
    }
  }

#### now KEGG stuff
  if ($d_flag == 3  || ($d_flag==4)) {
    print "Downloading relevant files for KEGG-annotation:\n\n";
    $flag=1;
    if (-e "$keggpathway_file") {$flag = &delete_file("$keggpathway_file"); }
    if ($flag ==1) {  
      print "Please wait, downloading $keggpathway_file now ...\n";
      system("$getdata7");   
      &filecheck($keggpathway_file);
      print "Done.\n";
    }
    $flag=1;
    if (-e "$kegguniprot_file") {$flag = &delete_file("$kegguniprot_file"); }
    if ($flag ==1) {  
      print "Please wait, downloading $kegguniprot_file now ...\n";
      system("$getdata8");   
      &filecheck($kegguniprot_file);
      print "Done.\n";
    }
    $flag=1;
    if (-e "$keggko_file") {$flag = &delete_file("$keggko_file"); }
    if ($flag ==1) {  
      print "Please wait, downloading $keggko_file now ...\n";
      system("$getdata9");   
      &filecheck($keggko_file);
      print "Done.\n";
    }
    $flag=1;
    if (-e "$keggmap_file") {$flag = &delete_file("$keggmap_file"); }
    if ($flag ==1) {  
      print "Please wait, downloading $keggmap_file now ...\n";
      system("$getdata0");   
      &filecheck($keggmap_file);
      print "Done.\n\n";
    }
  }
  chdir ("../");
  print "Back to main menu.\n\n";
  &options;
}   


##################################################################################
#################  2-annodb  ###########################
##################################################################################

sub annodb() {
  my $a_flag = $_[0]; my $flag = 1;
  unless (-e "downloads") { 
    print "Sub-directory 'downloads' does not exist. Please make\n";
    print "sure you are starting the program from the correct\n";
    print "directory. - Exiting now\n\n";
    exit;
  }

#### postgres sanity check
  print "\n\n";
  &postmaster_check();

#### file checks first
  chdir ("downloads");
  unless (-e "$sprot_file") {&warn_miss($sprot_file)}
  unless (-e "$trembl_file") {&warn_miss($trembl_file)}

  if ($a_flag==1 || $a_flag==4) {
    unless (-e "$goasso_file") {&warn_miss($goasso_file)}
    unless (-e "$gomap_file") {&warn_miss($gomap_file)} 
    unless (-e "$godef_file") {&warn_miss($godef_file)} 
  }
  if ($a_flag==2 || $a_flag==4) {
    unless (-e "$ecdata_file") {&warn_miss($ecdata_file)} 
  }
  if ($a_flag==3 || $a_flag==4) {
    unless (-e "$kegguniprot_file") {&warn_miss($kegguniprot_file)} 
    unless (-e "$keggko_file") {&warn_miss($keggko_file)} 
    unless (-e "$keggpathway_file") {&warn_miss($keggpathway_file)} 
    unless (-e "$keggmap_file") {&warn_miss($keggmap_file)} 
  }

#### if arrived here all files should at least exist - 
#### now start with go-stuff    
  if ($a_flag==1 || $a_flag==4) {  
    print "\nLooking up $go_db database.  Please wait ...\n";
    my $conn=DBI->connect("dbi:Pg:dbname=$go_db", "", "", {PrintError => 0}); #Last two values would be user/pass.
    if (! $conn)   { # Couldn't connect to the database  
      print "\nCouldn't connect to the database $go_db, creating it now.\n";
      &create_godb($go_db,1,1,1,1);
    }  
    else     { 
      print "\nDatabase $go_db does already exist, do you want to update it? ";
      my $answer = &yes_no();
      if ($answer == 1) {
        my @table = $conn->tables('','',undef,'TABLE');
        my $GO_table_flag=1; 
        my $GOslim_flag=1; 
        my $GOdef_flag=1;      
        for(my $n=0; $n < @table; $n++)    {
          $table[$n] =~ s/public\.//; #get rid of "public." which is present in some versions of DBD.Pg 
          if ($table[$n] eq "go") { $GO_table_flag=0;    my $del = $conn->do("DELETE from go;"); my $del2 = $conn->do("drop index upid_index;");} 
	  if($table[$n] eq "go_slim") { $GOslim_flag=0;  my $del = $conn->do("DELETE from go_slim;");} 
          if($table[$n] eq "go_def") { $GOdef_flag=1;  my $del = $conn->do("DROP TABLE go_def;");} #deals with update from old versions of a8r 
        }                           
        &create_godb($go_db,0,$GO_table_flag,$GOslim_flag,$GOdef_flag);  ### arguments: <db, does db exist, does table exist>     
      }
      else   {
        if ($conn) {$conn->disconnect or warn "Disconnection failed: $DBI::errstr\n";} 
	print "You have decided not to continue. Exiting now.\n\n";
	exit();
      } 
    }
    if ($conn) {$conn->disconnect or warn "Disconnection failed: $DBI::errstr\n";} 
     
#### IEA or not IEA 
    my $iea_flag;
    print "\n\nNow it is your turn to decide on the level of evidence that is considered to\n";
    print "be sufficient for including the respective GO-term in your reference set. The\n"; 
    print "main issue are GO-terms which are \"inferred from electronic annotation\" (IEA).\n\n";
    print "Including them will greatly increase the number of sequences in your reference\n";
    print "database. Therefore, the number of annotated GO-terms for your dataset is also\n";
    print "likely to increase.\n";
    print "Excluding them will make annot8r's GO-annotation faster and the annotated terms\n";
    print "will be more reliable (as IEA is less reliable than expert curator annotation).\n";
    print "Essentially your choice comes down to \"Coverage\" vs \"Quality\".\n\n";
    print "Do you want to include IEA in your reference database (if unsure opt for no)? ";
    $iea_flag=yes_no();

#### Ready to import data, delete old entries if any there, mere upgrading doesn't make sense
    my @entry; my $result; my $flatfile;   my $dbh; my %gocount;
    open(FH,"<$goasso_file") ||  die "Can't open $goasso_file\n";  
    $conn = DBI->connect("dbi:Pg:dbname=$go_db", "", ""); 
    
    $dbh = $conn->prepare_cached("INSERT INTO go values (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?);") or die "Cannot prepare: " . $dbh->errstr();
#### $conn->do('COPY go FROM STDIN'); faster, but breaks old DBD:pg    
    print "Now inserting GO entries.\n";   
    print "Note: depending on the size of $goasso_file this step may take some time\n"; 
    print "(up to a few hours) Please wait ...\n"; 

#### Processing file line by line  
    my $dummy_1 = "zero"; my $dummy_4 = "GO:zero"; 
    my $n=0;
    while (my $line=<FH>) {
      @entry = split (/\t/, $line); 
      chomp ($entry[14]);
####  now IEA or not 
      unless ($iea_flag ==1) {
        my $check = $entry[6];       
        chomp $check;
        unless (exists $go_hash{$check}) {next;}
      }
      $entry[9] =~ s/'/\{prime\}/g;        
#### can't get quote to work properly do it this way
    
#### Now filter for redundancy: ie file contains some terms more than once ie via interpro and swprot
#### We have to trust the annotations which are in the flat file therefore we 
#### don't care whether they are coming from interpro or sptr or wherever. We are taking advantage of the order of the flat
#### file. ie we only have to compare with neighboring entry to identify redundancy. 
#### Note this does not work for upgrades, but 
#### it is strongly recommended to delete the old stuff anyway.
      unless (("$dummy_1" eq "$entry[1]") && ("$dummy_4" eq "$entry[4]"))  { ### keeping protentry
        $dummy_1 = $entry[1];
        $dummy_4 = $entry[4];
        $n++; # if ($n==50) {exit;}
##  columns are go standard, even if we only need a few of them putting all in db makes it more flexible and transferable
##  and potentially useful for other stuff - column names more or less reminescent of respective go stuff  
## 0->db_text 1->dbo_id 2->dbo_syn 3->_not 4->go_id 5->dbref 6->evid 7->w_f 8->asp 9->dbo_name 10->dbo_syn 11->dbo_typ 12->taxon 12-> 13->date 14->as_by
#####        $result = $conn->do("INSERT INTO go values ('$entry[0]','$entry[1]','$entry[2]','$entry[3]','$entry[4]','$entry[5]','$entry[6]',
#####        '$entry[7]','$entry[8]','$entry[9]','$entry[10]','$entry[11]','$entry[12]','$entry[13]','$entry[14]');", {PrintError => 0});
        $dbh->execute ($entry[0],$entry[1],$entry[2],$entry[3],$entry[4],$entry[5],$entry[6],
           $entry[7],$entry[8],$entry[9],$entry[10],$entry[11],$entry[12],$entry[13],$entry[14]);
#$conn->pg_putline ("$entry[0]\t$entry[1]\t$entry[2]\t$entry[3]\t$entry[4]\t$entry[5]\t$entry[6]\t
#        $entry[7]\t$entry[8]\t$entry[9]\t$entry[10]\t$entry[11]\t$entry[12]\t$entry[13]\t$entry[14]\n");

#### counting instances of each go term in db   
        $gocount{$entry[4]}++;
        printf("\r%9d entries inserted",$n);
      }       
    }
    close(FH);
    print "\n$n GO-entries added to $go_db.\n";

#### A lot of look-ups to be done, therefore indexing ...    
    print "\nCreating indices for table 'go'...\n";
    $result = $conn->do("CREATE INDEX upid_index ON go (dbo_id);");

#### now same story for GO_slim
    open(FH,"$gomap_file") ||  die "Can't open $gomap_file\n";  
    $result = $conn->do("DELETE from go_slim;");
    print "Now adding GO_slim entries ...\n";
        
#### Processing file line by line  
    $n=0;
    while (my $line=<FH>) {
      unless ($line =~ /^!/) {
        @entry = split (/\t/, $line); 
        $n++; 
        chomp ($entry[1]);
#### 0->go_id 1->go_slim_id 
        $result = $conn->do("INSERT INTO go_slim values ('$entry[0]','$entry[1]');", {PrintError => 0});
        printf("\r%9d entries inserted",$n);
      }       
    }
    close(FH);
    print "\n$n GO_slim entries added to $go_db.\n\n";
 
#### and for GO_def
    open(FH,"$godef_file") ||  die "Can't open $godef_file\n";  
    print "Inserting GO_def entries ...\n";   
    $result = $conn->do("DELETE from go_def;");
 
#### Processing file line by line  
    $n=0;
    while (my $line=<FH>) {
      unless ($line =~ /^!/) {
        @entry = split (/\t/, $line); 
        $entry[1] =~ s/'/\{prime\}/g;   
        chomp ($entry[2]);
        $n++; 
## 0->go_id 1->description 2->pcf 3->goidnum 
        $result = $conn->do("INSERT INTO go_def values ('$entry[0]','$entry[1]','$entry[2]','$gocount{$entry[0]}');", {PrintError => 0});
        printf("\r%9d entries inserted",$n);
      }       
    }
    close(FH);
    print "\n$n GO_def entries added to $go_db. GO-preparations done.\n\n";
    if ($conn) {$conn->disconnect or warn "Disconnection failed: $DBI::errstr\n";} 
 }  


#### now start with ec-stuff    
  if ($a_flag==2 || $a_flag==4) {  
  
#### is $database available?  
    print "\nLooking up $ec_db  Please wait ...\n";
    my $conn=DBI->connect("dbi:Pg:dbname=$ec_db", "", "", {PrintError => 0}); #Last two values would be user/pass.
    if (! $conn)   { # Couldn't connect to the database  
      print "\nCouldn't connect to the database $ec_db, creating it now.\n";
      &create_ecdb($ec_db,1,1,1);
    }
  
#### Check tables are present - create them if they are not, maybe allow user to define table name in next vs (for multiple tables)  
    else     { 
      print "\nDatabase $ec_db does already exist, do you want to update it? ";
      my $answer = &yes_no();
      if ($answer == 1) {
        my @table = $conn->tables('','',undef,'TABLE');
        my $EC_table_flag=1;  my $ECcount_table_flag=1;       
        for(my $n=0; $n < @table; $n++)    {
          $table[$n] =~ s/public\.//; #get rid of "public." which is present in some versions of DBD.Pg 
          if($table[$n] eq "ec") { $EC_table_flag=0; } 
          if($table[$n] eq "ec_count") { $ECcount_table_flag=0; } 
	}                           
        &create_ecdb($ec_db,0,$EC_table_flag,$ECcount_table_flag);  ### arguments: <db, create db, create table, create counttable>     
      }
      else   {if ($conn) {$conn->disconnect or warn "Disconnection failed: $DBI::errstr\n";} exit();} 
    }
    if ($conn) {$conn->disconnect or warn "Disconnection failed: $DBI::errstr\n";}
     
#### Ready to import data, delete old entries if there are any
    my @entry; my $result;    
    open(FH,"$ecdata_file") ||  die "Can't open $ecdata_file\n";  
    $conn = DBI->connect("dbi:Pg:dbname=$ec_db", "", ""); 
    $result = $conn->selectall_arrayref("select sw_id from ec;");
    my $ntuples=@$result;  
    if($ntuples == 0)  { 
      print "Inserting EC entries...\n";   
    }
    else  {
      print "\nEC entries already exist for $ec_db - Update the db? ";
      print "\nNote: this will delete old entries (recommended)";
      my $answer=yes_no();
      if($answer==1)   {
        print "Deleting old EC entries. Please wait ...\n";
        $result = $conn->do("DELETE from ec;"); my $del = $conn->do("drop index swsyn_index;");
        $result = $conn->do("DELETE from ec_count;");my $del2 = $conn->do("drop index ecid_index;");
	print "Old entries deleted.\n\n";
        print "Now updating EC entries.\n";
      }
      else   {
        if ($conn) {$conn->disconnect or warn "Disconnection failed: $DBI::errstr\n";} 
	print "You have decided not to continue. Exiting now.\n\n";
	exit();
      } 
    }    
    print "Please wait ...\n"; 
 
#### Processing file line by line  
    my $dbh;
    $dbh = $conn->prepare_cached("INSERT INTO ec (sw_id,sw_syn,ec_id,ec_des) values (?,?,?,?);") or die "Cannot prepare: " . $dbh->errstr();
    my (%sw_id, $sw_syn, $ec_id, $ec_des, %eccount);
    my $n=0;
    while (my $line=<FH>) {
      if ($line =~ /^\/\//) {
        $ec_id= ''; $ec_des= ''; 
      }  
      if ($line =~ /^ID\s+([0-9\.]+)/) {
         $ec_id = $1; 
      }
      if ($line =~ /^DE\s+/) {
         $ec_des = $'; 
         $ec_des =~ s/\s+//g;
         $ec_des =~ s/'/\{prime\}/g;        
      }
      if ($line =~ /^DR\s+/) {
        my $entry = $';   
        $entry =~ s/;\s+$//; 
        @entry = split /;/,$entry; 
        foreach $entry (@entry) { 
          $entry =~ s/^\s+//;
          my @pair = split /,/,$entry;
	  $sw_id{$pair[0]} =  $pair[1];
          if (defined ($pair[0])) {
            #my $result = $conn->do("INSERT INTO ec (sw_id,sw_syn,ec_id,ec_des) values ('$pair[0]','$sw_id{$pair[0]}', '$ec_id', '$ec_des');");
             $dbh->execute ($pair[0],$sw_id{$pair[0]},$ec_id,$ec_des);
	     $n++;
	     $eccount{$ec_id}++;
	     printf("\r%9d entries inserted.",$n);
	  }   
        }
      }       
    }
    close(FH);
    print "\n\nDatabase $ec_db successfully filled with $n entries.\n";

#### Now count table
    my $dbh2;
    $dbh2 = $conn->prepare_cached("INSERT INTO ec_count (ec_id,ec_id_num) values (?,?);") or die "Cannot prepare: " . $dbh->errstr();

    foreach my $key (keys %eccount)  {
      $dbh2->execute ($key,$eccount{$key});    
    }
    print "\n\nEC-count table updated.\n";



#### A lot of look-ups to be done, therefore indexing ...    
    print "Creating indices for table 'ec'...\n";
    $result = $conn->do("CREATE INDEX swsyn_index ON ec (sw_id);");
    print "Creating indices for table 'ec_count'...\n";
    $result = $conn->do("CREATE INDEX ecid_index ON ec_count (ec_id);");

    
    if ($conn) {$conn->disconnect or warn "Disconnection failed: $DBI::errstr\n";} 
  }
  

#### now KEGG stuff  
  if ($a_flag==3 || $a_flag==4) {  
 
#### is $database available?  
    print "\nLooking up $kegg_db  Please wait ...\n";
    my $conn=DBI->connect("dbi:Pg:dbname=$kegg_db", "", "", {PrintError => 0}); #Last two values would be user/pass.
    if (! $conn)   { # Couldn't connect to the database  
      print "\nCouldn't connect to the database $kegg_db, creating it now.\n";
      &create_keggdb($kegg_db,1,1,1,1);
    }
  
#### Check tables are present - create them if they are not, maybe allow user to define table name in next vs (for multiple tables)  
    else     { 
      print "\nDatabase $kegg_db does already exist, do you want to update it? ";
      my $answer = &yes_no();
      if ($answer == 1) {
        my @table = $conn->tables('','',undef,'TABLE');
        my $kegg_table_flag=1;my $map_table_flag=1;my $def_table_flag=1;        
        for(my $n=0; $n < @table; $n++)    {
          $table[$n] =~ s/public\.//; #get rid of "public." which is present in some versions of DBD.Pg 
          if($table[$n] eq "kegg") { $kegg_table_flag=0; } 
          if($table[$n] eq "kegg_komap") { my $del = $conn->do("DROP TABLE kegg_komap;");} #deals with update from old versions of a8r } 
          if($table[$n] eq "kegg_def") { $def_table_flag=0; } 
	}                           
        &create_keggdb($kegg_db,0,$kegg_table_flag,$map_table_flag,$def_table_flag);  ### arguments: <db, does db exist, do tables exist>     
      }
      else   {if ($conn) {$conn->disconnect or warn "Disconnection failed: $DBI::errstr\n";} exit();} 
    }
    
#### now loop through files, create three hashes holding the data     
    my %komps; my @pathways; my %uniprots;   
    open(FH,"$kegguniprot_file") ||  die "Can't open $kegguniprot_file\n";  
    print "Now processing $kegguniprot_file ...\n"; 
    while (my $line=<FH>) { 
      if ($line =~ /\w+:(\S+)\s+up:(\w+)/) {
        $uniprots{$1}=$2;
      }	
      else {print "not dealt with:\n$line\n";}
    }
    print "Done. \n";
    close FH;
       
    open(FH,"$keggko_file") ||  die "Can't open $keggko_file\n";  
    print "Now processing $keggko_file ...\n";
    while (my $line=<FH>) {
      if ($line =~ /\w+:(\S+)\s+ko:(\w+)/) { 
        $komps{$1}=$2;
      }	
       else {print "not dealt with:\n$line\n";}
    }
    print "Done.\n";
    close FH;
    
    if ($conn) {$conn->disconnect or warn "Disconnection failed: $DBI::errstr\n";}
#### now getting everything into db, based on uniprot entries, not really interested in the rest     
    $conn = DBI->connect("dbi:Pg:dbname=$kegg_db", "", ""); 
    my $result = $conn->selectall_arrayref("select up_id from kegg;");
    my $ntuples=@$result;  
    if($ntuples == 0)  { 
      print "Inserting KEGG entries\n";   
    }
    else  {
      print "\nKEGG entries already exist for $kegg_db - Update the db? ";
      print "\nNote: this will delete old entries (recommended)";
      my $answer=yes_no();
      if($answer==1)   {
        print "Deleting old KEGG entries. Please wait ...\n";
        $result = $conn->do("DELETE from kegg;");
        print "Old entries deleted.\n\n";
        print "Now updating KEGG entries.\n";
      }
    }    
    print "Please wait ...\n"; 

#### adding everything that has an associated pathway
    my $unid;  my $koid; my $n; my $dbh; my %keggcount;
    $dbh = $conn->prepare_cached("INSERT INTO kegg (up_id,ko_id) values (?,?);") or die "Cannot prepare: " . $dbh->errstr();
    foreach  my $key (keys %uniprots) {
      $unid = $uniprots{$key};
      if (defined $komps{$key}) {
        $koid=$komps{$key};      
        #my $result = $conn->do("INSERT INTO kegg (up_id,ko_id) values ('$unid','$koid');");
        $dbh->execute ($unid,$koid);
	$n++;
	$keggcount{$koid}++;
        printf("\r%9d entries inserted",$n);
      }	
    }  

#### A lot of look-ups to be done, therefore indexing ...    
    print "\nCreating indices for table 'kegg'...\n";
    $result = $conn->do("CREATE INDEX upid_index ON kegg (up_id);");

    
#### ko path lookup table
    $result = $conn->selectall_arrayref("select ko_id from kegg_komap;");
    $ntuples=@$result;  # Ie size of array referenced by result
    if($ntuples == 0)  { 
      print "Inserting KEGG KO mappping entries\n";   
    }
    else  {
      print "\nKEGG KO mapping entries already exist for $kegg_db - Update the db? ";
      print "\nNote: this will delete old entries (recommended)";
      my $answer=yes_no();
      if($answer==1)   {
        print "Deleting old KEGG KO mapping entries. Please wait ...\n";
        $result = $conn->do("DELETE from kegg_komap;");
        print "Old entries deleted.\n\n";
        print "Now updating KEGG KO mapping entries.\n";
      }
    }    
    print "Please wait ...\n"; 

    my $ko; my $path;
    open(FH,"$keggpathway_file") ||  die "Can't open $keggpathway_file\n";  
    print "Now processing and databasing $keggpathway_file ...\n";
    $n=0; my $write=0;
    while (my $line=<FH>) {
      if ($line =~ /^ENTRY\s+(K\d\d\d\d\d)/) {# print "$line";
        $ko = $1;
        $write=1;
      }
      if (($write == 1) && ($line =~ /\[PATH:ko(\d\d\d\d\d)\]/)) {
        $path = $1;   
        $result = $conn->do("INSERT INTO kegg_komap (ko_id,path_id,ko_id_num) values ('$ko','$path','$keggcount{$koid}');");
        $n++;
        printf("\r%9d entries inserted",$n);
      }
      if  ($line =~ /\/\/\//) {  $write=0;}
    }
#### A lot of look-ups to be done, therefore indexing ...    
    print "\nCreating indices for table 'kegg_komap'...\n";
    $result = $conn->do("CREATE INDEX koid_index ON kegg_komap (ko_id);");

####  map now
    #$conn = DBI->connect("dbi:Pg:dbname=$kegg_db", "", ""); 
    $result = $conn->selectall_arrayref("select path_id from kegg_def;");
    $ntuples=@$result;  # Ie size of array referenced by result
    if($ntuples == 0)  { 
      print "Inserting KEGG pathway entries\n";   
    }
    else  {
      print "\nKEGG pathway entries already exist for $kegg_db - Update the db? ";
      print "\nNote: this will delete old entries (recommended)";
      my $answer=yes_no();
      if($answer==1)   {
        print "Deleting old KEGG pathway entries. Please wait ...\n";
        $result = $conn->do("DELETE from kegg_def;");
        print "Old entries deleted.\n\n";
        print "Now updating KEGG pathway entries.\n";
      }
    }    
    print "Please wait ...\n"; 
    my $descr;
    open(FH,"$keggmap_file") ||  die "Can't open $keggmap_file\n";  
    print "Now processing and databasing $keggmap_file ...\n";
    while (my $line=<FH>) {
      if ($line =~ /(\d+)\s+(.+)/) {
        $path = $1;
	$descr=$2;
	$descr =~ s/'//g;   # messes up sql commands otherwise   
	my $result = $conn->do("INSERT INTO kegg_def (path_id,descr) values ('$path','$descr');");
      }	
    }
      

    print "Done - Back to main menu.\n";
    close FH;
    if ($conn) {$conn->disconnect or warn "Disconnection failed: $DBI::errstr\n";}
    
  }
  chdir ("../");
  options();    
}



##################################################################################
#################  3-BLAST preparations ##########################################
##################################################################################

sub prepare_blast () {
#### tidy up and re-format blastdb
  my $p_flag = $_[0]; my $flag = 1;
  unless (-e "blast_dbs") { system ("mkdir blast_dbs");}
  unless (-e "downloads") { 
    print "Sub-directory 'downloads' does not exist. Please make\n";
    print "sure you are starting the program from the correct\n";
    print "directory. - Exiting now\n\n";
    exit;
  }
  chdir ("downloads");
#### postgres sanity check
  print "\n\n";
  &postmaster_check();

#### file checks first
  unless (-e "$sprot_file") {&warn_miss($sprot_file)}
  unless (-e "$trembl_file") {&warn_miss($trembl_file)}

#### check databases    
  if ($p_flag==1 || $p_flag==4) {  
    my $conn=DBI->connect("dbi:Pg:dbname=$go_db", "", "", {PrintError => 0}); #Last two values would be user/pass.
    if (! $conn)   { #### Couldn't connect to the database  
      print "\nCouldn't connect to the database $go_db\n";
      exit;
    }      
    if ($conn) {$conn->disconnect or warn "Disconnection failed: $DBI::errstr\n";} 
  }
  if ($p_flag==2 || $p_flag==4) {  
    my $conn=DBI->connect("dbi:Pg:dbname=$ec_db", "", "", {PrintError => 0}); #Last two values would be user/pass.
    if (! $conn)   { #### Couldn't connect to the database  
      print "\nCouldn't connect to the database $ec_db\n";
      exit;    
    }
    if ($conn) {$conn->disconnect or warn "Disconnection failed: $DBI::errstr\n";} 
  }  
  if ($p_flag==3 || $p_flag==4) {  
    my $conn=DBI->connect("dbi:Pg:dbname=$kegg_db", "", "", {PrintError => 0}); #Last two values would be user/pass.
    if (! $conn)   { #### Couldn't connect to the database  
      print "\nCouldn't connect to the database $kegg_db\n";
      exit;    
    }  
    if ($conn) {$conn->disconnect or warn "Disconnection failed: $DBI::errstr\n";} 
  }

#### all relevant db exist - let's start with go 
  if ($p_flag==1 || $p_flag==4) {  
    print "Preparing BLAST/GO database.\n";
    print "... loading distinct sequences from $go_db\n";
    print "Please wait, may take some time\n";
    #### could add option for more tables here
    my $conn=DBI->connect("dbi:Pg:dbname=$go_db", "", "", {PrintError => 0}); #Last two values would be user/pass.
    my $result = $conn->prepare("select distinct dbo_id from go;");
    $result->execute();
    my @list;
    while (my $array_ref = $result->fetchrow_arrayref) {
      push @list, @$array_ref;
      printf("\rSo far %9d sequences extracted ...",scalar @list);
    }         
    if ($conn) {$conn->disconnect or warn "Disconnection failed: $DBI::errstr\n";} 
    my $number=@list;
    print "\n$number distinct sequences extracted\n";  
     
####  now building hash out of list for half way fast look up
    my %entries;
    foreach my $element(@list) {
      $entries{$element}='';
    }   
    open(FH,"$sprot_file") ||  die "Can't open $sprot_file\n";    
    chdir ("../blast_dbs");  
    unless (-e "$GO_blastfile") {system "touch $GO_blastfile"}
    open (WRITE, ">$GO_blastfile") || die "Can't open $GO_blastfile\n"; 
    print "\nWriting $GO_blastfile  ... Please wait\n";  
    my $write_flag = 0;
    
#### going through .fsa, writing all entries that are in db to blastnew.fsa   
    while (my $line=<FH>) {
      if ($line =~ /^>\w+\|(\w+)\|/)  {
        my $id = $1; # print "$id\n";
        if (exists($entries{$id}))  {
          $write_flag=1;
        }     
        else {$write_flag = 0;}
      }        
      if ($write_flag == 1) {
        print WRITE "$line"; 
      }
    }
    close FH;
    chdir ("../downloads");
    open(FH,"$trembl_file") ||  die "Can't open $trembl_file\n"; 
    $write_flag = 0;
    while (my $line=<FH>) {
      if ($line =~ /^>(\w+)/)  {
        my $id = $1; # print "$id\n";
        if (exists($entries{$id}))  {
          $write_flag=1;
        }     
        else {$write_flag = 0;}
      }        
      if ($write_flag == 1) {
        print WRITE "$line"; 
      }
    }  
    close WRITE;
    chdir ("../blast_dbs");

#### Format new blastdb
    print "\nFormatting $GO_blastfile  ...\n";
    system "formatdb -i $GO_blastfile";
    print "\nFormatting $GO_blastfile done, now ready for running blast\n";
    chdir ("../");
  }

#### ec-stuff next
  if ($p_flag==2 || $p_flag==4) {  
    print "Preparing BLAST/EC database.\n";
    print "... loading distinct sequences from $ec_db\n";
    print "Please wait, may take some time\n";
    my $conn=DBI->connect("dbi:Pg:dbname=$ec_db", "", "", {PrintError => 0}); #Last two values would be user/pass.
    my $result = $conn->prepare("select distinct sw_id from ec;");
    $result->execute();
    my @list;
    while (my $array_ref = $result->fetchrow_arrayref) {
      push @list, @$array_ref; 
      printf("\rSo far %9d sequences extracted ...",scalar @list);
    }
         
    if ($conn) {$conn->disconnect or warn "Disconnection failed: $DBI::errstr\n";} 
    my $number=@list;
    print "\n$number distinct sequences extracted\n";  
    
####  now building hash out of list for half way fast look up
    my %entries;
    foreach my $element(@list) {
      $element =~ s/\s+//g;
      $entries{$element}='';
    }   
    chdir ("downloads");
    open(FH,"$sprot_file") ||  die "Can't open $sprot_file\n";    
    chdir ("../blast_dbs");
    unless (-e "$EC_blastfile") {system "touch $EC_blastfile"}
    open (WRITE, ">$EC_blastfile") || die "Can't open $EC_blastfile\n"; 
    my $write_flag = 0;
    while (my $line=<FH>) {
      if ($line =~ /^>\w+\|(\w+)\|/)  {
        my $id = $1;  #print "$id\n";
        if (exists($entries{$id}))  {
          $write_flag=1;
        }     
        else {$write_flag = 0;}
      }        
      if ($write_flag == 1) {
        print WRITE "$line"; 
      }
    }
    close FH;
    chdir ("../downloads");
    open(FH,"$trembl_file") ||  die "Can't open $trembl_file\n";    
    chdir ("../blast_dbs");
    $write_flag = 0;
    while (my $line=<FH>) {
      if ($line =~ /^>(\w+)/)  {
        my $id = $1;  #print "$id\n";
        if (exists($entries{$id}))  {
          $write_flag=1;
        }     
        else {$write_flag = 0;}
      }        
      if ($write_flag == 1) {
        print WRITE "$line"; 
      }
    }
    close FH;
    close WRITE;
  
#### Format new blastdb
    print "\nFormatting $EC_blastfile  ...\n";
    system "formatdb -i $EC_blastfile";
    print "\nFormatting $EC_blastfile done, now ready for running blast\n";
    chdir ("../");
  }  

#### KEGG stuff finally
  if ($p_flag==3 || $p_flag==4) {  
    print "Preparing BLAST/KEGG database.\n";
    print "... loading distinct sequences from $kegg_db\n";
    print "Please wait, may take some time\n";
    my $conn=DBI->connect("dbi:Pg:dbname=$kegg_db", "", "", {PrintError => 0}); #Last two values would be user/pass.
    my $result = $conn->prepare("select distinct up_id from kegg;");
    $result->execute();
    my @list;
    while (my $array_ref = $result->fetchrow_arrayref) {
      push @list, @$array_ref; print "      @$array_ref";
      printf("\rSo far %9d sequences extracted ...",scalar @list);
    }
         
    if ($conn) {$conn->disconnect or warn "Disconnection failed: $DBI::errstr\n";} 
    my $number=@list;
    print "\n$number distinct sequences extracted\n";  

    my %entries;
    foreach my $element(@list) {
      $element =~ s/\s+//g;
      $entries{$element}='';
    }   
 
    chdir ("./downloads");
    open(FH,"$sprot_file") ||  die "Can't open $sprot_file\n";    
    chdir ("../blast_dbs");
    unless (-e "$KEGG_blastfile") {system "touch $KEGG_blastfile"}
    open (WRITE, ">$KEGG_blastfile") || die "Can't open $KEGG_blastfile\n"; 
    my $write_flag = 0;
    while (my $line=<FH>) {
      if ($line =~ /^>\w+\|(\w+)\|/)  {
 	my $id = $1;
        if (exists($entries{$id}))  {
          $write_flag=1;
        }     
        else {$write_flag = 0;}
      }        
      if ($write_flag == 1) {
        print WRITE "$line"; 
      }
    }
    close FH;
    chdir ("../downloads");
    open(FH,"$trembl_file") ||  die "Can't open $trembl_file\n";    
    chdir ("../blast_dbs");
    $write_flag = 0;
    while (my $line=<FH>) {
      if ($line =~ /^>\w+ \((\w+)/)  {
        my $id = $1; 
      	$line = ">" . "$id" . " dummy description\n"; #KEGG uses synonym ids, modifying FASTA description line makes life much easier later on 
	if (exists($entries{$id}))  {
          $write_flag=1;
        }     
        else {$write_flag = 0;}
      }        
      if ($write_flag == 1) {
        print WRITE "$line";
      }
    }
    close FH;
    close WRITE;
  
#### Format new blastdb
    print "\nFormatting $KEGG_blastfile  ...\n";
    system "formatdb -i $KEGG_blastfile";
    print "\nFormatting $KEGG_blastfile done, now ready for running blast\n";
    chdir ("../");
  }  

#### done
  options();
}


##################################################################################
#################  4-BLAST run   #################################################
##################################################################################

sub blast () {
  my $b_flag = $_[0];
  my $blast_exec=find_program("blastall");
  my $blast_method;
  unless (-e "blast_out") { system ("mkdir blast_out");}

####  check blastdbs
  if ($b_flag==1 || $b_flag==4) { 
    unless (-e "blast_dbs/$GO_blastfile") {
      print "Couldn't find blast_dbs/$GO_blastfile\n";
      print "please make sure you are in the appropriate directory.\nExiting now\n\n";
      exit;
    }
  }
  if ($b_flag==2 || $b_flag==4) { 
    unless (-e "blast_dbs/$EC_blastfile") {
      print "Couldn't find blast_dbs/$EC_blastfile\n";
      print "please make sure you are in the appropriate directory.\nExiting now\n\n";
      exit;
    }
  } 
  if ($b_flag==3 || $b_flag==4) { 
    unless (-e "blast_dbs/$KEGG_blastfile") {
      print "Couldn't find blast_dbs/$KEGG_blastfile\n";
      print "please make sure you are in the appropriate directory.\nExiting now\n\n";
      exit;
    }
  }

#### get sequences, add option for getting sequences from db later
  print "\nPlease enter the name of the file containing the fasta sequences ";
  print "\nyou want to annotate. For example ";
  print colored ("translations.fsa","green bold");  
  print "\nPlease give either the full or relative location of the file.\n";  
  $seqs2blast = &get_file();
  

#### protein or nucleotide
  print "\nDoes $seqs2blast contain protein or nucleotide sequences ?\n";
  print "Enter '1' for protein or '2' for nucleotide sequences.\n";
  my $flag=0; my $answer;
  while($flag==0) {
    $answer=<>;
    if($answer=~/^[1|2]$/) { $flag=1; next; }
    else {print " You have entered $answer This is not an option. Please try again\n";}
  }
  if ($answer == 1) { $blast_method = "blastp"; }
  else {$blast_method = "blastx";}
  
#### get e-value - removed - user's choice now at parsing level ...

#### get blastenvironment right
  my $here = `pwd`;
  chomp $here; 
  $here .= "/blast_dbs";
  my $seqfile = $seqs2blast; 
  $seqfile =~ s/.+\///g;
#  $ENV{'BLASTDB'}=$here;
  print "\nBlasting now ...\n";
  print "Please note this step may take some time.\n";
#  chdir ("./blast_out"); 
  if ($b_flag==1 || $b_flag==4) { 
    my $blastbase = "$here/" . "$GO_blastfile";
    $blast_output = "blast_out/" ."$seqfile" . "_GO.out"; 
    &run_blast	($seqs2blast, $blastbase, $blast_method, $blast_output);
  }
  if ($b_flag==2 || $b_flag==4) { 
    my $blastbase = "$here/" . "$EC_blastfile";
    $blast_output = "blast_out/" . "$seqfile" . "_EC.out"; 
    &run_blast	($seqs2blast, $blastbase, $blast_method, $blast_output);
  } 
  if ($b_flag==3 || $b_flag==4) { 
    my $blastbase = "$here/" . "$KEGG_blastfile";
    $blast_output = "blast_out/" . "$seqfile" . "_KEGG.out"; 
    &run_blast	($seqs2blast, $blastbase, $blast_method, $blast_output);
  }
  print "\nBlasting done. Ready for annotation now\n";
#  chdir ("../");

  options();  
}


##################################################################################
#################  5- annotation  ################################################
##################################################################################

sub anno()  {
  my $a_flag = $_[0];
  my $database;
  my @golist;
  my @eclist;
  my @kegglist;
  my $flag;

  unless (-e "output") { system ("mkdir output");}


#### First some checks...pg running,  swtrgo db, partigene db, blastoutput   
#### Checking for BLAST results first, might save some hussle   
  chdir ("blast_out");
  if ($a_flag == 1 || $a_flag==4) {  
    print "Checking for BLAST results required for GO-annotation ...\n";
    sleep(1);
    @golist = glob ("*GO.out");
    my $gos = scalar @golist;
    my $dir = `pwd`;
    if ($gos == 1) {
      print "\n$golist[0] found - OK\n";
    }
    elsif ($gos > 1)   {
      print "Following files have been found in $dir:\n";
      foreach my $file (@golist)  {
        print "$file\n";
      }
      print "\nDo you want to use all of them? \nTyping \"y\" will use all of them,\n";
      print "Typing \"n\" will exit the program and allow you to remove some files.\n";
      $flag  = &yes_no();
      unless ($flag == 1) {exit;};
      print \"OK - using all files ...\n";
    }
    else {
      print "No files ending with \"GO.out\" found in $dir \n";
      print "Exiting now. \n";
      exit;
    }    
  }
  if ($a_flag == 2 || $a_flag==4) {  
    print "Checking for BLAST results required for EC-annotation ...\n";
    sleep(1);
    @eclist = glob ("*EC.out");
    my $ecs = scalar @eclist;
    my $dir = `pwd`;
    if ($ecs == 1) {
      print "\n$eclist[0] found - OK\n";
    }
    elsif ($ecs > 1)   {
      print "Following files have been found in $dir:\n";
      foreach my $file (@eclist)  {
        print "$file\n";
      }
      print "\nDo you want to use all of them? \nTyping \"y\" will use all of them,\n";
      print "Typing \"n\" will exit the program and allow you to remove some files.\n";
      $flag  = &yes_no();
      unless ($flag == 1) {exit;};
      print \"OK - using all files ...\n";
    }
    else {
      print "No files ending with \"EC.out\" found in $dir \n";
      print "Exiting now. \n";
      exit;
    }    
  }
  if ($a_flag == 3 || $a_flag==4) {  
    print "Checking for BLAST results required for KEGG-annotation ...\n";
    sleep(1);
    @kegglist = glob ("*KEGG.out");
    my $keggs = scalar @kegglist;
    my $dir = `pwd`;
    if ($keggs == 1) {
      print "\n$kegglist[0] found - OK\n";
    }
    elsif ($keggs > 1)   {
      print "Following files have been found in $dir:\n";
      foreach my $file (@kegglist)  {
        print "$file\n";
      }
      print "\nDo you want to use all of them? \nTyping \"y\" will use all of them,\n";
      print "Typing \"n\" will exit the program and allow you to remove some files.\n";
      $flag  = &yes_no();
      unless ($flag == 1) {exit;};
      print \"OK - using all files ...\n";
    }
    else {
      print "No files ending with \"KEGG.out\" found in $dir \n";
      print "Exiting now. \n";
      exit;
    }    
  }
  chdir ("../");

#### Now database checks      
  &postmaster_check;
    
#### is GO annotation database available and does it contain go-table ?  
  if ($a_flag == 1 || $a_flag==4) { 
    my $conn=DBI->connect("dbi:Pg:dbname=$go_db", "", "", {PrintError => 0}); #Last two values would be user/pass.  
    if (! $conn)   { 
      print "\nCouldn't connect to the database $go_db .";
      &query_for_exit();
    }
    else {
      my @table = $conn->tables('','',undef,'TABLE');
      my $GO_table_flag=0;        
      my $GOslim_table_flag=0;
      for(my $n=0; $n < @table; $n++)    {
        $table[$n] =~ s/public\.//; #get rid of "public." which is present in some versions of DBD.Pg 
        if($table[$n] eq "go") { $GO_table_flag=1; } 
        if($table[$n] eq "go_slim") { $GOslim_table_flag=1; }
      }
      if ($GO_table_flag==0) {
        print colored("\n$go_db does not contain GO table, please check $go_db .\n","red bold"); 
        &query_for_continue();     
      }
      if ($GOslim_table_flag==0) {
        print colored("\n$database does not contain GO_slim table, GO_slim annotation will not be available.","red bold"); 
      }  
      if ($conn) {$conn->disconnect or warn "Disconnection failed: $DBI::errstr\n";} 
    }   
  }

#### is EC annotation database available ?  
  if ($a_flag == 2 || $a_flag==4) { 
    my $conn=DBI->connect("dbi:Pg:dbname=$ec_db", "", "", {PrintError => 0}); #Last two values would be user/pass.  
    if (! $conn)   { 
      print "\nCouldn't connect to the database $ec_db .";
      &query_for_exit();
    }
    else {
      my @table = $conn->tables('','',undef,'TABLE');
      my $EC_table_flag=0;
      my $ECcount_table_flag=0;        
      for(my $n=0; $n < @table; $n++)    {
        $table[$n] =~ s/public\.//; #get rid of "public." which is present in some versions of DBD.Pg 
        if($table[$n] eq "ec") { $EC_table_flag=1; } 
	if($table[$n] eq "ec_count") { $ECcount_table_flag=1; } 
      }
      if ($EC_table_flag==0) {
        print colored("\n$ec_db does not contain ec table, please check $ec_db .\n","red bold"); 
        &query_for_continue();     
      }
      if ($ECcount_table_flag==0) {
        print colored("\n$ec_db does not contain ec_count table, please check $ec_db .\n","red bold"); 
        &query_for_continue();     
      }

      if ($conn) {$conn->disconnect or warn "Disconnection failed: $DBI::errstr\n";} 
    }   
  }
#### is KEGG annotation database available  ?  
  if ($a_flag == 3 || $a_flag==4) { 
    my $conn=DBI->connect("dbi:Pg:dbname=$kegg_db", "", "", {PrintError => 0}); #Last two values would be user/pass.  
    if (! $conn)   { 
      print "\nCouldn't connect to the database $kegg_db .";
      &query_for_exit();
    }
    else {
      my @table = $conn->tables('','',undef,'TABLE');
      my $KEGG_table_flag=0;        
      my $KEGGmap_table_flag=0;
      for(my $n=0; $n < @table; $n++)    {
        $table[$n] =~ s/public\.//; #get rid of "public." which is present in some versions of DBD.Pg 
        if($table[$n] eq "kegg") { $KEGG_table_flag=1; } 
        if($table[$n] eq "kegg_def") { $KEGGmap_table_flag=1; }
      }
      if ($KEGG_table_flag==0) {
        print colored("\n$database does not contain kegg table, please check $kegg_db .\n","red bold"); 
        &query_for_continue();     
      }
      if ($KEGGmap_table_flag==0) {
        print colored("\n$database does not contain kegg_def table, please check $kegg_db .\n .","red bold"); 

      }  
      if ($conn) {$conn->disconnect or warn "Disconnection failed: $DBI::errstr\n";} 
    }   
  }

#### Now get database ready to actually store annotations
#### check for partigene.conf file and partigene db
  my $filename = "~/.partigene.conf";
  my $pg_database; my $pg_flag = 0;
  if  (-e "$filename") {
    open (CONFILE,"$filename") ||  die "Can't open configuration file\n";
    while (my $line=<CONFILE>) {      
      if ($line=~/^DATABASE\=(.+)/i) { $pg_database=$1; }       
    }
    close (CONFILE);
    
#### user has a Partigene db defined in the config file
    print "\nDo you want to use PartiGene database: $pg_database to store your annotations?";
    $pg_flag  = &yes_no();
  }
   
#### don't use Pg database  
  if ($pg_flag != 1) {
    print "What database do you want to create or use?\n";
    $pg_database = <STDIN>; 
    chomp $pg_database; 
  }  
  my $conn=DBI->connect("dbi:Pg:dbname=$pg_database", "", "", {PrintError => 0}); #Last two values would be user/pass.
  my $create=0;
  if (! $conn)   { ### Couldn't connect to the database  
    print "\nCouldn't connect to the database $pg_database";
    print "\nDo you want to create it?";
    $create = &query_for_continue;
    &create_pgdb($pg_database,"$create","1","1","1");  ### arguments: <db, create db , create tables>  
  }

#### db exists already, create relevant tables if necessary 
  else {
    my @table = $conn->tables('','',undef,'TABLE');
    my $go_flag=1;
    my $ec_flag=1;
    my $kegg_flag=1;
    for(my $n=0; $n < @table; $n++)    {
      $table[$n] =~ s/public\.//; #get rid of "public." which is present in some versions of DBD.Pg 
      if($table[$n] eq "a8r_blastgo") { $go_flag=0; } 
      if($table[$n] eq "a8r_blastec") { $ec_flag=0; } 
      if($table[$n] eq "a8r_blastkegg") { $kegg_flag=0; } 
    }    

#### GO stuff
    if ($a_flag == 1 || $a_flag==4) {  
      if (($go_flag)) {
        print "\nCreating relevant go tables for $pg_database ..."; 
        &create_pgdb($pg_database,"0","$go_flag","0","0");
      }   
      else {
        print "\na8r_go entries already exist for $pg_database.";
        print "\nDo you want to remove all old entries? (recommended)";
        print "\n\n- Otherwise the new entries will be added to existing entries.";
        my $answer=yes_no();
        if($answer==1)   {#delete existing entries
          print "Deleting old entries. Please wait ...\n";
          my $result = $conn->do("DELETE from a8r_blastgo;");	  
          print "Old GO-entries have been deleted.\n\n";                 
        }
      }  
      if ($conn) {$conn->disconnect or warn "Disconnection failed: $DBI::errstr\n";}
    }   
  
#### EC stuff
    if ($a_flag == 2 || $a_flag==4) {  
      if ($ec_flag==1) {
        print "\nCreating relevant ec tables for $pg_database ..."; 
        &create_pgdb($pg_database,"0","0","$ec_flag","0");
      }   
      else {
        print "\na8r_ec entries already exist for $pg_database.";
        print "\nDo you want to remove all old entries? (recommended)";
        print "\n\n- Otherwise the new entries will be added to existing entries.";
        my $answer=yes_no();
        if($answer==1)   {#delete existing entries
          print "Deleting old entries. Please wait ...\n";
          my $result = $conn->do("DELETE from a8r_blastec;");
          print "Old EC-entries have been deleted.\n\n";                 
        }
      }  
      if ($conn) {$conn->disconnect or warn "Disconnection failed: $DBI::errstr\n";}
    }
  
#### Kegg stuff
    if ($a_flag == 3 || $a_flag==4) {  
      if (($kegg_flag==1) ) {
        print "\nCreating relevant KEGG tables for $pg_database ..."; 
        &create_pgdb($pg_database,"0","0","0","$kegg_flag");
      }   
      else {
        print "\na8r_kegg entries already exist for $pg_database.";
        print "\nDo you want to remove all old entries? (recommended)";
        print "\n\n- Otherwise the new entries will be added to existing entries.";
        my $answer=yes_no();
        if($answer==1)   {#delete existing entries
          print "Deleting old entries. Please wait ...\n";
	  my $result = $conn->do("DELETE from a8r_blastkegg;");
          print "Old KEGG-entries have been deleted.\n\n";                 
        }
      }  
      if ($conn) {$conn->disconnect or warn "Disconnection failed: $DBI::errstr\n";}
    }
  }

#### cutoffs
  print "Before parsing the BLAST results you need to define \n";
  print "an \"acceptable level of similarity\". \n";
  print "This is done by setting a cut-off. Only BLAST hits \n";
  print "better than this cut-off will be accepted for \n";
  print "annotation. You can either use a cut-off based on\n";
  print "e-values or based on BLAST scores.\n"; 
    
#### score or e-value
  my $e_cut = 10; my $e_flag = 0;
  my $s_cut = 1; my $s_flag =0;
  print "\nTo use a cutoff based on e-values enter '1', to\n";
  print "use a cut-off based on scores enter '2'\n";
  $flag=0; my $answer;
  while($flag==0) {
    $answer=<>;
    if($answer=~/^[1|2]$/) { $flag=1; next; }
    else {print "\nYou have entered $answer .\nThis is not an option. Please try again\n";}
  }
  if ($answer == 1) { $e_flag = 1; }
  else {$s_flag = 1;}
  
  
####  I know this is patronising, but it can't be emphasized enough - and at least we give a shut up option^^
  unless ($quiet == 1) {
    print "Before we continue a few notes on cut-offs ...\n\n";
    sleep(1);
    print "There are no fix rules how to set\n";
    print "meaningful cut-offs.\n\n";
    sleep(1);
    print "For protein families which are functionally very \n";
    print "conserved, but divergent in sequence, an e-value\n";
    print "cutoff of 1e-03 may be strict enough to obtain\n";
    print "meaningful annotation.\n\n";
    sleep(1);
    print "For protein families with multiple or divergent \n";
    print "functions you might get erronous annotations even\n";
    print "when using an e-value cut-off of 1e-50.\n\n";
    sleep (1);
    print "And as yet we have not even mentioned the multi-domain\n";
    print "problem, or the issue of erroneous annotations present\n";
    print "in the underlying databases.\n\n";
    sleep (3);
    print "If you do not know what to do now, try 1e-08 as an \n";
    print "e-value cut-off, or 55 as a BLAST bits-score cut-off. \n";
    print "But use your \"biological common sense\" when interpreting\n";
    print "your results. Be aware that the annotations will be \n";
    print "incomplete and that some of the annotations will \n";
    print "probably be wrong.\n\n\n";
    sleep(1);
  }
#### get cut-offs
  if ($e_flag ==1) {
    print "\nPlease enter the e-value you want to use as a cut-off for annotation.\n";
    my $e_check = 1;
    while ($e_check == 1) {
      $e_cut = <STDIN>;
      chomp $e_cut;    
      if ($e_cut =~ /^\d+e-\d+$/) {### matches e- type
        print "\n$e_cut accepted";
        $e_check = 2;
      }
      elsif ($e_cut =~ /^\d*\.?\d*$/) {### matches numerical type
        print "\n$e_cut accepted";
        $e_check = 2;
      }  
      else {print "\n$e_cut not accepted. Please try again.\n";}
    } 
  }
  elsif ($s_flag ==1) {
    print "\nPlease enter the bits score you want to use as a cut-off for annotation.\n";
    my $s_check = 1;
    while ($s_check == 1) {
      $s_cut = <STDIN>;
      chomp $s_cut;    
      if ($s_cut =~ /^\d*\.?\d*$/) {### matches numerical type
        print "\n$s_cut accepted";
        $s_check = 2;
      }  
      else {print "\n$s_cut not accepted. Please try again.\n";}
    }  
  }
  else {print "\n\nSomething strange happened ...\n"; sleep(2);exit}




#### number of hits to be used
  print "\n\nIn addition to the \"traditional\" BLAST-tophit annotation\n";
  print "annot8r allows you to take up to 50 hits per sequence into account.\n";
  print "Each annotation term will have the best hit and the number of\n";
  print "hits pointing to this particular term as additional information\n";
  print "available. annot8r will also list the number of sequences in the\n";
  print "reference database sharing this annotation term and in case of\n";
  print "conflicting information the fraction of matches for each term.\n\n"; 
  
  print "Please enter the number of hits you want to use for\n";
  print "annotation (min 1 - max 50).\n";
  my $hit_cut;
  my $n_check = 1;
  while ($n_check == 1) {
    $hit_cut = <STDIN>;
    chomp $hit_cut;    
    if ($hit_cut =~ /^\d+$/) {
      if (($hit_cut > 0) || ($hit_cut < 51)) {
        print "\n$hit_cut accepted";
        $n_check = 2;
      }	
    }
    unless ($n_check==2) {print "\n$hit_cut not accepted, Please try again.\n";}
  }  

#### finally we are all ready to do the work ...  
  chdir ("blast_out");

#### go parsing first
  if ($a_flag == 1 || $a_flag==4) {  

#### get file ready
    chdir ("../output");
    if (-e "$goout") {
      print colored("\nWARNING: $goout exists already\n","red bold");
      my $old = $goout;
      $old =~ s/csv/old/;
      system ("mv $goout $old");
      print "$goout moved to $old\n";
    }	
    open(FILE,">$goout");  
    chdir ("../blast_out");
#### db
    my $conn1=DBI->connect("dbi:Pg:dbname=$pg_database", "", "", {PrintError => 0}); #Last two values would be user/pass.
    my $conn2=DBI->connect("dbi:Pg:dbname=$go_db", "", "", {PrintError => 0}); #Last two values would be user/pass.  
    print "\nNow parsing GO-results. Please wait ...\n";  
    foreach my $file (@golist) {
      my $in = new Bio::SearchIO( -format => 'blast',  
        		    -file   => "$file");  
      my $slim=''; 
      my $score='';  
      my $sig='';
      my $db_name=''; 
      my $acc=''; 
      my $prot = ''; 
      my $desc='';
      my $name='';
      my $prog='';
      my $success='';  
      my $n = 0;
      my $m = 0;
      
      $success=$conn2->prepare_cached('SELECT go_id from go where dbo_id=?');

#### each query
      while( my $result = $in->next_result )  {
        my %hitnumber =();
        my %besthit =();
	my %besteval=();
	my %bestscore=();
	$n++; 

#### each hit
        my $hit_nr=0; my $go_id;
#        $success->bind_columns(undef, \$go_id); 
	while (my $hit = $result->next_hit) {
	  $hit_nr++;
	  $db_name=$result->database_name;
          $prot=$result->query_name;           
          $acc=$hit->accession; 
	  $score=$hit->bits;
          $sig=$hit->significance;
          if ($acc =~ /^(\w+)/) {$acc = $1}
	  else {print "\nWARNING, can not deal with $acc\n";}	  
	  #### e-value filter here
	  if ($sig > $e_cut) {last;}
	  if ($score < $s_cut) { last;} 
	  if ($hit_cut < $hit_nr) {last;}


#	  $success=$conn2->prepare("SELECT go_id, evid, asp from go where dbo_id='$acc'");      
          $success->execute($acc);
          my $array_ref = $success->fetchall_arrayref();		 
         foreach my $output(@$array_ref)  {        
            my ($go_id) =@$output;
	    $hitnumber{$go_id}++;
	    unless ($besthit{$go_id}) {$besthit{$go_id}=$acc;}
	    unless ( $besteval{$go_id}) {$besteval{$go_id}=$sig;}	
	    unless ( $bestscore{$go_id}) {$bestscore{$go_id}=$score;}
          }
        }
#### loop through keys, find go slims and database stuff

          my $insert;
	  $insert = $conn1->prepare_cached("INSERT INTO a8r_blastgo values (?,?,?,?,?,?,?,?,?,?,?);") or die "Cannot prepare: " . $insert->errstr();


        foreach my $goresult (keys %besthit) {
   	  my $success2=$conn2->prepare("SELECT slim from go_slim where go_id='$goresult'");
	  $success2->execute();
#	  print "goresult: $goresult\n";
	  my $go_ref = $success2->fetchrow_arrayref();
	  if ($go_ref) {
	    ($slim)=@$go_ref; 
          }
	  else {$slim = "NA"}
#	  print "slim: $slim\n";
	  my $success3=$conn2->prepare("SELECT descr,pcf,go_id_num from go_def where go_id='$goresult'");
          $success3->execute();
	  my $def_ref = $success3->fetchrow_arrayref();
#### catch inconcistencies here
          if ($def_ref) {	  
	    my($descr,$pcf,$goidnum)=@$def_ref;# print "slim: $slim\n";
	    my  $goround = $goidnum;
	    if ($goidnum> $hit_cut) {$goround=$hit_cut}
	    my $fraction = $hitnumber{$goresult}/$goround;
	    my $round = sprintf("%.2f", $fraction);
	    #print"$prot a $db_name b $acc c $sig d  $go_id plus $evid\n";    
	    my $check = $conn1->do("DELETE from a8r_blastgo WHERE pept_id='$prot' AND go_term='$goresult';");
#	    my $insert = $conn1->do("INSERT INTO a8r_blastgo values ('$prot','$goresult','$pcf','$descr','$slim','$besthit{$goresult}','$bestscore{$goresult}','$besteval{$goresult}','$hitnumber{$goresult}');", {PrintError => 0});
            $insert->execute ($prot, $goresult, $pcf, $descr, $slim, $besthit{$goresult}, $bestscore{$goresult}, $besteval{$goresult}, $hitnumber{$goresult}, $goidnum, $round);

            print FILE "$prot, $goresult, $pcf, \"$descr\", $slim, $besthit{$goresult}, $bestscore{$goresult}, $besteval{$goresult}, $hitnumber{$goresult}, $goidnum, $round\n";        

#print "INSERT INTO a8r_blastgo values ('$prot','$goresult','$pcf','$descr','$slim','$besthit{$goresult}','$besteval{$goresult}','$hitnumber{$goresult}');"
          }
	  else {print colored("\nWARNING: GO-Term $goresult can't be associated with entry in table go_def - IGNORED\n","red bold"); }
        }      
        printf("\r%9d sequences processed so far",$n);  
      }
    }
  

#### at the end of this round go-slim stuff      
    print "\nNow writing GO-slim output files.\n\n";  
#### preparing outputfiles for piecharts
    my @pie_flag = ("P", "C", "F");
    foreach (@pie_flag) {
      my $filename = "piedata_" . "$_";
      my $backup = "$filename" . ".old";
      if (-e "$filename") {system "mv $filename $backup";}
    }   
  
#### extracting data for pie charts
    my $success = $conn1->prepare("SELECT DISTINCT pept_id, slim FROM a8r_blastgo");
    $success->execute();
    my $go_ref = $success->fetchall_arrayref();
    my %go_slim; my $go_slim; my $prot_id;

#### hash: go_slim id => number of entries in db  
    foreach my $output(@$go_ref) {
      ($prot_id, $go_slim) = @$output;
      #print "$prot_id, $go_slim\n";
      if (exists $go_slim{$go_slim}) {
        $go_slim{$go_slim}++;
      }      
      else {$go_slim{$go_slim}=0}
    }

##### create 3 files, loop through keys, write info to files
    chdir ("../output");
    open (FHC, ">piedata_C");
    print FHC "### go_slim for $pg_database - cellular component ###\n";
    print FHC "### go_slim_id \t description\t occurences\n";
    open (FHP, ">piedata_P");
    print FHP "### go_slim for $pg_database - biological process ###\n";
    print FHP "### go_slim_id \t description\t occurences\n";
    open (FHF, ">piedata_F");
    print FHF "### go_slim for $pg_database - molecular function ###\n";
    print FHF "### go_slim_id \t description\t occurences\n";
        
    foreach my $key (keys %go_slim) { 
      $success = $conn2 ->prepare("SELECT go_id, descr, pcf FROM go_def WHERE go_id='$key'"); 
      $success->execute(); 
      my $pie_ref = $success->fetchrow_arrayref();
      if ($pie_ref) {
        my ($go_id, $desc, $pcf)=@$pie_ref;  
        if ($pcf eq "C") {print FHC "$key \t $desc \t $go_slim{$key}\n";}  
        elsif ($pcf eq "P") {print FHP "$key \t $desc \t $go_slim{$key}\n";}
        elsif ($pcf eq "F") {print FHF "$key \t $desc \t $go_slim{$key}\n";}
        else {print colored ("WARNING: no valid GO-aspect (C, P or F) found for $key\n", "red bold")}
      }
      else {
        print colored("\nWARNING: couldn't process GO-slim term $key\n","red bold");  
      }
    }
    close FHC;
    close FHP;
    close FHF;
    chdir ("../blast_out");
    $success->finish(  );
    if ($conn2) {$conn2->disconnect or warn "Disconnection failed: $DBI::errstr\n";}
    if ($conn1) {$conn1->disconnect or warn "Disconnection failed: $DBI::errstr\n";}
  }

#### EC-parsing
  if ($a_flag == 2 || $a_flag==4) {  

#### get file ready
    chdir ("../output");
    if (-e "$ecout") {
      print colored("\nWARNING: $ecout exists already\n","red bold");
      my $old = $ecout;
      $old =~ s/csv/old/;
      system ("mv $ecout $old");
      print "$ecout moved to $old\n";
    }	
    open(FILE,">$ecout");  
    chdir ("../blast_out");
#### db
    my $conn1=DBI->connect("dbi:Pg:dbname=$pg_database", "", "", {PrintError => 0}); #Last two values would be user/pass.
    my $conn2=DBI->connect("dbi:Pg:dbname=$ec_db", "", "", {PrintError => 0}); #Last two values would be user/pass.  
    print "\nNow parsing EC-results. Please wait ...\n";  
    foreach my $file (@eclist) {
      my $in = new Bio::SearchIO( -format => 'blast',  
        		    -file   => "$file");  
      my $slim=''; 
      my $score='';  
      my $sig='';
      my $db_name=''; 
      my $acc=''; 
      my $prot = ''; 
      my $desc='';
      my $name='';
      my $prog='';
      my $success='';my $success2='';  
      my $n = 0;
      my $m = 0;

#### each query
     $success=$conn2->prepare_cached('SELECT ec_id, ec_des from ec where sw_id~?');      
     $success2=$conn2->prepare_cached('SELECT ec_id_num from ec_count where ec_id=?');


      while( my $result = $in->next_result )  {
        my %hitnumber =();
        my %besthit =();
	my %besteval=();
        my %bestscore=();
	my %bestdesc=();
	$n++;
        printf("\r%9d sequences processed so far",$n);  
       
#### each hit
        my $hit_nr=0;
	while (my $hit = $result->next_hit) {
	  $hit_nr++;
	  $db_name=$result->database_name;
          $prot=$result->query_name;           
          $acc=$hit->accession; 
          $score=$hit->bits;
	  $sig=$hit->significance;
          if ($acc =~ /^(\w+)/) {$acc = $1}
	  else {print "\nWARNING, can not deal with $acc\n";}	  
	  #### e-value filter here
	  if ($sig > $e_cut) {last;}
	  if ($score < $s_cut) {last;}
  	  if ($hit_cut < $hit_nr) {last;}
#	  $success=$conn2->prepare_cached('SELECT ec_id, ec_des from ec where sw_id~?');      
          $success->execute($acc);
          my $array_ref = $success->fetchall_arrayref();		 
          foreach my $output(@$array_ref)  {        
            my ($ec_id, $ec_des) =@$output;
	    $hitnumber{$ec_id}++;
	    unless ($besthit{$ec_id}) {$besthit{$ec_id}=$acc;}
	    unless ($bestscore{$ec_id}) {$bestscore{$ec_id}=$score;}
	    unless ($besteval{$ec_id}) {$besteval{$ec_id}=$sig;}
	    unless ($bestdesc{$ec_id}) {$bestdesc{$ec_id}=$ec_des;}	
          }
        }
#### loop through keys and database stuff
        foreach my $ecresult (keys %besthit) {
	    $success2->execute($ecresult);
            my $array_ref = $success2->fetchrow_arrayref();
	    my $ecnum= @$array_ref[0];#print "$ecnum\n";
	    my $ecround = $ecnum;
	    if ($ecnum> $hit_cut) {$ecround=$hit_cut}
	    my $fraction = $hitnumber{$ecresult}/$ecround;
	    my $round = sprintf("%.2f", $fraction);
	    my $check = $conn1->do("DELETE from a8r_blastec WHERE pept_id='$prot' AND ec_id='$ecresult';");
	    my $insert = $conn1->do("INSERT INTO a8r_blastec values ('$prot','$ecresult','$bestdesc{$ecresult}','$besthit{$ecresult}','$bestscore{$ecresult}','$besteval{$ecresult}',$hitnumber{$ecresult},'$ecnum','$round');", {PrintError => 0});
            print FILE "$prot, $ecresult, \"$bestdesc{$ecresult}\", $besthit{$ecresult}, $bestscore{$ecresult}, $besteval{$ecresult},  $hitnumber{$ecresult}, $ecnum, $round\n";            	    
	}
	$success2->finish();      
      }
    }
    if ($conn2) {$conn2->disconnect or warn "Disconnection failed: $DBI::errstr\n";}
    if ($conn1) {$conn1->disconnect or warn "Disconnection failed: $DBI::errstr\n";}
  }


#### kegg stuff
  if ($a_flag == 3 || $a_flag==4) {  

#### get file ready
    chdir ("../output");
    if (-e "$keggout") {
      print colored("\nWARNING: $keggout exists already\n","red bold");
      my $old = $keggout;
      $old =~ s/csv/old/;
      system ("mv $keggout $old");
      print "$keggout moved to $old\n";
    }	
    open(FILE,">$keggout");  
    chdir ("../blast_out");
    
#### db
    my $conn1=DBI->connect("dbi:Pg:dbname=$pg_database", "", "", {PrintError => 0}); #Last two values would be user/pass.
    my $conn2=DBI->connect("dbi:Pg:dbname=$kegg_db", "", "", {PrintError => 0}); #Last two values would be user/pass.  
    print "\nNow parsing KEGG-results. Please wait ...\n";  
    foreach my $file (@kegglist) {
      my $in = new Bio::SearchIO( -format => 'blast',  
        		    -file   => "$file");  
      my $slim=''; 
      my $score='';  
      my $sig='';
      my $db_name=''; 
      my $acc=''; 
      my $prot = ''; 
      my $desc='';
      my $name='';
      my $prog='';
      my $success='';  
      my $n = 0;
      my $m = 0;

#### each query
      while( my $result = $in->next_result )  {
	my %hitnumber =();
        my %besthit =();
	my %besteval=();
	my %bestscore=();
	$n++;
        printf("\r%9d sequences processed so far",$n);  

#### each hit
        my $hit_nr=0;
        while (my $hit = $result->next_hit) {
	  $hit_nr++;
	  $db_name=$result->database_name;
          $prot=$result->query_name;           
          $acc=$hit->accession; 
	  $score=$hit->bits; 
	  $sig=$hit->significance;
          if ($acc =~ /^(\w+)/) {$acc = $1}
	  else {print "\nWARNING, can not deal with $acc\n";}	  
          #### e-value filter here
	  if ($sig > $e_cut) {last;}
	  if ($score < $s_cut) {last;}
	  if ($hit_cut < $hit_nr) {last;}
	  $success=$conn2->prepare("SELECT distinct ko_id from kegg where up_id='$acc'");      
          $success->execute();
          my $array_ref = $success->fetchall_arrayref();		 
          foreach my $output(@$array_ref)  {        
            my ($ko_id) =@$output;
	    $hitnumber{$ko_id}++;
	    unless ($besthit{$ko_id}) {$besthit{$ko_id}=$acc;}
	    unless ($besteval{$ko_id}) {$besteval{$ko_id}=$sig;}
	    unless ($bestscore{$ko_id}) {$bestscore{$ko_id}=$score;}
          }
        }

#### loop through keys, 
        foreach my $keggresult (keys %besthit) {	
 	  my $check = $conn1->do("DELETE from a8r_blastkegg WHERE pept_id='$prot' AND ko_id='$keggresult';");
	  unless ($keggresult =~ /NULL/) {
  	    my $success2=$conn2->prepare("SELECT distinct path_id,ko_id_num from kegg_komap where ko_id='$keggresult'");
	    $success2->execute();
	    my $kegg_ref = $success2->fetchall_arrayref();
	    my $nk = @$kegg_ref;
  	    for(my $i=0;$i<$nk;$i++)	{
              my $refs = @$kegg_ref[$i];
 	      my $path=@$refs[0]; 
	      my $keggnum=@$refs[1];
	      my $kegground = $keggnum;
	      if ($keggnum> $hit_cut) {$kegground=$hit_cut}
	      my $fraction = $hitnumber{$keggresult}/$kegground;
	      my $round = sprintf("%.2f", $fraction);
   	      my $success3=$conn2->prepare("SELECT descr from kegg_def where path_id~'$path'");
	      $success3->execute();
	      my $kegg_ref2 = $success3->fetchrow_arrayref();
	      my $descr;
	      if ($kegg_ref2) {
	  	($descr)=@$kegg_ref2;
	      }
	      else {
	  	$descr = "tba";
	  	print colored("\nWARNING: Now description available for $path\n","red bold");
	      }    
    	      my $insert = $conn1->do("INSERT INTO a8r_blastkegg values ('$prot','$keggresult','$path','$descr','$besthit{$keggresult}','$bestscore{$keggresult}','$besteval{$keggresult}','$hitnumber{$keggresult}','$keggnum','$fraction');", {PrintError => 0});	    
              print FILE "$prot, $keggresult, $path, \"$descr\", $besthit{$keggresult}, $bestscore{$keggresult}, $besteval{$keggresult}, $hitnumber{$keggresult},$keggnum, $fraction\n";
#             print "INSERT INTO a8r_blastkegg values ('$prot','$keggresult','$path','$descr','$besthit{$keggresult}','$besteval{$keggresult}','$hitnumber{$keggresult}');" ;		 
	    }
	  }      
        }
      }
    }
    if ($conn2) {$conn2->disconnect or warn "Disconnection failed: $DBI::errstr\n";}
    if ($conn1) {$conn1->disconnect or warn "Disconnection failed: $DBI::errstr\n";}   
  }
  chdir ("../");
  &query_for_exit();      
}


###################################################################################
###################################################################################
###                                                                             ###
###                               all the subsubs                               ###
###                                                                             ###
###################################################################################
###################################################################################


#########################################################################################################################
sub postmaster_check() {
#### check for postmaster/postgresql process
  my $postmaster=`ps -e|grep postmaster`; ### See if the process is running
#### also catch alternative installation
  unless ($postmaster)  { $postmaster=`ps -e|grep postgres`}; 


  if(!$postmaster)  {
    print colored("\n#### POSTGRESQL IS NOT RUNNING ####\n","red bold");
    print colored("Please ensure that postgreSQL is correctly installed and running\n","red bold");
    exit();
  }
  else {print "CHECK POSTGRESQL   OK => Postmaster running\n\n"} 

#### check whether user does exist   
  my $user_status = system ("psql -l > /dev/null"); ### command will fail unless (postgresql)user exists      
  unless ($user_status == 0) {
    my $username = `whoami`; chomp $username;
    print colored("\n\t#### CONNECTION TO POSTGRESQL FAILED ####\n","red bold");
    print colored("Most likely you have forgotten to run \"createuser $username\"\n","red bold");
    print colored("during the postgreSQL setup\n","red bold");
    exit();
  }
}
###########################################################################################################################


############################################################################################################################
sub query_for_exit() {
#exits program if 'n' entered back to main for y
   print "\nWould you like to continue? ";
   print colored(" [y/n] : ","green bold");
   my $input='';
   while($input!~/y|n/i)   {
     print "\b";
     $input=<>;
     chomp $input;
  }
  if ($input=~/^y/i) { print "Back to main menu\n"; &options;}
  if ($input=~/^n/i) { print "Exiting the program\n"; exit(); }
}
####################################################################################


###########################################################################################################################
sub query_for_continue() {
#exits program if 'n' entered carries on for 'y'
#   print "\nWould you like to continue? ";
   print colored(" [y/n] : ","green bold");
   my $input='';
   while($input!~/y|n/i)   {
     print "\b";
     $input=<>;
     chomp $input;
  }
  if ($input=~/^n/i) { print "Exiting the program\n"; exit(); }
  else {return 1}
}
####################################################################################


####################################################################################
sub yes_no() {
#### returns 1 for y
  my $yflag=0;
  print colored(" [y/n] : ","green bold");
  my $input='';
  while($input!~/y|n/i)   {
    print "\b";
    $input=<STDIN>;
    chomp $input;
  }
  if($input=~/^y/i) { $yflag=1; }
  return $yflag;
}
####################################################################################


####################################################################################
sub delete_file() {
#### returns 1 for y
  my $file = $_[0];
  my $yflag=0;
  print "\nDo you want to replace old file $file with a new version?\n"; 
  print colored(" [y/n] : ","green bold");
  my $input='';
  while($input!~/y|n/i)   {
    print "\b";
    $input=<STDIN>;
    chomp $input;
  }
  if($input=~/^y/i) {
    $yflag=1; 
    unlink ("$file");   
  }
  
  return $yflag;
}
####################################################################################


####################################################################################
sub warn_miss(){
  my $file = $_[0];
  print "$file does not exist in your 'downloads' \n";
  print "sub-directory. Please make sure you have downloaded \n"; 
  print "$file and started annot8r_blast2anno.pl \n";
  print "from the directory one level above 'downloads'.\n\nExiting now.\n\n";
  exit;
}   
####################################################################################


###################################################################################
sub get_file() {    
#### get file from user input, check existence
  my $flag = 0; my $file_name;
  while ($flag == 0) {
    $file_name = <STDIN>;
    chomp $file_name;
    $file_name =~ s/\s+//g;
    if (-e $file_name) {$flag = 1;}
    else {
      print "\nCan't find $file_name, do you want to try again?";
      my $answer = &yes_no();
      if ($answer != 1) {exit();} 
    }  
  }
  return "$file_name";
}  
###################################################################################


######################################################################################
sub filecheck {
  my $file = $_[0];   
  unless (-e $file && (-s $file > 0)) {
    print colored ("Download/unpacking error for $file - try again later\n\n","red bold");
    sleep (2);
    exit
  } 
}
########################################################################################


#######################################################################################
sub find_program()  {
#### search for argument in path, exit if not found
  my $prog=$_[0];
  my $pathflag=0;
  my $path;
  my $finalpath;
  foreach $path (@PATH) {
    if (-f "$path/$prog") {
      $pathflag=1; $finalpath=$path; last;
    }
  }
  if($pathflag==0)   { 
    print colored("\nCan't find the $prog utility on your system\n","red bold");
    print colored("Please ensure it is installed and in your path\n\n","red bold");
    exit();
  }
  else  {
    return "$finalpath/$prog";
  }
}
###############################################################################################################


################################################################################################################
sub run_blast	{
#### run blast, note alignment is needed for some of the info to be extracted downstream
	my($seqs2blast, $dbase_path, $blast_method, $o) = @_;
	
	my @params = ('database' => "$dbase_path" , 'program' => "$blast_method" , 'output' => "$o" , 'b' => '50' , 'v' => '50' ,
		'_READMETHOD' => "Blast");
#print "now here $seqs2blast and here"; system ("pwd"); 
	my $factory = Bio::Tools::Run::StandAloneBlast->new(@params);
	my $blast_report = $factory->blastall($seqs2blast);
	return ($blast_report);
}
#######################################################################################


###################################################################################
###  database subs
####################################################################################
sub create_godb()  { 
  my $database = shift;
  my $new_db = shift;
  my $go_flag = shift;
  my $go_slim_flag = shift;
  my $go_def_flag =shift;  
  my $createdb_exe= &find_program("createdb");
  if($new_db == 1) { system("$createdb_exe $database >& /dev/null"); } 
  my $conn=DBI->connect("dbi:Pg:dbname=$database", "", "");
  if($go_flag==1) {
    my $result=$conn->do("create table go (db text null, dbo_id varchar(12) not null, dbo_sym varchar(12) not null, _not text null, 
    go_id varchar(12) not null, dbref text null, evid text null, w_f text null, asp text null, dbo_name text null,
    dbo_syn text null, dbo_typ text null, taxon text null, date int null, as_by text null);");     
  }
  if ($go_slim_flag==1) {
    my $result=$conn->do("create table go_slim (go_id varchar(12) not null, slim varchar(50) not null);");     
  }  
  if ($go_def_flag==1) {
    my $result=$conn->do("create table go_def (go_id varchar(12) not null, descr text, pcf varchar(1) not null, go_id_num int null);");     
  }
  my $errorMessage = $conn->errstr;
  if ($errorMessage) {print "$errorMessage\n";}
  if ($conn) {$conn->disconnect or warn "Disconnection failed: $DBI::errstr\n";}
}
############################################################################################## 


##############################################################################################
sub create_ecdb()  { 
  my $database=shift;
  my $new_db=shift;
  my $table_flag =shift; 
  my $count_flag =shift; 
  my $createdb_exe= &find_program("createdb");
  if($new_db == 1) { system("$createdb_exe $database >& /dev/null"); } 
  my $conn=DBI->connect("dbi:Pg:dbname=$database", "", "");
  if($table_flag==1) {
    my $result=$conn->do("create table ec (sw_id varchar(14) not null, sw_syn varchar(14) not null,  
    ec_id varchar(20) not null, ec_des text null);");     
  }
  if($count_flag==1) {
    my $result=$conn->do("create table ec_count (ec_id varchar(20) not null, ec_id_num int null);");     
  }  
  my $errorMessage = $conn->errstr;
  if ($errorMessage) {print "$errorMessage\n";}
  if ($conn) {$conn->disconnect or warn "Disconnection failed: $DBI::errstr\n";}
}
############################################################################################## 


##############################################################################################
sub create_keggdb()  { 
  my $database=shift;
  my $new_db=shift;
  my $k_flag=shift;  
  my $m_flag=shift;
  my $d_flag=shift;
  my $createdb_exe= &find_program("createdb");
  if($new_db == 1) { system("$createdb_exe $database >& /dev/null"); } 
  my $conn=DBI->connect("dbi:Pg:dbname=$database", "", "");
  if($k_flag==1) {
    my $result=$conn->do("create table kegg (up_id varchar(16) not null, ko_id varchar(16) null);");     
  }
  if($m_flag==1) {
    my $result=$conn->do("create table kegg_komap (ko_id varchar(16) not null, path_id varchar(14) null, ko_id_num int null);");   
  }  
  if($d_flag==1) {
    my $result=$conn->do("create table kegg_def (path_id varchar(14) not null, descr text null);");     
  }
  my $errorMessage = $conn->errstr;
  if ($errorMessage) {print "$errorMessage\n";}
  if ($conn) {$conn->disconnect or warn "Disconnection failed: $DBI::errstr\n";}
}
############################################################################################## 


##############################################################################################
sub create_pgdb()  { 
  my $database=shift; 
  my $new_db=shift;
  my $go_flag=shift;
  my $ec_flag=shift;
  my $kegg_flag=shift;  
  my $createdb_exe= &find_program("createdb");
  if($new_db == 1) { system("$createdb_exe $database >& /dev/null"); } 
  my $conn=DBI->connect("dbi:Pg:dbname=$database", "", "");
  if($go_flag==1) {
    my $result=$conn->do("create table a8r_blastgo (pept_id varchar(14) not null, go_term varchar(16) not null,  
    pcf varchar(2) null, descr text null, slim varchar(16) null, besthit varchar(14) not null, bestscore float null, bestev float null, hitnum int null, maxhits int null, fraction float null);");     
  }
  if($ec_flag==1) {
    my $result=$conn->do("create table a8r_blastec (pept_id varchar(14) not null, ec_id varchar(16) not null,  
    descr text null, besthit varchar(14) not null, bestscore float null, bestev float null, hitnum int null, maxhits int null, fraction float null);");   
  }

  if($kegg_flag==1) {
    my $result=$conn->do("create table a8r_blastkegg (pept_id varchar(14) not null, ko_id varchar(16) not null,  
    path varchar(16) null, descr text null, besthit varchar(14) not null, bestscore float null, bestev float null, hitnum int null, maxhits int null, fraction float null);");     
  }

  my $errorMessage = $conn->errstr;
  if ($errorMessage) {print "$errorMessage\n";}
  if ($conn) {$conn->disconnect or warn "Disconnection failed: $DBI::errstr\n";}
}
############################################################################################## 

##############
###   fin  ###
##############
