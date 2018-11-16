#!/usr/bin/perl -w
use warnings;
use strict;
use Getopt::Std;

#####################################################################################################
# This script is made to extract sequences from sequence file provided by the user, based on the    #
# cordinates given in the blast file (tabular output). user can provide number of nuceotide to add  #
# on both sides This option is optional                                                             #
#                                                                                                   #
# Author : Ratnesh Singh                                                                            #
# version 1.0                                                                                       #
# contact for bugs: ratnesh@hawaii.edu                                                              #
#####################################################################################################


getopt('sbhtomq');
our ($opt_s,$opt_b,$opt_h,$opt_t,$opt_o,$opt_m,$opt_q, $opt_g);
our (%seq);


# like the shell getopt, "d:" means d takes an argument
print "-sequence file: $opt_s\n" if defined $opt_s;
print "-blast file/sequence: $opt_b\n" if defined $opt_b;
print "-add to head: $opt_h\n" if defined $opt_h;
print "-add to tail: $opt_t\n" if defined $opt_t;
print "-print output as: $opt_o\n" if defined $opt_o;
print "-input mode: $opt_m\n" if defined $opt_m;
print "-using names till first space while searching" if defined $opt_g;
#print "-reverse: $opt_r\n" if defined $opt_r;
print "Unprocessed by Getopt::Std:\n" if $ARGV[0];
foreach (@ARGV) {
  print "-- $_\n";
}

my $help="\n\nThis script will read sequence name and cordinates from 
blast file and extract sequences from sequence file provided

usage:\n use following options\n -m mode of input(auto or manual) [Defaulst is auto]
-s sequence file \n -b Blast out put file in table format 
-h (optional) number of extra nt to add at head [Default 0]
-t (optional) number of extra nt to add at tail [Default 0]
-o output file [Default is: output_seq_extract]
-g use names till first space while searching. [use full length names while]
\n\n\n usage for manual input mode:\n\n-s sequence file name containing all sequences
-b sequence/contig name to extract from
-h (optional) number of extra nt to add at head [Default 0]
-t (optional) number of extra nt to add at tail [Default 0]
-o output file [Default is: output_seq_extract]
-q query name for manual input[Default is: manual_query ]
\n";


die "\nThere is no sequence file specified with -s \n $help" if !defined ($opt_s);

# Check if blast file is provided with 'Auto' option.
if(!defined $opt_m ){
	die "\nThere is no blast file specified with -b \n $help" if !defined ($opt_b);
	print "\n option for -m flag is either missing or improper (not 'auto' or 'manual'). Program assuming default 'auto' option for file input\n"; 
	$opt_m = 'auto';
}
elsif(lc($opt_m) eq 'auto'){die "\nThere is no blast file specified with -b \n $help" if !defined ($opt_b);}
elsif(lc($opt_m) eq 'manual'){$opt_m = 'manual';}
else{print "Have problems with -m option. Check the options\n";}

$opt_h = 0 if !defined ($opt_h);
$opt_t = 0 if !defined ($opt_t);
$opt_o = 'output_seq_extract' if !defined ($opt_o);
#$opt_m = 'auto' if !defined ($opt_m); # moved to above section
$opt_q = 'manual_query' if !defined ($opt_q);

print "saving output in:$opt_o\n";
if($opt_m eq'auto'){open OUT,">$opt_o" or die "cannot open Output file\n";}
else {open OUT,">>$opt_o" or die "cannot open Output file\n";}

print "reading sequence file....Plz wait\n";
open FASTA,"$opt_s" or die "cannot open Sequence file\n";

$/="\n>";
while(<FASTA>){
	chomp;
	
	my($header,@sequence)=split(/\n/,$_);
	  $header=~s/>//;
	  
	  if(defined $opt_g){my @names = split(/\s/,$header); $header=$names[0];}

	my$sequence= join("",@sequence);
	$sequence=~s/\s+//g;
	$sequence=~s/\n+//g;
	
	$seq{$header}{'sequence'}=$sequence;
	$seq{$header}{'len'}= length($sequence);
		
}
close (FASTA);
print".............Done\n";
$/="\n";



###################################################################
#parse information from blast file and extract seq start and end  #
###################################################################
if(lc($opt_m) eq 'auto'){goto "AUTO";}
else{goto "MANUAL";}

close(OUT);
close(BLAST);
exit;


AUTO:{
	
		open BLAST,"$opt_b" or die "cannot read blast file \n";
	
		while(<BLAST>){
			my $line=$_;
			#	print "before line : $line\n";
			if ($line=~/^\s+$/){ next ;};
			if ($line=~/query/i or /match/i or /score/ or /gap/ or /mismatch/){ next; } ;
			#	print "line: $line\n";
			my @line_info= split(/\t/,$line);
			my $query=$line_info[0];
			#	print "query:$query\n";
			my $subject=$line_info[1];
			if(defined $opt_g){my @names = split(/\s/,$line_info[1]); $line_info[1]=$names[0];}

			#	print"subject:$subject\n";
			my $sstart=$line_info[8];
			#	print "Extracting --> $subject: sstart: $sstart\t";
			my $subend=$line_info[9];
			#	print " subend:$subend\n";
	
			#adjustment for increase in length at both ends
	
			my($start_s,$end_s,$strand);
	
			if($sstart > $subend){ $end_s=$sstart + $opt_h - 1; $start_s=$subend - $opt_t - 1;$strand='minus';}
			else{$start_s = $sstart - $opt_h - 1; $end_s = $subend + $opt_t - 1 ;$strand='plus'}
			
			if($start_s<0){$start_s=0;} 
			# added to avoid negative values of start_s.

			my $len = $end_s - $start_s + 1;
			
						
			#print "Extracting --> $subject: sstart: $start_s\t end: $end_s \t length : $len \n";

	
			#my ($new_header,$new_sequence)=extract_seq($query,$subject,$start_s,$end_s,$len,$opt_h,$opt_t,$strand);
			##print OUT">$new_header.$line_info[0].$line_info[1].$line_info[2].$line_info[3].$line_info[4].$line_info[5].$line_info[6].$line_info[7].$line_info[8].$line_info[9] \n$new_sequence\n" if defined $new_sequence;
			#print OUT">$new_header\n$new_sequence\n" if defined $new_sequence;
		
		
	
			if(my ($new_header,$new_sequence,$returnedlength)=extract_seq($query,$subject,$start_s,$end_s,$len,$opt_h,$opt_t,$strand)){
			
			#print OUT">$new_header.$line_info[0].$line_info[1].$line_info[2].$line_info[3].$line_info[4].$line_info[5].$line_info[6].$line_info[7].$line_info[8].$line_info[9] \n$new_sequence\n" if defined $new_sequence;
			print OUT">$new_header\n$new_sequence\n" if defined $new_sequence;
			#print "Header:$new_header\nSequence:$new_sequence\n" if defined $new_sequence;
			print "Extracted --> $subject: sstart: $start_s\t end: $end_s \t Expected length : $len\t Extracted length:$returnedlength \n";
			}
			else{
			
			print "\nError in extracting.\n";
			}
		
		
		
		
		
		
		
		}
	exit;
}

##################################################################################################
# for manual input method.
MANUAL:{	
#	my $subject=$opt_b if defined ($opt_b);
#	if(!defined($opt_b)){print "type the name of contig to look for:";$subject=<STDIN>;	chomp($subject);};
	print "type the name of contig to look for:";my$subject=<STDIN>;	chomp($subject);
	
	print "print start site to cut:"; my $sstart=<STDIN>; 
	if ($sstart<=0){print"Start site cannot be less than 1\n\n"; exit;}
	chomp($sstart);
	print"\nPrint subject end to cut till:";
	my $subend=<STDIN>;
	if ($subend<=0){print"End site cannot be less than 1\n\n"; exit;}
	my($start_s,$end_s);
	
	
	if($sstart > $subend){ $end_s=$sstart + $opt_h -1; $start_s=$subend - $opt_t - 1;}
	else{$start_s = $sstart - $opt_h -1; $end_s = $subend + $opt_t - 1 ; }
	my $len = $end_s - $start_s + 1;
	print "Extracting --> $subject: sstart: $start_s\t end: $end_s \t length : $len \n";
	my $query=$opt_q;
	
	
			
	
	
	####
	if(my ($new_header,$new_sequence,$returnedlength)=extract_seq($query,$subject,$start_s,$end_s,$len,$opt_h,$opt_t)){
			
			#print OUT">$new_header.$line_info[0].$line_info[1].$line_info[2].$line_info[3].$line_info[4].$line_info[5].$line_info[6].$line_info[7].$line_info[8].$line_info[9] \n$new_sequence\n" if defined $new_sequence;
			print OUT">$new_header\n$new_sequence\n" if defined $new_sequence;
			#print "Header:$new_header\nSequence:$new_sequence\n" if defined $new_sequence;
			print "Extracted --> $subject: sstart: $start_s\t end: $end_s \t Expected length : $len\t Extracted length:$returnedlength \n";
			}
			else{
			
			print "\nError in extracting.\n";
			}

	###
	print"type q  to quit or m to extract more:"; my $option=<STDIN>; chomp($option);
	if(lc($option) eq 'm'){ goto "MANUAL";}
	else{last;}
	}



################################################################################
# subroutine for extraction of sequences
################################################################################

sub extract_seq{
	my($query1,$subject1,$start_s1,$end_s1,$len1,$opt1_h,$opt1_t,$strand)=@_;
	
	my $new_sequence1= substr($seq{$subject1}{'sequence'},$start_s1,$len1) if defined $seq{$subject1};
	
	if($strand eq 'minus'){ 
	my$new_sequence2=reverse$new_sequence1;
	$new_sequence2=~tr/atgcATGC/tacgTACG/;
	$new_sequence1=uc$new_sequence2;
	print "seq reversed\n\n";}
	
#	my $length_seq = length($seq{$subject1}) if defined $seq{$subject1};
	my $length_seq = $seq{$subject1}{'len'} if defined $seq{$subject1};
	my$returnedlength=length($new_sequence1);
	print "$subject1 has : $length_seq nt \n" if defined $length_seq;
	my $new_header1=$query1.'_'.$subject1.'-'.'start-'.$start_s1.'_'.'end-'.$end_s1;
	
	return ($new_header1,$new_sequence1,$returnedlength) if defined $new_sequence1;
	print "Cannot find $subject1\n" if !defined $seq{$subject1}{'sequence'}
	}

