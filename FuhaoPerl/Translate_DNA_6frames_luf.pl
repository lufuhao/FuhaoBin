#!/usr/bin/env perl 
use strict;
use warnings;
use Getopt::Long;
use constant USAGE =><<END;

SYNOPSIS:

perl $0 -i input.fa -o outputfile [Options]
Version 20150310

Descriptions:
	Translate DNA/RNA sequences into protein in 6-frame 

Options:
	--help/-h
		Prints help/usage;
	--input/-i <FILE>
		sequences file in fasta format to translate;
	--frame/-f <INT>
		[Opt] Which frame to translate, 0 for all (1-6); 
		Or list (by comma): 1,3,5
		Default: 1
	--min_length/-l <INT>
		[Opt] Minimum number of AAs for one protein
		 selection, Only meaningful when --orfout was specified
		Default [5]
	--clean_stop_codon/-c <STR>
		[Opt] yes|no Clean sequence of stop codon in frame 1.
		Onlymeaningful when --dnaout/-do was specified
		Default: [no]
	--output/-o <FILE>
		Output file full /path/to/file_name;
		Default: Inputbasename.frame[0-6].fa
	--dnaout/-do <File>
		[Opt] DNA/RNA sequence output with codon replace
	--orfout/-ro <File>
		[Opt] Selected ORF output
	--longestorf/-lo <file>
		[Opt] longest ORF
	--verbose
		Detailed output for trouble-shooting;
	--version/-v
		Print current SCRIPT version;

Example:
	perl fasta_stat.pl -i contigs.fa -f 1 -o translated.fa

Author:
	Fu-Hao Lu
	Post-Doctoral Scientist in Micheal Bevan laboratory
	Cell and Developmental Department, John Innes Centre
	Norwich NR4 7UH, United Kingdom
	E-mail: Fu-Hao.Lu\@jic.ac.uk
END
###HELP ends##########################
die USAGE unless (@ARGV);



###Receving parameter#################
our ($help, $input, $output, $dnaout, $orfout, $file_longestORF, $frame, $minlength, $cleanstopcodon, $verbose, $version);

GetOptions(
	"help|h!" => \$help,
	"input|i:s" => \$input,
	"output|o:s" => \$output,
	"dnaout|do:s" => \$dnaout,
	"orfout|ro:s" => \$orfout,
	"longestorf|lo:s" => \$file_longestORF,
	"frame|f:i" => \$frame,
	"min_length|l:i" => \$minlength,
	"clean_stop_codon|c:s" => \$cleanstopcodon,
	"verbose!" => \$verbose,
	"version|v!" => \$version) or die USAGE;


($help or $version or scalar(@ARGV)) and die USAGE;



###Defaults, file and folder checks######
print "\n\n############### Translate starts ################## \n\n" if (defined $verbose);
$frame=0 unless (defined $frame);
our @frame_arr=();
if ($frame==0) {
	@frame_arr=(1, 2, 3, 4, 5, 6);
}
elsif ($frame=~/,/) {
	@frame_arr=split(/,/, $frame);
}
else {
	push (@frame_arr, $frame);
}
foreach (@frame_arr) {
	if ($_<1 or $_>6) {
		die "Please specify the frames in right format. Only numbers 1 to 6 are accepted\n";
	}
}
@frame_arr=sort(@frame_arr);
$frame=join(',', @frame_arr);
$minlength=5 unless (defined $minlength);
$cleanstopcodon='no' unless (defined $cleanstopcodon);
die "Can not find the --input file specified\n" unless (defined $input and (-e $input));
$output=&RetrvNoExt($input).".frame".join('', @frame_arr).".fa" unless (defined $output);
&FileExistsTest($output);
print "Running parameters:\nInput:$input,\nOutput: $output,\nFrame: $frame,\nMinLength: $minlength,\nCleanStopCodon: $cleanstopcodon\n" if (defined $verbose);
our %code=('TCA'=>'S','TCC'=>'S','TCG'=>'S','TCT'=>'S','TTC'=>'F','TTT'=>'F','TTA'=>'L','TTG'=>'L','TAC'=>'Y','TAT'=>'Y','TAA'=>'*','TAG'=>'*','TGC'=>'C','TGT'=>'C','TGA'=>'*','TGG'=>'W','CTA'=>'L','CTC'=>'L','CTG'=>'L','CTT'=>'L','CCA'=>'P','CCC'=>'P','CCG'=>'P','CCT'=>'P','CAC'=>'H','CAT'=>'H','CAA'=>'Q','CAG'=>'Q','CGA'=>'R','CGC'=>'R','CGG'=>'R','CGT'=>'R','ATA'=>'I','ATC'=>'I','ATT'=>'I','ATG'=>'M','ACA'=>'T','ACC'=>'T','ACG'=>'T','ACT'=>'T','AAC'=>'N','AAT'=>'N','AAA'=>'K','AAG'=>'K','AGC'=>'S','AGT'=>'S','AGA'=>'R','AGG'=>'R','GTA'=>'V','GTC'=>'V','GTG'=>'V','GTT'=>'V','GCA'=>'A','GCC'=>'A','GCG'=>'A','GCT'=>'A','GAC'=>'D','GAT'=>'D','GAA'=>'E','GAG'=>'E','GGA'=>'G','GGC'=>'G','GGG'=>'G','GGT'=>'G');



###Translate#############################
open(DNAFILE, $input) || die "Cannot open file $input specified by --input\n";
open (OUTPUT,">>$output") || die "Can not write to $output file\n";
if (defined $dnaout and $dnaout ne '') {
	&FileExistsTest($dnaout);
	open (DNAOUT,">>$dnaout") || die "Can not write to $output file\n";
}
&FileExistsTest($orfout) if (defined $orfout and $orfout ne '');
&FileExistsTest($file_longestORF) if (defined $file_longestORF and $file_longestORF ne '');
our %existedheader=();
$/="\n>";
while(<DNAFILE>){
	(my $DNAheader,my @DANseq)=split(/\n/,$_);
	chomp $DNAheader;
###operate DNA ID
	$DNAheader=~s/\s+$//g;
	$DNAheader=~s/^>//g;
	my $DNAheader_new='';
	my $DNAheader_desc='';
	if ($DNAheader=~/(\S+)\s+(.*)/) {
		$DNAheader_new=$1;
		$DNAheader_desc=$2;
	}
	else {
		$DNAheader_new=$DNAheader;
	}
	if (exists $existedheader{$DNAheader_new}) {
		unlink ($output);
		unlink ($dnaout) if (defined $dnaout and $dnaout ne '');
		unlink ($orfout) if (defined $orfout and $orfout ne '');
		die "The fasta sequence ID is not unique.\n";
	}
	else {
		$existedheader{$DNAheader_new}++;
	} 
###operate DNA sequences
	my $DNAseq = join('', @DANseq);
	$DNAseq =~ s/\s//g;
#	$DNAseq=~s/>//g;
	my $DNA_length=length($DNAseq);
	print "\nSeq:$DNAheader_new\t:$DNA_length nt\n\n" if (defined $verbose);
###reverse DNAseq
	my $DNArevSeq = reverse($DNAseq);
	$DNArevSeq=~tr/ATGCatgc/TACGtacg/;
	print "\nThe original DNA sequence is:\n$DNAseq \nThe reverse of DNA sequence is:\n$DNArevSeq\n" if (defined $verbose);
	
	my @dna=();
	my @protein=();
	my %longestORF=();
	@{$longestORF{$DNAheader_new}}=();
	foreach my $frame_indv (@frame_arr) {
		my $newframe;
		if (defined $dnaout and $dnaout ne '') {
			TEST1: {
				if($frame_indv==1){$newframe='+1'; ($dna[1], $protein[1])=&Frame1all($DNAseq); last TEST1;}
				if($frame_indv==2){$newframe='+2'; ($dna[2], $protein[2])=&Frame2all($DNAseq); last TEST1;}
				if($frame_indv==3){$newframe='+3'; ($dna[3], $protein[3])=&Frame3all($DNAseq); last TEST1;}
				if($frame_indv==4){$newframe='-1'; ($dna[4], $protein[4])=&Frame4all($DNArevSeq); last TEST1;}
				if($frame_indv==5){$newframe='-2'; ($dna[5], $protein[5])=&Frame5all($DNArevSeq); last TEST1;}
				if($frame_indv==6){$newframe='-3'; ($dna[6], $protein[6])=&Frame6all($DNArevSeq); last TEST1;}
			}
			print DNAOUT ">$DNAheader_new".'_'."$newframe\n$dna[$frame_indv]\n";
			print ">$DNAheader_new".'_'."$newframe\n$dna[$frame_indv]\n" if (defined $verbose);
		}
		else {
			TEST2: {
				if($frame_indv==1){$newframe='+1'; $protein[1]=&Frame1aa($DNAseq); last TEST2;}
				if($frame_indv==2){$newframe='+2'; $protein[2]=&Frame2aa($DNAseq); last TEST2;}
				if($frame_indv==3){$newframe='+3'; $protein[3]=&Frame3aa($DNAseq); last TEST2;}
				if($frame_indv==4){$newframe='-1'; $protein[4]=&Frame4aa($DNArevSeq); last TEST2;}
				if($frame_indv==5){$newframe='-2'; $protein[5]=&Frame5aa($DNArevSeq); last TEST2;}
				if($frame_indv==6){$newframe='-3'; $protein[6]=&Frame6aa($DNArevSeq); last TEST2;}
			}
		}
		print OUTPUT ">$DNAheader_new".'_'."$newframe\n$protein[$frame_indv]\n";
		print ">$DNAheader_new".'_'."$newframe\n$protein[$frame_indv]\n" if (defined $verbose);
		@{$longestORF{$DNAheader_new}}=(@{$longestORF{$DNAheader_new}}, split(/\s|[.*]/, $protein[$frame_indv]));
		if (defined $orfout and $orfout ne '') {
			&SelectORF($DNAheader_new, $protein[$frame_indv], $newframe);
		}
	}
	if (defined $file_longestORF and $file_longestORF ne ''){
		&SelectLongestORF($DNAheader_new, @{$longestORF{$DNAheader_new}});
	}
}
close DNAOUT;
close OUTPUT;
close DNAFILE;



###############################################
###sub functions###############################
###############################################

### Backup file
sub backup {
	my $BUori_name=shift @_;
	my $BU_new_name='';
	$BU_new_name=$BUori_name.".bak.".int(rand(10000));
	if (-e "$BU_new_name") {
		&backup($BUori_name);
	}
	rename($BUori_name, $BU_new_name);
	print "The Original file $BUori_name was renamed to $BU_new_name\n";
}

###Return file path
sub RetrvDir {
	my $RD_ori=shift @_;
	chomp $RD_ori;
	my $RD_new='';
	if ($RD_ori=~/\//) {
		($RD_new=$RD_ori)=~ s/(.*)\/.*$/$1/s;
	}
	else {
		$RD_new='.';
	}
	return ($RD_new);
	
	
}

###return file base name
sub RetrvBase {
	my $RB_ori=shift @_;
	chomp $RB_ori;
	my $RB_new='';
	($RB_new=$RB_ori)=~ s/.*\///s; 
	return $RB_new;
}

###Return file name without extension
sub RetrvNoExt {
	my $RNE_ori=shift @_;
	chomp $RNE_ori;
	my ($RNE_new, $RNE_base)=('', '');
	($RNE_base=$RNE_ori)=~ s/.*\///s; 
	($RNE_new=$RNE_base)=~s/^(.*)\.\w+$/$1/;
	return $RNE_new;
}

###determine how to process those files existed
sub FileExistsTest {
	my $FET_fileid=shift @_;
	chomp $FET_fileid;
	if (-e $FET_fileid) {
		if ($verbose) {
			print "The $FET_fileid file specified already existed\nNeed to overwrite [y] or backup [n] or others [AnyKey], default[y]:\n";
			my $FET_test=''; $FET_test=<STDIN>;
			chomp $FET_test;
			if (lc ($FET_test) eq ('y'||'yes')) {
				unlink($FET_fileid);
			}
			elsif (lc ($FET_test) eq ('n'||'no')) {
				&backup($FET_fileid);
			}
			else {
				die "Please specify a new name for output\n";
			}
		}
		else {
			unlink($FET_fileid);
#			&backup($FET_fileid);
		}
	}
}



sub Frame1all {
	my $origin_seq = shift @_;
	my $protein='';
	my $dna='';
	my $codon1='';
	for(my $i=0;$i<(length($origin_seq)-2);$i+=3){
		$codon1=substr($origin_seq,$i,3);
		$protein.= &codon2aa($codon1);
		$dna.=&codon2nt($codon1);
	}
	return ($dna, $protein);
}



sub Frame2all {
	my $origin_seq = shift @_;
	my $protein='';
	my $dna='';
	my $codon2='';
	for(my $i=1;$i<(length($origin_seq)-2);$i+=3){
		$codon2=substr($origin_seq,$i,3);
		$protein.= &codon2aa($codon2);
		$dna.=&codon2nt($codon2);
	}
	return ($dna, $protein);
}



sub Frame3all {
	my $origin_seq = shift @_;
	my $protein='';
	my $dna='';
	my $codon3='';
	for(my $i=2;$i<(length($origin_seq)-2);$i+=3){
		$codon3=substr($origin_seq,$i,3);
		$protein.= &codon2aa($codon3);
		$dna.=&codon2nt($codon3);
	}
	return ($dna, $protein);
}



sub Frame4all {
	my $origin_seq = shift @_;
	my $protein='';
	my $dna='';
	my $codon4='';
	for(my $i=0;$i<(length($origin_seq)-2);$i+=3){
		$codon4=substr($origin_seq,$i,3);
		$protein.= &codon2aa($codon4);
		$dna.=&codon2nt($codon4);
	}
	return ($dna, $protein);
}



sub Frame5all {
	my $origin_seq = shift @_;
	my $protein='';
	my $dna='';
	my $codon5='';
	for(my $i=1;$i<(length($origin_seq)-2);$i+=3){
		$codon5=substr($origin_seq,$i,3);
		$protein.= &codon2aa($codon5);
		$dna.=&codon2nt($codon5);
	}
	return ($dna, $protein);
}



sub Frame6all {
	my $origin_seq = shift @_;
	my $protein='';
	my $dna='';
	my $codon6='';
	for(my $i=2;$i<(length($origin_seq)-2);$i+=3){
		$codon6=substr($origin_seq,$i,3);
		$protein.= &codon2aa($codon6);
		$dna.=&codon2nt($codon6);
	}
	return ($dna, $protein);
}



sub Frame1aa {
	my $origin_seq = shift @_;
	my $protein='';
	my $dna='';
	my $codon1='';
	for(my $i=0;$i<(length($origin_seq)-2);$i+=3){
		$codon1=substr($origin_seq,$i,3);
		$protein.= &codon2aa($codon1);
	}
	return $protein;
}



sub Frame2aa {
	my $origin_seq = shift @_;
	my $protein='';
	my $dna='';
	my $codon2='';
	for(my $i=1;$i<(length($origin_seq)-2);$i+=3){
		$codon2=substr($origin_seq,$i,3);
		$protein.= &codon2aa($codon2);
	}
	return $protein;
}



sub Frame3aa {
	my $origin_seq = shift @_;
	my $protein='';
	my $dna='';
	my $codon3='';
	for(my $i=2;$i<(length($origin_seq)-2);$i+=3){
		$codon3=substr($origin_seq,$i,3);
		$protein.= &codon2aa($codon3);
	}
	return $protein;
}



sub Frame4aa {
	my $origin_seq = shift @_;
	my $protein='';
	my $dna='';
	my $codon4='';
	for(my $i=0;$i<(length($origin_seq)-2);$i+=3){
		$codon4=substr($origin_seq,$i,3);
		$protein.= &codon2aa($codon4);
	}
	return $protein;
}



sub Frame5aa {
	my $origin_seq = shift @_;
	my $protein='';
	my $dna='';
	my $codon5='';
	for(my $i=1;$i<(length($origin_seq)-2);$i+=3){
		$codon5=substr($origin_seq,$i,3);
		$protein.= &codon2aa($codon5);
	}
	return $protein;
}



sub Frame6aa {
	my $origin_seq = shift @_;
	my $protein='';
	my $dna='';
	my $codon6='';
	for(my $i=2;$i<(length($origin_seq)-2);$i+=3){
		$codon6=substr($origin_seq,$i,3);
		$protein.= &codon2aa($codon6);
	}
	return $protein;
}



sub codon2aa{
	my $codon=shift @_;
	$codon=uc($codon);
	IM: {
		if(exists $code{$codon}){
			return $code{$codon};
			last IM;
		}
		elsif($codon=~/GC./i){return 'A';last IM;}
		elsif($codon=~/GG./i){return 'G';last IM;}
		elsif($codon=~/CC./i){return 'P';last IM;}
		elsif($codon=~/AC./i){return 'T';last IM;}
		elsif($codon=~/GT./i){return 'V';last IM;}
		elsif($codon=~/CG./i){return 'R';last IM;}
		elsif($codon=~/TC./i){return 'S';last IM;}
		else{
			return('x');
			print "Bad codon \"$codon\"!!\n";
		}
	}
}



sub codon2nt{
	my $codon=shift @_;
	$codon=uc($codon);
	if(lc($cleanstopcodon) eq 'yes'){
		if($code{$codon} ne '*'){
			return $codon;
		}
		else{
			return('---');
		}
	}
	else{
		return $codon;
	}
}



# Selecting ORFs from translated protein
sub SelectORF {
	my ($DNAheader, $protein, $newframe)=@_;
	my (@protein_ORF)=split(/[\*]+/,$protein);
	my $protein_ORF=scalar(@protein_ORF);
	my $k=0;
	open (ORFOUT, ">>$orfout") || die "Can not write to file ORFout\n";
	for(my$j=0;$j<=$protein_ORF;$j++){
		next if !defined $protein_ORF[$j];
		$protein_ORF[$j]=~s/\s*//g;
		if ($protein_ORF[$j]=~/M/){$protein_ORF[$j]=~s/(\w*?M)/M/;} 
		else{next;}
		if(length($protein_ORF[$j])<$minlength){next;}
		$k++;
		my $prot_length=length($protein_ORF[$j]);
		print ORFOUT ">$DNAheader".'_'."$newframe\.$k\t$prot_length aa\n$protein_ORF[$j]\n" if defined $protein_ORF[$j];
		print ">$DNAheader".'_'."$newframe\.$k\t$prot_length aa\n$protein_ORF[$j]\n" if defined $protein_ORF[$j];
		}
	close ORFOUT;
}



###Select longest ORF for each DNAseq
#&SelectLongestORF($DNAheader_new, @{$longestORF{$DNAheader_new}});
sub SelectLongestORF {
	my ($SLOid, @SLOpro_arr)=@_;
	my $SLOmax_id=0;
	my $SLOlongest_length=0;
	for (my $i=0; $i<scalar(@SLOpro_arr);$i++) {
#		next unless ($SLOpro_arr[$i]=~/^M/i);
		if (length($SLOpro_arr[$i])>$SLOlongest_length) {
			$SLOlongest_length=length($SLOpro_arr[$i]);
			$SLOmax_id=$i;
		}
	}
	open (LONGESTORF, ">>$file_longestORF") || die "Error: Can not output longest ORF: $file_longestORF\n";
	print LONGESTORF '>'.$SLOid."\n".$SLOpro_arr[$SLOmax_id]."\n";
	close LONGESTORF;
}
