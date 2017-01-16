#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long;
use constant USAGE=><<EOH;

SYNOPSIS:

perl $0 --flashfile merge.fq -f R1.fq.gz -r R2.fq.gz -1 outR1.fa.gz -2 out.R2.fq.gz [Options]
Version: LUFUHAO20150224

Requirements:
	Programs: gzip, zcat
	Modiles: Getopt::Long

Descriptions:
	Trim overlapped parts for paired reads

Options:
	--help|-h
		Print this help/usage;
	--flashfile|-m	<fastq>
		[Msg] Merge reads from Paired reads by FLASH or fastq-join;
	--forward|-f	<fastqR1>
		[Msg] Original forward reads used to Merge, in fq or fq.gz format;
	--reverse|-r	<fastqR2>
		[Msg] Original reverse reads used to Merge, in fq or fq.gz format;
	--seedlength|-l	<INT>
		[Opt] Substr length to check if outward merge, default: 12;
	--min_seedlen|-x	<INT>
		[Opt] Minimum substr length to check if outward merge, default: 8;
	--outR1|-1	<OutR1>
		[Opt] Output forward reads for in fq or fq.gz format;
	--outR2|-2	<OutR1>
		[Opt] Output reverse reads for in fq or fq.gz format;
	--debug
		[Opt] Output more info for debugging;
	--verbose
		[Opt] Detailed output for trouble-shooting;
	--version|v!
		[Opt] Print current SCRIPT version;

Example:
	flash --min-overlap=20 --read-len=250 --fragment_len 600 --fragment_len_stddev 150 --output-prefix=PREFIX --compress --threads=10 ForwardR1 ForwardR2
		#PREFIX.extendedFrags.fastq.gzip
		#PREFIX.notCombined_1.fastq.gz
		#PREFIX.notCombined_2.fastq.gz
	zcat PREFIX.extendedFrags.fastq.gz | perl -ne 'chomp; if (\$_=~/^\@(\\S+)\\s*/) {\$ID=\$1;} else {die "IDmatcherror\\n";} \$seq=<>; chomp \$seq; \$length=length(\$seq); <>; <>; print \$ID."\\t".\$length."\\n";' > PREFIX.extendedFrags.fastq.length
	perl /usr/users/celldev/luf/lufuhao/SizeCollectBin_luf.pl PREFIX.extendedFrags.fastq.length 20 > PREFIX.extendedFrags.fasta.SizeBin
	cat PREFIX.extendedFrags.fasta.length | cut -f 1 > PREFIX.extendedFrags.fastq.ID
	seqtk subseq ForwardR1 PREFIX.extendedFrags.fastq.ID > PREFIX.seqtk.R1.fq
	seqtk subseq ForwardR2 PREFIX.extendedFrags.fastq.ID > PREFIX.seqtk.R2.fq
	perl $0 -m PREFIX.extendedFrags.fastq.gz -f PREFIX.seqtk.R1.fq/ForwardR1 -r PREFIX.seqtk.R2.fq/ForwardR2 -l 12 -1 Out.R1.fa.gz -2 Out.R2.fq.gz



Author:
	Fu-Hao Lu
	Post-Doctoral Scientist in Micheal Bevan laboratory
	Cell and Developmental Department, John Innes Centre
	Norwich NR4 7UH, United Kingdom
	E-mail: Fu-Hao.Lu\@jic.ac.uk

EOH
###HELP ends#########################################################
die USAGE unless @ARGV;



###Receving parameter################################################
our ($help, $verbose, $debug, $ver);
our ($flashfile, $forward, $reverse);
our ($substr_length, $min_substr_length, $readlengthF, $readlengthR);
our ($output_R1, $output_R2);
GetOptions(
	"help|h!" => \$help,
	"flashfile|m:s" => \$flashfile,
	"forward|f:s" => \$forward,
	"reverse|r:s" => \$reverse,
	"lengthF:s" => \$readlengthF,
	"lengthR:s" => \$readlengthR,
	"seedlength|l:s" => \$substr_length,
	"min_seedlen|x:s" => \$min_substr_length,
	"outR1|1:s" => \$output_R1,
	"outR2|2:s" => \$output_R2,
	"debug!" => \$debug,
	"verbose!" => \$verbose,
	"version|v!" => \$ver) or die USAGE;
($help or $ver) and die USAGE;



### Defaults ########################################################
$substr_length=12 unless (defined $substr_length);
$min_substr_length=8 unless (defined $min_substr_length);
$debug=0 unless (defined $debug);



### input and output ################################################
die "Invalid --flashfile file\n" unless (defined $flashfile and -s $flashfile);
die "Invalid --forward file\n" unless (defined $forward and -s $forward);
die "Invalid --reverse file\n" unless (defined $reverse and -s $reverse);
#Delete Output if already exist
unlink ("$output_R1") if (-s $output_R1);
unlink ("$output_R2") if (-s $output_R2);



### Main ############################################################
###Step1. read first $substr_length and length for each overlap
if ($flashfile=~m/(\.gz$)|(\.gzip$)/) {
	print "FLASH file in gz format\n";
	open (FLASH, "zcat $flashfile |") || die "Step1: Can not open --flashfile $flashfile file\n";
}
elsif ($flashfile=~m/(\.fq$)|(\.fastq$)/) {
	print "FLASH file in fastq format\n";
	open (FLASH, "cat $flashfile |") || die "Step1: Can not open --flashfile $flashfile file\n";
}
else {
	die "Step1: Can not guess format: --flashfile $flashfile\n";
}
our %subtring=();
our %flash_seq_length=();
my $num_seqID=0;
while (my $line1=<FLASH>) {
	chomp $line1;
	my $seqid='';
	if ($line1=~/^\@(\S+)\s*/) {
		$seqid=$1;
		die "Step1: invalid seqID: $line1\n" if ($seqid eq '');
	}
	else {
		print STDERR "Step1 Warning: Ignoring this line as not seqID: $line1\n";
		next;
	}
	my $seq=<FLASH>;
	chomp $seq;
	$subtring{$seqid}=substr($seq, 0, $substr_length);
	if (length ($seq) >0) {
		$flash_seq_length{$seqid}=length ($seq);
		$num_seqID++;
	}
	else {
		print "Step1: warning: empty sequence for merge ID: $seqid\n";
	}
	<FLASH>;
	<FLASH>;
}
close FLASH;
print "Step1: total number of IDs from flash: $num_seqID\n\n\n";



###Step2. read R1 length and check first $substr_length strings to confirm overlap pattern
###Pattern wanted
#                R1 ================>
#                           <==================R2
###Pattern unwanted
#                            R1 ================>
#                  <==================R2
our %lengthF=();
if ($forward=~m/(\.gz$)|(\.gzip$)/) {
	print "Forward file in gz format\n";
	open (FORWARD, "zcat $forward |") || die "Step2: Can not open --forward $forward file\n";
}
elsif ($forward=~m/(\.fq$)|(\.fastq$)/) {
	print "Forward file in fastq format\n";
	open (FORWARD, "cat $forward |") || die "Step2: Can not open --forward $forward file\n";
}
else {
	die "Step2: Can not guess format: --forward $forward\n";
}
$num_seqID=0;
my $abandon_num_seqID=0;
while (my $line2=<FORWARD>) {
	chomp $line2;
	my $seqid='';
	if ($line2=~/^\@(\S+)\s*/) {
		$seqid=$1;
		die "Step2: invalid seqID: $line2\n" if ($seqid eq '');
	}
	else {
		print STDERR "Step2 Warning: Ignoring this line as not seqID: $line2\n";
		next;
	}
	if (exists $flash_seq_length{$seqid}) {
		my $seq=<FORWARD>;
		chomp $seq;
		my $R1_substr=substr($seq, 0, $substr_length);
		if ($R1_substr eq $subtring{$seqid}) {
			if (length ($seq) >0) {
				$lengthF{$seqid}=(defined $readlengthF and $readlengthF>0) ? $readlengthF : length ($seq);
				$num_seqID++;
			}
			else {
				print "Step2: empty sequence for Forward ID: $seqid\n";
				delete $subtring{$seqid};
				delete $flash_seq_length{$seqid};
				$abandon_num_seqID++;
			}
		}
		else {
			delete $subtring{$seqid};
			delete $flash_seq_length{$seqid};
			$abandon_num_seqID++;
		}
		<FORWARD>;
		<FORWARD>;
	}
	else {
		<FORWARD>;
		<FORWARD>;
		<FORWARD>;
	}
}
close FORWARD;
print "Step2: read number of forward ID: $num_seqID\n";
print "Step2: Number of read ids in hash: ".scalar(keys %subtring)."\n";
%subtring=();###Empty useless hash 
print "Step2: abandon number of forward ID: $abandon_num_seqID\n\n\n";



###Step3: calculate trim length from 3'end
our %lengthTrim=();
$num_seqID=0;
if (defined $readlengthR and $readlengthR>0) {
	foreach my $mergeID (keys %flash_seq_length) {
		my $lengthR1=(defined $readlengthF and $readlengthF>0) ? $readlengthF : $lengthF{$mergeID};
		my $trim_length=$lengthR1+$readlengthR-$flash_seq_length{$mergeID};
		if ($trim_length > 0) {
			$lengthTrim{$mergeID}=$trim_length;
			$num_seqID++;
		}
		else {
			print STDERR "Step3: non-positive trim length $trim_length ($lengthR1+$readlengthR-$flash_seq_length{$mergeID}) for $mergeID\n";
		}
	}
	print "Step3: read number of reverse ID: $num_seqID\n";
	print "Step3: read number of hash for trim: ".scalar(keys %lengthTrim)."\n";
	goto STEP4;
}
STEP3: {
if ($reverse=~m/(\.gz$)|(\.gzip$)/) {
	print "Reverse file in gz format\n";
	open (REVERSE, "zcat $reverse |") || die "Step3: Can not open --reverse $reverse file\n";
}
elsif ($forward=~m/(\.fq$)|(\.fastq$)/) {
	print "Reverse file in fastq format\n";
	open (REVERSE, "cat $reverse |") || die "Step3: Can not open --reverse $reverse file\n";
}
else {
	die "Step3: Can not guess format: --reverse $reverse\n";
}

my $R2_num_exclude=0;
while (my $line3=<REVERSE>) {
	chomp $line3;
	my $seqid='';
	if ($line3=~/^\@(\S+)\s*/) {
		$seqid=$1;
		die "invalid seqID: $line3\n" if ($seqid eq '');
	}
	else {
		print STDERR "Step3 Warning: Ignoring this line as not seqID: $line3\n";
		next;
	}
	if (exists $flash_seq_length{$seqid}) {
		my $seq=<REVERSE>;
		chomp $seq;
		if (exists $flash_seq_length{$seqid} and exists $lengthF{$seqid}) {
			if ($flash_seq_length{$seqid} > 0 and $lengthF{$seqid} > 0) {
				my $R2_length=length($seq);
				if ($R2_length >0) {
					my $trimlength=$lengthF{$seqid}+$R2_length-$flash_seq_length{$seqid};
#					print "Step3: trim length: $trimlength\n";
					if ($trimlength >0) {
						$lengthTrim{$seqid}=$trimlength;
						$num_seqID++;
					}
					else {
						print STDERR "Step3: minus or zero value for trim\n";
					}
				}
				else {
					print STDERR "Step3: minus or zero value for read length of flash or R1\n";
					$R2_num_exclude++;
				}
			}
		}
		else {
			print STDERR "Step3: $seqid flash or R1 length error";
			$R2_num_exclude++;
		}
		<REVERSE>;
		<REVERSE>;
	}
	else {
		<REVERSE>;
		<REVERSE>;
		<REVERSE>;
	}
}
close REVERSE;
print "Step3: read number of reverse ID: $num_seqID\n";
print "Step3: read number of hash for trim: ".scalar(keys %lengthTrim)."\n";
print "Step3: number of excluded ID: $R2_num_exclude\n\n\n";
}###STEP3




###Step4: Trim R1 reads
STEP4: {
###Empty useless hash;
%flash_seq_length=();
%lengthF=();
if ($forward=~m/(\.gz$)|(\.gzip$)/) {
	print "Forward file in gz format\n" if ($debug);
	open (FORWARD, "zcat $forward |") || die "Step4: Can not open --forward $forward file\n";
}
elsif ($forward=~m/(\.fq$)|(\.fastq$)/) {
	print "Forward file in fastq format\n" if ($debug);
	open (FORWARD, "cat $forward |") || die "Step4: Can not open --forward $forward file\n";
}
else {
	die "Step4: Can not guess format: --forward $forward\n";
}
if ($output_R1=~m/(\.gz$)|(\.gzip$)/) {
	print "Output of forward reads in gz format\n";
	open (OUTR1, " | gzip -9 > $output_R1") || die "Step4: Can not write --outR1 $output_R1 file\n";
}
elsif ($output_R1=~m/(\.fq$)|(\.fastq$)/) {
	print "Output of forward reads in fastq format\n";
	open (OUTR1, " > $output_R1") || die "Step4: Can not write --outR1 $output_R1 file\n";
}
else {
	die "Step4: Can not guess format: --outR1 $output_R1\n";
}
$num_seqID=0;
while (my $line4=<FORWARD>) {
	chomp $line4;
	my $seqid='';
	if ($line4=~/^\@(\S+)\s*/) {
		$seqid=$1;
		die "Step4: invalid seqID: $line4\n" if ($seqid eq '');
	}
	else {
		print STDERR "Step4 Warning: Ignoring this line as not seqID: $line4\n";
		next;
	}
	if (exists $lengthTrim{$seqid}) {
		my $seq=<FORWARD>;
		chomp $seq;
		if (exists $lengthTrim{$seqid} and $lengthTrim{$seqid} >0) {
			my $substr_length=length($seq)-$lengthTrim{$seqid};
			if ($substr_length > 0) {
				my $trimmed_R1_seq=substr($seq, 0, $substr_length);
				my $R1_qual_header=<FORWARD>;
				chomp $R1_qual_header;
				my $R1_qual=<FORWARD>;
				chomp $R1_qual;
				my $trimmed_R1_qual=substr($R1_qual, 0, $substr_length);
				print OUTR1 "$line4\n$trimmed_R1_seq\n$R1_qual_header\n$trimmed_R1_qual\n";
				$num_seqID++;
			}
			elsif ($substr_length ==0) {
				print STDERR "STEP4: trim full length: $seqid\n";
				<FORWARD>;
				<FORWARD>;
			}
			else {
				print STDERR "Step4: trim length ($lengthTrim{$seqid}) >= read length (".length($seq).") for R2: seqID: $seqid\n";
				<FORWARD>;
				<FORWARD>;
			}
		}
		else {
			print STDERR "Step4: Output R1 seqid not exist in Trim hash: $seqid\n";
			<FORWARD>;
			<FORWARD>;
		}
	}
	else {
		print "Step4: Ignoring... dur to \n";
		<FORWARD>;
		<FORWARD>;
		<FORWARD>;
	}
}
close FORWARD;
close OUTR1;
print "Step 4: output number of forward IDs: $num_seqID\n\n\n";
goto STEP5;
}###STEP4



###STEP5: 
STEP5: {
if ($reverse=~m/(\.gz$)|(\.gzip$)/) {
	print "Reverse file in gz format\n" if ($debug);
	open (REVERSE, "zcat $reverse |") || die "Step5: Can not open --reverse $reverse file\n";
}
elsif ($forward=~m/(\.fq$)|(\.fastq$)/) {
	print "Reverse file in fastq format\n" if ($debug);
	open (REVERSE, "cat $reverse |") || die "Step5: Can not open --reverse $reverse file\n";
}
else {
	die "Step5: Can not guess format: --reverse $reverse\n";
}
if ($output_R2=~m/(\.gz$)|(\.gzip$)/) {
	print "Output of reverse reads in gz format\n";
	open (OUTR2, " | gzip -9 > $output_R2") || die "Step5: Can not write --outR2 $output_R2 file\n";
}
elsif ($output_R2=~m/(\.fq$)|(\.fastq$)/) {
	print "Output of reverse reads in fastq format\n";
	open (OUTR2, " > $output_R1") || die "Step5: Can not open --outR2 $output_R2 file\n";
}
else {
	die "Step5: Can not guess format: --outR2 $output_R2\n";
}
$num_seqID=0;
while (my $line5=<REVERSE>) {
	chomp $line5;
	my $seqid='';
	if ($line5=~/^\@(\S+)\s*/) {
		$seqid=$1;
		die "invalid seqID: $line5\n" if ($seqid eq '');
	}
	else {
		print STDERR "Step5 Warning: Ignoring this line as not seqID: $line5\n";
		next;
	}
	if (exists $lengthTrim{$seqid}) {
		my $seq=<REVERSE>;
		chomp $seq;
		if ($lengthTrim{$seqid} >0) {
			my $substr_length=length($seq)-$lengthTrim{$seqid};
			if ($substr_length > 0) {
				my $trimmed_R2_seq=substr($seq, 0, $substr_length);
				my $R2_qual_header=<REVERSE>;
				chomp $R2_qual_header;
				my $R2_qual=<REVERSE>;
				chomp $R2_qual;
				my $trimmed_R2_qual=substr($R2_qual, 0, $substr_length);
				print OUTR2 "$line5\n$trimmed_R2_seq\n$R2_qual_header\n$trimmed_R2_qual\n";
				$num_seqID++;
			}
			elsif ($substr_length ==0) {
				print STDERR "Step5: trim full length: seqID: $seqid\n";
				<REVERSE>;
				<REVERSE>;
			}
			else {
				print STDERR "Step5: trim length ($lengthTrim{$seqid}) >= read length (".length($seq).") for R2: seqID: $seqid\n";
				<REVERSE>;
				<REVERSE>;
			}
		}
		else {
			print STDERR "Step5: Output R2 sedid not exist in Trim hash: $seqid\n";
			<REVERSE>;
			<REVERSE>;
		}
	}
	else {
		<REVERSE>;
		<REVERSE>;
		<REVERSE>;
	}
}
close REVERSE;
close OUTR2;
print "Step5: output number of reverse IDs: $num_seqID\n\n\n";
}###STEP5
#####################################################################
###                         sub functions                         ###
#####################################################################
### ReadSam
###&ReadSam(sam,ref, 1/2/3)
###Global:
###Dependency:
###Note:

