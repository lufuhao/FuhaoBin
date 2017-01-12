#!/usr/bin/env perl
use strict;
use warnings;
use constant USAGE =><<EOH;

usage: $0 input.psl[.gz] query.fasta.fai out.psl[.gz]

REQUIREMENTS:

  gzip    if input ot output PSL in gz format

DESCRIPTIONS:

  This script is used to convert the choped sequence parts to its original coordinates
    Using samtools faidx to extract the sequence:
        eg. index fasta: samtools faidx your.fasta
        eg. extract 1-1000: samtools faidx your.fasta seqID:1-1000 > seqID.1-1000.fasta
        eg. extract 1001-2000: samtools faidx your.fasta seqID:1001-2000 > seqID.1001-2000.fasta

    Run Lastal for each sequence parts, and convert maf to PSL
        Note: PSL column 10 much be in this format:
              origin_query_seqID:start-end
    This script will correct the coordinates
        Each resulting PSL used for ChainNet

INPUT AND OUTPUT

	input.psl: output of lastal_in_parts.sh

    query.fasta.fai: file containing the full length for each queryID
        Note: fasta index should work
    eg:
        #seqID[tab]length
        seqID1	100000000
        seqID2	150000
        ..

v20170112

EOH
die USAGE if (scalar(@ARGV) !=3 or $ARGV[0] eq '-h' or $ARGV[0] eq '--help');


my ($inputpsl, $lengthfile, $outputpsl)=@ARGV;

die "Error: can not find input PSL\n" unless (defined $inputpsl and -s $inputpsl);
die "Error: can not find file for query.fasta.fai\n" unless (defined $lengthfile and -s $lengthfile);
die "Error: invalid output PSL\n" unless (defined $outputpsl and $outputpsl=~/^[-\w\.\/]+$/);
unlink $outputpsl if (-e $outputpsl);



my $linenum=0;
my %seqhash=();


open (SEQLENGTH, " < $lengthfile") || die "Error: can not open $lengthfile\n";
while (my $line=<SEQLENGTH>) {
	$linenum++;
	chomp $line;
	next if ($line=~/^#/);
	my @arr=split (/\t/, $line);
	unless (scalar(@arr)>=2 and defined $arr[0] and $arr[0]=~/^\S+$/ and defined $arr[1] and $arr[1]=~/^\d+$/) {
		die "Error: Invalid line$linenum: $line\n";
	}
	if (exists $seqhash{$arr[0]}) {
		if ($seqhash{$arr[0]}==$arr[1]) {
			print STDERR "Warnings: repeated seqID $arr[0] having the same length at line$linenum\n";
		}
		else {
			die "Error: repeated seqID $arr[0] having different length at line$linenum\n";
		}
	}
	else {
		$seqhash{$arr[0]}=$arr[1];
	}
}
close SEQLENGTH;


print "### Read sequence length:\n
    Lines:        $linenum\n
    Sequences:    ", scalar(keys %seqhash), "\n\n";

if ($inputpsl=~/\.psl\.gz$/i) {
	print "Info: input psl in GZ format\n";
	open (INPUTPSL, " gzip -dc $inputpsl | ") || die "Error: can not open input PSL: $inputpsl\n";
}
elsif ($inputpsl=~/\.psl$/i) {
	print "Info: input psl in text format\n";
	open (INPUTPSL, " < $inputpsl ") || die "Error: can not open input PSL: $inputpsl\n";
}
else {
	die "Info: can not guess input PSL format, must end with .psl or .psl.gz\n";
}
if ($outputpsl=~/\.psl\.gz$/i) {
	print "Info: output psl in GZ format\n";
	open (OUTPUTPSL, " | gzip -9 -c > $outputpsl") || die "Error: can not write output PSL: $outputpsl\n";
}
elsif ($inputpsl=~/\.psl$/i) {
	print "Info: output psl in text format\n";
	open (OUTPUTPSL, " > $outputpsl ") || die "Error: can not write output PSL: $outputpsl\n";
}
else {
	die "Error: can not guess output PSL format, must end with .psl or .psl.gz\n";
}
$linenum=0;
my $addnum='NaN';
my $reverseadd='NaN';
my $seqlength='NaN';
my $seqid="NANUNKNOWNSEQ";
while (my $line=<INPUTPSL>) {
	$linenum++;
	chomp $line;
	if ($line=~/^#/) {
		print OUTPUTPSL $line, "\n";
		next;
	}
	my @arr=split(/\t/, $line);
	unless (scalar(@arr)==21) {### checking 21 columns
		die "Error: non 21-columns in input PSL ar line$linenum: $line\n";
	}

	if ($arr[9]=~/^(\S+):(\d+)-(\d+)$/) {
		my $testprint=0;
		$testprint=1 unless ($addnum==($2-1));
		$addnum=$2-1;
		$testprint=1 unless ($seqid eq $1);
		$seqid=$1;
		unless (exists $seqhash{$seqid}) {
			die "Error: unknown sequence $seqid length at line: $line\n";
		}
		$testprint=1 unless ($seqlength == $seqhash{$seqid});
		$seqlength=$seqhash{$seqid};
		$testprint=1 unless ($reverseadd == ($seqlength-$3));
		$reverseadd=$seqlength-$3;
		
		unless (defined $addnum and $addnum=~/^\d+$/ and defined $reverseadd and $reverseadd=~/^\d+$/ and defined $seqlength and $seqlength=~/^\d+$/) {
			die "Error: invalid line$linenum: $line\n";
		}
		if ($testprint) {
			print "Info: Seq: $seqid\tLength: $seqlength\tAddnumber=$addnum\tREVERSE=$reverseadd\n";
		}
		$arr[9]=$seqid;
	}
	else {
		die "Error: invalid sequence name at line$linenum: $line\n";
	}
	$arr[10]=$seqlength;
	$arr[11]+=$addnum;
	$arr[12]+=$addnum;
	if ($arr[8] eq '+') {
		my @arr2=split(/,/, $arr[19]);
		for (my $i=0;$i<scalar(@arr2);$i++){
			$arr2[$i]+=$addnum if ($arr2[$i]=~/^\d+$/);
		}
		$arr[19]=join(',', @arr2);
	}
	elsif ($arr[8] eq '-') {
		my @arr2=split(/,/, $arr[19]);
		for (my $i=0;$i<scalar(@arr2);$i++){
			$arr2[$i]+=$reverseadd if ($arr2[$i]=~/^\d+$/);
		}
		$arr[19]=join(',', @arr2);
	}
	else {
		die "Error: invalid strand at line$linenum: $line\n";
	}
	print OUTPUTPSL join("\t", @arr), "\n"; 
}
close INPUTPSL;
close OUTPUTPSL;


print "### Process PSL: $linenum lines\n\n";

exit 0;
