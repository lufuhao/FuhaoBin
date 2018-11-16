#!/usr/bin/env perl
use warnings;
use strict;
use constant USAGE =><<EOH;

usage: $0 wcdout fasta.fai output minlen

To prepare fastaid:
	samtools faidx XX.fa
	cat XX.fa.fai |cut -f 1 > fastaid  ###ignore
	
v20160616

EOH
die USAGE if (scalar(@ARGV) !=4);

my $wcdout=$ARGV[0];
my $fileID=$ARGV[1];
my $output=$ARGV[2];
my $minlength=$ARGV[3];



###Default###
$minlength=0 unless (defined $minlength);


open (IDFILE, "<$fileID") || die "Error: can not open fileID: $fileID\n";
my $linenum=0;
my %idhash=();
my %notlength=();
while (my $line1=<IDFILE>) {
	chomp $line1;
	my @arr=split(/\s+/, $line1);
	if (scalar(@arr)<2 or $arr[1]!~/^\d+$/) {
		die "Error: invalid fasta index, Line: ", $linenum+1 ," File: $fileID\n";
	}
	$idhash{$linenum}=$arr[0];
	if ($arr[1]<$minlength) {
		$notlength{$linenum}++;
	}
	$linenum++;
}
close IDFILE;
print "Total sequence number: $linenum\n";



open (WCDOUT, "<$wcdout") || die "Error: can not open $wcdout\n";
open (OUTPUT, ">$output") || die "Error: can not write $output\n";
$linenum=0;
my $lines2keep=0;
my $lines2ignore=0;
while (my $line2=<WCDOUT>) {
	$linenum++;
	chomp $line2;
	if ($line2=~/\.$/) {
		$line2=~s/\.$//;
	}
	my @arr=();
	@arr=split(/\s+/, $line2);
	die "Error: empty wcdout at line $linenum\n" if (scalar(@arr)<1);
	my @arr2=();
	my $test2keep=0;
	foreach (@arr) {
		if (exists $idhash{$_}) {
			$test2keep=1 unless (exists $notlength{$_});
			push (@arr2, $idhash{$_});
		}
		else {
			die "Error: ID $_ not exists at line $linenum\n";
		}
	}
	die "Error: empty output at line $linenum\n" if (scalar(@arr2)<1);
	my $newline=join("\t", @arr2);
	if ($test2keep==1) {
		print OUTPUT $newline."\n";
		$lines2keep++;
	}
	else {
		print STDERR "LengthInfo: Ignored: Line $linenum, FastaID: $newline\n";
		$lines2ignore++
	}
}
close WCDOUT;
close OUTPUT;
print "Total  cluster lines: $linenum\n";
print "Kept   cluster lines: $lines2keep\n";
print "Ignore cluster lines: $lines2ignore\n";
exit 0;
