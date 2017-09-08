#!/usr/bin/env perl
use strict;
use warnings;
use constant USAGE =><<EOH;

usage: $0 Input.gff3[.gz] score Output.gff3[.gz]

Transfer score and note from old to new GFF3

score file: 3 columns
## geneid[tab]score[tab]note
gene1	5	Part5missing

v20170704

EOH
die USAGE if (scalar(@ARGV) !=3 or $ARGV[0] eq '-h' or $ARGV[0] eq '--help');


my $inputgff3=shift @ARGV;
die "Error: invalid gff3 input\n" unless (defined $inputgff3 and -s $inputgff3);
my $score=shift @ARGV;
die "Error: invalid score input\n" unless (defined $score and -s $score);
my $outputgff3=shift @ARGV;
die "Error: invalid gff3 output\n" unless (defined $outputgff3);
unlink $outputgff3 if (-e $outputgff3);

if ($inputgff3=~/(\.gff3\.gz$)|(\.gff\.gz$)/) {
	open (INPUTGFF3, " zcat $inputgff3 | ") || die "Error: can not open gzipped GFF3 file\n";
}
elsif ($inputgff3=~/(\.gff3$)|(\.gff$)/) {
	open (INPUTGFF3, " < $inputgff3") || die "Error: can not open gzipped GFF3 file\n";
}
else {
	die "Error: input suffix shoulf be .gff[3] or .gff[3].gz\n";
}



my %hash=();
open (SCORENOTE, "< $score ") || die "Error: can not open score file\n";
while (my $line=<SCORENOTE>) {
	chomp $line;
	my @arr=split(/\t/, $line);
	if (exists $hash{$arr[0]}) {
		die "Error: repeated gene id ($.): $line\n";
	}
	unless ($arr[1] =~/^\d+\.*\d*$/) {
		die "Error: invalid score (line $.): $line\n";
	}
	$hash{$arr[0]}{'score'}=$arr[1];
	if (defined $arr[2] and $arr[2]=~/\S+/) {
		$hash{$arr[0]}{'note'}=$arr[2];
	}
	else{
		$hash{$arr[0]}{'note'}='';
	}
}
close SCORENOTE;




if ($outputgff3=~/(\.gff3\.gz$)|(\.gff\.gz$)/) {
	open (OUTPUTGFF3, " ! bgzip -c > $inputgff3") || die "Error: can not write gzipped GFF3 file\n";
}
elsif ($outputgff3=~/(\.gff3$)|(\.gff$)/) {
	open (OUTPUTGFF3, " > $outputgff3") || die "Error: can not write gzipped GFF3 file\n";
}
else {
	die "Error: output suffix shoulf be .gff[3] or .gff[3].gz\n";
}
while (my $line=<INPUTGFF3>) {
	if ($line=~/^#/) {
		print OUTPUTGFF3 $line;
		next;
	}
	chomp $line;
	my @arr=split(/\t/, $line);
	unless ($arr[2] =~/^gene$/i) {
		print OUTPUTGFF3 $line, "\n";
		next;
	}
	my $geneid=$arr[8];
	$geneid=~s/.*ID=//; $geneid=~s/;.*//;
	unless (exists $hash{$geneid}) {
		print STDERR "Warings: geneid not exists in score: $geneid\n";
		print OUTPUTGFF3 $line, "\n";
		next;
	}
	$arr[5]=$hash{$geneid}{'score'};
	if (exists $hash{$geneid}{'note'} and $hash{$geneid}{'note'}=~/\S+/) {
		$arr[8]=$arr[8].';Note='.$hash{$geneid}{'note'};
	}
	print OUTPUTGFF3 join("\t", @arr), "\n";
}
close OUTPUTGFF3;
close INPUTGFF3;
