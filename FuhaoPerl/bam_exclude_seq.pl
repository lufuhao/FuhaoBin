#!/usr/bin/env perl
use strict;
use warnings;
use constant USAGE =><<EOH;

usage: $0 in.bam id1,id2 out.bam

EOH
die USAGE if (scalar(@ARGV)!=3);


my $input=$ARGV[0];
die "Error: invalid input\n" unless ( -s $input);
my $excludelist=quotemeta($ARGV[1]);
my $output=$ARGV[2];
die "Error: output exists\n" if ( -e $output);
my @idlist=split(/,/, $excludelist);
if ($input=~/\.sam$/i) {
	open (BAMIN, "samtools view -Sh $input |") || die "Can not open SAM input: $input\n";
}
elsif ($input=~/\.bam$/i) {
	open (BAMIN, "samtools view -h $input |") || die "Can not open BAM input: $input\n";
}
else {
	die "Error: input not ending with .sam or .bam\n";
}
if ($output=~/\.sam$/i) {
	open (BAMOUT, " | samtools view -Sh - > $output") || die "Can not open SAM input: $output\n";
}
elsif ($output=~/\.bam$/i) {
	open (BAMOUT, " | samtools view -bhS - > $output") || die "Can not write BAM ouput: $output\n";
}
else {
	die "Error: output not ending with .sam or .bam\n";
}
my %ids=();
map {$ids{$_}++} @idlist;


while (my $line=<BAMIN>) {
	chomp $line;
	if ($line=~/\@SQ.*SN.(\S+)\s*/) {
		unless (exists $ids{$1}) {
			print BAMOUT $line."\n";
		}
	}
	else {
		my @arr=split(/\t/, $line);
		unless (exists $ids{$arr[2]}) {
			print BAMOUT $line."\n";
		}
	}
}
close BAMIN;
close BAMOUT;
exit 0;
