#!/usr/bin/env perl
use strict;
use warnings;
use constant USAGE =><<EOH;

usage: $0 input filter.out output

    Only keep both mate mapped reads
    
#Input
	ReadName1	1	...
	ReadName1	2	...
	ReadName2	2	...
	ReadName3	1	...

#filter.out
	ReadName2	2	...
	ReadName3	1	...

#output
	ReadName1	1	...
	ReadName1	2	...

v20160921

EOH
die USAGE if (scalar(@ARGV) !=3 or $ARGV[0] eq '-h' or $ARGV[0] eq '--help');


die "Error: invalid input\n" unless (defined $ARGV[0] and -s $ARGV[0]);
die "Error: invalid filter.out name\n" unless (defined $ARGV[1]);
die "Error: filter.out existing: $ARGV[1]\n" if (-e $ARGV[1]);
die "Error: invalid output name\n" unless (defined $ARGV[2]);
die "Error: output existing: $ARGV[2]\n" if (-e $ARGV[2]);



my %idhash=();
my $linenum=0;
my $numfilterout=0;
my $numkeep=0;



open (INPUT, "< $ARGV[0]") || die "Error: Can not open input\n";
while (my $line=<INPUT>) {
	chomp $line;
	$linenum++;
	my @arr=split(/\t/, $line);
	die "Error: invalid line($linenum): $line\n" unless (scalar(@arr)>1);
	die "Error: invalid mate number at line($linenum): $line\n" unless ($arr[1]=~/^[1-2]{1}$/);
	$idhash{$arr[0]}{$arr[1]}++;
}
close INPUT;



my @allid=keys %idhash;
foreach my $indid (@allid) {
	unless (exists $idhash{$indid}{1} and exists $idhash{$indid}{2}) {
		delete $idhash{$indid};
	}
}

$linenum=0;
open (INPUT, "< $ARGV[0]") || die "Error: Can not open input\n";
open (FILTEROUT, "> $ARGV[1]") || die "Error: can not write filter.out\n";
open (OUTPUT, "> $ARGV[2]") || die "Error: can not write output\n";
while (my $line=<INPUT>) {
	chomp $line;
	$linenum++;
	my @arr=split(/\t/, $line);
	if (exists $idhash{$arr[0]}) {
		print OUTPUT $line, "\n";
		$numkeep++;
	}
	else {
		print FILTEROUT $line, "\n";
		$numfilterout++;
	}
}
close INPUT;
close FILTEROUT;
close OUTPUT;


print "### SUMMARY ###\nTotal lines: $linenum\nTotal Filterout: $numfilterout\nTotal Keep: $numkeep\n\n\n";
