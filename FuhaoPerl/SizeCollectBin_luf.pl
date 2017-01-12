#!/usr/bin/env perl
use strict;
use warnings;
use Data::Dumper qw /Dumper/;
use constant USAGE =><<END;

perl $0 input_file bin_size

version: 20161125

Requirements: 
    Data::Dumper

    *Reading the number column it first meet;
    *The next column would be count if it's number

END
die USAGE if (scalar(@ARGV) !=2 or $ARGV[0] eq '-h' or $ARGV[0] eq '--help');


my $input=$ARGV[0];
my $bin=$ARGV[1];
my %count=();
my $min;
my $max;
my $test=0;
my $test2=0;
open (INPUT, "$input") || die "Can not open input file\n";
while (my $line=<INPUT>) {
	next if ($line=~/^#|^\*/);
	chomp $line;
	my @arr=();
	@arr=split(/\s+|,|;/,$line);
	my $col=0;
	if ($test==0) {
		for (my $i=0; $i<scalar(@arr);$i++){
			if ($arr[$i]=~/^\d+\.{0,1}\d*$/) {
				$col=$i;
				print STDERR "Info: col ", ($col+1), " is counted\n";
				last;
			}
		}
		if ((scalar(@arr)-1)>$col and defined $arr[$col+1] and $arr[$col+1]=~/^\d+$/) {
			print STDERR "Info: col ", ($col+2), " used as count Number\n";
			$test2=1;
		}
		$test++;
	}
	unless (defined $arr[$col]) {
		die "Error: invalid lines: $line Array $col ", join(',', @arr), "\n";
	}
	my $index=int($arr[$col]/$bin);
	if ($test2==1) {
		$count{$index}+=$arr[$col+1];
	}
	else {
		$count{$index}++;
	}
	
	if (defined $min and $min=~/^\S+$/) {
		$min=$arr[$col] if ($arr[$col]<$min);
	}
	else {
		$min=$arr[$col];
	}
	if (defined $max and $max=~/^\S+$/) {
		$max=$arr[$col] if ($arr[$col]>$max);
	}
	else {
		$max=$arr[$col];
	}
	
}
close INPUT;
#print Dumper \%count; ### For test ###

my @countkey=sort {$a <=> $b} keys %count;
print "Min:$min\tMax:$max\tBin:$bin\n";
for (my $j = $countkey[0]; $j<=$countkey[-1]; $j++) {
	if (exists $count{$j}) {
		print $j*$bin, "\t", $count{$j}, "\n";
	}
	else {
		print $j*$bin."\tNaN\n";
	}
}
