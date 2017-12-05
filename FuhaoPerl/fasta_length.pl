#!/usr/bin/env perl
use strict;
use warnings;
use constant USAGE =><<'EOD';

usage: perl $0 input.fa output.list

v20171205

Description:
    Get seq length
    
Output
seq_id1[tab]length1
...

EOD
die USAGE if (scalar(@ARGV)==0);

our $input=$ARGV[0];

open (INPUT, "$input") ||die "Error: can not find input: $input\n";
my $seqid='';
my $seqlength='';
while (my $line=<INPUT>) {
	chomp $line;
	if ($line=~/^>(\S+)\s*/) {
		if ($seqid ne '') {
			print $seqid."\t".$seqlength."\n";
		}
		$seqid=$1;
		$seqlength=0;
	}
	else {
		$seqlength+=length($line);
	}
}
print $seqid."\t".$seqlength."\n";
close INPUT;
