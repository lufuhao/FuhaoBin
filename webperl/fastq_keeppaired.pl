#!/usr/bin/env perl
use warnings;
use strict;

die "perl $0 <IN1:PE1.fastq> <IN2:PE2.fastq> <OUT1:out.PE1.fastq> <OUT2:out.PE2.fastq>\n" if (@ARGV!=4);

#@ILLUMINA-57021F:7:1:1020:11315#0/1
my ($PE1,$PE2,$out1,$out2)=@ARGV;
my %hash;
open (IN,$PE1) || die $!;
while(my $line1=<IN>) {
	if ($line1=~/^\@(\S+)\/[12]\S*\s*/) {
		my $readid=quotemeta($1);
		$hash{$readid}=$line1;
		$hash{$readid}.=<IN>;
		$hash{$readid}.=<IN>;
		$hash{$readid}.=<IN>;
	}
}
close IN;

open (IN,$PE2) || die $!;
open (OUT1,">$out1") || die $!;
open (OUT2,">$out2") || die $!;
open (SE1,">$out1.single") || die $!;
open (SE2,">$out2.single") || die $!;
while(my $line2=<IN>) {
	if ($line2=~/^\@(\S+)\/[12]\S*\s*/){
		my $readid2=quotemeta($1);
		if (exists $hash{$readid2}){
			print OUT1 $hash{$readid2};
			undef $hash{$readid2};
			delete $hash{$readid2};
			print OUT2 $line2;
			$_=<IN>;
			print OUT2 $_;
			$_=<IN>;
			print OUT2 $_;
			$_=<IN>;
			print OUT2 $_;
		}
		else {
			print SE2 $line2;
			$_=<IN>;
			print SE2 $_;
			$_=<IN>;
			print SE2 $_;
			$_=<IN>;
			print SE2 $_;
		}
	}
}
foreach my $key(keys %hash) {
	if ((defined $hash{$key})&&($hash{$key} ne "")) {
		print SE1 $hash{$key};
	}
}

close IN;
close OUT1;
close OUT2;
close SE1;
close SE2;
