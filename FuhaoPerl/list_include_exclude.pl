#!/usr/bin/env perl
use warnings;
use strict;
use constant USAGE =><<EOH;

Usage: $0 include1 exclude2 out_include out_exclude out_notsure

Version: 20120410

EOH
die USAGE if (scalar(@ARGV) != 5);

my $in_include=$ARGV[0];
chomp $in_include;
die "MainError: invalid include1\n" unless (-s $in_include);
my $in_exclude=$ARGV[1];
chomp $in_exclude;
die "MainError: invalid exclude2\n" unless (-s $in_exclude);
my $out_include=$ARGV[2];
chomp $out_include;
unlink ($out_include) if (-s $out_include);
my $out_exclude=$ARGV[3];
chomp $out_exclude;
unlink ($out_exclude) if (-s $out_exclude);
my $out_notsure=$ARGV[4];
chomp $out_notsure;
unlink ($out_notsure) if (-s $out_notsure);
print "\n\n\n#####  SUMMARY  #####\n";
print "---Input:\n";
print "------Include: $in_include\n";
print "------Exclude: $in_exclude\n";
print "---Output\n";
print "------Include: $out_include\n";
print "------Exclude: $out_exclude\n";
print "------NotSure: $out_notsure\n";


our %include=();
our %exclude=();


open (IN1, "<$in_include") || die "MainError: open file $in_include\n";
open (IN2, "<$in_exclude") || die "mainError: open file2 $in_exclude\n";

while (my $line1=<IN1>) {
	chomp $line1;
	$line1=~s/\s+.*$//;
	$include{$line1}++;
}
close IN1;
while (my $line2=<IN2>) {
	chomp $line2;
	$line2=~s/\s+.*$//;
	$exclude{$line2}++;
}
close IN2;


print "##### SUMMARY ######\n";
print "---IN1 Include: \t".scalar(keys %include)."\t$in_include\n";
print "---IN2 Exclude: \t".scalar(keys %exclude)."\t$in_exclude\n";
print "\n";



my $num_include=0;
my $num_exclude=0;
my $num_notsure=0;
open (OUT3, ">$out_include") || die "MainError: write file $out_include\n";
open (OUT4, ">$out_exclude") || die "MainError: write file $out_exclude\n";
open (OUT5, ">$out_notsure") || die "MainError: write file $out_notsure\n";
foreach my $readid (keys %include) {
	if (exists $exclude{$readid}){
		print OUT5 "$readid\n";
		delete $exclude{$readid};
		$num_notsure++;
	}
	else {
		print OUT3 "$readid\n";
		$num_include++;
	}
}
close OUT3;
close OUT5;
foreach my $readid2 (keys %exclude) {
	print OUT4 "$readid2\n";
	$num_exclude++;
}
close OUT4;
print "---Out3 include:\t$num_include\t$out_include\n";
print "---Out4 exclude:\t$num_exclude\t$out_exclude\n";
print "---Out5 notsure:\t$num_notsure\t$out_notsure\n\n";

exit 0;
