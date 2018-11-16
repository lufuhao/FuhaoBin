#!/usr/bin/env perl

use strict;
use warnings;
my $usage = "Usage: $0 <bam_flag>\n";
my $flag = shift or die $usage;
die "Please enter a numerical value\n" if $flag =~ /\D+/;

if ($flag & 0x1){
	print "template having multiple segments in sequencing\n";
}
if ($flag & 0x2){
	print "each segment properly aligned according to the aligner\n";
}
if ($flag & 0x4){
	print "segment unmapped\n";
}
if ($flag & 0x8){
	print "next segment in the template unmapped\n";
}
if ($flag & 0x10){
	print "SEQ being reverse complemented\n";
}
if ($flag & 0x20){
	print "SEQ of the next segment in the template being reversed\n";
}
if ($flag & 0x40){
	print "the first segment in the template\n";
}
if ($flag & 0x80){
	print "the last segment in the template\n";
}
if ($flag & 0x100){
	print "secondary alignment\n";
}
if ($flag & 0x200){
	print "not passing quality controls\n";
}
if ($flag & 0x400){
	print "PCR or optical duplicate\n";
}
if ($flag & 0x800){
	print "supplementary alignment\n";
}

exit(0);
__END__
