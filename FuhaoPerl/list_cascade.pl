#!/usr/bin/env perl
use strict;
use warnings;
use FuhaoPerl5Lib::MiscKit qw/GetCascadeList/;
use constant USAGE =><<EOH;

usage: $0 list1in list2in list3out

###list1 input
ID1	subID1	subID2	...

###List2
subID1	subsubID1	subsubID2	...
subID2	subsubID3	subsubID4	...

### List3 output
ID1	subsubID1	subsubID2	subsubID3	subsubID4

v20171130

EOH
die USAGE if (scalar(@ARGV) !=3 or $ARGV[0] eq '-h' or $ARGV[0] eq '--help');




unless (GetCascadeList($ARGV[0], $ARGV[1], $ARGV[2])) {
	die "Error: failed to get list\n";
}


