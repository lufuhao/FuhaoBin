#!/usr/bin/perl
use strict;
use warnings;

#http://bioen-compbio.bioen.illinois.edu/RACA/
#Perl script for converting a SAM file to RACA input files (1.7 KB)
#The first argument : path to a SAM file (input)
#The second argument : path to a RACA mapping file (will be created)
#The third argument : path to a mapping stats. (mean and stdev) file (will be created)


my $f = shift;			# input SAM file
my $map_f = shift;		# output RACA mapping file
my $stat_f = shift;		# output mapping stats. file

open(O,">$map_f");
my @tlengths = ();
my %map_used = ();
open(F,$f);
while(<F>) {
	chomp;
	if (length($_) == 0 || $_ =~ /^#/ || $_ =~ /^@/) { next; }
	my @ar = split(/\s+/);
	my $qname = $ar[0];
	if (defined($map_used{$qname})) { next; }
	$map_used{$qname} = 1;	

	my ($rname, $pos, $rnext, $pnext) = ($ar[2],$ar[3],$ar[6],$ar[7]);
	my $flag = $ar[1];
	my $seq = $ar[9];
	my $len = length($seq);
	my $tlen = $ar[8];

	if ($rname eq "*" || $rnext eq "*") { next; }

	if ($rnext eq "=") { 
		$rnext = $rname; 
	}

	if ($flag & 0x0002) {
		push(@tlengths, abs($tlen));
		#printf STDERR "%d\n", abs($tlen);
	}

	my $d = "+";
	if ($flag & 0x0010) { $d = "-"; }
	my $dnext = "+";
	if ($flag & 0x0020) { $dnext = "-"; }	

	print O "1\t$len\t$d\t$rname\t$pos\n";
	print O "2\t$len\t$dnext\t$rnext\t$pnext\n";
}
close(F);
close(O);

open(O,">$stat_f");
my ($mean,$stdev) = get_stats(@tlengths);
print O "\t$mean\t$stdev";
close(O);

sub get_stats {
	my @nums = @_;

	my $total = 0;
	foreach my $num (@nums) {
		$total += $num;
	}
	my $mean = $total/(scalar @nums);

	my $total2 = 0;
	foreach my $num (@nums) {
		$total2 += ($mean - $num)**2;
	}
	my $mean2 = $total2/(scalar @nums);
	my $sdv = sqrt($mean2);
	
	return ($mean, $sdv);	
}
