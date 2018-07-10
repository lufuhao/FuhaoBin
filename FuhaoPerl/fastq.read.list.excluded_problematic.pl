#!/usr/bin/env perl
use strict;
use warnings;
use constant USAGE =><<EOH;

usage: $0 input fasta.fai MaxInsertSize filter.out output

    Filter Out those problematic reads
        * mapped to one ref_seq >= 2 times
        * mapped to the same strand of one ref_seq
        * the end mapped to '-' strand should > the start mapped to +

#input
readname1	1	seqID1	start	end	strand


v20180710

EOH
die USAGE if (scalar(@ARGV) !=5 or $ARGV[0] eq '-h' or $ARGV[0] eq '--help');


die "Error: invalid input\n" unless (defined $ARGV[0] and -s $ARGV[0]);
die "Error: invalid fasta.fai\n" unless (defined $ARGV[1] and -s $ARGV[1]);
die "Error: invalid MaxInsertSize\n" unless (defined $ARGV[2] and $ARGV[2]=~/^\d+$/);
die "Error: invalid filter.out name\n" unless (defined $ARGV[3]);
die "Error: filter.out existing: ".$ARGV[3]."\n" if (-e $ARGV[3]);
die "Error: invalid output name\n" unless (defined $ARGV[4]);
die "Error: output existing: ".$ARGV[4]."\n" if (-e $ARGV[4]);



my %idhash=();
my $linenum=0;
my %seqlength=();
my %problematic=();
my $numfilterout=0;
my $numkeep=0;



### Read fasta.fai to get seq length
open (FASTAINDEX, "< $ARGV[1]") || die "Error: can not open fasta.fai\n";
while (my $line=<FASTAINDEX>) {
	chomp $line;
	my @arr=split(/\t/, $line);
	die "Error: invalid fasta.fai format\n" unless (scalar(@arr)>=2);
	die "Error: duplicated seqID: $arr[0]\n" if (exists $seqlength{$arr[0]});
	$seqlength{$arr[0]}=$arr[1];
}
close FASTAINDEX;


### 
$linenum=0;
open (INPUT, "< $ARGV[0]") || die "Error: Can not open input\n";
while (my $line=<INPUT>) {
	chomp $line;
	$linenum++;
	my @arr=split(/\t/, $line);
	die "Error: invalid line($linenum): $line\n" unless (scalar(@arr)>1);
	die "Error: invalid mate number at line($linenum): $line\n" unless ($arr[1]=~/^[1-2]{1}$/);
	if (exists $problematic{$arr[0]}) {
		next;
	}
	if (exists $idhash{$arr[0]} and exists $idhash{$arr[0]}{$arr[2]} and exists $idhash{$arr[0]}{$arr[2]}{$arr[1]}) {
		$problematic{$arr[0]}++;
		next;
	}
	$idhash{$arr[0]}{$arr[2]}{$arr[1]}{'start'}=$arr[3];
	$idhash{$arr[0]}{$arr[2]}{$arr[1]}{'end'}=$arr[4];
	$idhash{$arr[0]}{$arr[2]}{$arr[1]}{'strand'}=$arr[5];
}
close INPUT;



my @allid=keys %idhash;
foreach my $indid (@allid) {
	if (exists $problematic{$indid}) {
		delete $idhash{$indid} if (exists $idhash{$indid});
		next;
	}
	foreach my $indseq (keys %{$idhash{$indid}}) {
		my @matenum=keys %{$idhash{$indid}{$indseq}};
		if (scalar(@matenum)==1) {
#			if ($idhash{$indid}{$indseq}{$matenum[0]}{'strand'} eq '+') {
#				if ($idhash{$indid}{$indseq}{$matenum[0]}{'start'}+$ARGV[2]<=$seqlength{$indseq}) {
#					$problematic{$indid}++;
#					last;
#				}
#			}
#			elsif ($idhash{$indid}{$indseq}{$matenum[0]}{'strand'} eq '-') {
#				if ($idhash{$indid}{$indseq}{$matenum[0]}{'end'}>=$ARGV[2]) {
#					$problematic{$indid}++;
#					last;
#				}
#			}
#			else {
#				die "Error: invalid strand for READ $indid MATE $matenum[0] REFERENCE $indseq\n";
#			}
		}
		elsif (scalar(@matenum)==2) {
			if ($idhash{$indid}{$indseq}{1}{'strand'} eq '+') {
#				if ($idhash{$indid}{$indseq}{2}{'strand'} eq '-') {
#					unless ($idhash{$indid}{$indseq}{2}{'end'}>$idhash{$indid}{$indseq}{1}{'start'}) {
#						$problematic{$indid}++;
#						last;
#					}
#				}
#				else {
#					$problematic{$indid}++;
#					last;
#				}
			}
			elsif ($idhash{$indid}{$indseq}{1}{'strand'} eq '-') {
				if ($idhash{$indid}{$indseq}{2}{'strand'} eq '+') {
					unless ($idhash{$indid}{$indseq}{1}{'end'}>$idhash{$indid}{$indseq}{2}{'start'}) {
						$problematic{$indid}++;
						last;
					}
				}
				else {
					$problematic{$indid}++;
					last;
				}
			}
		}
		else {
			die "Error: no mates for READ $indid MATE $matenum[0] REFERENCE $indseq\n";
		}
	}
}
%idhash=();



$linenum=0;
open (INPUT, "< $ARGV[0]") || die "Error: Can not open input\n";
open (FILTEROUT, "> $ARGV[3]") || die "Error: can not write filter.out\n";
open (OUTPUT, "> $ARGV[4]") || die "Error: can not write output\n";
while (my $line=<INPUT>) {
	chomp $line;
	$linenum++;
	my @arr=split(/\t/, $line);
	if (exists $problematic{$arr[0]}) {
		print FILTEROUT  $line, "\n";
		$numfilterout++;
	}
	else {
		print OUTPUT $line, "\n";
		$numkeep++;
	}
}
close INPUT;
close FILTEROUT;
close OUTPUT;

print "### SUMMARY ###\nTotal lines: $linenum\nTotal Filterout: $numfilterout\nTotal Keep: $numkeep\n\n\n";
