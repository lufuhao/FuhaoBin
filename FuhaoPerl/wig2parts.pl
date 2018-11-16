#!/usr/bin/env perl
use strict;
use warnings;
use constant USAGE =><<EOH;

usage: $0 1.file.wig 2.minimum_cov 3.output_GFF
Version: 20150529

EOH
die USAGE if (scalar(@ARGV) !=3 or $ARGV[0] eq '-h' or $ARGV[0] eq '--help' or $ARGV[0] eq 'help' or $ARGV[0] eq '-help');

my $filewig=$ARGV[0];
die "Wig2GffError: invalid WIG input\n" unless ( -s $filewig);
my $minicov=$ARGV[1];
my $outgff=$ARGV[2];
open (WIGIN, "< $filewig") || die "Wig2GffError: Can not open WIG file\n";
open (GFFOUT, "> $outgff") || die "Wig2GffError: Can not write GFF file\n";
my $test_start=0;
my $test_end=0;
my $chrom='';
my $start;
my $end;
my $previous_end;
my %poshash=();
my $linenum=0;
my $partnum=1;
while (my $line=<WIGIN>) {
	$linenum++;
	chomp $line;
	if ($line=~/chrom=(\S+)/) {
		if (scalar(keys %poshash) >0) {
#			print "First print\n"; ### For test ###
			$test_start=1;
			$partnum=1;
			foreach my $i (sort {$a <=> $b} keys %poshash) {
				if ($test_start==1) {
					$start=$i;
					$test_start=0;
					$previous_end=$i;
					next
				}
				else {
					if ($i != $previous_end +1 ) {
						print GFFOUT "$chrom\twig2gff\tCov$minicov\t$start\t$previous_end\t.\t.\t.\tID=$chrom.part$partnum\n";
#						print "$chrom\twig2gff\tCov$minicov\t$start\t$previous_end\t.\t.\t.\tID=$chrom.part$partnum\n";
						$partnum++;
						$start=$i;
					}
					$previous_end=$i;
				}
			}
			print GFFOUT "$chrom\twig2gff\tCov$minicov\t$start\t$previous_end\t.\t.\t.\tID=$chrom.part$partnum\n";
#			print "$chrom\twig2gff\tCov$minicov\t$start\t$previous_end\t.\t.\t.\tID=$chrom.part$partnum\n";
			%poshash=();
		}
		$chrom=$1;
	}
	elsif ($line=~/^(\d+)\s+(\d+)$/) {
		my $base=$1;
		my $cov=$2;
		if (exists $poshash{$base}) {
			die "Wig2GffError: existing $chrom:$base at $linenum\n";
		}
		else {
			if ($cov>=$minicov) {
				$poshash{$base}=$cov;
			}
		}
	}
	else {
		print STDERR "Wig2GffWarnings: unknown line at $linenum\n\t$line\n";
	}
}
$test_start=1;
$partnum=1;
#print "Second print\n"; ### For test ###
foreach my $i (sort {$a <=> $b} keys %poshash) {
	if ($test_start==1) {
		$start=$i;
		$test_start=0;
		$previous_end=$i;
		next
	}
	else {
		if ($i != $previous_end +1 ) {
			print GFFOUT "$chrom\twig2gff\tCov$minicov\t$start\t$previous_end\t.\t.\t.\tID=$chrom.part$partnum\n";
#			print "$chrom\twig2gff\tCov$minicov\t$start\t$previous_end\t.\t.\t.\tID=$chrom.part$partnum\n";
			$partnum++;
			$start=$i;
		}
		$previous_end=$i;
	}
}
print GFFOUT "$chrom\twig2gff\tCov$minicov\t$start\t$previous_end\t.\t.\t.\tID=$chrom.part$partnum\n";
#print "$chrom\twig2gff\tCov$minicov\t$start\t$previous_end\t.\t.\t.\tID=$chrom.part$partnum\n";
close WIGIN;
close GFFOUT;
