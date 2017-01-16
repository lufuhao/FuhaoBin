#!/usr/bin/env perl
use strict;
use warnings;
use constant USAGE =><<EOH;

usage: $0 sspace.final.evidence conversion.seqid output.sum
	
Links scaffold ID by SSPACE to its's previous original contigs IDs
	sspace.final.evidence
	conversion.seqid
		#Format: Col1:currentnameInFinal.fasta	Col2:previousInFormatted.fasta
	outputprefix


v20161031

EOH
die USAGE if (scalar(@ARGV) !=3 or $ARGV[0] eq '-h' or $ARGV[0] eq '--help');



die "Error: invalid evidence input\n" unless (defined "$ARGV[0]" and -s "$ARGV[0]");
my $inputfasta=$ARGV[0];
die "Error: invalid seqid conversion\n" unless (defined "$ARGV[1]" and -s "$ARGV[1]");
my $seqidconversion=$ARGV[1];
die "Error: invalid output\n" unless (defined $ARGV[2]);
my $outputsum=$ARGV[2];
unlink $outputsum;
unlink "$outputsum.scaffold";
unlink "$outputsum.contig";
my $linenum=0;
my $seqnum=0;
my %idhash=();



open (IDCONVERSION, "< $seqidconversion") || die "Error: invalid conversion.seqid\n";
while (my $line=<IDCONVERSION>) {
	chomp $line;
	$linenum++;
	my @arr=split(/\t/, $line);
	unless (scalar(@arr)==2 and $arr[0]=~/^\S+$/ and $arr[1]=~/^\S+$/) {
		die "Error: invalid conversion.seqid line($linenum): $line\n";
	}
	unless (exists $idhash{$arr[0]}) {
		$idhash{$arr[0]}=$arr[1];
	}
	else {
		die "Error: duplicated ID: $arr[0] in $seqidconversion\n";
	}
}
close IDCONVERSION;
print "\n\n\nSUMMARY: total conversion.seqid line: $linenum\nTotal hash: ", scalar(keys %idhash), "\n\n";
#map {print $_, "\n";} keys %idhash;



$linenum=0;
my $total_scaffolded=0;
my $outnum=0;
open (EVIDENCE, " < $inputfasta") || die "Error: invalid input fasta\n";
open (OUTPUT, " > $outputsum") || die "Error: can not write output\n";
open (SCAFFOLD, " > $outputsum.scaffold") || die "Error: can not write scaffold: $outputsum.scaffold\n";
open (CONTIG, " > $outputsum.contig") || die "Error: can not write scaffold: $outputsum.contig\n";
my @arr=();
my $scaffoldid='';
while (my $line=<EVIDENCE>) {
	chomp $line;
	$linenum++;
#>scaffold1|size4217569|tigs2
#r_tig1795|size43901|links295|gaps3825
#f_tig936|size4169843

#>scaffold2|size6264750|tigs3
#f_tig903|size3662090|links107|gaps-8536
#r_tig933|size600697|links84|gaps-3720
#r_tig461|size2001961
	if ($line=~/^>/) {
		if ($line=~/^>(scaffold\d+\|size\d+)\|/) {
			$scaffoldid=$1;
		}
		else {
			die "Error: invalid scaffold name at line($linenum): $line\n";
		}
	}
	elsif ($line=~/\S+/) {
		if ($line=~/^[fr]_tig(\d+)\|/) {
			push (@arr, $1);
		}
		else {
			die "Error: invalid contig name at line($linenum): $line\n";
		}
	}
	else {
		if (scalar(@arr)>1) {
			print OUTPUT "$scaffoldid";#, join ("\t", @arr), "\n";
			print SCAFFOLD "$scaffoldid\n";#, join ("\t", @arr), "\n";
			foreach my $id (@arr) {
				if (exists $idhash{"contig$id"}) {
					print OUTPUT "\t", $idhash{"contig$id"};
				}
				else {
					die "Error: seqID contig$id do not have converted name\n";
				}
				$total_scaffolded++;
			}
			print OUTPUT "\n";
			print CONTIG join ("\t", @arr), "\n";
		}
		$scaffoldid='';
		@arr=();
	}
}
close EVIDENCE;
close OUTPUT;
close SCAFFOLD;
close CONTIG;
print "\n\n\nSUMMARY: total scaffolds: $total_scaffolded\n";
