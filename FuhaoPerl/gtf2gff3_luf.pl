#!/usr/bin/env perl
use strict;
use warnings;
use constant USAGE =><<EOH;

usage: $0 in.gtf out.gff3

Convert GTF to GFF3
Chr3    Aet_PGSB_v2     transcript      28202   30768   .       +       .       gene_id "AET3G0000100"; transcript_id "AET3G0000100.2"; clst_id "TCONS_00107747";
Chr3    Aet_PGSB_v2     exon    28202   28356   .       +       .       gene_id "AET3G0000100"; transcript_id "AET3G0000100.2"; exon_number "1"
Chr3    Aet_PGSB_v2     5UTR    28202   28202   .       +       .       gene_id "AET3G0000100"; transcript_id "AET3G0000100.2"; 5UTR_number "1"
Chr3    Aet_PGSB_v2     3UTR    30767   30768   .       +       .       gene_id "AET3G0000100"; transcript_id "AET3G0000100.2"; 3UTR_number "1"
Chr3    Aet_PGSB_v2     transcript      29156   32580   .       +       .       gene_id "AET3G0000100"; transcript_id "AET3G0000100.9"; clst_id "TCONS_00155760";

v20160620

EOH
die USAGE if (scalar(@ARGV) !=2 or $ARGV[0] eq '-h' or $ARGV[0] eq '--help');


my $linenum=0;
my %genehash=();

open (GTFINPUT, "< $ARGV[0]") || die "Error: can not open TF input\n";
while (my $line=<GTFINPUT>) {
	chomp $line;
	$linenum++;
	next if ($line=~/^#/);
	next unless ($line=~/\ttranscript\t/);
	(my $geneid=$line)=~s/^.*gene_id\s+"//;
	$geneid=~s/".*$//;
	my @arr=split(/\t/, $line);
	push (@{$genehash{$geneid}}, $arr[3]);
	push (@{$genehash{$geneid}}, $arr[4]);
}
close GTFINPUT;

print STDERR "Info: read GTF lines: $linenum\n";
print STDERR "Info: total genes: ", scalar(keys %genehash), "\n\n";


$linenum=0;
open (GTFINPUT, "< $ARGV[0]") || die "Error: cannot open TF input\n";
open (GFF3OUT, "> $ARGV[1]") || die "Error: can not write GFF3 output\n";
print GFF3OUT "##gff-version 3\n";

while (my $line=<GTFINPUT>) {
	chomp $line;
	$linenum++;
	next if ($line=~/^#/);
	my @arr=split(/\t/, $line);
	(my $transcriptid=$line)=~s/^.*transcript_id\s+"//;
	$transcriptid=~s/".*$//;
	
	if ($arr[2]=~/^transcript$/i) {
		(my $geneid=$line)=~s/^.*gene_id\s+"//;
		$geneid=~s/".*$//;
		if (exists $genehash{$geneid}) {
			my @geneborder=sort {$a<=>$b} @{$genehash{$geneid}};
			my @newgenearr=@arr;
			$newgenearr[2]='gene';
			$newgenearr[3]=$geneborder[0];
			$newgenearr[4]=$geneborder[-1];
			$newgenearr[8]="ID=$geneid";
			print GFF3OUT join("\t", @newgenearr), "\n";
			delete $genehash{$geneid};
		}
		$arr[2]='mRNA';
		$arr[8]="ID=$transcriptid;Parent=$geneid";
		print GFF3OUT join("\t", @arr), "\n";
	}
	elsif ($arr[2]=~/^exon$/i) {
		(my $exonnum=$line)=~s/^.*exon_number\s+"//;
		$exonnum=~s/".*$//;
		$arr[8]="ID=$transcriptid.exon$exonnum;Parent=$transcriptid";
		print GFF3OUT join("\t", @arr), "\n";
	}
	elsif ($arr[2]=~/^CDS$/i) {
		$arr[8]="ID=$transcriptid.cds;Parent=$transcriptid";
		print GFF3OUT join("\t", @arr), "\n";
	}
	elsif ($arr[2]=~/^5UTR$/i) {
		$arr[2]='five_prime_UTR';
		$arr[8]="ID=$transcriptid.5UTR;Parent=$transcriptid";
		print GFF3OUT join("\t", @arr), "\n";
	}
	elsif ($arr[2]=~/^3UTR$/i) {
		$arr[2]='three_prime_UTR';
		$arr[8]="ID=$transcriptid.3UTR;Parent=$transcriptid";
		print GFF3OUT join("\t", @arr), "\n";
	}
	else {### For development
		print STDERR "Warnings: unknown line($linenum): $line\n";
	}
}
close GTFINPUT;
close GFF3OUT;
