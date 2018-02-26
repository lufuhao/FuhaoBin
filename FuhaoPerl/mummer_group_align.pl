#!/usr/bin/env perl
use strict;
use FuhaoPerl5Lib::FastaKit qw/GroupMummerShowcoords ExtractFastaSeqtk/;
use Cwd qw/abs_path getcwd/;
use Data::Dumper qw/Dumper/;
use FuhaoPerl5Lib::CmdKit qw/exec_cmd_return/;
use FuhaoPerl5Lib::FileKit qw /MergeFiles/;
use warnings;

use constant USAGE =><<EOH;

usage: $0 show-coords.in ref.fa query.fa [genome.fa]

show-coords -l -r -c -T nucmer.delta > show-coords.in

v20171221

EOH
die USAGE if (scalar(@ARGV) <3 or $ARGV[0] eq '-h' or $ARGV[0] eq '--help');


my $coordsfile=$ARGV[0];
my $reffasta=$ARGV[1];
my $queryfasta=$ARGV[2];
my $genomefasta=$ARGV[3];

$coordsfile=abs_path($coordsfile);
$reffasta=abs_path($reffasta);
$queryfasta=abs_path($queryfasta);
$genomefasta=abs_path($genomefasta);

print "Coords file: $coordsfile\n";
print "Reference: $reffasta\n";
print "Quert: $queryfasta\n";
print "Genome: $genomefasta\n" if(defined $genomefasta and -s $genomefasta);



my %querids=();
my %refids=();
my $curdir=getcwd;


open COORDSIN, "< $coordsfile" || die "Error: can not open coords file\n";
<COORDSIN>;<COORDSIN>;<COORDSIN>;<COORDSIN>;
while (my $line=<COORDSIN>) {
	chomp $line;
	my @arr=split(/\t/, $line);
	my $queryid=pop @arr;
	my $rfnid=pop @arr;
	$querids{$queryid}++;
	$refids{$rfnid}++;
}
close COORDSIN;


my ($test, $rethash)=GroupMummerShowcoords($coordsfile);
die "Error: GroupMummerShowcoords failed\n"unless ($test);

open (TEMPLIST, "> temp.list") || die "Error: can not write temp.list\n";

foreach my $num (sort {$a<=>$b} keys %{$rethash}) {
	mkdir "$curdir/group$num" || die "Error: can not mkdir: $curdir/group$num\n";
	chdir "$curdir/group$num" || die "Error: can not chdir: $curdir/group$num\n";
	print TEMPLIST "group", $num;
	my @queryies=();
	my @refseqs=();
	foreach my $id (sort keys %{$$rethash{$num}}) {
		if (exists $querids{$id}) {
			push (@queryies, $id);
		}
		if (exists $refids{$id}) {
			push (@refseqs, $id);
		}
	}
	unless (scalar(@queryies)>0) {
		print STDERR "Warnings: no query\n";
		print Dumper $$rethash{$num};
		print "\n";
	}
	unless (scalar(@refseqs)>0) {
		print STDERR "Warnings: no refs\n";
		print Dumper $$rethash{$num};
		print "\n";
	}
	
	open (RFNSEQS, " > $curdir/group$num/group$num.rfn.list") || die "Error: can not write ref list: $curdir/group$num/group$num.rfn.list\n";
	foreach (@refseqs) {
		print RFNSEQS $_, "\n";
		print TEMPLIST "\t", $_;
	}
	close RFNSEQS;
	unless (ExtractFastaSeqtk ($reffasta, "$curdir/group$num/group$num.rfn.list.fa", "$curdir/group$num/group$num.rfn.list", 'seqtk')) {
		print STDERR "Warnings: failed to extract rfn: $curdir/group$num/group$num.rfn.list.fa\n";
		next;
	}
	
	open (QRYSEQS, " > $curdir/group$num/group$num.qry.list") || die "Error: can not write query list: $curdir/group$num/group$num.qry.list\n";
	foreach (@queryies) {
		print QRYSEQS $_, "\n";
		print TEMPLIST "\t", $_;
	}
	close QRYSEQS;
	print TEMPLIST "\n";
	unless (ExtractFastaSeqtk ($queryfasta, "$curdir/group$num/group$num.qry.list.fa", "$curdir/group$num/group$num.qry.list", 'seqtk')) {
		print STDERR "Warnings: failed to extract qry: $curdir/group$num/group$num.qry.list.fa\n";
		next;
	}
	
	unless (exec_cmd_return("mum.stat -i 99 -r $curdir/group$num/group$num.rfn.list.fa -q $curdir/group$num/group$num.qry.list.fa -o group$num.test1")) {
		print STDERR "Warnings:mummerplot error 1: group$num\n";
		next;
	}
	if (defined $genomefasta and -s $genomefasta) {
		unless (MergeFiles("$curdir/group$num/group$num.test2.fa", "$curdir/group$num/group$num.rfn.list.fa", "$curdir/group$num/group$num.qry.list.fa")) {
			print STDERR "Error: merge files failed: group$num\n";
			next;
		}
		unless (exec_cmd_return("mum.stat -i 99 -r $genomefasta -q $curdir/group$num/group$num.test2.fa -o group$num.test2")) {
		print STDERR "Warnings:mummerplot error 2: group$num\n";
		next;
	}
	} 
}
close TEMPLIST;
