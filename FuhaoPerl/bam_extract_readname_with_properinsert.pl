#!/usr/bin/env perl
use strict;
use warnings;
use Data::Dumper qw /Dumper/;
use FuhaoPerl5Lib::BamKit qw /VerifyCigarLength CalCigarRefLength/;
use constant USAGE =><<EOH;

usage: $0 in.bam mininsert maxinsert out.reads

Requirements
	BAM are read name sorted using
		samtools sort -n original.bam namesorted

Descriptions
	Extract tose read names from BAM file based on proper insert
	Will
		verify cigar length and read sequence length
		verify if paired from BAM FLAG
		verify different strand for each pair from BAM FLAG
		verify proper insert

v20160611

EOH
die USAGE if (scalar(@ARGV) !=4 or $ARGV[0] eq '-h' or $ARGV[0] eq '--help');

my $path_samtools='samtools';

my $linenum=0;
my %readhash=();
my $mininsert=20000;
my $maxinsert=60000;
$mininsert=$ARGV[1] if (defined $ARGV[1] and $ARGV[1]=~/^\d+$/);
$maxinsert=$ARGV[2] if (defined $ARGV[2] and $ARGV[2]=~/^\d+$/);

#my ($test, $bamobj)=ReadSam()
open (BAMIN, "$path_samtools view -h $ARGV[0] | ") || die "Error: can not open BAM file\n";
open (READOUT, "> $ARGV[3]") || die "Error: can not write read output\n";
print "Start...\n";### For test ###
while (my $line=<BAMIN>) {
	$linenum++;
	next if ($line=~/^\@/);
	chomp $line;
	my @arr=split(/\t/, $line);
	unless (exists $readhash{$arr[0]}) {
		foreach my $readname (sort keys %readhash) {
#			print "Readhash:\n"; print Dumper \%readhash; print "\n";### For test ###
#			print "Read: $readname\t";### For test ###
			print READOUT "$readname\t";
			my $test2usethisread=0;
			foreach my $refer (sort keys %{$readhash{$readname}}) {
				if (exists $readhash{$readname}{$refer}{'nouse'}) {
					$test2usethisread=0;
					last;
				}
				if (! exists $readhash{$readname}{$refer}{1} or ! $readhash{$readname}{$refer}{2}) {
					next;
				}
				if ($readhash{$readname}{$refer}{1}{'strand'} eq $readhash{$readname}{$refer}{2}{'strand'}) {
					next;
				}
				my $firstend=$readhash{$readname}{$refer}{1}{'posit'} + CalCigarRefLength($readhash{$readname}{$refer}{1}{'cigar'});
				my $secondend=$readhash{$readname}{$refer}{2}{'posit'} + CalCigarRefLength($readhash{$readname}{$refer}{2}{'cigar'});
				my @arrregion=($readhash{$readname}{$refer}{1}{'posit'}, $firstend, $readhash{$readname}{$refer}{2}{'posit'}, $secondend);
				@arrregion=sort {$a<=>$b} @arrregion;
				my $insertlength=$arrregion[-1]-$arrregion[0] +1;
				if ($insertlength>=$mininsert and $insertlength<=$maxinsert) {
					$test2usethisread=1;
				}
#				print "Insert: $insertlength\n";### For test ###
				print READOUT $readhash{$readname}{$refer}{1}{'posit'}, "\t", $readhash{$readname}{$refer}{2}{'posit'}, "\tInsert: ", $insertlength, "\t";
			}
			if ($test2usethisread==1) {
#				print READOUT "$readname\n";
				print READOUT "1\n";
			}
			else {
#				print "Ignored";### For test ###
				print READOUT "0\n";
			}
			print "\n";
		}
		%readhash=();
	}
	if ($arr[1] & 0x0004) {##Not mapped
		next;
	}
	unless (VerifyCigarLength($arr[5], length($arr[9]))) {
		print STDERR "Warnings: CIGAR length problem at line($linenum): $line\n";
		next;
	}
	my $matepairnum=0;
	if ($arr[1] & 0x0040) {
		$matepairnum=1;
	}
	if ($arr[1] & 0x0080) {
		$matepairnum=2;
	}
	my $strand='';
	if ($arr[1] & 0x0010) {
		$strand='-';
	}
	else {
		$strand='+';
	}
	if (exists $readhash{$arr[0]} and exists $readhash{$arr[0]}{$arr[2]} and exists $readhash{$arr[0]}{$arr[2]}{$matepairnum}) {
		$readhash{$arr[0]}{$arr[2]}{'nouse'}=1;
		next;
	}
	$readhash{$arr[0]}{$arr[2]}{$matepairnum}{'cigar'}=$arr[5];
	$readhash{$arr[0]}{$arr[2]}{$matepairnum}{'strand'}=$strand;
	$readhash{$arr[0]}{$arr[2]}{$matepairnum}{'posit'}=$arr[3];
}
foreach my $readname (sort keys %readhash) {
	my $test2usethisread=0;
	foreach my $refer (sort keys %{$readhash{$readname}}) {
		if (exists $readhash{$readname}{$refer}{'nouse'}) {
			$test2usethisread=0;
			last;
		}
		if (! exists $readhash{$readname}{$refer}{1} or ! $readhash{$readname}{$refer}{2}) {
			next;
		}
		if ($readhash{$readname}{$refer}{1}{'strand'} ne $readhash{$readname}{$refer}{2}{'strand'}) {
			next;
		}
		my $firstend=$readhash{$readname}{$refer}{1}{'posit'} + CalCigarRefLength($readhash{$readname}{$refer}{1}{'cigar'});
		my $secondend=$readhash{$readname}{$refer}{2}{'posit'} + CalCigarRefLength($readhash{$readname}{$refer}{2}{'cigar'});
		my @arrregion=($readhash{$readname}{$refer}{1}{'posit'}, $firstend, $readhash{$readname}{$refer}{2}{'posit'}, $secondend);
		@arrregion=sort {$a<=>$b} @arrregion;
		my $insertlength=$arrregion[-1]-$arrregion[0];
		if ($insertlength>=$mininsert and $insertlength<=$maxinsert) {
			$test2usethisread=1;
		}
	}
	print READOUT "$readname\n" if ($test2usethisread==1);
}
close BAMIN;
close READOUT;
#D1SF08P1:48:h7gwvbcxx:1:1101:1045:75623	99	v443_0261	627374	255	205M	=	668992	41795	ATAATTTCTGGAGTCCGAATAGGGTCAAGCCGAACTATGGACCTTTTATATGTGTGCGGCCCCACACACCACACAAGGAGTTTTGAAGCCCTGTTTTACATATACAGGGAGACATGTACAGTCTTCTTGTGGTTTGATTTTGTAGGAGTATAGCCGCTTAATCAAGCCCAGATTTATCTGAGCTGCTACACCTAGGTGTGCGCCG	DDDDDIIIIIIHIIIIIHIIIIIIHIIIIGGIIIIHHIIHHIIIIIIIIIIIIIIIIIGHHHHHIIEHIHHECHH?DC1DFHHHIHIIIIIIHGHHHHIEHHFHHIIIHEHHIIIHIIIIIIIIIGHGHEFHEHHGIHIGIII<GCGHHIEHHICHIIHDHHHHHHHHHHCHG<DHGFHE?GEHIIH..9CHII@FH?GIHHIID	XA:i:0	MD:Z:205NM:i:0
#D1SF08P1:48:h7gwvbcxx:1:1101:1045:75623	147	v443_0261	668992	255	177M	=	627374	-41795	CCCACGACGACACATACAGAGAGGATTTCCAGCGGTTCCAATGGAGATCCCCGTGGCTCAAATGGCGCTCCCAGCGGCCGACGACGACAATGATGGGTAGTTCACCATGCATCACGTAGTGGAAACCCATGTGATGTGGAAGGTCTTCTCGGTGGTGTACACAAGTGAGCCGGTCTT	HHA7?-ECHHFAG.EGCG@HEG@HD:.C?C--EDEHE@EHHE@@F<<C<<HCGIIIHHEHH@CCHCHCFE/DDHDEHCGHHHFHGHFC@C@HIHIIIHIHHHHF@HFHHHFEHFCIHHHCHIHHGGEG@IIHGIHFIIHHIIHGHCCDCDHIHEHGHHHHHECGIIHGHCHCDDDDD	XA:i:0	MD:Z:177	NM:i:0
#D1SF08P1:48:h7gwvbcxx:1:1101:1045:85407	99	v443_1048	155523	255	213M	=	189073	33726	ATGTGGCTCCCACATTGTAAGTGTACACCAGAGTTATAGTCCTTGAAGATTTTGTTTAGCCGACAGAATACGGCCCCTAAAAATGTATTTCTTCAAGTGAAGAAATTGGTACCGCCTTTGAAGAAATTGATGTTAAGACTTTGAAGCTTGAAGACTTTTGTTTTCATAGTTTCTTTCTTCTTATTTGAGTCATAGGAAAAAACATAATGTTAA	DDDDDIIIIIIIIIIIIIIIIIIIHIIIIIIHIIIIIIFH?HHIIIIIHIIIIIIIIIIIIH/FHHIIIIIIIIGHIIIIHIEHIIIIEHHIIHHHGHHHHDGHHIGIIIIGFHHIIIIIIGHHIHHIHHHIIHHHIIIIIIHHEHICCHHHIIGIIHIEEHCGHEHG?GHIIGHI0FEHE?/DGHE/D@</DD@FHHGFHIHICGDC@D..:	XA:i:0	MD:Z:213	NM:i:0
#D1SF08P1:48:h7gwvbcxx:1:1101:1045:85407	147	v443_1048	189073	255	176M	=	155523	-33726	ATCCTCAGTCCGGCTTTTCTTTCTGAACACCGAACCAAGCATCAGGGGCTACTGTCTATGCGGTACTATTTTACATACATCAAAATTCTTACCTCAAAGCCTGTTAGAAGGCTGCTGCAGGCTTCGGTCAGCCCGCTCTTGGCGAGCTGAACCTTCTGAATCACCGCACTCATAAC	BBAF.HHCDECHHHHHFD/</FHEHC=HHD?@<@GGEHHIHEHGIHIIIIIIIHGHCIHDHGGHIIHHHHEHHG<C@<10F@GHF1IHGCCHIHHHHHHEHGEHHECHHHCHEHHEIIHHHCIIIIIIIHHHIIIIIGIGIIIHCIHEEIIIIHIIIIIIIHHHIIIHHHIDDDDD	XA:i:0	MD:Z:176	NM:i:0
#D1SF08P1:48:h7gwvbcxx:1:1101:1047:4540	99	v443_0801	1223723	255	226M	=	1260949	37403	GGAGAAGGGCGGAAAACGACAGCCCAACATGTAATGAATGATATCGATCAACCAACCATGCTCGAATATATCTTGGCCTCCGTGCGAGCCCGTCTGTGTGTTCTAGCTCCTCGTCTCTCTCAAGGAACGGCGGGTTGGATCTTTGCCGTCAATTGAGATCGGTTTGCTCACACTCCTGTCCCACATGTCAGTGATATGCCAAGAGGAGGTGCATTGACCGACTGGC	DDDDDIIIIIHIIIIIIIIIIIIIIGIIIIIIIIIIHHHIIIIIIIHIIIIIIIIIHHHHIIIIGHIIHIEHIHIIGHHHIHHHIIIIGIIIHIIIIGHIIDHHIHIIIIGIIIIIHIHIGIIIIHIIHHHHIIDHDHHHIIIGHHH?D?/BFHHIIIIGIIHHHHH.FFGHHHEHHHEECGHHG?A@AFHHHIEFGAHHG@-8-B@EG-B--8@-8-8?7>=C-@	XA:i:0	MD:Z:226	NM:i:0
