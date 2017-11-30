#!/usr/bin/env perl
use strict;
use warnings;
use FuhaoPerl5Lib::FastaKit qw/AnalyzeMummerShowcoords/;
use Bio::DB::Fasta;
use Bio::SeqIO;
use Bio::Seq;
use constant USAGE =><<EOH;

nucmer -c 100 -p ori rfn.fa query.fa
delta-filter -i 99 ori.delta > ori.delta.filter
show-coords -T -l -r ori.delta.filter > coordinates

usage: $0 fasta coordinates seqname output

v20171124

EOH
die USAGE if (scalar(@ARGV) !=4 or $ARGV[0] eq '-h' or $ARGV[0] eq '--help');



my $fastain=$ARGV[0];
my $coords=$ARGV[1];
my $seqname=$ARGV[2];
my $output=$ARGV[3];



unless (defined $fastain and -s $fastain) {
	die "Error: invalid fasta\n";
}
unless (defined $coords and -s $coords) {
	die "Error: invalid coords\n";
}
unless (defined $output) {
	die "Error: invalid output\n";
}
unlink $output if (-s $output);



my ($test, $hashindex)=AnalyzeMummerShowcoords($coords);
unless ($test) {
	die "Error: AnalyzeMummerShowcoords failed\n";
}
my $db=Bio::DB::Fasta->new($fastain);
my $seqioout= Bio::SeqIO->new( -format => 'Fasta' , -file => ">$output");
foreach (keys %{$hashindex}) {
	my $idvname=$seqname++;
	my $idvseq='';
	my $desc='';
	
	for (my $i=0; $i<scalar(@{${$hashindex}{$_}}); $i+=2) {
		unless (${$hashindex}{$_}[$i]=~/^(0)|(1)$/) {
			die "Error: unknown strand\n";
		}
		
		if (${$hashindex}{$_}[$i+1]=~/^(\S+):(\d+)\-(\d+)$/) {
			my $thispartseq='';
			my $thisseqname=$1;
			my $thisseqstart=$2;
			my $thisseqend=$3;
			if (${$hashindex}{$_}[$i]==0) {
				$thispartseq=$db->seq("$thisseqname", $thisseqstart => $thisseqend);
				$desc.=" $thisseqname:$thisseqstart-$thisseqend";
			}
			elsif (${$hashindex}{$_}[$i]==1) {
				$thispartseq=$db->seq("$thisseqname", $thisseqend => $thisseqstart);
				$desc.=" $thisseqname:$thisseqend-$thisseqstart";
			}
			else {
				die "Error3\n";
			}
			if ($thispartseq eq '') {
				die "Error: 4\n";
			}
			else {
				$idvseq.=$thispartseq;
			}
		}
	}
	my $seqobj=Bio::Seq->new( -id => "$idvname",
                              -desc => "$desc",
                             -seq => $idvseq);
                             
    $seqioout->write_seq($seqobj);
}

