#!/usr/bin/env perl
use strict;
use warnings;
use constant USAGE =><<EOH;

usage: perl $0 findorf.output.gtf output_prefix

v20140425

EOH
die USAGE if (scalar(@ARGV) !=2 or $ARGV[0] eq '-h' or $ARGV[0] eq '--help');

###Configureable
my $min_orf_length=100;


my $input=$ARGV[0];
chomp $input;
my $output_prefix=$ARGV[1];
chomp $output_prefix;
my $tabular_out=$output_prefix.".tab";
my $stat_out=$output_prefix.".stat";
unlink ($tabular_out) if (-e $tabular_out);
unlink ($stat_out) if (-e $stat_out);


my %diff_5prime_most_start_and_orf=();
my %contig_len=();
my %distant_start=();
my %num_5prime_ATG=();
my %num_orf_candidates=();
my %orf_type=();
my %orf100type=();
my %internal_stop=();
my %majority_frameshift=();
my %most_5prime_query_start=();
my %most_5prime_relative=();
my %no_prediction_reason=();
my %most_5prime_sbjct_start=();
my %pfam_extended_5prime=();
my %num_relatives=();
my %ORFid=();
#my $num_orf_less_100=0;
open (INPUT, $input) || die "can not input\n";
open (TABOUT, ">>$tabular_out") || die "can not tabout\n";
print TABOUT "SeqID\tstart\tend\tstrand\tframe\tdiff_5prime_most_start_and_orf\tcontig_len\tdistant_start\tnum_5prime_ATG\tnum_orf_candidates\torf_type\tinternal_stop\tmajority_frameshift\tmost_5prime_query_start\tmost_5prime_relative\tno_prediction_reason\tmost_5prime_sbjct_start\tpfam_extended_5prime\tnum_relatives\n";
while (my $input_lines=<INPUT>) {
	chomp $input_lines;
	next if (!defined $input_lines or $input_lines eq '');
	my @input_arr=split(/\t/, $input_lines);
	if ($input_arr[8]=~m/diff_5prime_most_start_and_orf\s+(.*);.*contig_len\s+(.*);.*distant_start\s+(.*);.*num_5prime_ATG\s+(.*);.*num_orf_candidates\s+(.*);.*orf_type\s+(.*);.*internal_stop\s+(.*);.*majority_frameshift\s+(.*);.*most_5prime_query_start\s+(.*);.*most_5prime_relative\s+(.*);.*no_prediction_reason\s+(.*);.*most_5prime_sbjct_start\s+(.*);.*pfam_extended_5prime\s+(.*);.*num_relatives\s+(.*)/) {
		print TABOUT "$input_arr[0]\t$input_arr[3]\t$input_arr[4]\t$input_arr[6]\t$input_arr[7]\t$1\t$2\t$3\t$4\t$5\t$6\t$7\t$8\t$9\t$10\t$11\t$12\t$13\t$14\n";

		$diff_5prime_most_start_and_orf{$1}++;
		$contig_len{$2}++;
		$distant_start{$3}++;
		$num_5prime_ATG{$4}++;
		$num_orf_candidates{$5}++;
		
		my $thisorflength=0;
		$thisorflength=($input_arr[4] ne '.' and $input_arr[3] ne '.')?(abs($input_arr[4]-$input_arr[3])+1):0;
		if ($thisorflength >= $min_orf_length) {
			$orf_type{$6}++;
		} else {
			$orf100type{$6}++;
		}		
		$internal_stop{$7}++;
		$majority_frameshift{$8}++;
		$most_5prime_query_start{$9}++;
		$most_5prime_relative{$10}++;
		$no_prediction_reason{$11}++;
		$most_5prime_sbjct_start{$12}++;
		$pfam_extended_5prime{$13}++;
		$num_relatives{$14}++;
		if ($6 eq "full_length") {
			$ORFid{$input_arr[0]}=$thisorflength;
		}
	}
	else {
		print "Error: at $input_arr[0]\n";
	}
}
close INPUT;
close TABOUT;



open (STATOUT, ">>$stat_out") || die "Can not statout\n";
#my @diff_5prime_most_start_and_orf=keys %diff_5prime_most_start_and_orf;
#print STATOUT "diff_5prime_most_start_and_orf:\n";
#foreach (@diff_5prime_most_start_and_orf) {
#	print STATOUT $_."\t".$diff_5prime_most_start_and_orf{$_}."\n";
#}
#print STATOUT "\n\n\n";

#my @distant_start=keys %distant_start;
#print STATOUT "distant_start:\n";
#foreach (@distant_start) {
#	print STATOUT $_."\t".$distant_start{$_}."\n";
#}
#print STATOUT "\n\n\n";

#my @num_5prime_ATG=keys %num_5prime_ATG;
#print STATOUT "num_5prime_ATG:\n";
#foreach (@num_5prime_ATG) {
#	print STATOUT $_."\t".$num_5prime_ATG{$_}."\n";
#}
#print STATOUT "\n\n\n";

my @orf_type=keys %orf_type;
print STATOUT "orf_type:\n";
foreach (@orf_type) {
	print STATOUT $_."\t".$orf_type{$_}."\n";
}
my @orf100type=keys %orf100type;
print STATOUT "orf100type:\n";
foreach (@orf100type) {
	print STATOUT $_."\t".$orf100type{$_}."\n";
}
print STATOUT "\n\n\n";

my @internal_stop=keys %internal_stop;
print STATOUT "internal_stop:\n";
foreach (@internal_stop) {
	print STATOUT $_."\t".$internal_stop{$_}."\n";
}
print STATOUT "\n\n\n";

my @majority_frameshift=keys %majority_frameshift;
print STATOUT "majority_frameshift:\n";
foreach (@majority_frameshift) {
	print STATOUT $_."\t".$majority_frameshift{$_}."\n";
}
print STATOUT "\n\n\n";

my @most_5prime_relative=keys %most_5prime_relative;
print STATOUT "most_5prime_relative:\n";
foreach (@most_5prime_relative) {
	print STATOUT $_."\t".$most_5prime_relative{$_}."\n";
}
print STATOUT "\n\n\n";

my @no_prediction_reason=keys %no_prediction_reason;
print STATOUT "no_prediction_reason:\n";
foreach (@no_prediction_reason) {
	print STATOUT $_."\t".$no_prediction_reason{$_}."\n";
}
print STATOUT "\n\n\n";

my @ORFid=keys %ORFid;
print STATOUT "ORFid:\n";
foreach (@ORFid) {
	print STATOUT $_."\t".$ORFid{$_}."\n";
}
close STATOUT;




=annotation;
diff_5prime_most_start_and_orf 0;
contig_len 325;
distant_start None;
num_5prime_ATG 0;
num_orf_candidates 2;
orf_type partial;
internal_stop False;
majority_frameshift False;
most_5prime_query_start 0;
most_5prime_relative bd;
no_prediction_reason None;
most_5prime_sbjct_start 331;
pfam_extended_5prime None;
num_relatives 9

0.seqID		1.sour	2.feature	3.start	4.end	5.score	6.std	7.frame	8.desc
TauRoot00039953	findorf	predicted_orf	1	324	.	1	0	diff_5prime_most_start_and_orf 0;contig_len 325;distant_start None;num_5prime_ATG 0;num_orf_candidates 2;orf_type partial;internal_stop False;majority_frameshift False;most_5prime_query_start 0;most_5prime_relative bd;no_prediction_reason None;most_5prime_sbjct_start 331;pfam_extended_5prime None;num_relatives 9
TauRoot00039952	findorf	predicted_orf	236	1597	.	1	1	diff_5prime_most_start_and_orf 0;contig_len 1867;distant_start None;num_5prime_ATG 0;num_orf_candidates 11;orf_type full_length;internal_stop False;majority_frameshift False;most_5prime_query_start 235;most_5prime_relative bd;no_prediction_reason None;most_5prime_sbjct_start 1;pfam_extended_5prime None;num_relatives 10
TauRoot00039951	findorf	predicted_orf	3	1604	.	-1	2	diff_5prime_most_start_and_orf 0;contig_len 2842;distant_start None;num_5prime_ATG 0;num_orf_candidates 18;orf_type partial_5prime;internal_stop True;majority_frameshift True;most_5prime_query_start 2;most_5prime_relative ae;no_prediction_reason None;most_5prime_sbjct_start 239;pfam_extended_5prime None;num_relatives 9
TauRoot00039950	findorf	predicted_orf	317	760	.	-1	1	diff_5prime_most_start_and_orf 0;contig_len 1236;distant_start None;num_5prime_ATG 0;num_orf_candidates 10;orf_type full_length;internal_stop False;majority_frameshift False;most_5prime_query_start 316;most_5prime_relative ae;no_prediction_reason None;most_5prime_sbjct_start 1;pfam_extended_5prime None;num_relatives 9
TauRoot00014989	findorf	predicted_orf	1	1050	.	1	0	diff_5prime_most_start_and_orf 0;contig_len 2342;distant_start None;num_5prime_ATG 0;num_orf_candidates 11;orf_type partial_5prime;internal_stop True;majority_frameshift True;most_5prime_query_start 0;most_5prime_relative bd;no_prediction_reason None;most_5prime_sbjct_start 122;pfam_extended_5prime None;num_relatives 10
TauRoot00014988	findorf	predicted_orf	2	1288	.	1	1	diff_5prime_most_start_and_orf 0;contig_len 1448;distant_start None;num_5prime_ATG 0;num_orf_candidates 14;orf_type partial_5prime;internal_stop False;majority_frameshift False;most_5prime_query_start 124;most_5prime_relative so;no_prediction_reason None;most_5prime_sbjct_start 53;pfam_extended_5prime None;num_relatives 10
TauRoot00039955	findorf	predicted_orf	83	559	.	1	1	diff_5prime_most_start_and_orf 0;contig_len 1720;distant_start None;num_5prime_ATG 0;num_orf_candidates 19;orf_type full_length;internal_stop True;majority_frameshift True;most_5prime_query_start 82;most_5prime_relative bd;no_prediction_reason None;most_5prime_sbjct_start 1;pfam_extended_5prime None;num_relatives 10
TauRoot00039954	findorf	predicted_orf	.	.	.	.	.	diff_5prime_most_start_and_orf None;contig_len 766;distant_start None;num_5prime_ATG None;num_orf_candidates None;orf_type None;internal_stop None;majority_frameshift None;most_5prime_query_start None;most_5prime_relative None;no_prediction_reason inconsistent_strand;most_5prime_sbjct_start None;pfam_extended_5prime None;num_relatives None
TauRoot00014985	findorf	predicted_orf	187	585	.	1	0	diff_5prime_most_start_and_orf 0;contig_len 867;distant_start None;num_5prime_ATG 0;num_orf_candidates 15;orf_type full_length;internal_stop False;majority_frameshift False;most_5prime_query_start 186;most_5prime_relative bd;no_prediction_reason None;most_5prime_sbjct_start 1;pfam_extended_5prime None;num_relatives 10
TauRoot00014984	findorf	predicted_orf	186	1325	.	1	2	diff_5prime_most_start_and_orf 0;contig_len 1502;distant_start None;num_5prime_ATG 0;num_orf_candidates 17;orf_type full_length;internal_stop False;majority_frameshift False;most_5prime_query_start 185;most_5prime_relative bd;no_prediction_reason None;most_5prime_sbjct_start 1;pfam_extended_5prime None;num_relatives 10
TauRoot00014987	findorf	predicted_orf	47	613	.	-1	1	diff_5prime_most_start_and_orf 0;contig_len 1542;distant_start None;num_5prime_ATG 0;num_orf_candidates 7;orf_type full_length;internal_stop True;majority_frameshift True;most_5prime_query_start 91;most_5prime_relative bd;no_prediction_reason None;most_5prime_sbjct_start 13;pfam_extended_5prime None;num_relatives 10
TauRoot00014986	findorf	predicted_orf	150	2198	.	-1	2	diff_5prime_most_start_and_orf 0;contig_len 2664;distant_start None;num_5prime_ATG 0;num_orf_candidates 25;orf_type full_length;internal_stop False;majority_frameshift False;most_5prime_query_start 248;most_5prime_relative bd;no_prediction_reason None;most_5prime_sbjct_start 33;pfam_extended_5prime None;num_relatives 10
TauRoot00014981	findorf	predicted_orf	299	2158	.	-1	1	diff_5prime_most_start_and_orf 0;contig_len 2415;distant_start None;num_5prime_ATG 0;num_orf_candidates 19;orf_type full_length;internal_stop False;majority_frameshift False;most_5prime_query_start 280;most_5prime_relative rp;no_prediction_reason None;most_5prime_sbjct_start 12;pfam_extended_5prime None;num_relatives 10
TauRoot00014980	findorf	predicted_orf	2	364	.	1	1	diff_5prime_most_start_and_orf 0;contig_len 635;distant_start None;num_5prime_ATG 0;num_orf_candidates 9;orf_type partial_5prime;internal_stop True;majority_frameshift False;most_5prime_query_start 1;most_5prime_relative bd;no_prediction_reason None;most_5prime_sbjct_start 433;pfam_extended_5prime None;num_relatives 10
TauRoot00014983	findorf	predicted_orf	558	1997	.	-1	2	diff_5prime_most_start_and_orf 0;contig_len 2461;distant_start None;num_5prime_ATG 0;num_orf_candidates 26;orf_type full_length;internal_stop False;majority_frameshift False;most_5prime_query_start 557;most_5prime_relative bd;no_prediction_reason None;most_5prime_sbjct_start 1;pfam_extended_5prime None;num_relatives 10
TauRoot00014982	findorf	predicted_orf	2	2731	.	1	1	diff_5prime_most_start_and_orf 0;contig_len 2893;distant_start None;num_5prime_ATG 0;num_orf_candidates 17;orf_type partial_5prime;internal_stop False;majority_frameshift False;most_5prime_query_start 205;most_5prime_relative ae;no_prediction_reason None;most_5prime_sbjct_start 87;pfam_extended_5prime None;num_relatives 9
=cut;
