#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long;
#use Bio::SeqIO;
#use Bio::Seq;
#use Data::Dumper;
use constant USAGE =><<EOM;


SYNOPSIS:
perl file.pl -i my_seqs.blast -d ref_db.fasta -o coverage.output -l int

Version 20150210

OPTIONS:
    --input|i <File>
	[Msg] The tab output by blasting my seq fasta file against db;
    --database|-d <File>
	[Msg] The blastdb used to define the full-length of reference
    --output|-o <File>
	[Msg] Output file (tab delimited)
    --alignlength|-l <Int>
	[Opt] Default: 0


EXAMPLES:
    perl script.pl -i xx -d xx -o xx

AUTHOR:
    Fu-Hao Lu
    Post-Doctoral Scientist in Micheal Bevan Lab
    Cell and Developmental Biology
    John Innes Centre
    Norwich NR4 7UH, UK
    Fu-Hao.Lu\@jic.ac.uk

EOM
die USAGE unless (@ARGV);


###Get options##################################################
our ($help, $input, $db, $output, $alignlength);

GetOptions (
	'help|h!'=>\$help,
	'input|i:s'=>\$input,
	"database|d:s"=>\$db,
	"output|o:s"=>\$output,
	"alignlength|l:i"=>\$alignlength,
)or die USAGE;

$help and die USAGE;



###Variables default############################################
our %covparts=();
our %covlength=();
$alignlength=0 unless defined $alignlength;
die USAGE unless (defined $input or defined $db or defined $output);
print "\n\n\n##### Summary #####\n";
print "Input file: $input\nDB file: $db\nOutput file: $output\n";         ###for test
die "Invalid blast input\n" unless (-s $input);
unlink ("$output") if (-s $output);
###Variables default####end#####################################



###Calculate the reference db seq length########################
#our ($seqio_obj, $seq_obj, %db_length);
#%db_length=();
#$seqio_obj=Bio::SeqIO->new(-file=>"$db", -format=>"fasta");
#while($seq_obj=$seqio_obj->next_seq){
#	if (exists $db_length{$seq_obj->id}) {
#		print "Error: check if the seq_id in my.fasta is unique\n" and die;
#	}
#	else {
#		$db_length{$seq_obj->id}=$seq_obj->length;
#	}
#}



our %db_length=();
open (DB, "$db") || die "Can not open fastaDB file\n";
my $temp_first=0;
my $seqid='';
my $seq_length=0;
while (my $line=<DB>) {
	chomp $line;
	if ($line=~/^>(\S+)\s*/) {
		unless ($temp_first>0) {
			$db_length{$seqid}=$seq_length;
			$seq_length=0;
		}
		$seqid=$1;
	}
	else {
		$seq_length+=length($line);
	}
}
$db_length{$seqid}=$seq_length;
###Calculate the reference db seq length########end#############



###sub functions################################################
sub swapinsert {###To Make interval arrays for each ref_seq
	my $contig_id="";
	my $arr_start=0;
	my $arr_end=0;
	($contig_id, $arr_start, $arr_end)=@_;
	my @contigse=($arr_start, $arr_end);
#	print "The parameter for swap: $contig_id\t@contigse\n";            ###for test###
	my $ext_start=0;
	my $ext_end=0;
	my $test_start=0;
	my $test_end=0;
	my $test_left_line=0;
	our @empty=(0,0);
	SW: {
	if (exists $covparts{$contig_id}) {
#		print "Insert:".$covparts{$contig_id}[0][0]."\t".$covparts{$contig_id}[0][1]."\n";  ###for test###
		if ($arr_end < $covparts{$contig_id}[0][0]) {
			unshift @{$covparts{$contig_id}},[@contigse];
			last SW;
		}
		elsif ($arr_end==$covparts{$contig_id}[0][0]) {
			$covparts{$contig_id}[0][0]=$arr_start;
			last SW;
		}
		elsif ($arr_start > $covparts{$contig_id}[scalar(@{$covparts{$contig_id}})-1][1]) {
			push @{$covparts{$contig_id}},[@contigse];
			last SW;
		}
		elsif ($arr_start== $covparts{$contig_id}[scalar(@{$covparts{$contig_id}})-1][1]) {
			$covparts{$contig_id}[scalar (@{$covparts{$contig_id}}-1)][1]=$arr_end;
			last SW;
		}
		else {
			for (my $i=0; $i<scalar @{$covparts{$contig_id}}; $i++) {
				if ($test_start==0) {
					if ($arr_start <=$covparts{$contig_id}[$i][0]){
						$ext_start=$arr_start;
						$test_start++;
						$test_left_line=$i;
						if ($arr_end<=$covparts{$contig_id}[$i][1]) {
							$covparts{$contig_id}[$i][0]=$ext_start;
							$test_end++;
							last SW;
						}
					}
					elsif ($arr_start >=$covparts{$contig_id}[$i][0] and $arr_start <=$covparts{$contig_id}[$i][1]) {
						$test_start++;
						$test_left_line=$i;
						if ($arr_end<=$covparts{$contig_id}[$i][1]){
							$test_end++;
							last SW;
						}
					}
				}
				elsif ($test_start>0 and $test_end ==0) {
					if ($arr_end<=$covparts{$contig_id}[$i][0]) {
						$covparts{$contig_id}[$test_left_line][1]=$arr_end;
						$test_end++;
						last SW;
					}
					elsif ($arr_end>=$covparts{$contig_id}[$i][0] and $arr_end<=$covparts{$contig_id}[$i][1]) {
						$covparts{$contig_id}[$test_left_line][1]=$covparts{$contig_id}[$i][1];
						@{$covparts{$contig_id}[$i]}=(0,0);
						last SW;
					}
					else {
						@{$covparts{$contig_id}[$i]}=(0,0);
						next;
					}
				}
			}

		}
	}
	else {
		@{$covparts{$contig_id}[0]}=@contigse;
#		print "FirstTime:".$covparts{$contig_id}[0][0]."\t".$covparts{$contig_id}[0][1]."\n";
	}
} ##SW_ends
my $rm_undef_i=0;
my @new_arr=();
for (my $rm_undef_j=0; $rm_undef_j<scalar @{$covparts{$contig_id}}; $rm_undef_j++){
#	if (defined @{$covparts{$contig_id}[$rm_undef_j]} and $covparts{$contig_id}[$rm_undef_j][1] != 0) {
	if ($covparts{$contig_id}[$rm_undef_j][1] != 0) {
#		print "rev.def1: @{$covparts{$contig_id}[$rm_undef_j]}\n";   ###important for test###
		@{$new_arr[$rm_undef_i]}=@{$covparts{$contig_id}[$rm_undef_j]};
#		print "rev.def2: @{$new_arr[$rm_undef_i]}\n";            ###important for test###
		$rm_undef_i++;
	}
}
@{$covparts{$contig_id}}=@new_arr;
#print $covparts{$contig_id}[0][0]."\t".$covparts{$contig_id}[0][1]."\n";     ###for test###
}



sub covlen {###Calculat the non-overlap length intervals
	my ($covlen_id) = @_;
#	$covlen_id."\n";        ###For test###
	chomp $covlen_id;
	my $covlen_total=0;
	for (my $covlen_i=0; $covlen_i < scalar @{$covparts{$covlen_id}}; $covlen_i++) {
		if (defined $covparts{$covlen_id}[$covlen_i]){
#			print $covparts{$covlen_id}[$covlen_i][0]."\t".$covparts{$covlen_id}[$covlen_i][1]."\n";   ###test###
			$covlen_total+=($covparts{$covlen_id}[$covlen_i][1]-$covparts{$covlen_id}[$covlen_i][0]+1);
		}
		else {
			die "Undefined combination in $covlen_id\n";
		}
	}
	return $covlen_total;
}
###sub functions#######################end######################


our %idcov=();
open (MYBLAST, "$input") || die "Failed to open the MyFasta.file\n";
while (our $input_line=<MYBLAST>) {
	chomp $input_line;
	our @input_arr=();
	@input_arr=split(/\t/, $input_line);
#	print "@input_arr\n";             ###test###
#format
#0	1	2	3	4	5	6	7	8	9	10	11
#qid	sid	ident	length	mism	gap	qs	qe	ss	se	e	bits
	if (defined $input_arr[3] and $input_arr[3] >= $alignlength){
#		print "Send to swap: $input_arr[3]\n$input_arr[8]\n$input_arr[9]\n";
		if ($input_arr[8] < $input_arr[9]) {
			&swapinsert($input_arr[1], $input_arr[8], $input_arr[9]);
		}
		else {
			&swapinsert($input_arr[1], $input_arr[9], $input_arr[8]);
		}
	}
	else {
#		print "Please make sure your input is tab delimited blast/blast++ output\n";
#		die USAGE;
		next;
	}
}
close MYBLAST;


###calcaulate###################################################
open (OUTPUT, ">>$output") || die " can not write the output\n";
print OUTPUT "Contig_id\tRef_length\tCov_Length\tCov_Perc\tNum_Frag\n";
our @contigs_ids=keys %db_length;
foreach my $ind_contig_id (@contigs_ids) {
#	print $ind_contig_id."\n";             ###for test###
	if (exists $covparts{$ind_contig_id}) {
		print OUTPUT $ind_contig_id."\t".$db_length{$ind_contig_id}."\t".covlen($ind_contig_id)."\t".covlen($ind_contig_id)/$db_length{$ind_contig_id}."\t".scalar(@{$covparts{$ind_contig_id}})."\n";
	}
	else {
		print OUTPUT $ind_contig_id."\t".$db_length{$ind_contig_id}."\t0\t0\t0\n";
	}
}
###calcaulate#########################ends######################
