#!/usr/bin/perl
###usage: script.pl input output
our $input=$ARGV[0];
chomp $input;
our $output=$ARGV[1];
chomp $output;

open (INPUT, "$input") ||  die "Can not input\n";
open (OUTPUT, ">>$output") || die "Can not output\n";

our $temp01=0;
our (@pre_blast_line, $temp02, $temp03);

while (our $blast_line=<INPUT>) {
	chomp $blast_line;
	our @blast_comp=split(/\t/, $blast_line);
#0	1	2	3	4		5	6	7	8	9	10	11		12	13
#qseqid	sseqid	pident	length	mismatch	gapopen	qstart	qend	sstart	send	evalue	bitscore	qcovs	qcovhsp
	if ($temp01==0) {
		$pre_blast_line[1] = $blast_line;
		our $pre_query=$blast_comp[0];
		our $pre_subject=$blast_comp[1];
		our $pre_qcovhsp=$blast_comp[13];
		our $total_qcovhsp=$pre_qcovhsp;
		$temp01++;
		our $rep=1;
	}
	elsif ($blast_comp[0] eq $pre_query){
		if ($blast_comp[1] eq $pre_subject) {
			$rep++;
			$pre_blast_line[$rep]=$blast_line;
			if ($rep ==2) {
				$total_qcovhsp=$pre_qcovhsp + $blast_comp[13];
			}
			elsif ($rep>2){
				$total_qcovhsp=$total_qcovhsp + $blast_comp[13];
			}
		}
		else {
			if ($rep>1) {
				if ($total_qcovhsp >= 80) {
					for ($temp02=1; $temp02<=$rep; $temp02++) {
						if ($total_qcovhsp <=102) {
							print OUTPUT $pre_blast_line[$temp02]."\n";
						}
						else {
							my $rep_line=$pre_blast_line[$temp02];
							my @rep_lines=split(/\t/, $rep_line);
							if ($rep_lines[13] >=80){
								print OUTPUT $pre_blast_line[$temp02]."\n";
							}
						}
					}
				}
				$temp03=1;
				$rep =1;
				$total_qcovhsp=0;
			}
			elsif ($pre_qcovhsp >= 80) {
				print OUTPUT $pre_blast_line[1]."\n";
				$temp03=0;
			}
			$pre_blast_line[1] = $blast_line;
			$pre_subject=$blast_comp[1];
			$pre_qcovhsp=$blast_comp[13];
		}
	}
	
	else {
		if ($rep>1) {
			if ($total_qcovhsp >= 80) {
				for ($temp02=1; $temp02<=$rep; $temp02++) {
					if ($total_qcovhsp <=102) {
						print OUTPUT $pre_blast_line[$temp02]."\n";
					}
					else {
						my $rep_line=$pre_blast_line[$temp02];
						my @rep_lines=split(/\t/, $rep_line);
						if ($rep_lines[13] >=80){
							print OUTPUT $pre_blast_line[$temp02]."\n";
						}
					}
				}
			}
			$temp03=1;
			$rep =1;
			$total_qcovhsp=0;
			}
		elsif ($temp03 !=1 && $pre_qcovhsp >= 80) {
			print OUTPUT $pre_blast_line[1]."\n";
		}		
		$pre_blast_line[1] = $blast_line;
		$pre_query=$blast_comp[0];
		$pre_subject=$blast_comp[1];
		$pre_qcovhsp=$blast_comp[13];
		$temp03=0;
	}
}

#lastline
if ($rep>1) {
	if ($total_qcovhsp >= 80) {
		for ($temp02=1; $temp02<=$rep; $temp02++) {
			if ($total_qcovhsp <=102) {
				print OUTPUT $pre_blast_line[$temp02]."\n";
			}
			else {
				my $rep_line=$pre_blast_line[$temp02];
				my @rep_lines=split(/\t/, $rep_line);
				if ($rep_lines[13] >=80){
					print OUTPUT $pre_blast_line[$temp02];
				}
			}
		}
	}
	$temp03=1;
	$rep =1;
	$total_qcovhsp=0;
	}
elsif ($temp03 !=1 && $pre_qcovhsp >= 80) {
	print OUTPUT $blast_line[1];
}

close OUTPUT;
close INPUT;
