#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long;
use constant USAGE=><<EOH;

SYNOPSIS:

perl $0 pcrdup.output max_sd
Version: LUFUHAO20140401

Descriptions:
	Collect insertsize after 'bamutils pcrdup' (ngsutils)
	bamutils pcrdup -counts pcrdup.count 1094_LIB9374_LDI7704_L0_bwaalnQ5.sort.bam
	

Example:
	perl $0 pcrdup.count 20 20 > pcrdup.cluster

Author:
	Fu-Hao Lu
	Post-Doctoral Scientist in Micheal Bevan laboratory
	Cell and Developmental Department, John Innes Centre
	Norwich NR4 7UH, United Kingdom
	E-mail: Fu-Hao.Lu\@jic.ac.uk
EOH
###HELP ends#########################################################
die USAGE unless @ARGV;



###Default and initalization#########################################
our $input=$ARGV[0];
our $max_sd=(defined $ARGV[1]) ? $ARGV[1] : 10;
our $max_is_sd=(defined $ARGV[2]) ? $ARGV[2] : 10;
print "##############\nInput:\t$input\nMax_SD:\t$max_sd\nMax_insert_sd:\t$max_is_sd\n##################\n";

my $temp_01=0;
our $temp_02=0;
our @groups=();
our @prev_arr=();
open (INPUT, "$input") || die "Can not find input file: $input\n";
while (my $line=<INPUT>) {
	chomp $line;
	my @arr=();
	@arr=split(/\t/, $line);
	next if (scalar(@arr) != 4);
	if ($temp_01 == 0) {
#		print "\n\n\narr\t@arr\n";###Test###
		@{$groups[$temp_02++]}=@arr;
		@prev_arr=@arr;
		$temp_01++;
		next;
	}
	if ($arr[0] eq $prev_arr[0]) {
#		print "equal\n";###Test###
#		print "$arr[1]-$prev_arr[1]=",abs($arr[1]-$prev_arr[1]),"\n";###Test###
		if (abs($arr[1]-$prev_arr[1]) <= $max_sd){
			@{$groups[$temp_02++]}=@arr;
		}
		else {
			&ps_groups() unless (scalar(@groups)==0);
			$temp_02=0;
			@{$groups[$temp_02++]}=@arr;
			$temp_01=1;
		}
	}
	else {
		&ps_groups() unless (scalar(@groups)==0);
		$temp_02=0;
		@{$groups[$temp_02++]}=@arr;
		$temp_01=1;
	}
	@prev_arr=@arr;
#	print "\n\n\narr\t@arr\n";###Test###
#	print "prev_arr\t@prev_arr\n";###Test###
}
&ps_groups() unless (scalar(@groups)==0);


#####################################################################
###                         sub functions                         ###
#####################################################################

###regroup @groups and print
###&ps_groups()
###global variable: @groups
sub ps_groups {
#	print "Number: ",scalar(@groups),"\n";			###Test###
#	for (my $i=0; $i<scalar(@groups);$i++) {		###Test###
#		print "@{$groups[$i]}\n";			###Test###
#	}							###Test###
	if (scalar(@groups)==1) {
		print join("\t", @{$groups[0]}), "\n";
		@groups=();
		return 0;
	}
	else {
		my @new_groups=();
		my $temp_03=0;
#		print "@{$groups[0]}\n";			###Test###
		@{$new_groups[0]}=@{$groups[0]};
#		print "\@new_groups: @{$new_groups[0]}\n";	###Test###
		$temp_03++;
		splice(@groups,0,1);
#		print "Number2: ",scalar(@groups),"\n";		###Test###
		for (my $i=0; $i<scalar(@groups);$i++) {
#			print "***\n";				###Test###
			for (my $j=0; $j<scalar(@new_groups);$j++){
#				print "*** ***\n";					###Test###
#				print "\@new_groups: @{$new_groups[$j]}\n";		###Test###
#				print "${$groups[$i]}[2] - ${$new_groups[$j]}[2]=", ${$groups[$i]}[2]-${$new_groups[$j]}[2],"\n";###Test###
				if (abs(${$groups[$i]}[2] - ${$new_groups[$j]}[2]) <= $max_is_sd) {
					@{$new_groups[$temp_03++]}=@{$groups[$i]};
					$groups[$i]=0;
					last;
				}
			}
		}
		my $temp04=0;
		my @new2_groups=();
		for (my $i=0; $i<scalar(@groups);$i++) {
			if ($groups[$i] != 0) {
				@{$new2_groups[$temp04++]}=@{$groups[$i]};
			}
		}
		@groups=@new2_groups;

		my $sum=0;
		for (my $a=0; $a<scalar(@new_groups);$a++){
			$sum+=${$new_groups[$a]}[3];
		}
		pop @{$new_groups[0]};
		print join("\t",@{$new_groups[0]}),"\t$sum\n";
		if (scalar(@groups)>0){
			&ps_groups();
		}
		else{
			@groups=();
			return 0;
		}
	}
}
