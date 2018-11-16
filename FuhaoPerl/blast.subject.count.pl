#!/usr/bin/perl
#usage: blast subject count in an order provided
#usage: perl blast.target.count.pl input list output

our $input = $ARGV[0];
our $target_list=$ARGV[1];
our $output=$ARGV[2];

open (INPUT, "$input") || die "Can not open input file\n";
open (OUTPUT, ">>$output") || die "Can not output\n";


our ($temp, %count);
$temp=0;

while (our $blasts = <INPUT>) {
	chomp $blasts;
	@blast=split(/\t/, $blasts);
	if ($temp==0){
		$pre_query=$blast[0];
		$pre_subject=$blast[1];
		$count{$blast[1]}++;
		$temp++;
	}
	elsif ($pre_query eq $blast[0]) {
		if ($pre_subject ne $blast[1]){
		$count{$blast[1]}++;
		}
	}
	else {
		$count{$blast[1]}++;
	}
}
close INPUT;

open (LIST, "$target_list") || die "can not open list file\n";
while (our $list_line=<LIST>) {
	chomp $list_line;
	if (exists $count{$list_line}) {
	}
	else {
		$count{$list_line}=0;
	}
	print OUTPUT $list_line."\t".$count{$list_line}."\n";
}
close LIST;
close OUTPUT;
