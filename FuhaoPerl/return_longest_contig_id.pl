#!/usr/bin/env perl

sub USAGE {
	print "perl $0 fasta.file\n\nPrint the longest contig ID on STDOUT\n";
}


$input=$ARGV[0];
chomp $input;
open (INPUT, "$input") || die "Can not open input file\n";
our $line=<INPUT>;
chomp $line;
our @lines = split(/\t/, $line);
our $id=$lines[0];
our $length=$lines[1];
@lines = ();
while ($line=<INPUT>) {
	chomp $line;
	@lines = split(/\t/, $line);
	if ($lines[1] > $length) {
		$id=$lines[0];
		$length=$lines[1];}
	@lines=();
}
print $id."\n";
