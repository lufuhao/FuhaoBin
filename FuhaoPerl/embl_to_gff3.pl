#!/usr/bin/env perl
use warnings;
use strict;
use Bio::SeqIO;

#!/usr/bin/env perl
use strict;
use warnings;
use constant USAGE =><<EOH;

usage: embl_to_gff3.pl in.embl out.gff3

v20170321

EOH
die USAGE if (scalar(@ARGV) !=2 or $ARGV[0] eq '-h' or $ARGV[0] eq '--help');


### Parameter
my $emblin=$ARGV[0];
my $outgff=$ARGV[1];



### Input and output
unless (defined $emblin and -s $emblin) {
	die "Error: invalid EMBL input\n";
}
unlink $outgff if (-e $outgff);



### Main
my $in = Bio::SeqIO->new(-file=>$emblin,-format=>'EMBL');
open (OUTGFF3, ">$outgff") || die "Error: can not write GFF3 file: $outgff\n";
while (my $seq = $in->next_seq) {
	for my $feat ($seq->top_SeqFeatures) {
		print OUTGFF3 $feat->gff_string,"\n";
	}
}
close OUTGFF3;
