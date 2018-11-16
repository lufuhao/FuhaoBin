#!/usr/bin/env perl
#http://code.google.com/p/amyottecode/w/list
# parse through a signalp results list and return any protein that is listed 
# as "Y" for the predictive test in at least N/7 tests, which is set by the cutoff value
use strict;
use seq::tools;

open(REPORT, "$ARGV[0]");
open(OUT, '>', "$ARGV[1]");

my $multifasta = $ARGV[2];
my %seq_hash = %{ tools::FASTA2HASH(\$multifasta) };

my $cutoff = $ARGV[3];

my $count = "";
while(<REPORT>) 
{
        chomp($_);
    unless($_ =~ /^#/) 
        {
                my @list = split(/\s+|\t+/, $_);
                my %hash = ();
                my $internal_count = "";
                foreach my $list_items (@list)
                {
                        $hash{$list_items} ++
                }
                if($hash{"Y"} >= $cutoff) 
                {
                        print OUT "\>$list[14]\n";
                        print OUT $seq_hash{$list[14]} . "\n";
                        $count++;
                }
        }
}
print "$count\n";


