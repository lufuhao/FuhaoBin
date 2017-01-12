#!/usr/bin/perl -w

#-------------------------------------------------------------#
#         kmer_counting_in_FASTQ.pl - By Jacob Shreve         #
#           Bioinformatics Core, Purdue University            #
#                          v1.0.0                             #
#-------------------------------------------------------------#

#  The purpose of this script is to read-in a FASTQ file line-
#  by-line and record every kmer found as defined by the user.
#  All kmers are recorded and quantified, with several metrics
#  reported after the FASTQ has been fully transversed.

#  Usage:     perl kmer_counting_in_FASTQ.pl /
#             sample_file.FASTQ
#             kmer

#  Example:   perl kmer_counting_in_FASTQ.pl 01928_1.fastq 7


#---- STEP1: build a hash with all kmers from the reads ----#
$position = 0;
$reads = 0;
open(FASTQ,$ARGV[0]) || die("Can't open the FASTQ file.");
while(<FASTQ>) {
        $position++;
        if($position == 2) {
                $reads++;
                $_ =~ s/\n//;
                for($i=0; $i<=(length($_)-$ARGV[1]); $i++) {
                        $kmer = substr($_, $i, $ARGV[1]);
                        if(not defined $HASH{$kmer}) {
                                $HASH{$kmer} = 1;
                        }
                        else {
                                $HASH{$kmer} += 1;
                        }
                }
        }
        if($position == 4) {
                $position = 0;
        }
}
close(FASTQ);

#---- STEP2: scan through the hash to calculate metrics ----#
$count = 0;
$running_total = 0;
$top_5 = "";
$bottom_5 = "";
foreach (sort { $HASH{$b} <=> $HASH{$a} } keys(%HASH) ) {
        $count++;
        $running_total += $HASH{$_};
        if($count <= 5) {
                $top_5 .= "$HASH{$_}\t$_\n";
        }
        if($count == 10) {
                $at_10 = $running_total;
        }
        if($count == 100) {
                $at_100 = $running_total;
        }
        if($count == 1000) {
                $at_1000 = $running_total;
        }
        if($count == 10000) {
                $at_10000 = $running_total;
        }
}
for($q=-1; $q>=-5; $q--) {
        $bottom_key = (sort { $HASH{$b} <=> $HASH{$a} } keys(%HASH))[$q];
        $bottom_5 .= "$HASH{$bottom_key}\t$bottom_key\n";
}

#---- STEP3: report the metrics ----#
print "Used $0 to find kmer size $ARGV[1] metrics\n";
print "Reads file: $ARGV[0]\n\n";
print "Total number of reads: $reads\n";
print "Total number of kmers: ".(keys %HASH)."\n\n";
print "Top 10 kmers combined:       $at_10/$running_total:  ".sprintf("%.2f", (($at_10/$running_total)*100))."% of all kmers\n";
print "Top 100 kmers combined:      $at_100/$running_total:  ".sprintf("%.2f", (($at_100/$running_total)*100))."% of all kmers\n";
print "Top 1000 kmers combined:     $at_1000/$running_total:  ".sprintf("%.2f", (($at_1000/$running_total)*100))."% of all kmers\n";
print "Top 10000 kmers combined:    $at_10000/$running_total:  ".sprintf("%.2f", (($at_10000/$running_total)*100))."% of all kmers\n";
print "\nMost common 5 kmers:\n$top_5\nLeast common 5 kmers:\n$bottom_5\n";
