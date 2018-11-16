#!/usr/bin/env perl
#http://code.google.com/p/amyottecode/w/list
use warnings;
use strict;

my $fasta  = $ARGV[0];
my $In = $ARGV[1];      # Path to blast intput file
my $Out = $ARGV[2];     # Path to gff output file
my $Prog = $ARGV[3];    # descriptive term to descrive program used

# begin with adding info on each contig
open (GFFOUT, ">".$Out) || die "Can not open GFF ouput file.$Out.\n";
print GFFOUT "##gff-version 3\n";
open(FASTA, $fasta);
my $temp;
while(<FASTA>) {
        $temp .= $_;
}

my @contigs = split(/>/, $temp);
shift @contigs;
foreach my $contig(@contigs) {
#       print $contig . "\n";
        my @id_seq = split(/\n/, $contig, 2);
        my $id = $id_seq[0];
        $id =~ s/ .*$//;
        my $seq = $id_seq[1];
        $seq =~ s/\n//g;
        my $length = length $seq;
        print OUT "$id\t.\tcontig\t1\t$length\t.\t.\t.\tID=$id;Name=$id;\n";
}

# open and convert blast ouput to gff3
my $GStart;             
my $GEnd;               
my $HitNum = "0";
open (BLASTIN, "<".$In) || die "Can not open BLAST input file.$In.\n";
print GFFOUT "###\n";
while (<BLASTIN>) {

        next if (/^\#/ or /^\s*$/); # filter comments and empty lines

        $HitNum++;

        my ($QryId, $SubId, $PID, $Len, $MisMatch, $GapOpen, $QStart,$QEnd, $SStart, $SEnd, $EVal, $BitScore) = split(/\t/);
        
        my $Strand;
        my $Frame = ".";
        
        # Set the start to be less then the end
        # This info can be used to deduce the strand
        if ($QStart < $QEnd) {
                $GStart = $QStart;
                $GEnd = $QEnd;
                $Strand = "+";
        }
        elsif ($QStart > $QEnd) {
                $GStart = $QStart;
                $GEnd = $QEnd;
                $Strand = "-";
        } 
        else {
                die "Unexpected Query Start and End:\n\tS:$QStart\n\tE:$QEnd";
        }

        # Trim leading white space from Bit score
        $BitScore =~ s/^\s*(.*?)\s*$/$1/;

        print GFFOUT "$QryId\t$Prog\tmatch\t$GStart\t$GEnd\t$BitScore\t$Strand\t$Frame\tID=$QryId:hsp:$HitNum;Name=$SubId;score=$BitScore;\n";
        print GFFOUT "$QryId\t$Prog\tmatch_part\t$GStart\t$GEnd\t$BitScore\t$Strand\t$Frame\tID=$QryId:hsp:$HitNum;Parent=$QryId:hsp:$HitNum;Name=$SubId;Target=$SubId $SStart $SEnd +;\n";
}
