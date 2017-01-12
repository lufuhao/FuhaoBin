#!/usr/bin/perl -w
use strict;
use warnings;
use constant USAGE =><<EOH;
Descriptions:
	For each mRNA declared in -gff, write to STDOUT its spliced CDS sequence taken from -fasta

Requirements: 
	Bio::SeqIO, Bio::DB::SeqFeature::Store, Bio::SeqFeatureI

Usage: 
	$0 -fasta genome.fa -gff xx.gff

AUTHOR
	malcolm.cook\@gmail.com
EOH
use File::Basename;
use Bio::SeqIO;
use Bio::Seq;
use Bio::DB::SeqFeature::Store;
use Bio::SeqFeatureI;

my $fasta = $ARGV[0];
my $outfile = basename($fasta, ".fa") . "_cds.fa"; #formats name of infile for outfile

my $gff = $ARGV[1];

open (OUT, ">$outfile") or die "Can't open file for writing: $!"; 

my $db = Bio::DB::SeqFeature::Store->new(    -adaptor => 'memory',
                                            -fasta => $fasta,
                                            -gff => $gff) or die $!; # This is where gff and fasta are loaded

my $ss = $db->get_seq_stream(-type => 'gene'); #looks at gff, find genes and their associated information

while (my $gene = $ss->next_seq) #while looking at each gene
{
    print OUT ">" , $gene -> display_name , "\n"; #fasta id for CDS

    #print OUT $gene -> start , "\t", $gene -> end , "\t", $gene -> strand, "\n"; #test to see if script was reading gff correctly

    my @segments = $gene -> segments(-type => 'CDS');

    my %START;
    my %END;
    foreach my $segment(@segments)
    {
        $START{$segment} = $segment -> start; #create hash of CDS id and start nucleotide
        $END{$segment} = $segment -> end; #create hash of CDS id and end nucleotide
    }

    my @CDS;
    foreach my $segment (sort{$START{$a} <=> $START{$b}} keys %START) #sorts CDS segments by order of starting nucleotide
    {
        my $cds_for_seq = $db->fetch_sequence($gene -> seq_id , $START{$segment}, $END{$segment}); #extracts CDS segment sequence
        push(@CDS, $cds_for_seq); #puts CDS segment sequence to end of array

        #print OUT $segment , "\t" , $START{$segment} , "\t", $END{$segment} ,"\n"; #test to see if CDS segments are ordered correctly
        #print OUT $cds_for_seq , "\n";
    }

    my $CDS_seq = join('', @CDS); #concatenate array of CDS segment sequences

    my $for_obj = Bio::Seq -> new(-seq => $CDS_seq , -alphabet => 'dna'); #load sequence string as sequence object
    my $rev_obj = $for_obj -> revcom; #reverse complement sequence object
    my $rev_seq = $rev_obj -> seq; #convert sequence object to sequence string

    if ($gene -> strand == +1) #check if gene was on plus strand
    {
        print OUT $CDS_seq , "\n";
    }
    elsif ($gene -> strand == -1) #check if gene was on minus strand
    {
        print OUT $rev_seq , "\n";
    }
}

exit;
