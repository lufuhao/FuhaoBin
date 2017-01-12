#!/usr/bin/perl
use strict;
use warnings;
use Bio::Seq;
use Bio::SeqIO;
use Bio::DB::Fasta;

$| = 1;    # Flush output
my $outfile_cds = Bio::SeqIO->new( -format => 'fasta', -file => ">$ARGV[2].cds.fasta" );
my $outfile_pep = Bio::SeqIO->new( -format => 'fasta', -file => ">$ARGV[2].pep.fasta" );



### First, index the genome
my $file_fasta = $ARGV[0];
my $db                   = Bio::DB::Fasta->new($file_fasta);
print ("Genome fasta parsed\n");



### Second, parse the GFF3
my @mRNA;
my $mRNA_name;
my $frame;
open GFF, "<$ARGV[1]" or die $!;
while ( my $line = <GFF> ) {
    chomp $line;
    my @array = split( "\t", $line );
    my $type = $array[2];
    next if ( $type eq 'gene' or $type eq 'UTR' );

    if ( ( $type eq 'mRNA' ) and ( $. > 2 ) ) {
        # Collect CDSs and extract sequence of the previous mRNA
        my $mRNA_seq;
        foreach my $coord (@mRNA) {
            my @cds_coord = split( " ", $coord );
            my $cds_seq = $db->seq( $cds_coord[0], $cds_coord[1], $cds_coord[2] );
            $mRNA_seq .= $cds_seq;
        }

        my $output_nucleotide = Bio::Seq->new(
            -seq        => $mRNA_seq,
            -id         => $mRNA_name,
            -display_id => $mRNA_name,
            -alphabet   => 'dna',
        );
        if ($frame eq '-') {
            $output_nucleotide = $output_nucleotide->revcom();
        }
        my $output_protein = $output_nucleotide->translate();
        $outfile_cds->write_seq($output_nucleotide);
        $outfile_pep->write_seq($output_protein);

        # Now initialize the next mRNA
        my @attrs = split( ";", $array[8] );
        $attrs[0] =~ s/ID=//;
        $mRNA_name = $attrs[0];
        $frame=$array[6];
        @mRNA = (); # Empty the mRNA
    }
    elsif ( $type eq 'mRNA' ) {    # First mRNA
        my @attrs = split( ";", $array[8] );
        $attrs[0] =~ s/ID=//;
        $mRNA_name = $attrs[0];
        $frame=$array[6];
    }
    elsif ( $type eq 'CDS' ) {
        my $cds_coord = $array[0] . " " . $array[3] . " " . $array[4];
        push( @mRNA, $cds_coord );
    }
}

close GFF;
