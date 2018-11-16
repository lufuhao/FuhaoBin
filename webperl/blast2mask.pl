#!/usr/local/env perl
#http://code.google.com/p/amyottecode/wiki/blast2mask
# code is not completely modified yet for use


# Blast repeat against masked fasta blast db and mask all sequences of hits
if( system("blastn","-query","$repeat_fasta","-db","$multifasta","-outfmt","6","-evalue","0.00001","-out","temp_rep_vs_genome.BLASTN" ) !=0 )
{ die "Error executing BLASTN of repeat sequence"; }
 
open(REPEAT_FASTA, '>>', $repeat_fasta);
open(TEMP, '<', "temp_rep_vs_genome.BLASTN");
while (<TEMP>)
{
        my $line = $_;
        my @temp = split(/\t/, $line);
        my $contig = $temp[1];
        my $temp_substr = $temp[8];
        my $temp_substp = $temp[9];
        my $mask_start = "";
        my $mask_stop = "";
        my $inverted = "";
        if ($temp_substp < $temp_substr )   { $mask_start = $temp_substp; $mask_stop = $temp_substr; $inverted = "y" }
        else                                { $mask_start = $temp_substr; $mask_stop = $temp_substp; $inverted = "n" }
        my $dna = $seqhash{$contig};
        my $length = ($mask_stop - $mask_start + 1 );
        my $hit = substr($dna, $mask_start -1, $mask_stop - $mask_start + 1);
        my $mask = "";
        my @masked_length = (1 .. $length);
        foreach my $position(@masked_length)
        {
                $position = "n";
                $mask .= $position 
        }
        $dna =~ s/$hit/$mask/;
        $seqhash{$contig} = $dna;
        unless(($contig eq $seed_contig) && ($mask_start == $final_start) && ($mask_stop == $final_stop))
        {   
                if( $inverted eq "n")
                {
                        print REPEAT_FASTA "\n>$contig:$mask_start..$mask_stop\n$hit";
                }
                else
                {
                        $hit =~ tr/atcg/ATCG/;
                        $hit =~ tr/ATCG/TAGC/;
                        my $revComp = reverse($hit);
                        print REPEAT_FASTA "\n>$contig:$mask_stop..$mask_start\n$revComp";
                }
        }
}
 
# print masked genome
$multifasta = "temp_masked_genome.fasta";
open(FASTA, '>', $multifasta);
foreach my $masked_contig(keys %seqhash)
{
        print FASTA ">$masked_contig\n$seqhash{$masked_contig}\n";
}
