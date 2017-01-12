#!/usr/bin/perl
use strict;
use warnings;
use File::Basename;
use Getopt::Long;

# set the usage text
my @usage = ("
USAGE: perl $0 --read_file <R1 read file> --vector <vector FASTA> --tags <output FASTA>

Description: Finds vector sequences in reads and attempts to identify 
the number of unique sequences attached to the vector, ie. for BAC splash screening
Output is a count of the number of unique seqeunces and the actual sequences in FASTA format

Mandatory arguments:
--read_file                     STRING          path to read 1 file
--vector                        STRING		vector seqeunce in FASTA format
--tags                          STRING		FASTA file to append sequence tags to

Jon Wright (TGAC) <Jon.Wright\@tgac.ac.uk>

\n");

# define and set command line arguments
my ($read_file, $vector, $tags);

GetOptions ("read_file=s" => \$read_file, "vector=s" => \$vector, "tags=s" => \$tags);

# exit if args are not set
die (@usage) if !($read_file && $vector && $tags);

##########################

# check read file exists
unless (-e $read_file) {
 print "Read file $read_file doesn't exist!\n";
 exit;
} 

# check vector file exists
unless (-e $vector) {
 print "Vector file $vector doesn't exist!\n";
 exit;
}

# get the filename
my $base_fastq = basename($read_file,".gz");
my $base = basename($base_fastq,".fastq");
#print $base_fastq;
#print $base;

# if file is gzipped then unzip it
my $rmv_fastq = 0;		# flag to remove fastq - we only want to remove the fastq if we have created it by unzipping
my $fastq = $read_file;
if (substr ($read_file, -3) eq ".gz") {
  #print "Unzipping read file $read_file...\n";
  my $status = qx(gunzip -c $read_file > $base_fastq);
  $rmv_fastq = 1;
  if ($status) {
    print "Unzip failed!\n";
    exit;
  }
  $fastq=$base_fastq;
}

#print "Operating on fastq file $fastq\n";

# convert the fastq to fasta
my $i=0;
open FASTQ, "<$fastq";
open FASTA, ">$base.fasta";
while(my $line = <FASTQ>){
  if ($line =~ /^\@/ && $i==0){
    $line =~ s/^\@/\>/;
    print FASTA $line;
  } elsif ($i==1) {
    $line =~ s/\./N/g;
    print FASTA $line;
    $i=-3
  }
  $i++;
}
close FASTQ;
close FASTA;

# run the alignments and record the sequences adjacent to the vector in a hash
my %hash_seqs;
#my @exonerate = qx(source exonerate-2.2.0; exonerate $base.fasta $vector --bestn 1 --ryo "%qi:%qab-%qae:%qS:%qal-%ql\t%ti:%tab-%tae:%tS:%tal-%tl\n");
my @exonerate = qx(source exonerate-2.2.0; exonerate $base.fasta $vector --bestn 1 --ryo \"ryo,%qi,%qab,%qae,%ql,%qS,%ti,%tab,%tae,%tl,%tS,%pi\n\");
# ryo format is ryo,read_id,read_aln_start,read_aln_end,read_length,read_strand,vector_id,vector_aln_start,vector_aln_end,vector_length,vector_strand,pid 
foreach my $aln (@exonerate) {
  if ($aln =~ /^ryo,.*/) {
    chomp $aln;		# remove newline
    $aln =~ s/ryo,//;	# remove ryo bit
    my ($read_id,$r_aln_start,$r_aln_end,$r_length,$r_strand,$vector_id,$v_aln_start,$v_aln_end,$v_length,$v_strand,$pid) = split /,/, $aln;
    #print "$read_id,$r_aln_start,$r_aln_end,$r_length,$r_strand,$vector_id,$v_aln_start,$v_aln_end,$v_length,$v_strand,$pid\n";
    my $trimmed_read;
    if ($v_strand eq "+") { # read aligns to forward strand of vector
      # if we can get 20 bases or more sticking off the end of the vector we are ok
      if ($r_aln_start > 20) {
        # grep the read from the fasta (%gs in exonerate only retrieves 70bp)
        my @seq = `grep -A 1 "$read_id" $base.fasta`;
        chomp $seq[1];
        # get the portion of read that is not vector
        #print "bac seq is pos 0 - $r_aln_start in $seq[1]\n";
        # take the first 20 bases next to teh vector and record in a hash
        $trimmed_read = substr($seq[1], $r_aln_start - 20, 20);
        #print $trimmed_read,"\n";
      }
    } else {	# read aligns to reverse strand of vector
      if ($r_aln_end < ($r_length - 20)) {
        my @seq = `grep -A 1 "$read_id" $base.fasta`;
        chomp $seq[1];
        #print "bac seq is pos $r_aln_end to $r_length in $seq[1]\n";
        $trimmed_read = substr($seq[1], $r_aln_end, 20);	# take 20 bases from the end of the read adjacent to teh vector
        # this also needs reverse-complementing as match is on teh reverse strand
        $trimmed_read = reverse_complement($trimmed_read);
        #print $trimmed_read," - revcomp\n";
      }
    }
    # now if we have a trimmed read seq - store it in teh hash 
    if ($trimmed_read) {
      if (exists $hash_seqs{$trimmed_read}) {
        $hash_seqs{$trimmed_read} = $hash_seqs{$trimmed_read} + 1;
      } else {
        $hash_seqs{$trimmed_read} = 1;
      }
    }
  }
}

# print the trimmed reads
#print "Trimmed reads\n";
#while (my ($k,$v)=each %hash_seqs){
#  print "$k $v\n"
#}

# we need to account for sequencing errors so assuming that the BAC sequence that occures teh highest number of times is the 
# correct one we can count all seqeunces that show <3 bases mismatches to this one as the same seqeunce
# ie. the hamming distance is less than 3

if (%hash_seqs) {	
  # open the file to append the unique seqeunces
  open UNIQUE, ">>$tags";  

  # first find the sequence that occurs with the highest frequency
  my @top_hit = sort {$hash_seqs{$b} <=> $hash_seqs{$a}} keys %hash_seqs;
  #print "$top_hit[0] $hash_seqs{$top_hit[0]}\n";

  # declare a new hash and add the top hit
  # this is because we want to first collapse all teh seqeunces that are similar to the sequence that occurs the most frequently, ie. the correct one
  my %it1_hash;
  $it1_hash{$top_hit[0]} = $hash_seqs{$top_hit[0]};
  # now iterate through the hash and refine the counts
  while (my ($k,$v)=each %hash_seqs){
    my $hd = get_hamming_dist($top_hit[0],$k);
    if ($hd == 0) {	# this is the top hit so will already be in the new hash
      #print "top hit so don't do anything\n";
      #$new_hash{$top_hit[0]} = $hash_seqs{$top_hit[0]};
    } elsif ($hd < 3) {	# less than 3 mismatches so increment the top hit entry
      #print "similar so add $v to $top_hit[0] entry ($it1_hash{$top_hit[0]})\n";
      $it1_hash{$top_hit[0]} = $it1_hash{$top_hit[0]} + $v;
    } else {		# hd is larger so this must be another sequence, add to teh hash
      #print "new one $k of $v\n";
      $it1_hash{$k} = $v;
    }
  } 

  # print the hash
  #print "it1 hash\n";
  #while (my ($k,$v)=each %it1_hash){
  #  print "$k $v\n"
  #}

  # now try to collapse the hash again, this time to collapse other sequences together
  my %it2_hash;
  foreach my $k1 ( keys %it1_hash ) {
    if (!%it2_hash) {	# hash has no elements so add one
      $it2_hash{$k1} = $it1_hash{$k1};
    } else {		# hash has elements so compare to what is already saved
      foreach my $k2 ( keys %it2_hash ) {
        my $hd = get_hamming_dist($k1,$k2);
        #print "hd between $k1 and $k2 is $hd\n";
        if ($hd < 3) {   # less than 3 mismatches so add the it1_hash value to the entry in it2_hash
          $it2_hash{$k2} = $it2_hash{$k2} + $it1_hash{$k1};
        } else {              # hd is larger so this must be another sequence, add to it2_hash
          $it2_hash{$k1} = $it1_hash{$k1};
        }
      }
    }
  }

  # print the final hash
  #print "it2 hash\n";
  #print "$fastq\n";
  #while (my ($k,$v)=each %it2_hash){
  #  print "$k $v\n"
  #}

  # now remove any seqeunces we see 2 times or fewer
  foreach my $key (keys %it2_hash) {
    #if ($it2_hash{$key} == 1) {
    if ($it2_hash{$key} < 3) {
      delete($it2_hash{$key});
    }
  }

  # print the seqeuences
  #while (my ($k,$v)=each %it2_hash){
  #  print "$k $v\n";
  #}

  # print teh seqs to a fasta file
  while (my ($k,$v)=each %it2_hash){
    print UNIQUE ">",$fastq,"_",$v,"\n$k\n";
  }

  # print the number of unique seqeunces found bordering the vector for this BAC
  my $count = keys %it2_hash;
  print "$fastq\t$count\n";

  close UNIQUE;
} else {
  print "$fastq\t0\n";
}

# clean up
unlink "$base.fasta" or warn "Could not remove file: $base.fasta\n";
if ($rmv_fastq) {
  unlink "$fastq" or warn "Could not remove file: $fastq\n";
}

#########################
# sub to reverse complement a DNA sequence
sub reverse_complement {
  my $dna = shift;

  # reverse the DNA sequence
  my $revcomp = reverse($dna);

  # complement the reversed DNA sequence
  $revcomp =~ tr/ACGTacgt/TGCAtgca/;
  return $revcomp;
}

# sub to get the hamming distance (edit distance allowing only substitutions) assuming equal length strings
sub get_hamming_dist {

  my ($k,$l) = @_;
  return 0 if $k eq $l; 
                      
  my $len = length ($k);
  my $num_mismatch = 0;
  for (my $i=0; $i<$len; $i++) {
    ++$num_mismatch if substr($k, $i, 1) ne substr($l, $i, 1);
  }
  
  return $num_mismatch;
}

