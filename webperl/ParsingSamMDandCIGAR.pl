#!/usr/bin/env perl
use Bio::DB::Sam;
use warnings;
use strict;
#http://bioinfomative.blogspot.co.uk/2012/07/parsing-cigar-and-md-fields-of-sambam.html
# Prototypes
sub prettyPrint($);
sub getAuxValue($$);
sub splitCigar($);
sub returnVar($);

# globals
my $bamFileName = $ARGV[1];
my $refFileName = $ARGV[0];

# high level API
my $sam = Bio::DB::Sam->new(-bam  =>$bamFileName, -fasta=>$refFileName);

#my @alignments = $sam->features();
my @alignments = $sam->get_features_by_location(-seq_id => 'ParMerge_MblContig65',
                                                 -start  => 1550,
                                                 -end    => 1650,);



for my $a (@alignments) {
  if(my $variants = returnVar($a) ) {
    print $variants;
    prettyPrint($a);
  }
}

1;

# Given a Bio::DB::Bam::Alignment, print it pretty like
sub prettyPrint($) {
  my $a = shift;
  my $query_start = $a->query->start;     
  my $query_end   = $a->query->end;
  my $ref_dna   = $a->dna;        # reference sequence bases
    my $query_dna = $a->query->dna; # query sequence bases
    my $cigar  = $a->cigar_str;
  my $MD_value = getAuxValue($a, "MD");
  my $query_stream;
  my $ref_stream;
  my $align_stream;
  my @cigarOperations = @{splitCigar($cigar)};
  my $query_pos = 0;
  my $ref_pos = 0;
  my $charPerLine = 50;



# print meta information
  print "Cigar: $cigar\n";
  print "MD: $MD_value\n";

  foreach my $operation (@cigarOperations) {
    my $op_length = $operation->[0];
    my $op_type = $operation->[1];

    if($op_type =~ /^S$/) {
      for(my $op_pos = 0; $op_pos < $op_length; $op_pos++) {
        $ref_stream .= "-";
        $query_stream .= substr ($query_dna, $query_pos++, 1);
        $align_stream .= " ";
      } 
    }
    elsif($op_type =~ /^M$/) {
      for(my $op_pos = 0; $op_pos < $op_length; $op_pos++) {
        my $ref_char = substr ($ref_dna, $ref_pos++, 1);
        my $query_char = substr ($query_dna, $query_pos++, 1);
        $ref_stream .= $ref_char;
        $query_stream .= $query_char;
        if ($ref_char eq $query_char)  {
          $align_stream.="|" 
        }
        else {
          $align_stream .=".";
        }
      }
    }
    elsif($op_type =~ /^D$/) {
      for(my $op_pos = 0; $op_pos < $op_length; $op_pos++) {
        $ref_stream .= substr ($ref_dna, $ref_pos++, 1);
        $query_stream .= "-";
        $align_stream .= " ";
      }
    }
    elsif($op_type =~ /^I$/) {
      for(my $op_pos = 0; $op_pos < $op_length; $op_pos++) {
        $ref_stream .= ".";
        $query_stream .= substr ($query_dna, $query_pos++, 1);
        $align_stream .= " ";
      }
    }
  }

  for(my $stringPos = 0; $stringPos < (length($ref_stream) - length($ref_stream) % $charPerLine); $stringPos += $charPerLine) {
    print join "\n", substr($ref_stream,$stringPos,$charPerLine) . " " . ($stringPos+$charPerLine), 
          substr($align_stream, $stringPos,$charPerLine),
          substr($query_stream, $stringPos,$charPerLine),
          "\n";
  }

  print join "\n", substr($ref_stream,length($ref_stream) - length($ref_stream) % $charPerLine),
        substr($align_stream,length($ref_stream) - length($ref_stream) % $charPerLine),
        substr($query_stream,length($ref_stream) - length($ref_stream) % $charPerLine),
        "\n";
}



# Given a Bio::DB::Bam::Alignment, return the value associated with the given auxillary string
sub getAuxValue($$) {
  my $a = shift;
  my $requestedField = shift;
  my $aux_values = $a->aux;
  my @auxStrings = split /\t/, $aux_values;
  my %hash;
  foreach my $auxString (@auxStrings) {
    my ($field, $type, $value) = split /:/, $auxString;
    $hash{$field} = $value
  }
  (exists $hash{$requestedField})? return $hash{$requestedField}: 0;
}

# Given a Cigar string return a double array
# such that each sub array contiains the [ length_of_operation, Cigar_operation]
sub splitCigar($) {
  my $cigar_string = shift;
  my @returnable;
  my (@matches) = ($cigar_string =~ /(\d+\w)/g);
  foreach (@matches) {
    my @operation = ($_ =~ /(\d+)(\w)/);
    push @returnable, \@operation;
  }
  return \@returnable;
}

# Given a Bio::DB::Bam object, print out the variants
sub returnVar($) {
  my $a = shift;
  my $query_seq = $a->query->dna;
  my $cigarOperations = splitCigar($a->cigar_str);
  my $mdString = getAuxValue($a, "MD");
  my $refPos = $a->start;
  my $seqPos = 0;
  my @variants;
  my $returnString;


  foreach my $operation (@$cigarOperations) {
    my $cig_length = $operation->[0];
    my $cig_op = $operation->[1];

    if ($cig_op =~ /^D$/) {
      $returnString .= "Deletion, $refPos, $cig_length\n";
      $refPos+=$cig_length;
    }
    elsif($cig_op =~ /^I$/) {
      my $insertedBases = substr($query_seq, $seqPos, $cig_length);
      $returnString .= "Insertion, $refPos, $insertedBases\n";
      $seqPos+=$cig_length;

    }
    elsif($cig_op =~ /^M$/) {
      $refPos += $cig_length;
      $seqPos += $cig_length;
    }
    elsif($cig_op =~ /^S$/) {
      $seqPos+=$cig_length;
# Don't increment refPos
    }
    else {
      die ("Unrecognized Cigar State: $cig_op\n");
    }

  }

  $refPos = $a->start;

# Doing this like a stream would be helpful because it would mirror the c string case

# create a buffer variable
  my $matchedBaseBuffer;
  my $mdLength = length($mdString);
# Read until you hit a special character
  for(my $mdIndex = 0;$mdIndex < $mdLength; $mdIndex++) {
# Check if this character is a special character
    if(substr($mdString,$mdIndex,1) =~/[ACTG]/) {
# increment refPos by the number of matched bases
      $refPos += int($matchedBaseBuffer);
      $matchedBaseBuffer = '';
      $returnString .= "Mismatch, $refPos, ".substr($mdString,$mdIndex,1)."\n";
      $refPos++;
    }
# Check for deletion case
    elsif(substr($mdString,$mdIndex,1) =~/[\^]/) {
      $refPos += int($matchedBaseBuffer);
      $matchedBaseBuffer = '';
      while (substr($mdString,($mdIndex+1),1) =~/[ACTG]/ && $mdIndex < $mdLength) {
        $refPos++;
        $mdIndex++;
      }
    }
# Last character
    elsif($mdIndex == ($mdLength -1)) {
# May not actually have to do anything here..
    }
    else {
      $matchedBaseBuffer.=substr($mdString,$mdIndex,1);
    } 
  }


  return $returnString; 
}
