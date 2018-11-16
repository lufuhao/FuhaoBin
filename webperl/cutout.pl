#! /usr/bin/perl

use warnings;
use strict;

use DBI;
use Bio::Tools::Run::StandAloneBlast;
use Bio::SeqIO;
use Bio::SearchIO;
use Getopt::Long;
use Bio::DB::Fasta;
use SeqFunctions;
use constant USAGE =><<END;

SYNOPSIS:

 cutout.pl [OPTIONS] <file_to_cut_from> [gff_file]

DESCRIPTION:

Cuts out Fasta entries from a sequence file based on GFF lines of features in the sequence entries.

OPTIONS:

       --help
             Prints a help.
       --format <format of input sequence file>
             Can be GenBank, EMBL, Fasta ... Defaults to Fasta
       --outfile <file>
             Specifies output file. STDOUT is default.
       --index
             Creates an index for extraction of subseqs (good for genomes!)
       --reindex
             Forces new creation of index.
       --flanking <int>
             Makes cutout include <int> bp flanking sequence
       --suffix
             Default adds a suffix '.M:<start>:<end>:<strand> to the sequenceid if anything else than
             than ID is used as sequcence id (eg. --idbyseqname or --cutbyparent)
             and no suffix if an ID is given. You can counter this by giving --nosuffix
       --revcompl
             Converts the sequence to its reverse complement if it is on the minus strand
             and puts 'REV_COMPL' in the description line. Default does NOT!
       --coordonly|--sparse
             Allows input as lines with seqname, start, end, strand, [id] seperated either ',' space or tabs.
       --nout <file>
             A gff file specifying features on the mother sequence that are to be replaced by Ns on the
             cut out sequence. A strand specified as '.' is interpreted here as both + and -.
       --forceplus
             Forces cutting out from the plus strand.
       --feature '<featuretag>,<featuretag>...'
             Specifies a list of features to cut out.
       --cutbyregex <regex>
             Cuts out based on id start and end [and strand] returned by a (perl) regex (\$1, \$2, \$3, [\$4])
             Eg. --cutbyparent is equivalent to --cutbyregex 'Sequence (\\S+) (\\d+) (\\d+) ;'
       --idbyseqname
             Generates the IDs on the cut Fastas from the seqname (field one) id
             and not from the id given by the ID tag that is default.
       --idbyregexp '<regexp>'
             Generates the IDs on the cut Fastas from a regexp on the gff line.
       --cutbyparent
             Cut out the corresonding parent seqence based on the info in the attribute field.

EXAMPLES:

 cutout.pl file.fa file.gff > output.fa

 cutout.pl --format EMBL file.embl file.gff > output.fa

 cat file.gff | cutout.pl --index file.fa > output.fa

 cutout.pl --feature 'intron,5utr' --format GenBank file.gb > file.fa

AUTHOR:

Kasper Munch

COPYRIGHT:

This program is free software. You may copy and redistribute it under
the same terms as Perl itself.

END

my $help = 0;
my $index = 0;
my $reindex = 0;
my $flanking = 0;
my $cutbyparent = 0;
my $cutbyregex = 0;
my $idbyseqname = 0;
my $idbyregexp = '';
my $suffix = 1;
my $revcompl = 0;
my $coordonly = 0;
my $format = 'Fasta';
my $feature = '';
my $outfile = '';
my $nout = '';
my $forceplus = '';

GetOptions("help" => \$help,
           "cutbyparent|maskbyparent" => \$cutbyparent,
           "cutbyregex|maskbyregex=s" => \$cutbyregex,
           "index" => \$index,
           "reindex" => \$reindex,
           "idbyregexp=s" => \$idbyregexp,
           "flanking=s" => \$flanking,
           "idbyseqname" => \$idbyseqname,
           "suffix|masksuffix!" => \$suffix,
           "revcompl!" => \$revcompl,
           "sparse|coordonly!" => \$coordonly,
           "format=s" => \$format,
           "nout=s" => \$nout,
           "forceplus" => \$forceplus,
           "feature=s" => \$feature,
           "outfile|output=s" => \$outfile) or die USAGE;

$help and die USAGE;

my $infile = shift @ARGV or die USAGE;
# Unless format is specified, we do some guessing based on the suffix of the infile:
unless ($format) {
  $format = guess_format($infile) or die USAGE;
}

my $in = Bio::SeqIO->newFh(-file => "$infile" , '-format' => $format );

my $out;
if ($outfile = @ARGV ? shift @ARGV : '') {
  open $out, ">$outfile" or die "Couldn't open $outfile: $1\n";
} else {
  open $out, ">&STDOUT" or die "Couldn't copy STDOUT filehandle\n";
}

my $incrementid;
my %cutout;

# If --nout is specified, get all the gff lines into a hash of arrays:
my %nout = ();
if ($nout) {
  open my $fh, "$nout" or die "$nout: $!\n";;
  while (my $gff = <$fh>) {
    my $seqname = (split /\t/, $gff)[0];
    push @{$nout{$seqname}}, $gff;
  }
}

if ($feature) {
  my $regex =  join("|", split(/,/, $feature));
  while (my $seq = <$in>) {
    my $id = $seq->id();

    # Get the features included in the list from the seq:
    foreach my $f ($seq->get_SeqFeatures()) {
      # Check if the prim tag matches regex:
      next unless $f->primary_tag() =~ $regex;
      if (ref($f->location()) eq 'Bio::Location::Split') {

        # Get the array of Bio::Location::Simple objs.:
        my @sl = $f->location->sub_Location();
        for (my $j = 0; $j < @sl; $j++) {
#          my $cutid = "$sl[$j]->{seqname}.M:$sl[$j]->{start}:$sl[$j]->{end}:$sl[$j]->{strand}";
          my $cutid = sprintf "$id.M:%s:%s:%s", $sl[$j]->start, $sl[$j]->end, $sl[$j]->strand;
          my $gff = join "\t", $seq->id, '.', '.', $sl[$j]->start, $sl[$j]->end, '.', $sl[$j]->strand, '.';
          push @{$cutout{$id}}, {start => $sl[$j]->start, end => $sl[$j]->end,
                                  strand => $sl[$j]->strand, cutid => $cutid, gffstr => $gff};
        }
      } else {
        my $cutid = $seq->id . ".M:" . $f->location->start. ":" . $f->location->end . ":" . $f->location->strand;
        my $gff = join "\t", $seq->id, '.', '.', $f->location->start, $f->location->end, '.', $f->location->strand, '.';
        push @{$cutout{$id}}, {start => $f->location->start, end => $f->location->end,
                                strand => $f->location->strand, cutid => $cutid, gffstr => $gff};
      }
    }
    for my $s (@{$cutout{$id}}) {
      my $start = $s->{start} - $flanking;
      my $end = $s->{end} + $flanking;
      my $substr = $seq->subseq($start,$end);

      nout($id, $start, $end, \$substr) if $nout;
      printFasta($substr, $s, $out);
      }
  }
}
else {
  # Make a hash of arrays of hashes each holding start, end, string and the GFF line.
  my $cut;
  if (my $cutfile = @ARGV ? shift @ARGV : '') {
    open $cut, "$cutfile" or die "$cutfile: $!\n";
  } else {
    open $cut, ">&STDIN" or die "Couldn't copy STDIN filehandle\n";
  }

  while (my $gff = <$cut>) {

    next if $gff =~ /^#/;

    # If there are no digits in the first row then we probably forgot
    # to remove the header from the gff input:
    if ($. == 1 && $gff !~ /\d+/) {
      warn "Removing table header\n";
      next;
    }

    next unless my $f = hashgff($gff);

    my $cutid = 0;
    my $id = 0;
    my $start = 0;
    my $end = 0;
    my $strand = 0;

    if ($cutbyparent) {
      ($id, $start, $end, $strand) = ($f->{parent}, $f->{parent_start}, $f->{parent_end}, '+');
      $cutid = "$f->{parent}.M:$f->{parent_start}:$f->{parent_end}:1";
    }
    if ($cutbyregex) {
      $gff =~ /$cutbyregex/;
      ($id, $start, $end, $strand) = ($1, $2, $3, $4);
      my $strand = $strand eq '+' ? '1' : $strand eq '-' ? '-1' : '?';
      $cutid = "$id.M:$start:$end:$strand";
    }
    if ($idbyseqname) {
      my $strand = $f->{strand} eq '+' ? '1' : $f->{strand} eq '-' ? '-1' : '?';
      $cutid = "$f->{seqname}.M:$f->{start}:$f->{end}:$strand";
    }
    if ($idbyregexp) {
      $gff =~ /$idbyregexp/;
      $cutid = $1 ? $1 : $&;
    }

    $id ||= $f->{seqname};
    $start ||= $f->{start};
    $end ||= $f->{end};
    $strand ||= $f->{strand};

    $cutid ||= "$f->{ID}" if defined($f->{ID});

    $cutid = ++$incrementid unless $cutid;

    $strand = '+' if $forceplus;

    $cutout{$id} ||= ();
    push @{$cutout{$id}}, {start => $start, end => $end,
                            strand => $strand, cutid => $cutid, gffstr => $gff};
  }

  # Write the cut seqs as Fasta to a file:
  if ($index or $reindex) {

    $feature && die "option --feature incompatible with --index and --reindex\n";

    unless ( -e "$infile.index") {
      $reindex = 1;
    }

    sub indexid {
      my $description = shift;
      if (/^>gi\|[^\|]+\|\w+\|([^\|]+)\|/) {
        $description = $1;
      } else {
        $description =~ /^>(\S+)/;
        $description = $1;
      }
      return $description;
    }

    my $fileindex = Bio::DB::Fasta->new($infile, -makeid => \&indexid, -reindex => $reindex);

    for my $seqid (%cutout) {
      for my $s (@{$cutout{$seqid}}) {
        my $start = $s->{start} - $flanking;
        my $end = $s->{end} + $flanking;
        my $substr = $fileindex->seq($seqid, $start, $end)
          or die "Can't find sequence id: '$seqid' in index. Consider the --reindex option.\n";

        nout($seqid, $start, $end, \$substr) if $nout;

        printFasta($substr, $s, $out);
      }
    }
  } else {
    while (my $seq = <$in>) {
      my $seqid = $index ? $seq : $seq->id();
      # If a cut file of gff strings is given, then find the proper substring(s) and write tmpfile(s):
      if (defined( $cutout{$seqid}) ) {

        for my $s (@{$cutout{$seqid}}) {
          my $start = $s->{start} - $flanking;
          my $end = $s->{end} + $flanking;
          my $substr = $seq->subseq($start,$end);

          nout($seqid, $start, $end, \$substr) if $nout;
          printFasta($substr, $s, $out);
        }
      }
    }
  }
}

sub printFasta {
  (my $substr, my $s, my $out) = @_;

  my $revcompl_suffix = '';
  if ($s->{strand} eq '-' && $revcompl) {
    $substr = rev_compl($substr);
    $revcompl_suffix = ' REV_COMPL';
  }
  $substr =~ tr/a-z/A-Z/;
  if (!$suffix) {
    $s->{cutid} =~ s/\.M.*$//;
  }
  print $out ">$s->{cutid}$revcompl_suffix\n";
  my $start = 0;
  my $length = length($substr);
  my $p = int($length / 70);
  while ($p--) {
    print $out substr($substr, $start, 70),"\n";;
    $start += 70;
  }
  print $out substr($substr, $start),"\n" unless $start == $length;

}

# Returns the reverse complementary strand
# sub rev_compl {
#   my $string = shift;
#   my @compl;
#   # Complement it:
#   $string =~ tr/AGCTN/TCGAN/;
#   # Reverse it:
#   my @string = split //, $string;
#   while (@string) {
#     push @compl, (pop @string);
#   }
#   return join "", @compl;
# }
sub rev_compl {
  my $string = shift;
  my @compl;
  my @string = split //, $string;

  while (my $base = pop @string) {
  SWITCH: {
      if ($base =~ /A/i) {
        push @compl, ('T'); last SWITCH;
      }
      if ($base =~ /T/i) {
        push @compl, ('A'); last SWITCH;
      }
      if ($base =~ /C/i) {
        push @compl, ('G'); last SWITCH;
      }
      if ($base =~ /G/i) {
        push @compl, ('C'); last SWITCH;
      }
      if ($base =~ /N/i) {
        push @compl, ('N'); last SWITCH;
      }
      print STDERR "Sequence letter not recognized\n";
    }
  }
  return join "", @compl;
}

sub guess_format {
  my $s = shift @_;
  $s =~ s/\.(.+)$/$1/;
  my $failed = 0;
  my $format;
 SW: {
    if ($s =~ /(^fasta)|(^fast)|(^fst)|(^fsa)|(^ft)|(^fs)/i) {$format = 'Fasta'; last SW};
    if ($s =~ /(lfasta)|(lfast)|(lfst)|(lfsa)|(lft)|(lfs)/i) {$format = 'LabeledFasta'; last SW};
    if ($s =~ /(embl)|(emb)|(em)|(eml)/i) {$format = 'EMBL'; last SW};
    if ($s =~ /(genebank)|(genbank)|(genb)|(geneb)|(gbank)|(gb)/i) {$format = 'GenBank'; last SW};
    if ($s =~ /(swissprot)|(sprt)|(swissp)|(sprot)|(sp)|(spr)/i) {$format = 'Swissprot'; last SW};
    if ($s =~ /pir/i) {$format = 'PIR'; last SW};
    if ($s =~ /gcg/i) {$format = 'GCG'; last SW};
    if ($s =~ /scf/i) {$format = 'SCF'; last SW};
    if ($s =~ /ace/i) {$format = 'Ace'; last SW};
    if ($s =~ /phd/i) {$format = 'phd'; last SW};
    if ($s =~ /phred/i) {$format = 'phred'; last SW};
    if ($s =~ /raw/i) {$format = 'raw'; last SW};
    $failed++;
  }
  return eval{$failed ? 0 : $format};
}

sub hashgff {
  my $l = shift @_;
  chomp $l;
  return 0 if $l =~ /(^\#)|(^\s*$)/;

  $l = coord2gff($l) if $coordonly;

  my @f = split(/\t/, $l);
  my %gff = (seqname => $f[0],
             source => $f[1],
             feature => $f[2],
             start => $f[3],
             end => $f[4],
             score => $f[5],
             strand => $f[6],
             frame => $f[7]);
  my @a =  split(/ *; */, $f[8]) if $f[8];
  if (@a) {
    my $parentseq = shift @a;
    $parentseq =~ /^\S+\s+(\S+)/;
    $gff{parent} = $1;
    if ($parentseq =~ /^\S+\s+\S+\s+(\d+)\s+(\d+)/) {
      $gff{parent_start} = $1;
      $gff{parent_end} = $2;
    }
    for my $a (@a) {
      (my $key, my $value) = split /\s+/, $a, 2;
      $gff{$key} = $value;
    }
  }
  if (
      $gff{seqname} and
      $gff{source} and
      $gff{feature} and
      $gff{start} =~ /\d+/ and
      $gff{end} =~ /\d+/ and
      $gff{score} and
      $gff{strand} =~ /-|\+|\./ and
      defined $gff{frame}
     ) {
    return \%gff;
  } else {
    die "Couldn't make sense of the GFF input.\n";;
  }
}

sub coord2gff {
  my $l = shift;
  my @r;
  if ($l =~ /,/) {
    @r = split /,/, $l;
  } elsif ($l =~ /\t/) {
    @r = split /\t/, $l;
  } elsif ($l =~ / /) {
    @r = split / /, $l;
  }
  @r < 4 and die "Coordinate entry contains less than the required four entries.\n";

  my $strand = $r[3] =~ /(-1)|(-)/ ? '-' : $r[3] =~ /(1)|(\+)/ ? '+' : die "Couldn't esablish strand.\n";
  my $id = defined($r[4]) ? " ; ID $r[4]" : '';

  return "$r[0]\t.\t.\t$r[1]\t$r[2]\t.\t$strand\t.\tSequence .$id";
}

sub nout {
  my ($seqid, $start, $end, $substr) = @_;

  foreach my $line (@{$nout{$seqid}}) {
    my ($nout_start, $nout_end) = (split /\t/, $line)[3,4];
    # Unless the nout feature doen't overlap with the cute out feature:
    my $offset = undef;
    my $length = undef;
    if ($nout_end < $start || $nout_start > $end) {
      # Not overlapping.
      next;
    } elsif ($nout_start >= $start && $nout_end <= $end) {
      # Both start end end within segment.
      $offset = $nout_start - $start;
      $length = $nout_end - $nout_start + 1;
    } elsif ($nout_end > $end) {
      # Only start is within segment.
      $offset = $nout_start - $start;
      $length = $end - $start + 1 - $offset;
    } elsif ($nout_start < $start) {
      # Only end is within segment.
      $offset = 0;
      $length = $nout_end - $start + 1;
    } else {
      die "Something is wrong at line ", __LINE__, ".\n";
    }

    $end - $start + 1 >= $offset + $length or die "substring was is passed coordinates outside the substring. exiting...\n";

    # TEST LIGE AT DET HER ER RIGTIGT (AT DER IKKE ER NOGLE +1/-1 FEJL ET STED, DET ER DER NOK.)

    substr($$substr,$offset,$length) =~ tr/ATGCatgc/NNNNnnnn/;
  }
}

=head1 SYNOPSIS:

 cutout.pl [OPTIONS] <file_to_cut_from> [gff_file]

=head1 DESCRIPTION:

Cuts out Fasta entries from a sequence file based on GFF lines of features in the sequence entries.

=head1 OPTIONS:

=over 4

=item --help

Prints a help.

=item --format <format of input sequence file>

Can be GenBank, EMBL, Fasta ... Defaults to Fasta

=item --outfile <file>

Specifies output file. STDOUT is default.

=item --index

Creates an index for extraction of subseqs (good for genomes!)

=item --reindex

Forces new creation of index.

=item --flanking <int>

Makes cutout include <int> bp flanking sequence

=item --suffix

Default adds a suffix '.M:<start>:<end>:<strand> to the sequenceid if anything else than
than ID is used as sequcence id (eg. --idbyseqname or --cutbyparent)
and no suffix if an ID is given. You can counter this by giving --nosuffix

=item --revcompl

Converts the sequence to its reverse complement if it is on the minus strand
and puts 'REV_COMPL' in the description line. Default does NOT!

=item --coordonly|--sparse

Allows input as lines with seqname, start, end, strand, [id] seperated either ',' space or tabs.

=item --nout <file>

A gff file specifying features on the mother sequence that are to be replaced by Ns on the
cut out sequence. A strand specified as '.' is interpreted here as both + and -.

=item --forceplus

Forces cutting out from the plus strand.

=item --feature '<featuretag>,<featuretag>...'

Specifies a list of features to cut out.

=item --cutbyregex <regex>

Cuts out based on id start and end [and strand] returned by a (perl) regex ($1, $2, $3, [$4])
Eg. --cutbyparent is equivalent to --cutbyregex 'Sequence (\\S+) (\\d+) (\\d+) ;'

=item --idbyseqname

Generates the IDs on the cut Fastas from the seqname (field one) id
and not from the id given by the ID tag that is default.

=item --idbyregexp '<regexp>'

Generates the IDs on the cut Fastas from a regexp on the gff line.

=item --cutbyparent

Cut out the corresonding parent seqence based on the info in the attribute field.

=back

=head1 EXAMPLES:

 cutout.pl file.fa file.gff > output.fa

 cutout.pl --format EMBL file.embl file.gff > output.fa

 cat file.gff | cutout.pl --index file.fa > output.fa

 cutout.pl --feature 'intron,5utr' --format GenBank file.gb > file.fa

=head1 AUTHOR:

Kasper Munch

=head1 COPYRIGHT:

This program is free software. You may copy and redistribute it under
the same terms as Perl itself.


=cut
