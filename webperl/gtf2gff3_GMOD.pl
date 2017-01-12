#!/usr/bin/perl
# gtf2gff.pl
#http://iubio.bio.indiana.edu/gmod/GMODTools/GMODTools/bin/gtf2gff.pl
# $modeltype="gene.CDS"; or "gene.mRNA.exon" or ...

my $suf=".gff";
my %renameft = (
'5UTR' => 'five_prime_UTR',
'3UTR' => 'three_prime_UTR',
'transcript' => 'mRNA',
);
my %dropft = (
'start_codon' => 1,'stop_codon' => 1, # these are all subsumed by CDS/UTR ?
);
my %renamea = (
# 'gene_id' => 'Parent',  # for mRNA only
'transcript_id' => 'Parent', # for mRNA kids
);
my %dropa = (
'transcript_id' => 1,
# gene_id => 1 if mrna.kid
);

warn "FIXME: stop_codon add to last CDS (also change to pipe)\n";

foreach my $in (@ARGV) {
  my $out= $in;
  unless( $out =~ s/\.\w+$/$suf/ && $out ne $in) { $out.= $suf; }
  open(IN,$in) or die "open $in";
  open(OUT,">$out") or die "write $out";
 
  print OUT "##gff-version 3\n";
  while(<IN>){
    if(/^#/ && !/##gff/) { print OUT; next; }
    next unless(/^\w/);
    chomp;
    @v=split"\t";
    next if ($dropft{$v[2]});
    $v[2]=$renameft{$v[2]} || $v[2];
    @a=split( /\s*;\s*/, $v[8]);
    foreach (@a) {
      ($k,$v)= split " ",$_,2;
      $k= $renamea{$k} || $k;
      $v=~ s/"//g;
      $_= $dropa{$k} ? "" : "$k=$v";
      }
    $v[8]= join(";",@a);
    print OUT join("\t",@v),"\n";
  }
  close(OUT); close(IN);
}

## this isn't helpful; use plain perl
__END__
use strict;
use Bio::Tools::GFF;
use Bio::FeatureIO;
use Getopt::Long;

my ($output,$input,$format,$type,$help,$cutoff,$sourcetag,$comp,
    $gffver,$match,$quiet);
$format = 'gtf'; # by default
$gffver = 3;
# GTF, is also known as GFF v2.5

GetOptions(
           'i|input:s'  => \$input,
           'o|output:s' => \$output,
           'f|format:s' => \$format,
           'v|version:i'=> \$gffver,
           'q|quiet'    => \$quiet,
           'h|help'     => sub{ exec('perldoc',$0); exit(0) },
           );

# if no input is provided STDIN will be used
my $parser  = new Bio::Tools::GFF(-gff_version => 2.5, -file =>  $input);
#my $parser  = new Bio::FeatureIO(-format => $format, -file   => $input);

my $out;
if( defined $output ) {
  $out = new Bio::Tools::GFF(-gff_version => $gffver, -file => ">$output");
  #$out = new Bio::FeatureIO(-format => 'gff', version => $gffver, -file => ">$output");
} else { 
  $out = new Bio::Tools::GFF(-gff_version => $gffver); # STDOUT
  #$out = new Bio::FeatureIO(-format => 'gff', version => $gffver); # STDOUT
}

while( my $result = $parser->next_feature ) {
  $out->write_feature($result);
  }
__END__
