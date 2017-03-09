#!/usr/bin/perl -w
use strict;
use Bio::LITE::Taxonomy;
use Bio::LITE::Taxonomy::NCBI;
use Getopt::Long;
use Pod::Usage;
use Data::Dumper;

# local modules
use PhyloStratiphytUtils;

our ( $help, $man, $tax_folder, $db, $blastdbcmd, $user_provided_query_taxon_id, $tax_info,$out );

GetOptions(
    'help' => \$help,
    'man' => \$man,
    'tax_folder=s' => \$tax_folder,
    'db=s' => \$db,
    'query_taxon=s' => \$user_provided_query_taxon_id,
    'tax_info=s' => \$tax_info,
    'out|o=s' => \$out,
) or pod2usage(0);

pod2usage(0) if (defined $help);
pod2usage(-exitstatus => 2, -verbose => 2) if (defined $man);

#### LOAD TAXONOMY DB
print_OUT("Reading taxonomy information");

# load tree of life information, both node' connections and names of nodes.
print_OUT(" '-> Reading phylogenetic tree and species information");
my $nodesfile = $tax_folder . "nodes.dmp";
my $namefile = $tax_folder . "names.dmp";
my $taxNCBI = Bio::LITE::Taxonomy::NCBI->new( db=>"NCBI", names=> $namefile, nodes=>$nodesfile, dict=>"$tax_folder/gi_taxid_prot.bin");
print_OUT(" ... done ...");

my @ql = $taxNCBI->get_taxonomy( $user_provided_query_taxon_id);

my %taxons = ();


if (defined $db){
    
    defined $blastdbcmd or $blastdbcmd = 'blastdbcmd';
    
    my $tax_info = shift;
    print_OUT("Getting sequence ids from sequence db [ $db ]");
    my $ids = `fastacmd -D 2 -d $db`;
    my @seq_ids = split(/\n/,$ids);
    print_OUT(" '-> got [ " . scalar @seq_ids . " ] sequence ids.");
    my $batch_size = 1000;
    for (my $i = 0; $i < scalar @seq_ids; $i += $batch_size) {
open(IDS,">tmp.$$.txt") or die $!;
print IDS join "\n", $seq_ids[ $i .. $i + $batch_size ];
close(IDS);
`$blastdbcmd -outfmt \"%a,%g,%T,%L,%S\" -entry_batch tmp.$$.txt -db $db -out tmp.$$.tax_info_file`;
         unlink("tmp.$$.txt");

print "check now\n";
getc;
    }



# print_OUT("Getting taxons ids for hit sequences from sequence db [ $seq_db ]");
# my $tmp_file = "tmp.seq_tax_ids.$$";
# open(IDS,">$tmp_file.txt") or die $!;
# foreach my $id (@$seq_ids){ print IDS $id,"\n"; }
# close(IDS);
# $tax_info_file = "$out.gi_tax_id.csv";
# print_OUT(" '-> Running blastdbcmd to get sequences information");
# `$blastdbcmd -outfmt \"%a,%g,%T,%L,%S\" -entry_batch $tmp_file.txt -db $seq_db -out $tax_info_file`;
# unlink("$tmp_file.txt");
#
#
# print_OUT(" '-> Parsing sequence information from [ $tax_info_file ]");
# open (TAX_IDS,$tax_info_file ) or die $!;
# my %back_gi_to_taxinfo = ();
# my %target_taxons = ();
# my $taxon_counter = 0;
# my %seen_taxon = ();
# while (my $line = <TAX_IDS>){
# chomp($line);
# my @data = split(/,/,$line);
# $back_gi_to_taxinfo{$data[0]} = {'accession' => $data[0], 'gi' => $data[0],'taxid' => $data[2], 'common_tax_name' => $data[3],'scientific_name' => $data[4]};
# push @{ $target_taxons{ $back_gi_to_taxinfo{$data[0]}->{'taxid'} }->{'seqs'} }, $back_gi_to_taxinfo{$data[0]}->{'accession' };
# if (not exists $seen_taxon{ $back_gi_to_taxinfo{$data[0]}->{'taxid'} } ){
# $seen_taxon{$back_gi_to_taxinfo{$data[0]}->{'taxid'}} = $taxon_counter;
# $taxon_counter++;
# }
# $target_taxons{ $back_gi_to_taxinfo{$data[0]}->{'taxid'} }->{'matrix_number'} = $seen_taxon{$back_gi_to_taxinfo{$data[0]}->{'taxid'}};
# }
# return(\%back_gi_to_taxinfo,\%target_taxons);
}

if (defined $tax_info){
    my %back_gi_to_taxinfo = ();
    my %target_taxons = ();
    my $taxon_counter = 0;
    my %seen_taxon = ();
    open (TAX_IDS,$tax_info ) or die $!;
    while (my $line = <TAX_IDS>){
        chomp($line);
        my ($accn,$gi,$taxid,$common_name,$scientific_name) = split(/,/,$line);
        if (not defined $taxons{ $taxid }){
            $taxons{ $taxid } = {
                'taxid' => $taxid,
                'common_tax_name' => $common_name,
                'scientific_name' => $scientific_name,
                'lineage' => [],
                'sequence_count' => 0,
            };
        }
        $taxons{ $taxid }->{sequence_count}++;
    }
}

my %ancestors = ();
foreach my $taxid (keys %taxons){
    next if ($taxons{ $taxid}->{sequence_count} == 0);
    my @lineage = $taxNCBI->get_taxonomy( $taxid);
    next if ($ql[0] ne $lineage[0]);
    $taxons{ $taxid}->{ lineage } = \@lineage;
    foreach my $level (@lineage){
        $ancestors{ $level } += $taxons{ $taxid}->{sequence_count};
    }
}
open (OUT,">$out") or die $!;
while ( my ($level,$count) = each %ancestors){
    print "$level\t$count\n";
}

print_OUT("Finished");
exit;
