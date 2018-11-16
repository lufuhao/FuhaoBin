#!/usr/bin/env perl
#http://code.google.com/p/amyottecode/wiki/blast2kingdom_taxonomy
use strict;
use Getopt::Long;
use Bio::DB::Taxonomy;
use Bio::SearchIO;
use DBI;
use Env;
use File::Spec;
#use vars qw($SEP);
#$SEP = '_';
my $DEBUG = 0;

# set default for variables
my $evalue_filter = 1e-3;
my $file;
my $zcat = 'zcat';                              # or gunzip -c
my $prefix = undef;                             # path to nodes.dmp and names.dmp files
my $gi2taxidfile = undef;               # needs to be full path and file name
my $top_level_ids_string = '';  # comma separated list of extra IDs to look for, other than kingdom and superkingdom.

GetOptions
(
        'i|in:s' => \$file,
        't|taxonomy:s'=> \$prefix,
        'g|gi|gi2taxid:s' => \$gi2taxidfile,
        'l|top_level_ids:s' => \$top_level_ids_string,
        'e|evalue:f' => \$evalue_filter,
        'z|zca:s' => \$zcat,
        'v|verbose|debug' => \$DEBUG,
        'h|help' => sub { system('perldoc', $0); exit },
);

# comma separated to list of IDs
my @top_level_ids_array = split ',', $top_level_ids_string;
my %top_level_ids_hash = map { $_ => 1 } @top_level_ids_array;

# Ensure idx location is created
mkdir(File::Spec->catfile($prefix,'idx')) unless -d File::Spec->catfile($prefix,'idx');

# these files came from ftp://ftp.ncbi.nih.gov/pub/taxonomy
my $taxdb = Bio::DB::Taxonomy->new
(
        -source    => 'flatfile',
        -directory => File::Spec->catfile($prefix,'idx'),
        -nodesfile => File::Spec->catfile($prefix,'nodes.dmp'),
        -namesfile => File::Spec->catfile($prefix,'names.dmp')
);

# Create an sqlite database to hold taxonomy information in a more useful format.
my $giidxfile = File::Spec->catfile($prefix,'idx','gi2taxid');
my $done = -e $giidxfile;
my $dbh = DBI->connect("dbi:SQLite:dbname=$giidxfile","","");
# if the db does not exist, create it
if( ! $done ) {
        $dbh->do("CREATE TABLE gi2taxid ( gi integer PRIMARY KEY,taxid integer NOT NULL)");
        $dbh->{AutoCommit} = 0;
        my $fh;
        # this file came from ftp://ftp.ncbi.nih.gov/pub/taxonomy
        # I'm interested in protein hits therefore the _prot file.
        warn "Using gi to taxonomy ID database file `$gi2taxidfile'\n";
        if( $gi2taxidfile =~ /\.gz$/ ) {
                open($fh, "$zcat $gi2taxidfile |" ) || die "$zcat $gi2taxidfile: $!";
        } else {
                open($fh, $gi2taxidfile ) || die $!;
        }
        my $i = 0;
        my $sth = $dbh->prepare("INSERT INTO gi2taxid (gi,taxid) VALUES (?,?)");
        while(<$fh>) {
                my ($gi,$taxid) = split;
                $sth->execute($gi,$taxid);
                $i++;
                if( $i % 500000 == 0 ) {
                        $dbh->commit;
                }
        }
        $dbh->commit;
        $sth->finish;
}
warn "Finished SQlite database.\n";

my %query;
my %gi2node;

my $sth = $dbh->prepare("SELECT taxid FROM gi2taxid WHERE gi=?");

my $searchio = Bio::SearchIO->new(
        -file => $file,
        -format => 'blastxml'
) or die "parse failed";

# Run through all of the hits from the BLAST output file provided.
while( my $result = $searchio->next_result()) {
        # go through all hits to each query
        while( my $hit = $result->next_hit ) {
                # gather information about hsp for every hit, then retrieve taxonomy information
            my $hsp    = $hit->next_hsp; 
                my $qname  = $result->query_name;
                my $hname  = $hit->name;
                my $evalue = $hsp->evalue;                       
        
                # only go forward with hsp that are below the cutoff value
                next if( $evalue > $evalue_filter );

                # check if subject is following a gi naming scheme
                if( $hname =~ /gi\|(\d+)/) {
                        my $gi = $1;    
                        my $gi_genus;
                        # see if we cached the results from before
                        if( ! $gi2node{$gi} ){ 
                                # query to SQLite db to find taxid 
                                $sth->execute($gi);
                                my $taxid;
                                $sth->bind_columns(\$taxid);
                                if( ! $sth->fetch ) {
                                        warn("no taxid for $gi\n");
                                        next;
                                }
                                my $node = $taxdb->get_Taxonomy_Node($taxid);
                                if( ! $node ) {
                                        warn("cannot find node for gi=$gi ($hname) (taxid=$taxid)\n");
                                        next;
                                }
                                my $parent = $taxdb->get_Taxonomy_Node($node->parent_id);

                                # THIS IS WHERE THE KINGDOM DECISION IS MADE
                                # loops through all node_names until it finds the taxonomy ceiling desired, then finishes loop
                                while( defined $parent && $parent->node_name ne 'root' ) {
                                        # this is walking up the taxonomy hierarchy
                                        # can be a little slow, but works...
                                        # deal with Eubacteria, Archea separate from
                                        # Metazoa, Fungi, Viriplantae differently
                                        # (everything else Eukaryotic goes in Eukaryota)
                                        if( $parent->rank eq 'genus') {         
                                                $gi_genus = $parent->node_name;
                                        }
                                        if( $parent->rank eq 'kingdom') {
                                                $gi2node{$gi} = $parent->node_name;
                                                last;
                                        } elsif( $parent->rank eq 'superkingdom' ) {
                                                $gi2node{$gi} = $parent->node_name;
                                                $gi2node{$gi} =~ s/ \<(bacteria|archaea)\>//g;
                                                last;
                                        } elsif (exists($top_level_ids_hash{$parent->id})) { #user-specified taxonomy ceiling found?
                                                $gi2node{$gi} = $parent->node_name;
                                                last;
                                        }
                                        $parent = $taxdb->get_Taxonomy_Node($parent->parent_id);
                                }
                        }

                print($qname."\t".$hname."\t".$gi2node{$gi}."\t".$gi_genus."\n");

                        # keep track of counts for each kingdom type
                        my $kingdom = $gi2node{$gi};
                        unless( defined $kingdom && length($kingdom) ) {
                                warn("no kingdom for $hname\n");
                        } else {
                                $query{$kingdom} ++;
                        }

                } else { 
                        # if query does not have GI name, simply go to the next hsp
                        warn "No GI in $hname\n"; 
                }
        }
}

# finished querying the taxon sqlite db
$sth->finish;

# print out total number of hits to each taxonomic category
foreach my $kingdom_keys (keys %query) {
        my $total = scalar $query{$kingdom_keys};
        warn "$kingdom_keys\ttotal=$total\n";
}


=head1 NAME

retrieve taxonomic information for blast hits with gi_codes

=head2 USAGE

classify_hits_kingdom [-i tab_file] [-i second_BLAST_file] [-e evalue_cutoff]
[-t dir_where_TAXONOMY_files_are] [-g gi2taxid]
[-z PATH_TO_zcat] [-v]

=head2 DESCRIPTION

Will print out the taxonomic distribution (at the kingdom level) for a
set of hits against the NR database. This script assumes you've done
a search against the protein database, you'll have to make minor
changes in the gi_taxid part to point to the gi_taxid_nuc.dump file.
The gi_taxid files and nodes.dmp (part of taxdump.tar.gz) can be
downloaded from ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/.

Input values:
-t/--taxonomy directory where the taxonomy .dmp files are (from NCBI)
-g/--gi Location of gi_taxid_prot.dmp (or gi_taxid_nucl.dmp if
the search was against a NT db)
-i/--in The name of the tab delimited -m8/-m9 output files to
process.

-e/--evalue Provide an E-value cutoff for hits to be considered
-z/--zcat Path to the 'zcat' executable, can also be 'gunzip -c'
if no zcat on your system.
Flags
-v/--verbose To turn on verbose messages
-h/--help Display this helpful information

This is intended to be useful starting script, but users may want to
customize the output and parameters. Note that I am summarizing the
kingdoms here and Eukaryota not falling into Metazoa, Viridiplantae,
or Fungi gets grouped into the general superkingdom Eukaryota. for
simplicity. There are comments in the code directing you to where
changes can be made if you wanted to display hits by phylum for
example. Note that you must wipe out the cache file 'gi2class' that
is craeed in your directory after making these changes.

=head2 AUTHOR

original by Jason Stajich jason_at_bioperl_dot_org,
modified by Stefan Amyotte

=cut


