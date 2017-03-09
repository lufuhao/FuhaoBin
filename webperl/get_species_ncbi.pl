#!/usr/bin/perl -w

=pod

=head1 NAME

  get_species_ncbi.pl

=head1 USAGE

 'species:s' => Species or common name
 'genus:s'  => Genus
 'ncbi:i'   => NCBI Taxonomy ID
 'all'	    => Return entire phylogeny
 'd|flatfile_dir:s'  => Bio::DB::Taxonomy flatfile (optional)
 GI	NCBI GI

=head1 DESCRIPTION

 An alternative is to use a flatfile acquired from ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz

=cut

use strict;
use Data::Dumper;
use Getopt::Long;
use Pod::Usage;
use Bio::Taxon;
use Bio::DB::Taxonomy;

my ($genus,$species,$ncbi_taxid,$all,$verbose,$gi,$entrez);
my $flatfile_dir='/databases/ncbi_taxonomy/';
GetOptions(
 'species:s' => \$species,
 'genus:s'  => \$genus,
 'ncbi:i'   => \$ncbi_taxid,
 'all'	   => \$all,
 'd|flatfile_dir:s' => \$flatfile_dir,
 'verbose' => \$verbose,
 'gi:i' => \$gi,
 'entrez:s' =>\$entrez
);
if ($ncbi_taxid && $ncbi_taxid!~/^\d+$/){$ncbi_taxid='';}

unless ($entrez || $ncbi_taxid || $genus || $species || $gi){
	$ncbi_taxid=shift;
	if ($ncbi_taxid!~/^\d+$/){$genus=$ncbi_taxid;$species=shift;$ncbi_taxid='';}
	warn ("Need an NCBI id or both a genus and species\n") unless ($ncbi_taxid || ($genus && $species));
	pod2usage unless ($ncbi_taxid || ($genus && $species));
}
my $taxondb = Bio::DB::Taxonomy->new(-source => 'entrez');
my $ncbi_dictionary;
if ($flatfile_dir && (-d $flatfile_dir) && (-s $flatfile_dir.'/nodes.dmp') && (-s $flatfile_dir.'/names.dmp')){
	warn "Using $flatfile_dir as NCBI TaxDB directory\n" if $verbose;
	if (!-s "$flatfile_dir/gi_taxid_nucl.bin" && -s "$flatfile_dir/gi_taxid_nucl.dmp"){
		use Bio::LITE::Taxonomy::NCBI::Gi2taxid qw/new_dict/;
		new_dict (in => "$flatfile_dir/gi_taxid_nucl.dmp", out => "$flatfile_dir/gi_taxid_nucl.bin");
	}
	$taxondb = Bio::DB::Taxonomy->new(-source => 'flatfile',-nodesfile => $flatfile_dir.'/nodes.dmp',-namesfile => $flatfile_dir.'/names.dmp',-directory => $flatfile_dir);
	$ncbi_dictionary = Bio::LITE::Taxonomy::NCBI::Gi2taxid->new(dict=>"$flatfile_dir/gi_taxid_nucl.bin") if (-s "$flatfile_dir/gi_taxid_nucl.bin");
}
else{
	warn "Using remote Entrez query\n" if $verbose;
}
if ($entrez){
	use Bio::DB::Query::GenBank;
	use Bio::DB::GenBank;
	my $query = Bio::DB::Query::GenBank->new(-db => 'nucleotide',-query => $entrez,-verbose=>-1);
	my @ids= $query->ids if $query->count && $query->count==1;
	$gi = $ids[0];
}

if ($gi && !$ncbi_taxid){
	die "No NCBI dictionary\n" unless $ncbi_dictionary;
	$ncbi_taxid = $ncbi_dictionary->get_taxid($gi);
}
elsif ($genus && $species && !$ncbi_taxid){
	$ncbi_taxid = $taxondb->get_taxonid($genus.' '.$species);
	if (!$ncbi_taxid){
		die "Could not find NCBI taxid for $genus $species.\n";
}
}elsif($genus && !$species && !$ncbi_taxid){
	$ncbi_taxid = $taxondb->get_taxonid($genus);
	if (!$ncbi_taxid){
		my $str='';
		$str.=$genus if $genus;
		$str.=$species if $species;
		die "Could not find NCBI taxid for $str.\n";
	}
}
my (%all_taxa,%all_taxa_sorted);
my $sorted_ref=\%all_taxa_sorted;
my $taxon = $taxondb->get_taxon(-taxonid => $ncbi_taxid)  || die "Cannot find your query\n";;
# we want class, order, family, genus, species, common name
my ($class,$order,$family,$i);
my $species_rank=$taxon->rank();
my $species_name=$taxon->scientific_name();
if ($species_rank eq 'genus'){
	$species = $species_name.' sp.';
	$genus=$species_name;
}elsif($species_rank eq 'species'){
	$species = $species_name;
}if (!$species){
	my $t=$taxon;
	while ($t=$t->ancestor()){
		$species_rank=$t->rank();
		next if !$species_rank;
		$species = $taxon->scientific_name() if ($species_rank eq 'species'); 
	}
	if (!$species){
		my $t=$taxon;
		while ($t=$t->ancestor()){
			$species_rank=$t->rank();
			next if !$species_rank;
			$species = $taxon->scientific_name().' sp.' if ($species_rank eq 'genus');
			$genus=$taxon->scientific_name();
		}
		if (!$species){
			$species = 'Unknown';
			$genus = 'Unknown';
			$i++;$sorted_ref->{$i}={
				'rank' => $species_rank,
				'name' => $species_name,
			};
		}
	}
}
$species =~s/^(\w+)\s//;
$i++;$sorted_ref->{$i}={
	'rank' => 'species',
	'name' => $species,
} if $species ne 'Unknown';
$all_taxa{$species_rank}=$species;
my $t=$taxon;
unless ($genus && $genus eq 'Unknown'){
while ($t=$t->ancestor()){
	my $rank=$t->rank();
	next if !$rank;
	my $name = $t->scientific_name();
	next if !$name;
	$i++;$sorted_ref->{$i}={
		'rank' => $rank,
		'name' => $name,
	};
	$all_taxa{$rank}=$name;
	if ($rank eq 'genus'){
		$genus=$name;
		last;
	}
 }
}
if (!$genus || !$species){die "Cannot find genus ($genus) or species ($species)\n";}
my $common="No common name";
$common=$taxon->common_names() ? $taxon->common_names() : $genus.' '.$species if $genus ne 'Unknown';
$t=$taxon;
while ($t=$t->ancestor()){
	my $rank=$t->rank();
	next if !$rank;
	my $name = $t->scientific_name() ;
	next if !$name;
	$i++;
	$sorted_ref->{$i}={
		'rank' => $rank,
		'name' => $name,
	};
	$all_taxa{$rank}=$name;
	if ($rank eq 'class'){
		$class=$name;
		last if (!$all && $genus ne 'Unknown');
	}elsif($rank eq 'order'){
		$order=$name;
	}elsif($rank eq 'family'){
		$family=$name;
	}
}
if (!$all && $genus ne 'Unknown'){
	if (!$class){
		# go up
		if ($all_taxa{'phylum'}){
			$class=$all_taxa{'phylum'};
		}elsif ($all_taxa{'kingdom'}){
			$class=$all_taxa{'kingdom'};
		}else{
			$class='unknown';
		}
	}if (!$order){
		#go down one, then up one
		if ($all_taxa{'suborder'}){
			$order=$all_taxa{'suborder'};
		}
		elsif ($all_taxa{'subclass'}){
			$order=$all_taxa{'subclass'};
		}
		else{
			$order='unknown';
		}
	}if (!$family){
		#go down one, then up one
		if ($all_taxa{'subfamily'}){
	                $family=$all_taxa{'subfamily'};
        	}elsif ($all_taxa{'superfamily'}){
	                $family=$all_taxa{'superfamily'};
        	}else{
	                $family='unknown';
        	}
	}
		print "class:$class;order:$order;family:$family;genus:$genus;species:$species;common:$common;ncbi:$ncbi_taxid\n";
}elsif($genus eq 'Unknown' && !$all){
        my $print="";
        foreach my $sort (sort {$b<=>$a} keys %all_taxa_sorted){
                $print.=$all_taxa_sorted{$sort}{'rank'}.':'.$all_taxa_sorted{$sort}{'name'}.';' unless $all_taxa_sorted{$sort}{'name'} eq 'unknown';
        }
        print $print."common:$common;ncbi:$ncbi_taxid\n";


}else{
	my ($print);
	foreach my $sort (sort {$b<=>$a} keys %all_taxa_sorted){
		$print.=$all_taxa_sorted{$sort}{'rank'}.':'.$all_taxa_sorted{$sort}{'name'}.';';
	}
	print $print."common:$common;ncbi:$ncbi_taxid\n";
}
