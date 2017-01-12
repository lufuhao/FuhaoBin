#!/usr/bin/env perl
use strict;
use warnings;
use Bio::Graphics;
use Bio::SeqIO;
use Bio::SeqFeature::Generic;
use constant USAGE =><<EOH;

usage:  perl graphics.pl NM_172587.gb >1.png

Dexcription:
	把Genbank格式的序列转换为png图片-基因结构图，显示序列的长度，CDS区，exon区，STS区等。简单地讲就是把该序列的Genbank格式里的信息用图片表示。
Require: 
	BioPerl (	Bio::Graphics
				Bio::SeqIO
				Bio::SeqFeature::Generic
			)
Source: http://liucheng.name/531/

v20151102

EOH
die USAGE if (scalar(@ARGV) !=1 or [0] eq '-h' or [0] eq '--help');



my $file = shift or die "provide a sequence file as the argument";
my $io = Bio::SeqIO->new(-file=>$file) or die "couldn't create Bio::SeqIO";
my $seq = $io->next_seq or die "couldn't find a sequence in the file";
my @features = $seq->all_SeqFeatures;
# sort features by their primary tags
my %sorted_features;
for my $f (@features) {
    my $tag = $f->primary_tag;
    push @{$sorted_features{$tag}},$f;
}
my $panel = Bio::Graphics::Panel->new(
-length => $seq->length,
-key_style => 'between',
-width => 800,
-pad_left => 10,
-pad_right => 10);
$panel->add_track(arrow =>
	Bio::SeqFeature::Generic->new(-start => 1,
	-end => $seq->length),
	-bump => 0,
	-double=>1,
	-tick => 2
);
$panel->add_track(generic =>
	Bio::SeqFeature::Generic->new(-start => 1,
	-end => $seq->length,
	-bgcolor => 'blue',
	-label => 1)
);
# general case
my @colors = qw(cyan orange blue purple green chartreuse magenta yellow aqua);
    my $idx = 0;
    for my $tag (sort keys %sorted_features) {
        my $features = $sorted_features{$tag};
		$panel->add_track($features,
			-glyph => 'generic',
			-bgcolor => $colors[$idx++ % @colors],
			-fgcolor => 'black',
			-font2color => 'red',
			-key => "${tag}s",
			-bump => +1,
			-height => 8,
			-label => 1,
			-description => 1);
	}
print $panel->png;
