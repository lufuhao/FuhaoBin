#!/usr/local/bin/perl
#http://ib.berkeley.edu/labs/slatkin/munch/scriptlist/filtergff.txt
# Time-stamp: <2004-02-23 11:49:08 kasper>

use warnings;
use strict;
use Getopt::Long;
use constant USAGE =><<END;

SYNOPSIS:

 filtergff.pl [OPTIONS] [inputfile] [outputfile]

DESCRIPTION:

Filters a GFF stream based on given options.

OPTIONS:

       --help
             Prints a help.
       --seqname <regex>
             Prints all lines with a seqname (first entry) that matches a perl regex.
       --source <regex>
             Prints all lines with a source field that matches a perl regex.
       --feature <regex>
             Prints all lines with a feature field that matches a perl regex.
       --maxstart <int>
             Prints all lines with a start at at most <int>.
       --minstart <int>
             Prints all lines with a start at at least <int>.
       --maxend <int>
             Prints all lines with a end at at most <int>.
       --minend <int>
             Prints all lines with a end at at least <int>.
       --maxlength <int>
             Prints all lines with at most length <int>.
       --minlength <int>
             Prints all lines with at least length <int>.
       --maxscore <float>
             Prints all lines with a score of at most <float>.
       --minscore <float>
             Prints all lines with a score of at least <float>.
       --strand <+|->
             Prints all lines with strand corresponding to <-1|0|1>.
       --frame <regex>
             Prints all lines with a frame field that matches a perl regex.
       --attribute <regex>
             Prints all lines with a attribute field that matches a perl regex.

       The following options are only valid if the GFF lines have a perent and a id tag
       in the attribute field (Eg. Sequence ds342_4 ; ID kj2344_d)

       --id <regex>
             Prints all lines where the ID matches a perl regex.
       --parent <regex>
             Prints all lines where the parent matches a perl regex.
       --parent_start <int>

       --parent_end <int>

       --noself
             Prints all lines that is not a match to itself.

       The following options are only valid if the GFF lines are blast_hsps
       with the propor additional info in the atribute field. This is for use
       with blast.pl which produces such GFF lines.

       --query_desc <regex>
             Prints all lines with a query description matching <regex>.
       --minhit_frac_ident <float>
             Prints all lines with min hit fractional identity of <float>.
       --minquery_frac_ident <float>
             Prints all lines with min query fractional identity of <float>.
       --minhit_frac_cons <float>
             Prints all lines with min hit fractional conservation of <float>.
       --minquery_frac_cons <float>
             Prints all lines with min query fractional conservation of <float>.
       --mintotal_frac_cons <float>
             Prints all lines with min total fractional conservation of <float>.
       --minpersent_ident <float>
             Prints all lines with min percent indentity of <int>.
       --maxhit_gaps <int>
             Prints all lines with at most <int> gaps on the hit.
       --maxquery_gaps <int>
             Prints all lines with at most <int> gaps on the query.
       --maxtotal_gaps <int>
             Prints all lines with at most <int> gaps in the HSP.
       --nogaps
             Prints all lines with no gaps at all .

       --not
             Negates the entire query. Like the grep option -v
       --strict
             Returns with an error unless ALL GFF lines have the information used by the query.


EXAMPLES:

 cat file.gff | filtergff.pl > filtered.gff

 filtergff.pl file.gff filtered.gff

AUTHOR:

Kasper Munch

COPYRIGHT:

This program is free software. You may copy and redistribute it under
the same terms as Perl itself.

END

my $help;
my $not = 0;
my $strict = 0;
my $maxlength;
my $minlength;
my $minstart;
my $maxstart;
my $minend;
my $maxend;
my $strand;
my $seqname;
my $attribute;
my $frame;
my $maxscore;
my $minscore;
my $source;
my $feature;
my $regex;
my $id;
my $parent;
my $parent_start;
my $parent_end;
my $noself;
my $query_desc;
my $minhit_frac_ident;
my $minquery_frac_ident;
my $minhit_frac_cons;
my $minquery_frac_cons;
my $mintotal_frac_cons;
my $minpercent_ident;
my $maxhit_gaps;
my $maxquery_gaps;
my $maxtotal_gaps;
my $nogaps;

GetOptions("help" => \$help, 
           "not" => \$not, 
           "strict" => \$strict, 
           "seqname=s" => \$seqname, 
           "source=s" => \$source, 
           "feature=s" => \$feature, 
           "maxstart=s" => \$maxstart, 
           "minstart=s" => \$minstart, 
           "maxend=s" => \$maxend, 
           "minend=s" => \$minend, 
           "maxlength=s" => \$maxlength, 
           "minlength=s" => \$minlength, 
           "maxscore=s" => \$maxscore, 
           "minscore=s" => \$minscore, 
           "strand=s" => \$strand, 
           "frame=s" => \$frame, 
           "attribute=s" => \$attribute, 
           "regex=s" => \$regex, 
           "id=s" => \$id,
           "parent=s" => \$parent,
           "parent_start=s" => \$parent_start,
           "parent_end=s" => \$parent_end,
           "noself" => \$noself,
           "query_desc=s" => \$query_desc,
           "minhit_frac_ident=f" => \$minhit_frac_ident, 
           "minquery_frac_ident=f" => \$minquery_frac_ident, 
           "minhit_frac_cons=f" => \$minhit_frac_cons, 
           "minquery_frac_cons=f" => \$minquery_frac_cons, 
           "mintotal_frac_cons=f" => \$mintotal_frac_cons,
           "minpercent_ident=f" => \$minpercent_ident, 
           "maxhit_gaps=s" => \$maxhit_gaps, 
           "maxquery_gaps=s" => \$maxquery_gaps, 
           "maxtotal_gaps=s" => \$maxtotal_gaps, 
           "nogaps" => \$nogaps) or die USAGE;

$help and die USAGE;


@ARGV = ('-') unless @ARGV;
my $input = shift @ARGV;
open my $in, "$input" or die "$input: $!\n";
@ARGV = ('>&STDOUT') unless @ARGV;
my $output = shift @ARGV;
open my $out, ">$output" or die "$output: $!\n";

# 0 seqname
# 1 source
# 2 feature
# 3 start
# 4 end
# 5 score 
# 6 strand
# 7 frame
# 8 [attribute]

while (my $line = <$in>) {

    next if $line =~ /^(\#)/;
    my @f = split(/\t/, $line);
    my @a =  split(/ ; /, $f[8]);
    my %extrainfo;
    my $parentseq = shift @a;
    $parentseq =~ /^\S+\s+(\S+)/;
    $extrainfo{parent} = $1;
    if ($parentseq =~ /^\S+\s+\S+\s+(\d+)\s+(\d+)/) {
        $extrainfo{parent_start} = $1;
        $extrainfo{parent_end} = $2;    
    }
    for my $a (@a) {
        $a =~ /^(\S+)\s+(\S+.*)/;
        (my $key, my $value) = ($1, $2);
        if ($key eq 'ID') {
            $key = 'id';
        } 
        $extrainfo{$key} = $value;
    }
    if ($seqname) { 
        if ($f[0] =~ /$seqname/) { next if $not } else { next unless $not }
    }
    if ($source) {
        if ($f[1] =~ /$source/) { next if $not } else { next unless $not }
    }
    if ($feature) {
        if ($f[2] =~ /$feature/) { next if $not } else { next unless $not }
    }
    if ($maxstart) {
        if ($f[3] <= $maxstart) { next if $not } else { next unless $not }
    }
    if ($minstart) {
        if ($f[3] >= $minstart ) { next if $not } else {next unless $not }
    }
    if ($maxend) {
        if ($f[4] <= $maxend ) { next if $not } else { next unless $not }
    }
    if ($minend) {
        if ($f[4] >= $minend) { next if $not } else { next unless $not }
    }
    if ($maxlength) {
        if ($f[4] - $f[3] + 1 <= $maxlength ) { next if $not } else { next unless $not }
    }
    if ($minlength) {
        if ($f[4] - $f[3] + 1 >= $minlength) { next if $not } else { next unless $not }
    }
    if ($maxscore) {
        if ($f[5] <= $maxscore) { next if $not } else { next unless $not }
    }
    if ($minscore) {
        if ($f[5] >= $minscore) { next if $not } else { next unless $not }
    }
    if ($strand) {
        if ($f[6] eq $strand) { next if $not } else { next unless $not }
    }
    if ($frame) {
        if ($f[7] =~ /$frame/) { next if $not } else { next unless $not }
    }
    if ($attribute) {
        if ($f[8] =~ /$attribute/) { next if $not } else { next unless $not }
    }
    if ($regex) {
        if ($line =~ /$regex/) { next if $not } else { next unless $not }
    }
    if ($id) {
        !$strict or defined($extrainfo{id}) or die "Option --id not surported by a GFF string.\n";
        if ($extrainfo{id} =~ /$id/) { next if $not } else { next unless $not }
    }
    if ($parent) {
        !$strict or defined($extrainfo{parent}) or die "Option --parent not surported by a GFF string.\n";
        if ($extrainfo{parent} =~ /$parent/) { next if $not } else { next unless $not }
    }
    if ($parent_start) {
        !$strict or defined($extrainfo{parent_start}) or die "Option --parent_start not surported by a GFF string.\n";
        if ($extrainfo{parent_start} =~ /$parent_start/) { next if $not } else { next unless $not }
    }
    if ($parent_end) {
        !$strict or defined($extrainfo{parent_end}) or die "Option --parent_end not surported by a GFF string.\n";
        if (defined($extrainfo{parent_end}) && $extrainfo{parent_end} =~ /$parent_end/) { next if $not } else { next unless $not }
    }
    if ($noself) {
        !$strict or defined($extrainfo{parent}) or die "Option --noself not surported by a GFF string.\n";
        if (defined($extrainfo{parent}) && $extrainfo{parent} ne $f[0]) { next if $not } else { next unless $not }
    }
    # Additional attribute information in blast GFFs:
    if ($query_desc) {
        !$strict or defined($extrainfo{query_desc}) or die "Option --query_desc not surported by a GFF string.\n";
        if (defined($extrainfo{query_desc}) && $extrainfo{query_desc} =~ /$query_desc/) { next if $not } else { next unless $not }
    }
    if ($minhit_frac_ident) {
        !$strict or defined($extrainfo{hit_frac_ident}) or die "Option --minhit_frac_ident not surported by a GFF string.\n";
        if (defined($extrainfo{hit_frac_ident}) && $extrainfo{hit_frac_ident} >= $minhit_frac_ident) { next if $not } else { next unless $not }
    }
    if ($minquery_frac_ident) {
        !$strict or defined($extrainfo{query_frac_ident}) or die "Option --minquery_frac_ident not surported by a GFF string.\n";
        if (defined($extrainfo{query_frac_ident}) && $extrainfo{query_frac_ident} >= $minquery_frac_ident) { next if $not } else { next unless $not }
    }
    if ($minhit_frac_cons) {
        !$strict or defined($extrainfo{hit_frac_cons}) or die "Option --minhit_frac_cons not surported by a GFF string.\n";
        if (defined($extrainfo{hit_frac_cons}) && $extrainfo{hit_frac_cons} >= $minhit_frac_cons) { next if $not } else { next unless $not }
    }
    if ($minquery_frac_cons) {
        !$strict or defined($extrainfo{query_frac_cons}) or die "Option --minquery_frac_cons not surported by a GFF string.\n";
        if (defined($extrainfo{query_frac_cons}) && $extrainfo{query_frac_cons} >= $minquery_frac_cons) { next if $not } else { next unless $not }
    }
    if ($mintotal_frac_cons) {
        !$strict or defined($extrainfo{total_frac_cons}) or die "Option --mintotal_frac_cons not surported by a GFF string.\n";
        if (defined($extrainfo{total_frac_cons}) && $extrainfo{total_frac_cons} >= $mintotal_frac_cons) { next if $not } else { next unless $not }
    }
    if ($minpercent_ident) {
        !$strict or defined($extrainfo{percent_ident}) or die "Option --minpercent_ident not surported by a GFF string.\n";
        if (defined($extrainfo{percent_ident}) && $extrainfo{percent_ident} >= $minpercent_ident) { next if $not } else { next unless $not }
    }
    if ($maxhit_gaps) {
        !$strict or defined($extrainfo{hit_gaps}) or die "Option --maxhit_gaps not surported by a GFF string.\n";
        if (defined($extrainfo{hit_gaps}) && $extrainfo{hit_gaps} <= $maxhit_gaps) { next if $not } else { next unless $not }
    }
    if ($maxquery_gaps) {
        !$strict or defined($extrainfo{query_gaps}) or die "Option --maxquery_gaps not surported by a GFF string.\n";
        if (defined($extrainfo{query_gaps}) && $extrainfo{query_gaps} <= $maxquery_gaps) { next if $not } else { next unless $not }
    }
    if ($maxtotal_gaps) {
        !$strict or defined($extrainfo{total_gaps}) or die "Option --maxtotal_gaps not surported by a GFF string.\n";
        if (defined($extrainfo{total_gaps}) && $extrainfo{total_gaps} <= $maxtotal_gaps) { next if $not } else { next unless $not }
    }
    if ($nogaps) {
        !$strict or defined($extrainfo{total_gaps}) or die "Option --nogaps not surported by a GFF string.\n";
        if (defined($extrainfo{total_gaps}) && $extrainfo{total_gaps} == 0) { next if $not } else { next unless $not }
    }
    print $out $line;
}

=head1 SYNOPSIS:

 filtergff.pl [OPTIONS] [inputfile] [outputfile]

=head1 DESCRIPTION:

Filters a GFF stream based on given options.

=head1 OPTIONS:

=over 4

=item --help

Prints a help.

=item --seqname <regex>

Prints all lines with a seqname (first entry) that matches a perl regex.

=item --source <regex>

Prints all lines with a source field that matches a perl regex.

=item --feature <regex>

Prints all lines with a feature field that matches a perl regex.

=item --maxstart <int>

Prints all lines with a start at at most <int>.

=item --minstart <int>

Prints all lines with a start at at least <int>.

=item --maxend <int>

Prints all lines with a end at at most <int>.

=item --minend <int>

Prints all lines with a end at at least <int>.

=item --maxlength <int>

Prints all lines with at most length <int>.

=item --minlength <int>

Prints all lines with at least length <int>.

=item --maxscore <float>

Prints all lines with a score of at most <float>.

=item --minscore <float>

Prints all lines with a score of at least <float>.

=item --strand <+|->

Prints all lines with strand corresponding to <-1|0|1>.

=item --frame <regex>

Prints all lines with a frame field that matches a perl regex.

=item --attribute <regex>

Prints all lines with a attribute field that matches a perl regex.

The following options are only valid if the GFF lines have a perent and a id tag
in the attribute field (Eg. Sequence ds342_4 ; ID kj2344_d)

=item --id <regex>

Prints all lines where the ID matches a perl regex.

=item --parent <regex>

Prints all lines where the parent matches a perl regex.

=item --parent_start <int>


=item --parent_end <int>


=item --noself

Prints all lines that is not a match to itself.

The following options are only valid if the GFF lines are blast_hsps
with the propor additional info in the atribute field. This is for use
with blast.pl which produces such GFF lines.

=item --query_desc <regex>

Prints all lines with a query description matching <regex>.

=item --minhit_frac_ident <float>

Prints all lines with min hit fractional identity of <float>.

=item --minquery_frac_ident <float>

Prints all lines with min query fractional identity of <float>.

=item --minhit_frac_cons <float>

Prints all lines with min hit fractional conservation of <float>.

=item --minquery_frac_cons <float>

Prints all lines with min query fractional conservation of <float>.

=item --mintotal_frac_cons <float>

Prints all lines with min total fractional conservation of <float>.

=item --minpersent_ident <float>

Prints all lines with min percent indentity of <int>.

=item --maxhit_gaps <int>

Prints all lines with at most <int> gaps on the hit.

=item --maxquery_gaps <int>

Prints all lines with at most <int> gaps on the query.

=item --maxtotal_gaps <int>

Prints all lines with at most <int> gaps in the HSP.

=item --nogaps

Prints all lines with no gaps at all .

=item --not

Negates the entire query. Like the grep option -v

=item --strict

Returns with an error unless ALL GFF lines have the information used by the query.


=back

=head1 EXAMPLES:

 cat file.gff | filtergff.pl > filtered.gff

 filtergff.pl file.gff filtered.gff

=head1 AUTHOR:

Kasper Munch

=head1 COPYRIGHT:

This program is free software. You may copy and redistribute it under
the same terms as Perl itself.


=cut
