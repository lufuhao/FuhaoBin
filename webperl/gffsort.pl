#!/usr/bin/env perl
use strict;
use warnings;
# Time-stamp: <2004-02-29 23:07:28 kasper>
#http://ib.berkeley.edu/labs/slatkin/munch/scriptlist/gffsort.txt
use warnings;
use strict;
use Getopt::Long;
use constant USAGE =><<END;

SYNOPSIS:

 gffsort.pl [OPTIONS]

DESCRIPTION:

This script sorts gff lines from standard input by first sequence
name, then start and then end.

OPTIONS:

       --help
             Prints this help.
       --decending
             Sorts decending and not ascending.

EXAMPLES:

 cat file.gff | gffsort.pl > sorted.gff

AUTHOR:

Kasper Munch

COPYRIGHT:

This program is free software. You may copy and redistribute it under
the same terms as Perl itself.

END

my $help = 0;
my $decending = 0;
GetOptions(
           "help|h" => \$help,
          ) or die USAGE;

$help and die USAGE;

my @f=();
my @a=();
my @s=();
my @n=();
my @e=();

while (<>) {
    s/#.*//;
    next unless /\S/;
    @f = split /\t/;
    push @a, $_;
    push @n, $f[0];
    push @s, $f[3];
	push @e, $f[4];
}

if ($decending) {
	foreach my $i (sort { $n[$a] cmp $n[$b] or $s[$a] <=> $s[$b] or $e[$a] <=> $e[$b] } 0..$#a) {
		print $a[$i];
	}
}
else {
	foreach my $i (sort { $n[$b] cmp $n[$a] or $s[$b] <=> $s[$a] or $e[$b] <=> $e[$a] } 0..$#a) {
		print $a[$i];
	}
}

=head1 SYNOPSIS:

 gffsort.pl [OPTIONS]

=head1 DESCRIPTION:

This script sorts gff lines from standard input by first sequence
name, then start and then end.

=head1 OPTIONS:

=over 4

=item --help

Prints this help.

=item --decending

Sorts decending and not ascending.

=back

=head1 EXAMPLES:

 cat file.gff | gffsort.pl > sorted.gff

=head1 AUTHOR:

Kasper Munch

=head1 COPYRIGHT:

This program is free software. You may copy and redistribute it under
the same terms as Perl itself.


=cut
