#!/usr/bin/env perl
### https://www.jianshu.com/p/948b90815dca
use strict;
use warnings;
use File::Basename;
use Getopt::Long;

sub usage{
    print <<USAGE;
Name:
    $0
Description:
    fetch fasta sequence form gtf and output gene2transcript information
Update:
    2016-12-14  edit by Alipe
Usage:
    perl $0 <gtf> <fa> [options]
Options:
    -g, --gtf       <string>*   refernece gtf file
    -f, --fasta     <string>*   reference fasta file
    -p, --prefix    <string>    prefix of output file       [default: gtfname]
    -o, --outdir    <string>    output directory            [default: ./ ]
    help                        print this help information
e.g
    perl $0 -g unreference.gtf -f hg38.fa -p fetched -o ./
    perl $0 --gtf unreference.gtf --fasta hg38.fa --prefix fetched --outdir ./
USAGE
    exit 0;
}

my ($help,$outdir,$prefix,$fasta,$gtf);
GetOptions(
    "gtf|g=s"    => \$gtf,
    "fasta|f=s"  => \$fasta,
    "prefix|p:s" => \$prefix,
    "outdir|o:s" => \$outdir,
    "help|?"     => \$help
);
die &usage if (!defined $gtf || !defined $fasta || defined $help);

$prefix ||= basename($gtf);
$outdir ||= ".";
mkdir $outdir if (not -e $outdir);

## read gtf file
my %gtf;
my %exon;
my %transcript;
my %fasta;

open GTF,$gtf or die "Can't open the gtf file $!";
while (<GTF>){
    chomp;
    my @a = split /\t/,$_;
    if ($a[2] eq "exon"){
        (my $gen) = $_ =~ /\s+gene_id\s+"(\S+)";/;
        (my $tra) = $_ =~ /\s+transcript_id\s+"(\S+)";/;
        (my $num) = $_ =~ /\s+exon_number\s+"(\S+)";/;
#        push @{$exon{$id}},$_;
#        my $n = @{$exon{$id}};
        $transcript{$tra}{"strand"} = $a[6];
        $transcript{$tra}{"geneid"} = $gen ;
        $gtf{$a[0]}{$tra}{$num}{"start"}  = $a[3] - 1;
        $gtf{$a[0]}{$tra}{$num}{"length"} = $a[4] - $a[3] + 1;
    } else {
        next;
    }

}
close GTF;

open FAO,">$outdir/$prefix.fa"  || die "Can't write the file $outdir/$prefix.fa $! \n";
open G2T,">$outdir/gene2tr.txt" || die "Can't write the file $outdir/gene2tr.txt $!\n";
open FA,$fasta || die "Can't open the reference fasta file $!\n";
$/ = ">";
<FA>;
while (<FA>){
    my ($info,$seq) = split(/\n/,$_,2);
    $seq =~ s/\n+//g;
    my @b = split(/\s/,$info);
    my $z = $gtf{$b[0]};
    foreach my $k (sort keys %$z){

        my $strand = $transcript{$k}{"strand"};
        my $name   = $transcript{$k}{"geneid"};
        my $x = $z->{$k};

        foreach my $j (keys %$x ) {
            my $start  = $x->{$j}->{"start"};
            my $length = $x->{$j}->{"length"};
            $fasta{$k} .= substr($seq,$start,$length);
        }

        my $fa = $fasta{$k};

        if ($strand eq "-" ) {
            $fa = reverse $fa;
            $fa =~ tr/AaGgCcTt/TtCcGgAa/;
        }
        $fa = out_fasta($fa,60);

        print FAO ">$k\t$name\n$fa\n";
        print G2T "$name\t$k\n";
    }
}
close FA;
close FAO;
close G2T;

sub out_fasta{
    my ($seq,$num) = @_;
    my $len = length $seq;
    $seq =~ s/([A-Za-z]{$num})/$1\n/g;
    chop($seq) unless $len % $num;
    return $seq;
}

__END__
