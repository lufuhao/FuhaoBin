#!/usr/bin/env perl
use strict;
use warnings;
use Bio::Seq;
use Bio::SeqIO;
use Bio::DB::Fasta;
use Cwd;
use Data::Dumper qw /Dumper/;
use FuhaoPerl5Lib::FastaKit qw /Codon2AA/;
use FuhaoPerl5Lib::GffKit qw /ReadGff3 WriteGff3/;
use Getopt::Long;
use constant USAGE =><<EOH;

usage: $0 -f genome.fa -g infile.gff3 -p output.prefix -l length_bp_updownstream -s Outoption

Desccriptions:
    Extract CDS/cDNA/genomic sequence sequences from genome by GFF3
    Option 10-11 are designed for RATT
        10 output EMBL for each old seq
        11 output EMBL for each gene; -l option available

Options:
    --help | -h
        Print this help/usage;
    --fasta | -f    [STR]
        Fasta file
    --gff | -g      [STR]
        GFF3 file
    --prefix | -p   [STR]
        output file Prefix for OUTOPT 1-9; [MyOutput]
    --length | -l   [INT]
        Length up/down - stream; [2000]
    --options | -s [INT,INT...]
             1=genomic seq
             2=CDS
             3=cDNA
             4=protein
             5=upstream.only
             6=genome.withupstream
             7=downstream.only
             8=genome.with.downstream
             9=genome.with.up/downstream
             10=embl per seq
             11=embl per gene
             12=genome.with.up/downstream + corresponding GFF
                Will output 9 and GFF3
             13=embl in one file
          default:1,2,3,4
    --outdir | -d
        Output dir for OUTOPT 10-11; [./MyEMBL]
    --verbose
        Detailed output for trouble-shooting;
    --version|v!
        Print current SCRIPT version;


v20180213

EOH
die USAGE unless @ARGV;



###Receving parameter################################################
my ($help, $verbose, $debug, $version);
my ($input, $output);
my ($file_fasta, $gffin, $outprefix, $myoptions, $outdir);
my $stream=2000;

GetOptions(
	"help|h!" => \$help,
	"fasta|f:s" => \$file_fasta,
	"gff|g:s" => \$gffin,
	"prefix|p:s" => \$outprefix,
	"length|l:i" => \$stream,
	"options|s:s" => \$myoptions,
	"outdir|d:s" => \$outdir,
#	"output|o:s" => \$output,
#	!:s:i
	"debug!" => \$debug,
	"verbose!" => \$verbose,
	"version|v!" => \$version) or die USAGE;
($help or $version) and die USAGE;




### Defaults ########################################################
$debug=0 unless (defined $debug);
$verbose=0 unless (defined $verbose);
$outprefix='MyOutput' unless (defined $outprefix);
unless (defined $outdir) {
	$outdir=getcwd;
	$outdir.='/MyEMBL';
}

my %outcode=();
my @outarr=(1,2,3,4);
my $cds={};
my $referenceids={};
my $mrnas={};
my $gene={};
my $gene2mrna=();
my %mergedregion=();
my $exons={};
my $printID=0;


### input and output ################################################
print "\n\n\n";
die "Error: invalid fasta input\n" unless (defined $file_fasta and -s $file_fasta);
die "Error: invalid gff input\n" unless (defined $gffin and -s $gffin);
unless (defined $stream and $stream=~/^\d+$/ and $stream >= 0) {
	die "Error: invalid Up/down - stream length\n"
}
print "INFO: UpDoanStream options: $stream bp\n";
### Decide output option
if (defined $myoptions) {
	@outarr=split(/,/, $myoptions);
}
foreach my $codeind (@outarr) {
	$outcode{$codeind}++ if (defined $codeind and $codeind=~/^\d+$/);
}
print "INFO: OUTPUT options: ", join(",", sort {$a<=>$b} keys %outcode), "\n";

###1. Genome
unlink "$outprefix.genome.fasta" if (exists $outcode{'1'} and -e "$outprefix.genome.fasta");
my $outfile_gene = Bio::SeqIO->new( -format => 'fasta', -file => ">$outprefix.genome.fasta" ) if (exists $outcode{'1'});

###2. CDS
unlink "$outprefix.cds.fasta" if (exists $outcode{'2'} and -e "$outprefix.cds.fasta");
my $outfile_cds = Bio::SeqIO->new( -format => 'fasta', -file => ">$outprefix.cds.fasta" ) if (exists $outcode{'2'});

###3. cDNA
unlink "$outprefix.cdna.fasta" if (exists $outcode{'3'} and -e "$outprefix.cdna.fasta");
my $outfile_cdna = Bio::SeqIO->new( -format => 'fasta', -file => ">$outprefix.cdna.fasta" ) if (exists $outcode{'3'});

### 4. Protein
unlink "$outprefix.pep.fasta" if (exists $outcode{'4'} and -e "$outprefix.pep.fasta");
my $outfile_pep = Bio::SeqIO->new( -format => 'fasta', -file => ">$outprefix.pep.fasta" ) if (exists $outcode{'4'});

###5. upstream.only
unlink "$outprefix.upstream$stream.only.fasta" if (exists $outcode{'5'} and -e "$outprefix.upstream$stream.only.fasta");
my $outfile_upstream3000only = Bio::SeqIO->new( -format => 'fasta', -file => ">$outprefix.upstream$stream.only.fasta" ) if (exists $outcode{'5'});

###6. genome.withupstream
unlink "$outprefix.genome.with.upstream$stream.fasta" if (exists $outcode{'6'} and -e "$outprefix.genome.with.upstream$stream.fasta");
my $outfile_upstream3000with = Bio::SeqIO->new( -format => 'fasta', -file => ">$outprefix.genome.with.upstream$stream.fasta" ) if (exists $outcode{'6'});

###7. downstream.only
unlink "$outprefix.downstream$stream.only.fasta" if (exists $outcode{'7'} and -e "$outprefix.downstream$stream.only.fasta");
my $outfile_downstream3000only = Bio::SeqIO->new( -format => 'fasta', -file => ">$outprefix.downstream$stream.only.fasta" ) if (exists $outcode{'7'});

###8. genome.with.downstream
unlink "$outprefix.genome.with.downstream$stream.fasta" if (exists $outcode{'8'} and -e "$outprefix.genome.with.downstream$stream.fasta");
my $outfile_downstream3000with = Bio::SeqIO->new( -format => 'fasta', -file => ">$outprefix.genome.with.downstream$stream.fasta" ) if (exists $outcode{'8'});

###9. genome.with.up/downstream
unlink "$outprefix.genome.with.updownstream$stream.fasta" if ((exists $outcode{'9'} or exists $outcode{'12'}) and -e "$outprefix.genome.with.updownstream$stream.fasta");
my $outfile_updownstream3000with = Bio::SeqIO->new( -format => 'fasta', -file => ">$outprefix.genome.with.updownstream$stream.fasta") if (exists $outcode{'9'} or exists $outcode{'12'});

### 10. 11.
if (exists $outcode{'10'} or exists $outcode{'11'} or exists $outcode{'13'}) {
	my $numembl=0;
	$numembl++ if (exists $outcode{'10'});
	$numembl++ if (exists $outcode{'11'});
	$numembl++ if (exists $outcode{'13'});
	unless ($numembl==1) {
		die "Error: can not specify options 10/11/13 at the same time\nSeparated running is recommended\n";
	}
	if (exists $outcode{'10'} or exists $outcode{'11'}) {
		if ( -d $outdir) {
			my @existingfiles=glob "$outdir/*";
			if (scalar(@existingfiles)>0) {
				die "Error: Outdir not empty\n";
			}
		}
		else {
			mkdir $outdir, 0766 || die "Error: can not make dir $outdir\n";
		}
	}
	if (exists $outcode{'13'}) {
		unlink "$outprefix.embl" if (-e "$outprefix.embl");
	}
}



### 12 . OUTGFF
my $outputgff3="$outprefix.genome.with.updownstream$stream.gff3" if (exists $outcode{'12'});
unlink $outputgff3 if (exists $outcode{'12'} and -e $outputgff3);
### 13 EMBL in one file






### Main ############################################################
(my $test_success, $referenceids, $gene, $gene2mrna, $mrnas, $exons, $cds)=ReadGff3($gffin, $file_fasta);
unless ($test_success) {
	die "Error: failed to read GFF3\n";
}


### Reading GFF3 ####################################################
#print "\n\n\n### Reading GFF3\n";
#open (GFFIN, "< $gffin") || die "Error: can not open GFF input\n";
#while (my $line=<GFFIN>) {
#	chomp $line;
#	next if ($line=~/^#/);
#	my @arr=split(/\t/, $line);
#	unless (scalar(@arr)==9 and $arr[3]=~/^\d+$/ and $arr[4]=~/^\d+$/ and $arr[3]<=$arr[4]) {
#		print STDERR "Error: border error ($.): $line\n";
#		exit 1;
#	}
#	unless (defined $arr[6] and $arr[6]=~/(^\+$)|(^-$)/) {
#		unless ($arr[2]=~/chrom/i) {
#			print STDERR "Error: unknown strand line($.): $line\n" ;
#			exit 1;
#		}
#	}
#	if ($arr[2] =~ /^gene$/i) {
#		next unless (exists $outcode{'1'} or exists $outcode{'5'} or exists $outcode{'6'} or exists $outcode{'7'} or exists $outcode{'8'} or exists $outcode{'9'} or exists $outcode{'10'} or exists $outcode{'11'});
#		my $thisgeneid=$arr[8];
#		$thisgeneid=~s/^.*ID=//;$thisgeneid=~s/;.*$//;
#		if (exists ${$gene}{$thisgeneid}) {
#			print STDERR "Error: duplicated gene ID: $thisgeneid\n";
#			exit 1;
#		}
#		${$referenceids}{$arr[0]}{$arr[3]}{$thisgeneid}++;
#		${$gene}{$thisgeneid}{'reference'}=$arr[0];
#		${$gene}{$thisgeneid}{'start'}=$arr[3];
#		${$gene}{$thisgeneid}{'end'}=$arr[4];
#		${$gene}{$thisgeneid}{'strand'}=$arr[6];
#		${$gene}{$thisgeneid}{'score'}{$arr[5]}++;
#		if ($arr[8]=~/Note=(\S+)/) {
#			my $notes=$1;
#			$notes=~s/;.*$//;
#			$notes=~s/\s+$//;
#			$notes=~s/^\s+//;
#			${$gene}{$thisgeneid}{'note'}=$notes;
#		}
#		print "Test: geneid: $thisgeneid\tgenenote: ", ${$gene}{$thisgeneid}{'note'}, "\n"; ### For test ###
#	}
#	elsif ($arr[2] =~ /^mRNA$/i) {
#		next unless (exists $outcode{'10'} or exists $outcode{'11'});
#		my $thismrnaid=$arr[8];
#		$thismrnaid=~s/^.*ID=//; $thismrnaid=~s/;.*$//;
#		my $parent=$arr[8];
#		$parent=~s/^.*Parent=//; $parent=~s/;.*$//;
#		my @temparr=split(/,/, $parent);
#		${$mrnas}{$thismrnaid}{'start'}=$arr[3];
#		${$mrnas}{$thismrnaid}{'end'}=$arr[4];
#		${$mrnas}{$thismrnaid}{'strand'}{$arr[6]}++;
#		${$mrnas}{$thismrnaid}{'score'}{$arr[5]}++;
#		${$mrnas}{$thismrnaid}{'reference'}{$arr[0]}++;
#		foreach my $indpar (@temparr) {
#			${$gene2mrna}{$indpar}{$thismrnaid}++;
#		}
#	}
#	elsif ($arr[2] =~ /^exon$/i) {
#		next unless (exists $outcode{'3'} or exists $outcode{'10'} or exists $outcode{'11'});
#		my $parent=$arr[8];
#		$parent=~s/^.*Parent=//; $parent=~s/;.*$//;
###		print "Test: exon $arr[3]-$arr[4] \tParent: $parent\n"; ### For test ###
#		my @temparr=split(/,/, $parent);
#		my $exonid=$arr[8];
#		$exonid=~s/^.*ID=//; $exonid=~s/;.*$//;
#		foreach my $indpar (@temparr) {
#			if (exists $exons{$indpar} and exists $exons{$indpar}{'exon'} and exists $exons{$indpar}{'exon'}{$arr[3]} and exists $exons{$indpar}{'exon'}{$arr[3]}{$arr[4]}) {
#				die "Error: repeated exons ".$arr[3].'-'.$arr[4]. "for mRNA $indpar\n"; 
#			}
#			$exons{$indpar}{'exon'}{$arr[3]}{$arr[4]}=$exonid;
#			$exons{$indpar}{'strand'}{$arr[6]}++;
#			$exons{$indpar}{'score'}{$arr[5]}++;
#			$exons{$indpar}{'reference'}{$arr[0]}++;
#		}
#	}
#	elsif ($arr[2] =~ /^CDS$/i) {
#		next unless (exists $outcode{'2'} or exists $outcode{'4'} or exists $outcode{'10'} or exists $outcode{'11'});
#		my $parent=$arr[8];
#		$parent=~s/^.*Parent=//; $parent=~s/;.*$//;
###		print "Test: CDS $arr[3]-$arr[4] \tParent: $parent\n";### For test ###
#		my @temparr=split(/,/, $parent);
#		foreach my $indpar (@temparr) {
#			${$cds}{$indpar}{'exon'}{$arr[3]}{$arr[4]}++;
#			${$cds}{$indpar}{'strand'}{$arr[6]}++;
#			${$cds}{$indpar}{'score'}{$arr[5]}++;
#			${$cds}{$indpar}{'reference'}{$arr[0]}++;
#			${$cds}{$indpar}{'phase'}{$arr[3]}{$arr[4]}=$arr[7];
#		}
#	}
#}
#close GFFIN;



### Test
#print "Test: \$gene\n"; print Dumper $gene; print "\n"; ### For test ###
#print "Test: \$cds\n"; print Dumper $cds; print "\n"; ### For test ###
#print "Test: \%exons\n"; print Dumper \%exons; print "\n"; ### For test ###
#print "Test: \$referenceids\n"; print Dumper $referenceids; print "\n"; ### For test ###
#print "Test: \$mrnas\n"; print Dumper $mrnas; print "\n"; ### For test ###
#print "Test: \$gene2mrna\n"; print Dumper $gene2mrna; print "\n"; ### For test ###
#exit 0;


print "\n\n### Reading DB: $file_fasta\n";
my $db = Bio::DB::Fasta->new($file_fasta);



#Write gene sequences
if (exists $outcode{'1'} or exists $outcode{'5'} or exists $outcode{'6'} or exists $outcode{'7'} or exists $outcode{'8'} or exists $outcode{'9'} or exists $outcode{'12'}) {
	print "### 1. Info: output genome/UpDownstream sequences\n";
	foreach my $idvgene (sort keys %{$gene}) {
	#	print "Info: gene: $idvgene\n"; ### For test ###
		my $geneseq='';
		my $upstream3000seq='';
		my ($upstart, $upend);
		my $downstream3000seq='';
		my ($downstart,$downend);
		my $seqlength=$db->length(${$gene}{$idvgene}{'reference'});
		unless (defined $seqlength and $seqlength=~/^\d+$/ and $seqlength>0) {
			die "Error: invalid seq length: SeqID ", ${$gene}{$idvgene}{'reference'}, "\n";
		}
		my $test_upstream=0;
		my $test_downstream=0;
		my $rangestart=${$gene}{$idvgene}{'start'};
		my $rangeend=${$gene}{$idvgene}{'end'};
		if (${$gene}{$idvgene}{'strand'} eq '+') {
			$geneseq=$db->seq(${$gene}{$idvgene}{'reference'}, ${$gene}{$idvgene}{'start'} => ${$gene}{$idvgene}{'end'});
			if (exists $outcode{'5'} or exists $outcode{'6'} or exists $outcode{'9'} or exists $outcode{'12'}) {
				if (${$gene}{$idvgene}{'start'}>=1) {
					$test_upstream=1;
					$upstart=${$gene}{$idvgene}{'start'}-$stream;
					if ($upstart<1) {
						$upstart=1;
					}
					$upend=${$gene}{$idvgene}{'start'}-1;
				}
				else {
					die "Error: invalid gene start: GENE $idvgene START: ", ${$gene}{$idvgene}{'start'}, "\n";
				}
				if ($test_upstream==1) {
					$upstream3000seq=$db->seq(${$gene}{$idvgene}{'reference'}, $upstart => $upend);
					$rangestart=$upstart;
				}
			}
			if (exists $outcode{'7'} or exists $outcode{'8'} or exists $outcode{'9'} or exists $outcode{'12'}) {
				if (${$gene}{$idvgene}{'end'}<=$seqlength) {
					$test_downstream=1;
					$downstart=${$gene}{$idvgene}{'end'}+1;
					$downend=${$gene}{$idvgene}{'end'}+$stream;
					if ($downend>$seqlength) {
						$downend=$seqlength;
					}
				}
				else {
					print Dumper ${$gene}{$idvgene};
					die "Error: invalid gene end: GENE $idvgene end > seqlength : ", ${$gene}{$idvgene}{'end'}, " > ", $seqlength, "\n";
					
				}
			
				if ($test_downstream==1) {
					$downstream3000seq=$db->seq(${$gene}{$idvgene}{'reference'}, $downstart => $downend);
					$rangeend=$downend;
				}
			}
		}
		elsif (${$gene}{$idvgene}{'strand'} eq '-') {
			$geneseq=$db->seq(${$gene}{$idvgene}{'reference'}, ${$gene}{$idvgene}{'end'} => ${$gene}{$idvgene}{'start'});
			if (exists $outcode{'7'} or exists $outcode{'8'} or exists $outcode{'9'} or exists $outcode{'12'}) {
				if (${$gene}{$idvgene}{'start'}>=1) {
					$test_downstream=1;
					$downstart=${$gene}{$idvgene}{'start'}-$stream;
					if ($downstart<1) {
						$downstart=1;
					}
					$downend=${$gene}{$idvgene}{'start'}-1;
				}
				else {
					die "Error2: invalid gene start: GENE $idvgene START: ", ${$gene}{$idvgene}{'start'}, "\n";
				}
				if ($test_downstream==1) {
					$downstream3000seq=$db->seq(${$gene}{$idvgene}{'reference'}, $downend => $downstart);
					$rangestart=$downstart;
				}
			}
			if (exists $outcode{'5'} or exists $outcode{'6'} or exists $outcode{'9'} or exists $outcode{'12'}) {
				if (${$gene}{$idvgene}{'end'}<=$seqlength) {
					$test_upstream=1;
					$upstart=${$gene}{$idvgene}{'end'}+1;
					$upend=${$gene}{$idvgene}{'end'}+$stream;
					if ($upend>$seqlength) {
						$upend=$seqlength;
					}
				}
				else {
					print Dumper ${$gene}{$idvgene};
					die "Error: invalid gene end: GENE $idvgene end > seqlength : ", ${$gene}{$idvgene}{'end'}, " > ", $seqlength, "\n";
				}
				if ($test_upstream==1) {
					$upstream3000seq=$db->seq(${$gene}{$idvgene}{'reference'}, $upend => $upstart);
					$rangeend=$upend;
				}
			}
		}
		
		my $genescore='.';
		$genescore=${$gene}{$idvgene}{'score'} if (exists ${$gene}{$idvgene}{'score'});
	### Gene
	#	print "Info: Write gene: $idvgene\n"; ### For test ###
		my $seqnote="Score $genescore Strand ".${$gene}{$idvgene}{'strand'};
		if (exists ${$gene}{$idvgene}{'note'} and ${$gene}{$idvgene}{'note'}=~/\S+/) {
			$seqnote.=" Note ";
			$seqnote.=${$gene}{$idvgene}{'note'};
		}
		
		if (exists $outcode{'1'}) {
			my $geneseqobj = Bio::Seq->new(
					-seq        => $geneseq,
					-id         => $idvgene,
					-display_id => $idvgene,
					-alphabet   => 'dna',
					-desc       => $seqnote
			);
			$outfile_gene->write_seq($geneseqobj) || print STDERR "Warnings: write genome $idvgene failed\n";
		}
	### UpStream
	#	print "Info: write upstream: $idvgene\n"; ### For test ###
		if (exists $outcode{'5'} or exists $outcode{'6'} and $test_upstream==1) {
	#		print "Info: write upstream only: $idvgene\n"; ### For test ###
			my $usonlyseqobj = Bio::Seq->new(
				-seq        => $upstream3000seq,
				-id         => $idvgene,
				-display_id => $idvgene,
				-alphabet   => 'dna',
				-desc       => "UpstreamOnly $stream ".$seqnote
			);
			if (exists $outcode{'5'}) {
				$outfile_upstream3000only->write_seq($usonlyseqobj) || print STDERR "Warnings: upstream.only $idvgene failed\n";
			}
	#		print "Info: write upstream with: $idvgene\n"; ### For test ###
			my $uswithseqobj = Bio::Seq->new(
				-seq        => $upstream3000seq.$geneseq,
				-id         => $idvgene,
				-display_id => $idvgene,
				-alphabet   => 'dna',
				-desc       => "WithUpstream $stream ".$seqnote
			);
			if (exists $outcode{'6'}) {
				$outfile_upstream3000with->write_seq($uswithseqobj) || print STDERR "Warnings: upstream.with $idvgene failed\n";
			}
		}
	### DownStream
	#	print "Info: write downstream: $idvgene\n"; ### For test ###
		if (exists $outcode{'7'} or exists $outcode{'8'} and $test_downstream==1) {
			my $dsonlyseqobj = Bio::Seq->new(
				-seq        => $downstream3000seq,
				-id         => $idvgene,
				-display_id => $idvgene,
				-alphabet   => 'dna',
				-desc       => "DownstreamOnly $stream".$seqnote
			);
			if (exists $outcode{'7'}) {
				$outfile_downstream3000only->write_seq($dsonlyseqobj) || print STDERR "Warnings: downstream.only $idvgene failed\n";
			}
			my $dswithseqobj = Bio::Seq->new(
				-seq        => $geneseq.$downstream3000seq,
				-id         => $idvgene,
				-display_id => $idvgene,
				-alphabet   => 'dna',
				-desc       => "WithUpstream $stream".$seqnote,
			);
			if (exists $outcode{'8'}) {
				$outfile_downstream3000with->write_seq($dswithseqobj) || print STDERR "Warnings: downstream.with $idvgene failed\n";
			}
		}
	### UpDownStream
		if (exists $outcode{'9'} and $test_upstream==1 and $test_downstream==1) {
	#		print "Info: write updownstream: $idvgene\n"; ### For test ###
			my $udsonlyseqobj = Bio::Seq->new(
				-seq        => $upstream3000seq.$geneseq.$downstream3000seq,
				-id         => $idvgene,
				-display_id => $idvgene,
				-alphabet   => 'dna',
				-desc       => "WithUpDownstream $stream".$seqnote
			);
			$outfile_updownstream3000with->write_seq($udsonlyseqobj) || print STDERR "Warnings: up/downstream.with $idvgene failed\n";
		}
		if (exists $outcode{'12'} and ($test_upstream==1) and ($test_downstream==1) ) {
			unless ($rangestart=~/^\d+$/ and  $rangeend=~/^\d+$/ and $rangeend>=$rangestart) {
				die "Error: Seq end $rangeend <= start $rangestart\n";
			}
			unless (exists $mergedregion{${$gene}{$idvgene}{'reference'}} and exists $mergedregion{${$gene}{$idvgene}{'reference'}}{$rangestart} and $mergedregion{${$gene}{$idvgene}{'reference'}}{$rangestart}=~/^\d+$/ and $mergedregion{${$gene}{$idvgene}{'reference'}}{$rangestart} >= $rangeend) {
				$mergedregion{${$gene}{$idvgene}{'reference'}}{$rangestart}= $rangeend;
			}
		}
	}
}



#Write exons sequences
if (exists $outcode{'3'}) {
	print "### 2. Info: output cDNA sequences\n";
	foreach my $indexon (sort keys %{$exons}) {
		my $exonrefer=${$exons}{$indexon}{'reference'};
		my $exonstrand=${$exons}{$indexon}{'strand'};
		my $exonscore='.';
		$exonscore=${$exons}{$indexon}{'score'} if (exists ${$exons}{$indexon}{'score'});
		#check overlaps
		my ($testoverlap, $exonarr)=&CheckOverlap(${$exons}{$indexon}{'exon'});
		unless ($testoverlap) {
			print STDERR "Error: transcript $indexon have exon overlaps\n";
			next;
		}
		my $exonseq='';
		foreach my $exonnum (@{$exonarr}) {
			$exonseq.=$db->seq($exonrefer, ${$exonnum}[0] => ${$exonnum}[1]);
		}
	
		my $genename=$indexon;
		$genename=~s/\.trans\d+.*//;
		my $transnote='';
		if (exists ${$gene}{$genename} and ${$gene}{$genename}{'note'}) {
			$transnote=" Note ".${$gene}{$genename}{'note'};
		}
		my $exonseqobj = Bio::Seq->new(
				-seq        => $exonseq,
				-id         => $indexon,
				-display_id => $indexon,
				-alphabet   => 'dna',
				-desc       => "Score $exonscore Strand $exonstrand GeneID $genename$transnote",
		);
		$exonseqobj=$exonseqobj->revcom() if ($exonstrand eq '-');
		if (exists $outcode{'3'}) {
			$outfile_cdna->write_seq($exonseqobj) || print STDERR "Warnings: write exons $indexon failed\n";
		}
	}
}



### CDS
if (exists $outcode{'2'} or exists $outcode{'4'}) {
	print "### 3. Info: output CDS sequences\n" if (exists $outcode{'2'});
	print "### 3. Info: output Protein sequences\n" if (exists $outcode{'4'});
	foreach my $indcds (sort keys %{$cds}) {
		my ($cdsseq, $translateseq)=&getProtein($indcds);

		if (exists $outcode{'2'}) {
		
		
			$outfile_cds->write_seq($cdsseq) || print STDERR "Warnings: write cds $indcds failed\n";
		}

		if (exists $outcode{'4'}) {
			$outfile_pep->write_seq($translateseq)  || print STDERR "Warnings: write CDS $indcds protein failed\n";
		}
	}
}



### 3. gff2embl

if (exists $outcode{'10'} or exists $outcode{'11'} or exists $outcode{'13'}) {
	open (EMBLOUT, " > $outprefix.temp.embl") if (exists $outcode{'13'});
	print "\n\n### Info: output EMBL\n";
	my @all_seqids=$db->get_all_primary_ids;
	@all_seqids=sort @all_seqids;
	foreach my $ind_ref (@all_seqids) {
		if (exists $outcode{'10'} or exists $outcode{'11'}) {
			next unless (exists ${$referenceids}{$ind_ref});
		}
		my $seqlength=$db->length("$ind_ref");
		unless (defined $seqlength and $seqlength=~/^\d+$/ and $seqlength>0) {
			die "Error: invalid seq length2: SeqID ", $ind_ref, "\n";
		}
		my $outfile;
		if (exists $outcode{'10'} or exists $outcode{'13'}) {
			if (exists $outcode{'10'}) {### write new file
				($outfile=$ind_ref)=~s/\|/_/g;
				close EMBLOUT if (defined fileno(EMBLOUT));
				open (EMBLOUT, " > $outdir/$outfile.temp.embl ") 
				  || die "Error: Can not write EMBL out: $outdir/$outfile.temp.embl\n";
			}
			
#			print "Test: output filename: $outdir/$outfile.embl\n"; ### For test ###
			
			&PrintEmbl('ID', "$ind_ref; SV1; linear; genomic DNA; STD; PLN; $seqlength BP.");
			&PrintEmbl('AC', "$ind_ref;");
			&PrintEmbl('FH');
			&PrintEmbl('FT', 'source', "1..".$seqlength);
			&PrintEmbl('FT', '', "/mol_type=\"genomic DNA\"");
		}
		
		if (exists ${$referenceids}{$ind_ref}) {
		foreach my $ind_num (sort {$a <=> $b} keys %{${$referenceids}{$ind_ref}}) {
### write gene
			foreach my $ind_geneid (sort keys %{${$referenceids}{$ind_ref}{$ind_num}}) {
				my ($cutstart, $cutend, $tobecut);
				###Strand: ${$gene}{$ind_geneid}{'strand'}
				unless (exists ${$gene}{$ind_geneid} and exists ${$gene}{$ind_geneid}{'strand'} and ${$gene}{$ind_geneid}{'strand'} =~/^(\+)|(\-)$/) {
					die "Error: unknown strand GENE $ind_geneid\n";
				}
				if (exists $outcode{'11'}) {
					$cutstart=${$gene}{$ind_geneid}{'start'}-$stream;
					$cutstart=1 unless ($cutstart>0);
					$cutend=${$gene}{$ind_geneid}{'end'}+$stream;
					$cutend=$seqlength if ($cutend>$seqlength);
					$tobecut=$cutstart-1;
					($outfile=$ind_geneid)=~s/\|/_/g;($outfile=$ind_geneid)=~s/:/_/g;
#					print "Test: output filename: $outdir/$outfile.embl\n"; ### For test ###
					close EMBLOUT if (defined fileno(EMBLOUT));
					open (EMBLOUT, " > $outdir/$outfile.temp.embl ") || die "Error: Can not write EMBL out: $outdir/$outfile.temp.embl\n";
					&PrintEmbl('ID', "$outfile; SV1; linear; genomic DNA; STD; PLN; $seqlength BP.");
					&PrintEmbl('AC', "$outfile;");
					&PrintEmbl('FH');
					&PrintEmbl('FT', 'source', "1..".($cutend-$cutstart+1));
					&PrintEmbl('FT', '', "/mol_type=\"genomic DNA\"");
					if (${$gene}{$ind_geneid}{'strand'} eq '+') {
						&PrintEmbl('FT', 'gene', (${$gene}{$ind_geneid}{'start'}-$tobecut)."..".(${$gene}{$ind_geneid}{'end'}-$tobecut));
					}
					elsif (${$gene}{$ind_geneid}{'strand'} eq '-') {
						&PrintEmbl('FT', 'gene', 'complement('.(${$gene}{$ind_geneid}{'start'}-$tobecut)."..".(${$gene}{$ind_geneid}{'end'}-$tobecut).')');
					}
					else {
						print STDERR "Warnings1: strand error: GENE $ind_geneid\n";
						next;
					}
					
				}
				if (exists $outcode{'10'} or exists $outcode{'13'}) {
					if (${$gene}{$ind_geneid}{'strand'} eq '+') {
						&PrintEmbl('FT', 'gene', ${$gene}{$ind_geneid}{'start'}."..".${$gene}{$ind_geneid}{'end'});
					}
					elsif (${$gene}{$ind_geneid}{'strand'} eq '-') {
						&PrintEmbl('FT', 'gene', 'complement('.${$gene}{$ind_geneid}{'start'}."..".${$gene}{$ind_geneid}{'end'}.')');
					}
					else {
						print STDERR "Warnings2:  strand error: GENE $ind_geneid\n";
						next;
					}
				}
				&PrintEmbl('FT', '', "/gene=\"$ind_geneid\"");
				if (exists ${$gene}{$ind_geneid}{'pseudogene'}) {
					my $info_pseudo=${$gene}{$ind_geneid}{'pseudogene'};
					$info_pseudo=~s/^"//;$info_pseudo=~s/"$//;
					&PrintEmbl('FT', '', "/pseudogene=\"$info_pseudo\"");
				}
				if (exists ${$gene}{$ind_geneid}{'Ontology_term'}) {
					my $info_ontology=${$gene}{$ind_geneid}{'Ontology_term'};
					my @arr_ontology=split(/,/, $info_ontology);
					foreach my $ind_ontology (sort @arr_ontology) {
						$ind_ontology=~s/^"//;$ind_ontology=~s/"$//;
						&PrintEmbl('FT', '', "/db_xref=\"$ind_ontology\"");
					}
				}
### write mRNA
				unless (exists ${$gene2mrna}{$ind_geneid}) {### check if gene got mrna
					print STDERR "Warnings: GENE $ind_geneid no mRNAs\n";
				}
				foreach my $ind_mrnaid (sort keys %{${$gene2mrna}{$ind_geneid}}) { 
					###mRNA
					my %uniqueexons=();
					my $exonidname=$ind_geneid."exon00000001";
					my ($testoverlap, $exonarr)=&CheckOverlap(${$exons}{$ind_mrnaid}{'exon'});
					unless ($testoverlap) {
						print STDERR "Error2: transcript $ind_mrnaid have exon overlaps\n";
						print Dumper ${$exons}{$ind_mrnaid}{'exon'};
						next;
					}
					my @mrnaexons=();
					foreach my $exonnum (@{$exonarr}) {
						if (exists $outcode{'10'} or exists $outcode{'13'}) {
							push (@mrnaexons, ${$exonnum}[0].'..'.${$exonnum}[1]);
						}
						if (exists $outcode{'11'}) {
							push (@mrnaexons, (${$exonnum}[0]-$tobecut).'..'.(${$exonnum}[1]-$tobecut));
						}
					}
					unless (scalar(@mrnaexons)>0) {
						print STDERR "Warnings: GENE $ind_geneid mRNA $ind_mrnaid no exons\n";
					}
					my $mrnastr='';
					if (${$gene}{$ind_geneid}{'strand'} eq '+') {
						$mrnastr='join('. join(',', @mrnaexons).')';
					}
					elsif (${$gene}{$ind_geneid}{'strand'} eq '-') {
						$mrnastr='complement(join('. join(',', @mrnaexons).'))';
					}
					else {
						print STDERR "Warnings: GENE $ind_geneid mRNA $ind_mrnaid no mRNAs\n";
						next;
					}
					&PrintEmbl('FT', 'mRNA', "$mrnastr");
					&PrintEmbl('FT', '', "/mrna=\"$ind_mrnaid\"");
#					&PrintEmbl('FT', '', "/note=\"transcript_id=$ind_mrnaid\"");
					&PrintEmbl('FT', '', "/gene=\"$ind_geneid\"");
					@mrnaexons=();
###exon
					foreach my $exonnum (@{$exonarr}) { 
						if (exists $outcode{'10'} or exists $outcode{'13'}) {
							if (${$gene}{$ind_geneid}{'strand'} eq '+') {
								&PrintEmbl('FT', 'exon', ${$exonnum}[0].'..'.${$exonnum}[1]);
							}
							else {
								&PrintEmbl('FT', 'exon', 'complement(' . ${$exonnum}[0] . '..' . ${$exonnum}[1].')');
							}
						}
						if (exists $outcode{'11'}) {
							if (${$gene}{$ind_geneid}{'strand'} eq '+') {
								&PrintEmbl('FT', 'exon', (${$exonnum}[0]-$tobecut).'..'.(${$exonnum}[1]-$tobecut));
							}
							else {
								&PrintEmbl('FT', 'exon', 'complement('.(${$exonnum}[0]-$tobecut).'..'.(${$exonnum}[1]-$tobecut).')');
							}
						}
						my $thisexonid='';
						my $testunique=0;
						
						if (exists ${$exons}{$ind_mrnaid} and exists ${$exons}{$ind_mrnaid}{'exon'} and exists ${$exons}{$ind_mrnaid}{'exon'}{${$exonnum}[0]} and exists ${$exons}{$ind_mrnaid}{'exon'}{${$exonnum}[0]}{${$exonnum}[1]}) {
							$thisexonid=${$exons}{$ind_mrnaid}{'exon'}{${$exonnum}[0]}{${$exonnum}[1]};
						}
						else {
							while ($testunique==0) {
								$thisexonid=$exonidname++;
								unless (exists $uniqueexons{$thisexonid}) {
									$testunique=1;
								}
							}
						}
						&PrintEmbl('FT', '', "/gene=\"$ind_geneid\"");
						&PrintEmbl('FT', '', "/mrna=\"$ind_mrnaid\"");
						&PrintEmbl('FT', '', "/note=\"$thisexonid\"");
					}
					###CDS
					@mrnaexons=();
					my ($testoverlap2, $cdsarr)=&CheckOverlap(${$cds}{$ind_mrnaid}{'cds'});
					unless ($testoverlap2) {
						print Dumper ${$cds}{$ind_mrnaid}{'cds'};
						die "Error: transcript $ind_mrnaid have CDS overlaps\n";
					}
					if (scalar(@{$cdsarr})>0) {
						foreach my $cdsnum (@{$cdsarr}) {
							if (exists $outcode{'10'} or exists $outcode{'13'}) {
								push (@mrnaexons, ${$cdsnum}[0].'..'.${$cdsnum}[1]);
							}
							if (exists $outcode{'11'}) {
								push (@mrnaexons, (${$cdsnum}[0]-$tobecut).'..'.(${$cdsnum}[1]-$tobecut));
							}
						}
						my $cdsstr='';
						if (${$gene}{$ind_geneid}{'strand'} eq '+') {
							$cdsstr='join('. join(',', @mrnaexons).')';
						}
						elsif (${$gene}{$ind_geneid}{'strand'} eq '-') {
							$cdsstr='complement(join('. join(',', @mrnaexons).'))';
						}
						else {
							print STDERR "Warnings: GENE $ind_geneid mRNA $ind_mrnaid no CDS\n";
							next;
						}
						&PrintEmbl('FT', 'CDS', "$cdsstr");
						&PrintEmbl('FT', '', "/gene=\"$ind_geneid\"");
						&PrintEmbl('FT', '', "/mrna=\"$ind_mrnaid\"");
#						&PrintEmbl('FT', '', "/protein_id=\"$ind_mrnaid\"");
#						&PrintEmbl('FT', '', "/note=\"transcript_id=$ind_mrnaid\"");

#						my $protseqobj={};
#						my $cdsseqobj={};
#						($cdsseqobj, $protseqobj)=&getProtein($ind_mrnaid);
#						&PrintEmbl('FT', '', '/translation="'. $protseqobj->seq .'"');
					}
				}
				if (exists $outcode{'11'}) {
					&PrintEmbl('SQ', $db->seq("$ind_ref", $cutstart, $cutend));
					close EMBLOUT;
					unless (-s "$outdir/$outfile.temp.embl") {
						die "Error: EMBL temp failed: $outdir/$outfile.temp.embl\n";
					}
					unless (&ReformatEmbl("$outdir/$outfile.temp.embl", "$outdir/$outfile.embl")) {
						die "Error: ReformatEmbl failed: $outdir/$outfile.temp.embl\n";
					}
				}
			}
		}}
		if (exists $outcode{'10'} or exists $outcode{'13'}) {
			&PrintEmbl('SQ', $db->seq("$ind_ref"));
		}
		if (exists $outcode{'10'}) {
			close EMBLOUT;
			unless (-s "$outdir/$outfile.temp.embl") {
				print STDERR "Error: EMBL temp failed: $outdir/$outfile.temp.embl\n";
				next;
			}
			unless (&ReformatEmbl("$outdir/$outfile.temp.embl", "$outdir/$outfile.embl")) {
				print STDERR "Error: ReformatEmbl failed: $outdir/$outfile.temp.embl\n";
				next;
			}
		}
	}
	if (exists $outcode{'13'}) {
		close EMBLOUT;
		unless (-s "$outprefix.temp.embl") {
			die "Error: EMBL temp failed: $outprefix.temp.embl\n";
		}
		unless (&ReformatEmbl("$outprefix.temp.embl", "$outprefix.embl")) {
			die "Error: ReformatEmbl failed: $outprefix.temp.embl\n";
		}
	}
}




if (exists $outcode{'12'}) {
	my %sequence_writen=();
	foreach my $idvref (sort keys %mergedregion) {
		my @regionunique=();
		my $region_start=0;
		my $region_end=0;
		my $test1=0;
		foreach my $num1 (sort {$a<=>$b} keys %{$mergedregion{$idvref}}) {
			if ($test1==0) {
				$region_start=$num1;
				$region_end=$mergedregion{$idvref}{$num1};
				$test1++;
			}
			else {
				if ($num1>=$region_start and $num1<=($region_end+1)) {
					if ($region_end<$mergedregion{$idvref}{$num1}) {
						$region_end=$mergedregion{$idvref}{$num1};
					}
				}
				elsif ($num1> ($region_end+1)) {
					push (@regionunique, [$region_start, $region_end]);
					$region_start=$num1;
					$region_end=$mergedregion{$idvref}{$num1};
				}
				else {
					print Dumper $mergedregion{$idvref};
					die "Error: unknown error: REF $idvref\n";
				}
			}
		}
		push (@regionunique, [$region_start, $region_end]);
		unless (exists ${$referenceids}{$idvref}) {
			die "Error: no gene on $idvref\n";
		}
		my @genenames=();
		my $test2=0;
		foreach (sort {$a<=>$b} keys %{${$referenceids}{$idvref}}) {
			foreach my $idvgenename (sort keys %{${$referenceids}{$idvref}{$_}}) {
				$genenames[$test2++]=$idvgenename;
			}
		}
		delete ${$referenceids}{$idvref};
		foreach my $idvgenename (@genenames) {
			unless (exists ${$gene}{$idvgenename}) {
				die "Error: gene id not exists: REF $idvref GENE $idvgenename\n";
			}
			unless (exists ${$gene}{$idvgenename}{'start'} and ${$gene}{$idvgenename}{'start'}=~/^\d+$/ and ${$gene}{$idvgenename}{'start'}>0) {
				die "Error: gene id start problem: REF $idvref GENE $idvgenename\n";
			}
			unless (exists ${$gene}{$idvgenename}{'end'} and ${$gene}{$idvgenename}{'end'}=~/^\d+$/ and ${$gene}{$idvgenename}{'end'}>0) {
				die "Error: gene id end problem: REF $idvref GENE $idvgenename\n";
			}
			my $genecutoff=0;
			my $test3=0;
			my @regioncoord=();
			foreach (@regionunique) {
				if (${$gene}{$idvgenename}{'start'}>=$_->[0] and ${$gene}{$idvgenename}{'end'}<=$_->[1]) {
					$test3++;
					$genecutoff=$_->[0] - 1;
					@regioncoord=($_->[0], $_->[1]);
				}
			}
			unless ($test3==1) {
				die "Error: unique region problem: REF $idvref GENE $idvgenename\n";
			}
			if ($genecutoff<0) {
				die "Error: unique region start: REF $idvref GENE $idvgenename\n";
			}
			### write sequence
			my $idvrefname=${$gene}{$idvgenename}{'reference'}.'_'.$regioncoord[0].'_'.$regioncoord[1];
			
			unless (exists $sequence_writen{$idvrefname}) {
				my $mergedgeneseq=$db->seq(${$gene}{$idvgenename}{'reference'}, $regioncoord[0] => $regioncoord[1]);
				my $mergedgeneseqobj = Bio::Seq->new(
				-seq        => $mergedgeneseq,
				-id         => $idvrefname,
				-display_id => $idvrefname,
				-alphabet   => 'dna',
				-desc       => "Source:".${$gene}{$idvgenename}{'reference'}.":".$regioncoord[0]."-".$regioncoord[1],
												);
				$outfile_updownstream3000with->write_seq($mergedgeneseqobj) || print STDERR "Warnings: failed to write sequence region: $idvrefname\n";
				$sequence_writen{$idvrefname}++;
			}
			${$gene}{$idvgenename}{'start'}=${$gene}{$idvgenename}{'start'}-$genecutoff;
			${$gene}{$idvgenename}{'end'}=${$gene}{$idvgenename}{'end'}-$genecutoff;
			${$gene}{$idvgenename}{'reference'}=${$gene}{$idvgenename}{'reference'}.'_'.$regioncoord[0].'_'.$regioncoord[1];
			${$referenceids}{$idvrefname}{${$gene}{$idvgenename}{'start'}}{$idvgenename}++;
			unless (exists ${$gene2mrna}{$idvgenename}) {
				die "Error: gene do not have mRNA: $idvgenename\n"
			}
			foreach my $idvmrnaname (sort keys %{${$gene2mrna}{$idvgenename}}) {
				### mRNA
				${$mrnas}{$idvmrnaname}{'start'}=${$mrnas}{$idvmrnaname}{'start'}-$genecutoff;
				${$mrnas}{$idvmrnaname}{'end'}=${$mrnas}{$idvmrnaname}{'end'}-$genecutoff;
				${$mrnas}{$idvmrnaname}{'reference'}=$idvrefname;
				### exon
				if (exists ${$exons}{$idvmrnaname}) {
					${$exons}{$idvmrnaname}{'reference'}=$idvrefname;
					if (exists ${$exons}{$idvmrnaname}{'exon'}) {
						my @temparr1=sort {$a<=>$b} keys %{${$exons}{$idvmrnaname}{'exon'}};
						foreach my $x (@temparr1) {
							my @temparr2=sort {$a<=>$b} keys %{${$exons}{$idvmrnaname}{'exon'}{$x}};
							foreach my $y (@temparr2) {
								my $exonid=${$exons}{$idvmrnaname}{'exon'}{$x}{$y};
								delete ${$exons}{$idvmrnaname}{'exon'}{$x}{$y};
								${$exons}{$idvmrnaname}{'exon'}{$x-$genecutoff}{$y-$genecutoff}=$exonid;
							}
							delete ${$exons}{$idvmrnaname}{'exon'}{$x} if (scalar(keys %{${$exons}{$idvmrnaname}{'exon'}{$x}})==0);
						}
					}
				}
				### CDS
				if (exists ${$cds}{$idvmrnaname}) {
					${$cds}{$idvmrnaname}{'reference'}=$idvrefname;
					if (exists ${$cds}{$idvmrnaname}) {
						my @temparr1=sort {$a<=>$b} keys %{${$cds}{$idvmrnaname}{'cds'}};
						foreach my $x (@temparr1) {
							my @temparr2=sort {$a<=>$b} keys %{${$cds}{$idvmrnaname}{'cds'}{$x}};
							foreach my $y (@temparr2) {
								my $exonid=${$cds}{$idvmrnaname}{'cds'}{$x}{$y};
								delete ${$cds}{$idvmrnaname}{'cds'}{$x}{$y};
								${$cds}{$idvmrnaname}{'cds'}{$x-$genecutoff}{$y-$genecutoff}=$exonid;
							}
							delete ${$cds}{$idvmrnaname}{'cds'}{$x} if (scalar(keys %{${$cds}{$idvmrnaname}{'cds'}{$x}})==0);
						}
					}
					if (exists ${$cds}{$idvmrnaname}{'phase'}) {
						my @temparr1=sort {$a<=>$b} keys %{${$cds}{$idvmrnaname}{'phase'}};
						foreach my $x (@temparr1) {
							my @temparr2=sort {$a<=>$b} keys %{${$cds}{$idvmrnaname}{'phase'}{$x}};
							foreach my $y (@temparr2) {
								my $phase=${$cds}{$idvmrnaname}{'phase'}{$x}{$y};
								delete ${$cds}{$idvmrnaname}{'phase'}{$x}{$y};
								${$cds}{$idvmrnaname}{'phase'}{$x-$genecutoff}{$y-$genecutoff}=$phase;
							}
							delete ${$cds}{$idvmrnaname}{'phase'}{$x} if (scalar(keys %{${$cds}{$idvmrnaname}{'phase'}{$x}})==0);
						}
					}
				}
			}
		}
	}
#	print Dumper $referenceids; ### For test ###
	unless (WriteGff3($outputgff3, $referenceids, $gene2mrna, $gene, $mrnas, $exons, $cds)) {
		die "Error: failed to write output GFF3\n";
	}
}



#####################################################################
####################### sub function ################################
#####################################################################



sub ReformatEmbl {
	my ($REemblin, $REemblout)=@_;
	
	unless (-s $REemblin) {
		print STDERR "Error: EMBL temp file not found: $REemblin\n";
		return 0;
	}
	unlink $REemblout if (-e $REemblout);

	my $REemblinobj=Bio::SeqIO->new( -format => 'EMBL', -file => $REemblin );
	my $REembloutobj=Bio::SeqIO->new( -format => 'EMBL', -file => ">$REemblout" );
	while (my $REseqobject=$REemblinobj->next_seq()) {
#			print "Test: SeqIO object format\n"; print Dumper $seqobject; print "\n"; ### For test ###
		$REembloutobj->write_seq($REseqobject);
	}
	if ( -s $REemblout) {
		unlink $REemblin if (-e $REemblin);
		return 1;
	}
	else {
		print STDERR "Error: ReformatEmbl failed: $REemblin";
		return 0;
	}
}



### Print EMBL
sub PrintEmbl {
	my ($PElinestart, $PEstring, $PEnote)=@_;
	if ($PElinestart=~/^(ID)|(AC)$/) {
		print EMBLOUT "\n" if ($printID>0);
		print EMBLOUT $PElinestart, ' ' x 3 , $PEstring, "\n";
		print EMBLOUT "XX\n";
		$printID++;
	}
	elsif ($PElinestart=~/^(DE)$/) {
		print EMBLOUT $PElinestart, ' ' x 3 , $PEstring, "\n";
		print EMBLOUT "XX\n";
	}
	elsif ($PElinestart=~/^(FH)$/) {
		print EMBLOUT "FH   Key             Location/Qualifiers\n";
		print EMBLOUT "FH\n";
	}
	elsif ($PElinestart=~/^(FT)$/) {
		print EMBLOUT 'FT   ', $PEstring, ' ' x (16-length($PEstring)), $PEnote, "\n";
	}
	elsif ($PElinestart=~/^(SQ)$/) {
		print EMBLOUT "XX\n";
		my $seqlength=length($PEstring);
		my $maxlength=length($seqlength) + 5 + 66;
		my $numA=$PEstring=~tr/Aa/Aa/;
		my $numT=$PEstring=~tr/Tt/Tt/;
		my $numC=$PEstring=~tr/Cc/Cc/;
		my $numG=$PEstring=~tr/Gg/Gg/;
		my $numO=$seqlength-$numA-$numC-$numG-$numT;
		print EMBLOUT "SQ   Sequence $seqlength BP; $numA A; $numC C; $numG G; $numT T; $numO other;\n";

		for (my $i=0; $i<$seqlength; $i+=60) {
			print EMBLOUT '     '; my $thislineused=5; ### 5 spaces here
			for (my $j=0; $j<60; $j+=10) {
				my $x=$i+$j;
				my $numremainbases=$seqlength-$x;
				if ($numremainbases>=10) {
					print EMBLOUT substr($PEstring, $x, 10), ' ';
					$thislineused+=11;
				}
				elsif ($numremainbases<10 and $numremainbases>0) {
					print EMBLOUT substr($PEstring, $x, $numremainbases);
					$thislineused+=$numremainbases;
				}
#				else {
#					die "Error501\n;"
#				}
			}
			
			if (($seqlength-$i)<=60) {
				print EMBLOUT ' ' x ($maxlength-$thislineused-length($seqlength));
				print EMBLOUT $seqlength, "\n";
			}
			else {
				my $marknum=$i+60;
				print EMBLOUT ' ' x ($maxlength-$thislineused-length($marknum));
				print EMBLOUT $marknum, "\n";
			}
		}
		print EMBLOUT '//';
	}
	else {
		die "Error510\n"
	}
}



### Get Protein sequence using %CDS
### getProtein($mrna_id)
### Return: CDS Protein seqobj for that mRNA ID, and desc
### Dependency: 
### Global: $cds, $gene
### Note:
sub getProtein {
	my $GPmrnaid=shift;

	my $GPsubinfo="SUB(getProtein)";
	my $GPcdsseq='';
	my $GPtranslateseq='';
	my @GPtemparrcds=();

	my $GPcdsrefer=${$cds}{$GPmrnaid}{'reference'};
	#check strand
	my $GPcdsstrand=${$cds}{$GPmrnaid}{'strand'};
	#check scores
	my $GPcdsscore='.';
	$GPcdsscore=${$cds}{$GPmrnaid}{'score'} if (exists ${$cds}{$GPmrnaid} and exists ${$cds}{$GPmrnaid}{'score'});

	#check overlaps
	my ($GPtestoverlap, $GPcdsarr)=&CheckOverlap(${$cds}{$GPmrnaid}{'cds'});
	unless ($GPtestoverlap) {
		die "Error: transcript $GPmrnaid have cds overlaps\n";
	}
	
	my @GPphases=();
	foreach my $GPcdsnum (@{$GPcdsarr}) {
		$GPcdsseq.=$db->seq($GPcdsrefer, ${$GPcdsnum}[0] => ${$GPcdsnum}[1]);
		if (exists ${$cds}{$GPmrnaid}{'phase'} and exists ${$cds}{$GPmrnaid}{'phase'}{${$GPcdsnum}[0]} and exists ${$cds}{$GPmrnaid}{'phase'}{${$GPcdsnum}[0]}{${$GPcdsnum}[1]} and defined ${$cds}{$GPmrnaid}{'phase'}{${$GPcdsnum}[0]}{${$GPcdsnum}[1]} and ${$cds}{$GPmrnaid}{'phase'}{${$GPcdsnum}[0]}{${$GPcdsnum}[1]}=~/^\d+$/ and ${$cds}{$GPmrnaid}{'phase'}{${$GPcdsnum}[0]}{${$GPcdsnum}[1]}>=0 and ${$cds}{$GPmrnaid}{'phase'}{${$GPcdsnum}[0]}{${$GPcdsnum}[1]}<=2) {
			push (@GPphases, ${$cds}{$GPmrnaid}{'phase'}{${$GPcdsnum}[0]}{${$GPcdsnum}[1]});
		}
		else {
			print STDERR "Warnings: unknown mRNA ($GPmrnaid) CDS (${$GPcdsnum}[0]-${$GPcdsnum}[1]) phase\n";
		}
	}
	
	my $GPgenename=$GPmrnaid;
	$GPgenename=~s/\..*//;
	my $GPtransnote='';
	if (exists ${$gene}{$GPgenename} and ${$gene}{$GPgenename}{'note'}) {
		$GPtransnote=" Note ".${$gene}{$GPgenename}{'note'};
	}
	my $GPcdsseqobj = Bio::Seq->new(
			-seq        => $GPcdsseq,
			-id         => $GPmrnaid,
			-display_id => $GPmrnaid,
			-alphabet   => 'dna',
			-desc       => "Score $GPcdsscore Strand $GPcdsstrand GeneID $GPgenename$GPtransnote",
	);
	$GPcdsseqobj=$GPcdsseqobj->revcom() if ($GPcdsstrand eq '-');

	my $GPfirstphase=0;
	if (scalar(@{$GPcdsarr}) == scalar(@GPphases) and scalar(@GPphases)>0) {
		if ($GPcdsstrand eq '+') {
			$GPfirstphase=$GPphases[0];
		}
		elsif ($GPcdsstrand eq '-') {
			$GPfirstphase=$GPphases[-1];
		}
		if (1) {
			&GetCdsPhase($GPcdsarr, $GPmrnaid, $GPcdsstrand, \@GPphases);
		}
	}
	else {
		print STDERR "Warnings: mRNA ($GPmrnaid) 5end CDS phase was set to 0\n";
	}
	
	if (0) {###Test intact CDS
		unless ((length($GPcdsseq)-$GPfirstphase)%3 ==0) {
			print STDERR "Warnings: mRNA ($GPmrnaid) are not in 3* $GPtransnote\n";
		}
	}
	
	for (my $GPi=($GPfirstphase+1); $GPi<=($GPcdsseqobj->length() - 2); $GPi+=3) {
		my $GPcodon='';
		$GPcodon=$GPcdsseqobj->subseq($GPi, $GPi+2);
		$GPtranslateseq.=Codon2AA($GPcodon);
	}
	my $GPpeptideobj=Bio::Seq->new(
			-seq        => $GPtranslateseq,
			-id         => $GPmrnaid,
			-display_id => $GPmrnaid,
			-alphabet   => 'protein',
			-desc       => "Score $GPcdsscore Strand $GPcdsstrand GeneID $GPgenename$GPtransnote",
	);
	
	if ($verbose) {###Test intact protein
		if ($GPtranslateseq=~/^M/) { print STDERR "CDS\t$GPmrnaid:\tM";}else{print STDERR "CDS\t$GPmrnaid:\tnoM";}
		if ($GPtranslateseq=~/\*$/) { print STDERR "\t*";}else{print STDERR "\tno*";}
		my $GPnumstop=$GPtranslateseq=~s/\*/\*/g;
		print STDERR "\t$GPnumstop\t$GPtransnote\n";
	}
	
	$GPcdsseq='';
	$GPtranslateseq='';
	return ($GPcdsseqobj, $GPpeptideobj);
}


sub CheckOverlap {
	my $COborderhash=shift;
	my @COarr=();
	if (scalar(keys %{$COborderhash})<1) {
		return (1, \@COarr);
	}
	my $COx=0;
	foreach my $COstart (sort {$a<=>$b} keys %{$COborderhash}) {
		my @COarrend=sort {$a<=>$b} keys %{${$COborderhash}{$COstart}};
		unless (scalar(@COarrend)==1 and defined $COarrend[0] and $COarrend[0]=~/^\d+/ and $COarrend[0] >0) {
			return 0;
		}
		@{$COarr[$COx++]}=($COstart, $COarrend[0]);
	}
	for (my $COi=0; $COi<scalar(@COarr); $COi++) {
		for (my $COj=$COi+1; $COj<scalar(@COarr); $COj++) {
			if (($COarr[$COi][0]>=$COarr[$COj][0] and $COarr[$COi][0]<=$COarr[$COj][1]) or ($COarr[$COi][1]>=$COarr[$COj][0] and $COarr[$COj][1]<=$COarr[$COi][1]) or ($COarr[$COi][0]<=$COarr[$COj][0] and $COarr[$COi][1]>=$COarr[$COj][1]) or ($COarr[$COi][0]>=$COarr[$COj][0] and $COarr[$COi][1]<=$COarr[$COj][1])) {
				print STDERR "Error: feature overlap ($COarr[$COi][0] $COarr[$COi][1]) and ($COarr[$COj][0] $COarr[$COj][1])\n";
				return 0;
			}
		}
	}
	return (1, \@COarr);
}



### GetCDS phase
### GetCdsPhase($CDSarr, 5/3)
### Global:
### Dependency:
### Note: 
sub GetCdsPhase {
	my ($GCPcdsarr, $GCPmrnaid, $GCPstrand, $GCPexistingphase)=@_;
	
	my $GCPsubinfo='SUB(GetCdsPhase)';
	my @GCPphase=();
	my $GCPwhichend=5;
	$GCPmrnaid='unknown' unless (defined $GCPmrnaid);
	
#	print Dumper $GCPcdsarr; ### For test ###
	my $GCPphase=0;
	if ($GCPstrand eq '+') {
		$GCPphase=${$GCPexistingphase}[0];
		if ($GCPwhichend==5) {
			foreach my $GCPindcds (@{$GCPcdsarr}) {
				push (@GCPphase,  $GCPphase);
				$GCPphase=3-($GCPindcds->[1]-$GCPindcds->[0]+1-$GCPphase)%3;
				$GCPphase=0 if ($GCPphase==3);
			}
		}
		$GCPphase=$GCPphase[0];
	}
	elsif ($GCPstrand eq '-') {
		$GCPphase=${$GCPexistingphase}[-1];
		if ($GCPwhichend==5) {
			for (my $GCPi=(scalar(@{$GCPcdsarr})-1); $GCPi>=0; $GCPi--) {
				unshift (@GCPphase,  $GCPphase);
	#			print "Start: ${$GCPcdsarr}[$GCPi][0]\n"; ### For test ###
	#			print "End: ${$GCPcdsarr}[$GCPi][1]\n"; ### For test ###
				$GCPphase=3- ((${$GCPcdsarr}[$GCPi][1]-${$GCPcdsarr}[$GCPi][0] + 1 - $GCPphase)%3);
				$GCPphase=0 if ($GCPphase==3);
			}
		}
		$GCPphase=$GCPphase[-1];
	}
	
	for (my $GCPj=0; $GCPj<scalar(@{$GCPcdsarr}); $GCPj++) {
		unless (${$GCPexistingphase}[$GCPj] == ${$GCPexistingphase}[$GCPj]) {
			print STDERR $GCPsubinfo, "Warnings: mRNA $GCPmrnaid CDS ${$GCPcdsarr}{$GCPj}[0] - ${$GCPcdsarr}{$GCPj}[1] original phase ${$GCPexistingphase}[$GCPj] => ${$GCPexistingphase}[$GCPj]\n";
		}
	}
}
