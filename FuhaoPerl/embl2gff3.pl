#!/usr/bin/env perl
use strict;
use warnings;
use Data::Dumper qw /Dumper/;
use FuhaoPerl5Lib::GffKit qw/WriteGff3/;
use FuhaoPerl5Lib::FastaKit qw/Frame3Translation SeqRevComp/;
use Bio::DB::Fasta;
use List::Util qw /min max/;
use constant USAGE =><<EOH;

usage: $0 input.embl output.gff3

v20170407

EOH
die USAGE if (scalar(@ARGV) !=2 or $ARGV[0] eq '-h' or $ARGV[0] eq '--help');



my $emblinput=$ARGV[0];
my $outgff3=$ARGV[1];
my $debug=0;



my ($test, $refs, $gene2mrna, $gene, $mrnas, $exons, $cds)=&ReadEMBL($ARGV[0], $emblinput);
unless ($test) {
	die "Error: read EMBL file failed\n";
}
if ($debug or 0) {### For test ###
	print "Test: \$refs\n"; print Dumper $refs; print "\n";
}
if ($debug or 0) {### For test ###
	print "Test: \$gene2mrna\n"; print Dumper $gene2mrna; print "\n";
}
if ($debug or 0) {### For test ###
	print "Test: \$gene\n"; print Dumper $gene; print "\n";
}
if ($debug or 0) {### For test ###
	print "Test: \$mrnas\n"; print Dumper $mrnas; print "\n";
}
if ($debug or 0) {### For test ###
	print "Test: \$exons\n"; print Dumper $exons; print "\n";
}
if ($debug or 0) {### For test ###
	print "Test: \$cds\n"; print Dumper $cds; print "\n";
}

unless (WriteGff3($outgff3, $refs, $gene2mrna, $gene, $mrnas, $exons, $cds)) {
	die "Error: write GFF3 file failed\n";
}



### Read GFF3 into hash
###
###
### %referenceids  => ( $reference_id => $gene_start_pos => $gene_id => num++)
### %gene2mrna     => ( $gene_id => $mrna_id => num++ )
### %gene=($geneid => ('reference' => $arr[0],
###                    'start'     => $arr[3], 
###                    'end'       => $arr[4],
###                    'strand'    => $arr[6],
###                    'score'     => $arr[5],
###                    'note'      => $note ###Not necessarily exist
###                   )
###       )
### %mrnas=($geneid => ('reference' => $arr[0],
###                     'start'      => $arr[3], 
###                     'end'        => $arr[4],
###                     'strand'     => $arr[6],
###                     'score'      => $arr[5],
###                     'note'       => $note ###Not necessarily exist
###                    )
###       )
### %exon=($mrnaid => ('reference' => $arr[0],
###                    'exon' => ({$arr[3]} => ($arr[4] => $exonid)),
###                    'strand'    => $arr[6],
###                    'score'     => $arr[5]
###                    )
###       )
### %cds=($mrnaid => ('reference' => $arr[0],
###                   'cds'       => ({$arr[3]} => ($arr[4] => num++)),
###                   'strand'    => $arr[6],
###                   'score'     => $arr[5],
###                   'phase'     => ({$arr[3]} => ($arr[4] => $arr[7]))
###                   )
###       )
sub ReadEMBL {
	my $REemblfile=shift;
	
	my $REsubinfo='SUB(ReadEMBL)';
	my $REref2posgene=();
	my $REref2gene={};
	my $REgene2mrna={};
	my $REgene={};
	my $REmrnas={};
	my %REgenedirection=(); ### trust gene and exon strand NOT CDS or mRNA
	my $REexons={};
	my $REcds={};
	my $RElinenum=0;
	my $REseqstart=0;
	my $REtotalIDlines=0;
	my $REi=0;
	my %REarray=();
	my $REtemp_fasta=$REemblfile.'.fa';
	my $REdb;
	my %REorfinfo=();
	my $REgetcdsbyarray=0;###set to 0 to use exon coordinates and 1 to use CDS coordinates, 0 is better and 1 for find errors)
	my $REseqterminal=0;
	my %RECDSphases=();
	local *REEMBLFILE; local *RETEMPFASTA;
	
	
	unless (defined $REemblfile and -s $REemblfile) {
		print STDERR $REsubinfo, "Error: invalid EMBL file: $REemblfile\n";
		return 0;
	}
	unlink $REtemp_fasta if (-e $REtemp_fasta);
	
	close REEMBLFILE if (defined fileno(REEMBLFILE));
	unless (open (REEMBLFILE, " < $REemblfile ")) {
		print STDERR $REsubinfo, "Error: can not open EMBL file: $REemblfile\n";
		return 0;
	}
	close RETEMPFASTA if (defined fileno(RETEMPFASTA));
	unless (open(RETEMPFASTA, " > $REtemp_fasta")) {
		print STDERR $REsubinfo, "Error: can not write file: $REemblfile\n";
		return 0;
	}
	
	while (my $REline=<REEMBLFILE>) {
		chomp $REline;
		$RElinenum++;
		
		my $RErefname;
		my $REreflent;
		
		unless ($REline=~/(^ID)|(^FT)|(^\/\/)|(^SQ)|(^XX)|(^\s+\S+.*\d+\s*)/) {
			next;
		}
		if ($REline=~/^\/\//) {
			$REseqterminal=1;
			$REtotalIDlines=0;
			$REseqstart=0;
			#process array
			if ($debug) {### For test ###
#				print $REsubinfo, "test: \%REarray\n"; print Dumper \%REarray; print "\n" ;
			}
			$RErefname='NaN';
			$REreflent=0;
			foreach my $REy (sort {$a<=>$b} keys %REarray) {
				unless (scalar(@{$REarray{$REy}})>0) {
					### do sth
					next;
				}
				if ($REy==0) { ### SEQID
					if ($RErefname eq 'NaN') {
						$RErefname=$REarray{$REy}[0];
					}
					else {
						print STDERR $REsubinfo, "Error: different seq ID: $RErefname\n";
						return 0;
					}
					if ($REreflent==0) {
						$REreflent=$REarray{$REy}[1];
					}
					else {
						print STDERR $REsubinfo, "Error: different seq length: SEQ $RErefname LEN $REreflent\n";
						return 0;
					}
					next;
				}
				if (($RErefname eq 'NaN') or ($REreflent==0)) {
					
					print STDERR $REsubinfo, "Error: different seq length: SEQ $RErefname LEN $REreflent\n";
					return 0;
				}
				my $REfeatures=&GetGeneHash($REarray{$REy});
				if ($debug) {### For test ###
					print $REsubinfo, "Test: \$REfeatures\n"; print Dumper $REfeatures; print "\n" ; 
				}
				if ($REarray{$REy}[0] eq 'source') {
					next;
				}
				elsif ($REarray{$REy}[0] eq 'gene') {### GENE
					my $REgenestrand; my $REgenestart; my $REgeneend; my $REgeneid; 
					my $REgenescore='.'; my $REgenenote='';
					if ($REarray{$REy}[1] =~/^(\d+)\.\.(\d+)$/) {
						$REgenestrand='+';
						$REgenestart=$1;
						$REgeneend=$2;
					}
					elsif ($REarray{$REy}[1] =~/^complement\((\d+)\.\.(\d+)\)$/) {
						$REgenestrand='-';
						$REgenestart=$1;
						$REgeneend=$2;
					}
					else {
						print Dumper $REarray{$REy}, "\n";
						print STDERR $REsubinfo, "Error: unknown gene feature\n";
						return 0;
					}
					unless ($REgenestart<=$REgeneend) {
						print Dumper $REarray{$REy}, "\n";
						print STDERR $REsubinfo, "Error: GENE end<start\n";
						return 0;
					}
					if (exists ${$REfeatures}{'gene'}) {
						$REgeneid=${$REfeatures}{'gene'};
						$REgenedirection{$REgeneid}{$REgenestrand}++;
					}
					else {
						print Dumper $REarray{$REy}, "\n";
						print STDERR $REsubinfo, "Error: unknown gene ID\n";
						return 0;
					}
					if (exists ${$REgene}{$REgeneid}) {
						print Dumper $REarray{$REy}, "\n";
						print STDERR $REsubinfo, "Error: duplicated gene ID\n";
						return 0;
					}
					${$REref2gene}{$RErefname}{$REgeneid}++;
					${$REgene}{$REgeneid}{'reference'}=$RErefname;
					${$REgene}{$REgeneid}{'start'}=$REgenestart;
					${$REgene}{$REgeneid}{'end'}=$REgeneend;
					${$REgene}{$REgeneid}{'strand'}=$REgenestrand;
					
					if (exists ${$REfeatures}{'score'} and ${$REfeatures}{'score'}=~/^\d+$/) {
						${$REgene}{$REgeneid}{'score'}=${$REfeatures}{'score'};
					}
					if (exists ${$REfeatures}{'note'}) {
						${$REgene}{$REgeneid}{'note'}=${$REfeatures}{'note'}
					}
				}
				elsif ($REarray{$REy}[0] eq 'mRNA') {### mRNA
					my $REmrnastrand; my $REmrnastart; my $REmrnaend; my $REmrnaid; my $REgeneid;
					my $REmrnascore='.'; my $REmrnanote='';
					if ($REarray{$REy}[1] =~/^(\d+)\.\.(\d+)$/) {
						$REmrnastrand='+';
						$REmrnastart=$1;
						$REmrnaend=$2;
					}
					elsif ($REarray{$REy}[1] =~/^join\((\d+)\.\..*\.\.(\d+)\)$/) {
						$REmrnastrand='+';
						$REmrnastart=$1;
						$REmrnaend=$2;
					}
					elsif ($REarray{$REy}[1] =~/^complement\((\d+)\.\.(\d+)\)$/) {
						$REmrnastrand='-';
						$REmrnastart=$1;
						$REmrnaend=$2;
					}
					elsif ($REarray{$REy}[1] =~/^complement\(join\((\d+)\.\..*\.\.(\d+)\)\)$/) {
						$REmrnastrand='-';
						$REmrnastart=$1;
						$REmrnaend=$2;
					}
					else {
						print Dumper $REarray{$REy}, "\n";
						print STDERR $REsubinfo, "Error: unknown gene feature\n";
						return 0;
					}
					unless ($REmrnastart<=$REmrnaend) {
						print Dumper $REarray{$REy}, "\n";
						print STDERR $REsubinfo, "Error: mRNA end<start\n";
						return 0;
					}
					if (exists ${$REfeatures}{'mrna'}) {
						$REmrnaid=${$REfeatures}{'mrna'};
					}
					else {
						print Dumper $REarray{$REy}, "\n";
						print STDERR $REsubinfo, "Error: unknown mRNA ID\n";
						return 0;
					}
					if (exists ${$REmrnas}{$REmrnaid}) {
						print Dumper $REarray{$REy}, "\n";
						print STDERR $REsubinfo, "Error: duplicated mRNA ID\n";
						return 0;
					}
					if (exists ${$REfeatures}{'gene'}) {
						$REgeneid=${$REfeatures}{'gene'};
						${$REgene2mrna}{$REgeneid}{$REmrnaid}++;
						${$REref2gene}{$RErefname}{$REgeneid}++;
					}
					else {
						print Dumper $REarray{$REy}, "\n";
						print STDERR $REsubinfo, "Error: unknown mRNA parent gene ID\n";
						return 0;
					}
					${$REmrnas}{$REmrnaid}{'reference'}=$RErefname;
					${$REmrnas}{$REmrnaid}{'start'}=$REmrnastart;
					${$REmrnas}{$REmrnaid}{'end'}=$REmrnaend;
					${$REmrnas}{$REmrnaid}{'strand'}=$REmrnastrand;
					
					if (exists ${$REfeatures}{'score'} and ${$REfeatures}{'score'}=~/^\d+$/) {
						${$REmrnas}{$REmrnaid}{'score'}=${$REfeatures}{'score'};
					}
					if (exists ${$REfeatures}{'note'}) {
						${$REmrnas}{$REmrnaid}{'note'}=${$REfeatures}{'note'}
					}
				}
				elsif ($REarray{$REy}[0] eq 'exon') {### exon
					my $REexonstrand; my $REexonstart; my $REexonend; my $REexonid; 
					my $REexonscore='.'; my $REexonnote='';
					my $REgeneid;
					my $REmrnaid; 
					if ($REarray{$REy}[1] =~/^(\d+)\.\.(\d+)$/) {
						$REexonstrand='+';
						$REexonstart=$1;
						$REexonend=$2;
					}
					elsif ($REarray{$REy}[1] =~/^complement\((\d+)\.\.(\d+)\)$/) {
						$REexonstrand='-';
						$REexonstart=$1;
						$REexonend=$2;
					}
					else {
						print Dumper $REarray{$REy}, "\n";
						print STDERR $REsubinfo, "Error: unknown exon feature\n";
						return 0;
					}
					unless ($REexonstart<=$REexonend) {
						print Dumper $REarray{$REy}, "\n";
						print STDERR $REsubinfo, "Error: EXON end<start\n";
						return 0;
					}
					if (exists ${$REfeatures}{'gene'}) {
						$REgeneid=${$REfeatures}{'gene'};
						$REgenedirection{$REgeneid}{$REexonstrand}++;
						${$REref2gene}{$RErefname}{$REgeneid}++;
						$REgenedirection{$REgeneid}{$REexonstrand}++;
					}
					else {
						print Dumper $REarray{$REy}, "\n";
						print STDERR $REsubinfo, "Error: unknown EXON gene ID\n";
						return 0;
					}
					if (exists ${$REfeatures}{'mrna'}) {
						$REmrnaid=${$REfeatures}{'mrna'};
					}
					else {
						print Dumper $REarray{$REy}, "\n";
						print STDERR $REsubinfo, "Error: unknown EXON mRNA ID\n";
						return 0;
					}
					${$REgene2mrna}{$REgeneid}{$REmrnaid}++;
					if (exists ${$REfeatures}{'note'}) {#EXON_ID
						$REexonid=${$REfeatures}{'note'};
					}
					else {
						print Dumper $REarray{$REy}, "\n";
						print STDERR $REsubinfo, "Error: unknown EXON ID\n";
						return 0;
					}
					if (exists ${$REfeatures}{'phase'}) {###phase
						unless (${$REfeatures}{'phase'}=~/^(0)|(1)|(2)$/) {
							print Dumper $REarray{$REy}, "\n";
							print STDERR $REsubinfo, "Error: invalid phase\n";
							return 0;
						}
						if (exists $RECDSphases{$REmrnaid} and (exists $RECDSphases{$REmrnaid}{$REexonstart} or exists $RECDSphases{$REmrnaid}{$REexonend})) {
							print Dumper $REarray{$REy}, "\n";
							print STDERR $REsubinfo, "Error: existing phase\n";
							return 0;
						}
						$RECDSphases{$REmrnaid}{$REexonstart}=${$REfeatures}{'phase'};
						$RECDSphases{$REmrnaid}{$REexonend}=${$REfeatures}{'phase'};
					}
					${$REexons}{$REmrnaid}{'reference'}=$RErefname;
					${$REexons}{$REmrnaid}{'exon'}{$REexonstart}{$REexonend}=$REexonid;
					${$REexons}{$REmrnaid}{'strand'}=$REexonstrand;
					if (exists ${$REfeatures}{'score'} and ${$REfeatures}{'score'}=~/^\d+$/) {
						${$REexons}{$REmrnaid}{'score'}=${$REfeatures}{'score'};
					}
				}
				elsif ($REarray{$REy}[0] eq 'CDS') {### CDS
					my $REcdsstrand; my $REcdsstart; my $REcdsend; my $REcdsid; 
					my $REcdsscore='.'; my $REcdsnote='';
					my $REgeneid;
					my $REmrnaid;
					my @REcdsarr=();
					if ($REarray{$REy}[1] =~/^(\d+\.\.\d+)$/) {
						$REcdsstrand='+';
						push (@REcdsarr, $1);
					}
					elsif ($REarray{$REy}[1] =~/^join\(\d+.*,.*\d+\)$/) {
						$REcdsstrand='+';
						my $REstr=$REarray{$REy}[1];
						$REstr=~s/^join\(//; $REstr=~s/\)$//;
						@REcdsarr=split(/,/, $REstr);
					}
					elsif ($REarray{$REy}[1] =~/^complement\(\d+\.\.\d+\)$/) {
						$REcdsstrand='-';
						my $REstr=$REarray{$REy}[1];
						$REstr=~s/^complement\(//; $REstr=~s/\)$//;
						@REcdsarr=split(/,/, $REstr);
					}
					elsif ($REarray{$REy}[1] =~/^complement\(join\(\d+.*,.*\d+\)\)$/) {
						$REcdsstrand='-';
						my $REstr=$REarray{$REy}[1];
						$REstr=~s/^complement\(join\(//; $REstr=~s/\)\)$//;
						@REcdsarr=split(/,/, $REstr);
					}
					else {
						print Dumper $REarray{$REy}, "\n";
						print STDERR $REsubinfo, "Error: unknown CDS feature\n";
						return 0;
					}
					if (exists ${$REfeatures}{'gene'}) {
						$REgeneid=${$REfeatures}{'gene'};
						${$REref2gene}{$RErefname}{$REgeneid}++;
					}
					else {
						print Dumper $REarray{$REy}, "\n";
						print STDERR $REsubinfo, "Error: unknown CDS gene ID\n";
						return 0;
					}
					if (exists ${$REfeatures}{'mrna'}) {
						$REmrnaid=${$REfeatures}{'mrna'};
					}
					else {
						print Dumper $REarray{$REy}, "\n";
						print STDERR $REsubinfo, "Error: unknown CDS mRNA ID\n";
						return 0;
					}
					${$REgene2mrna}{$REgeneid}{$REmrnaid}++;
					${$REcds}{$REmrnaid}{'reference'}=$RErefname;
					${$REcds}{$REmrnaid}{'strand'}=$REcdsstrand;
					foreach my $REcdsrange (@REcdsarr) {
						if ($REcdsrange =~ /^(\d+)\.\.(\d+)$/) {
							$REcdsstart=$1;
							$REcdsend=$2;
							push (@{$REorfinfo{$REmrnaid}{'array'}}, $REcdsstart);
							push (@{$REorfinfo{$REmrnaid}{'array'}}, $REcdsend);
						}
						elsif ($REcdsrange =~ /^(\d+)$/) {
							$REcdsstart=$1;
							$REcdsend=$1;
							push (@{$REorfinfo{$REmrnaid}{'array'}}, $REcdsstart);
							push (@{$REorfinfo{$REmrnaid}{'array'}}, $REcdsend);
						}
						else {
							print Dumper $REarray{$REy}, "\n";
							print STDERR $REsubinfo, "Error: unknown CDS range: $REcdsrange\n";
							return 0;
						}
						unless ($REcdsstart<=$REcdsend) {
							print Dumper $REarray{$REy}, "\n";
							print STDERR $REsubinfo, "Error: CDS end<start\n";
							return 0;
						}
						if ($REgetcdsbyarray==1) {
							${$REcds}{$REmrnaid}{'cds'}{$REcdsstart}{$REcdsend}++;
						}
					}
					if (exists ${$REfeatures}{'score'} and ${$REfeatures}{'score'}=~/^\d+$/) {
						${$REcds}{$REmrnaid}{'score'}=${$REfeatures}{'score'};
					}
				}
				else {
					print STDERR $REsubinfo, 'Warnings: unknown feature: ', $REarray{$REy}[0], "\n";
				}
			}
			%REarray=();
			next;
		}
		if ($REline=~/^ID/) {
			unless ($REtotalIDlines<2) {
				print STDERR $REsubinfo, "Error: invalid ID lines at line ($RElinenum): $REline\n";
				return 0;
			}
			if ($REline=~/^ID\s*([\w+\-\.]+).*\s+(\d+)\s+BP\.$/) {
#ID                   wheat3DLnew.chr3DL_BAC_00000001.final ; ; ; ; ; 1633443 BP.
				$REarray{0}=[$1, $2];
			}
			else {
				print STDERR $REsubinfo,  "Error: invalid ID line ($RElinenum): $REline\n";
				return 0;
			}
			$REtotalIDlines++;
		}
		elsif ($REline=~/^FT/) {
			my $REfeat;
			if ($REline=~/^FT   (\w+)\s+(\S+)/) {
				$REi++;
				push (@{$REarray{$REi}}, $1);
				push (@{$REarray{$REi}}, $2);
			}
			elsif ($REline=~/^FT   [ ]{15,23}(\/\S+)/) {
				push (@{$REarray{$REi}}, $1);
			}
			elsif ($REline=~/^FT   [ ]{15,23}(\S+.*)$/) {
				$REarray{$REi}[-1]=$REarray{$REi}[-1].$1;
			}
			else {
				print STDERR $REsubinfo, "Error: invalid FT line ($RElinenum): $REline\n";
				return 0;
			}
		}
		elsif ($REline=~/^SQ/) {
			# Sequence
			$REseqstart=1;
			if (exists $REarray{0} and defined $REarray{0}[0]) {
				print RETEMPFASTA '>', $REarray{0}[0];
				if (defined $REarray{0}[1]) {
					print RETEMPFASTA '    ', $REarray{0}[1], "BP";
				}
				print RETEMPFASTA "\n";
			}
			else {
				die $REsubinfo, "Error: no ID line\n";
			}
			next;
		}
		elsif ($REseqstart==1 and $REline=~/^\s+/) {
			# Sequence
			$REline=~s/\s+\d+\s*$//;
			$REline=~s/\s+//g;
			print RETEMPFASTA $REline, "\n";
			next;
		}
		elsif ($REline=~/^XX/) {
			$REi++;
		}
		else {
			print STDERR $REsubinfo, "Error: unknown line ($RElinenum): $REline\n";
			return 0;
		}
	}
	close REEMBLFILE;
	close RETEMPFASTA;
	
	my $REthereisfasta=0;
	if (-s $REtemp_fasta) {
		$REdb=Bio::DB::Fasta->new($REtemp_fasta);
		$REthereisfasta=1;
	}
	else {
		print STDERR $REsubinfo, "Warnings: can not find fasta sequences, CDS phase would be ignored\n";
	}
	
	foreach my $RErefseq (sort keys %{$REref2gene}) {
		foreach my $REgeneid (sort keys %{${$REref2gene}{$RErefseq}}) {
			unless (exists ${$REgene2mrna}{$REgeneid} and scalar(keys %{${$REgene2mrna}{$REgeneid}})>0) {
				print STDERR $REsubinfo, "Warnings: GENE got no mRNAs: $REgeneid\n";
				next;
			}
			my $REgenestrand;
			if (exists ${$REgene}{$REgeneid} and exists ${$REgene}{$REgeneid}{'strand'} and ${$REgene}{$REgeneid}{'strand'} =~ /(^\+$)|(^\-$)/) {
				$REgenestrand=${$REgene}{$REgeneid}{'strand'};
			}
			else {
				if (exists $REgenedirection{$REgeneid} and scalar(keys %{$REgenedirection{$REgeneid}}) ==1 and (exists $REgenedirection{$REgeneid}{'+'} or exists $REgenedirection{$REgeneid}{'-'})) {
					if (exists $REgenedirection{$REgeneid}{'+'}) {
						$REgenestrand='+';
						${$REgene}{$REgeneid}{'strand'}='+';
					}
					elsif (exists $REgenedirection{$REgeneid}{'-'}) {
						$REgenestrand='-';
						${$REgene}{$REgeneid}{'strand'}='-';
					}
					if (exists ${$REgene}{$REgeneid}{'reference'}) {
						if (${$REgene}{$REgeneid}{'reference'} ne $RErefseq) {
							print STDERR $REsubinfo, "Error: duplicated GENE : $REgeneid on $RErefseq and ",${$REgene}{$REgeneid}{'reference'}, "\n";
						}
					}
					else {
						${$REgene}{$REgeneid}{'reference'}=$RErefseq;
					}
				}
				else {
					print STDERR $REsubinfo, "Error: no strand, REF $RErefseq GENE $REgeneid\n";
				}
			}
			my $REgeneminstart=0;
			my $REgenemaxend=0;
			foreach my $REmrnaid (sort keys %{${$REgene2mrna}{$REgeneid}}) {
				my $REmrnaminstart=0;
				my $REmrnamaxend=0;
				if (exists ${$REmrnas}{$REmrnaid}) {
					if (exists ${$REmrnas}{$REmrnaid}{'strand'} and ${$REmrnas}{$REmrnaid}{'strand'}  =~ /(^\+$)|(^\-$)/) {
						unless (${$REmrnas}{$REmrnaid}{'strand'} eq $REgenestrand) {
							print STDERR $REsubinfo, "Warnings: correct mRNA strand: $REmrnaid\n";
							${$REmrnas}{$REmrnaid}{'strand'}=$REgenestrand;
						}
					}
					else {
						die $REsubinfo, "Error: mRNA no strand: $REmrnaid\n";
					}
				}
				###EXON
				if (exists ${$REexons}{$REmrnaid}) {
					if (exists ${$REexons}{$REmrnaid}{'strand'} and ${$REexons}{$REmrnaid}{'strand'}  =~ /(^\+$)|(^\-$)/) {
						unless (${$REexons}{$REmrnaid}{'strand'} eq $REgenestrand) {
							print STDERR $REsubinfo, "Warnings: correct EXON strand: $REmrnaid\n";
							${$REexons}{$REmrnaid}{'strand'}=$REgenestrand;
						}
					}
					else {
						die $REsubinfo, "Error: mRNA no strand: $REmrnaid\n";
					}
					if (exists ${$REexons}{$REmrnaid}{'exon'}) {
						my @REexonarr=sort {$a<=>$b} keys %{${$REexons}{$REmrnaid}{'exon'}};
						if (scalar(@REexonarr)>0) {
							my $REgenebordleft=$REexonarr[0];
							my @REends=sort {$a<=>$b} keys %{${$REexons}{$REmrnaid}{'exon'}{$REexonarr[-1]}};
							my $REgenebordright=$REends[-1];
							
							if ($REmrnaminstart==0) {
								$REmrnaminstart=$REgenebordleft;
							}
							elsif ($REmrnaminstart>$REgenebordleft) {
								$REmrnaminstart=$REgenebordleft;
							}
							if ($REmrnamaxend==0) {
								$REmrnamaxend=$REgenebordright;
							}
							elsif ($REmrnamaxend<$REgenebordright) {
								$REmrnamaxend=$REgenebordright;
							}
							if ($REgeneminstart==0) {
								$REgeneminstart=$REgenebordleft;
							}
							elsif ($REgeneminstart>$REgenebordleft) {
								$REgeneminstart=$REgenebordleft;
							}
							if ($REgenemaxend==0) {
								$REgenemaxend=$REgenebordright;
							}
							elsif ($REgenemaxend<$REgenebordright) {
								$REgenemaxend=$REgenebordright;
							}
							if ($REgetcdsbyarray==0) {
								if (exists $REorfinfo{$REmrnaid} and exists $REorfinfo{$REmrnaid}{'array'} and scalar(@{$REorfinfo{$REmrnaid}{'array'}})>0) {
									my @REallcoordinates=sort {$a<=>$b} @{$REorfinfo{$REmrnaid}{'array'}};
									my $REmincoord=$REallcoordinates[0];
									my $REmaxcoord=$REallcoordinates[-1];
#									print $REsubinfo, "Test: MIN $REmincoord MAX $REmaxcoord: mRNA $REmrnaid\n"; ### For test ###
									foreach my $REexoncoord1 (@REexonarr) {
										foreach my $REexoncoord2 (sort {$a<=>$b} keys %{${$REexons}{$REmrnaid}{'exon'}{$REexoncoord1}}) {
#											print $REsubinfo, "Test: CO1 $REexoncoord1 CO2 $REexoncoord2\n"; ### For test ###
											if ($REmincoord>=$REexoncoord1) {
												if ($REmaxcoord<=$REexoncoord2) {
													${$REcds}{$REmrnaid}{'cds'}{$REmincoord}{$REmaxcoord}++;
												}
												elsif ($REmaxcoord>$REexoncoord2 and $REmincoord<=$REexoncoord2) {
													${$REcds}{$REmrnaid}{'cds'}{$REmincoord}{$REexoncoord2}++;
												}
											}
											else {
												if ($REmaxcoord<=$REexoncoord2 and $REmaxcoord>=$REexoncoord1) {
													${$REcds}{$REmrnaid}{'cds'}{$REexoncoord1}{$REmaxcoord}++;
												}
												elsif ($REmaxcoord>$REexoncoord2) {
													${$REcds}{$REmrnaid}{'cds'}{$REexoncoord1}{$REexoncoord2}++;
												}
											}
										}
#										print Dumper ${$REcds}{$REmrnaid}{'cds'}; ### For test ###
									}
								}
								else {
									print STDERR $REsubinfo, "Warnings: unknown CDS min max\n";
								}
							}
						}
						else {
							print STDERR $REsubinfo, "Warnings: mRNA got no EXONs: $REmrnaid\n";
						}
					}
				}
				else {
					print STDERR $REsubinfo, "Warnings: mRNA got no EXONs: $REmrnaid\n";
				}
				
				### CDS
				if (exists ${$REcds}{$REmrnaid}) {
					if (exists ${$REcds}{$REmrnaid}{'strand'}) {
						if (${$REcds}{$REmrnaid}{'strand'}  =~ /(^\+$)|(^\-$)/) {
							unless (${$REcds}{$REmrnaid}{'strand'} eq $REgenestrand) {
								print STDERR $REsubinfo, "Warnings: correct CDS strand: $REmrnaid\n";
								${$REcds}{$REmrnaid}{'strand'}=$REgenestrand;
							}
						}
						else {
							print STDERR $REsubinfo, "Warnings: guess CDS strand: $REmrnaid\n";
							${$REcds}{$REmrnaid}{'strand'}=$REgenestrand;
						}
					}
					else {
						die $REsubinfo, "Error: mRNA no strand: $REmrnaid\n";
					}
					if (exists ${$REcds}{$REmrnaid}{'cds'}) {
						my @REcdsarr=sort {$a<=>$b} keys %{${$REcds}{$REmrnaid}{'cds'}};
						if (scalar(@REcdsarr)>0) {
							my $REgenebordleft=$REcdsarr[0];
							my @REends=sort {$a<=>$b} keys %{${$REcds}{$REmrnaid}{'cds'}{$REcdsarr[-1]}};
							my $REgenebordright=$REends[-1];
							
							if ($REmrnaminstart==0) {
								$REmrnaminstart=$REgenebordleft;
							}
							elsif ($REmrnaminstart>$REgenebordleft) {
								$REmrnaminstart=$REgenebordleft;
							}
							if ($REmrnamaxend==0) {
								$REmrnamaxend=$REgenebordright;
							}
							elsif ($REmrnamaxend < $REgenebordright) {
								$REmrnamaxend=$REgenebordright;
							}
							if ($REgeneminstart==0) {
								$REgeneminstart=$REgenebordleft;
							}
							elsif ($REgeneminstart>$REgenebordleft) {
								$REgeneminstart=$REgenebordleft;
							}
							if ($REgenemaxend==0) {
								$REgenemaxend=$REgenebordright;
							}
							elsif ($REgenemaxend<$REgenebordright) {
								$REgenemaxend=$REgenebordright;
							}
							if ($REthereisfasta==1) {
								my ($REcotest, $REcdsarr)=&CheckOverlap(${$REcds}{$REmrnaid}{'cds'});
								unless ($REcotest) {
									print $REsubinfo, "Test: mRNA $REmrnaid CDS coord: \n";print Dumper ${$REcds}{$REmrnaid}{'cds'}; print "\n";
									die $REsubinfo, "Error: transcript $REmrnaid have cds overlaps\n";
								}
								my $REcdsseq='';
								my $REtestcdsseq=0;
								my $REtest_completephase=0;
								my $REtest_partialphase=0;
								foreach my $REcdsnum (@{$REcdsarr}) {
									if (${$REcdsnum}[0] =~/^\d+$/ and ${$REcdsnum}[1]=~/^\d+$/) {
										$REcdsseq.=$REdb->seq($RErefseq, ${$REcdsnum}[0] => ${$REcdsnum}[1]);
										if (exists $RECDSphases{$REmrnaid}) {
											if (exists $RECDSphases{$REmrnaid}{${$REcdsnum}[0]} and $RECDSphases{$REmrnaid}{${$REcdsnum}[0]}=~/^[0-2]{1}$/) {
												${$REcds}{$REmrnaid}{'phase'}{${$REcdsnum}[0]}{${$REcdsnum}[1]}=$RECDSphases{$REmrnaid}{${$REcdsnum}[0]};
											}
											elsif (exists $RECDSphases{$REmrnaid}{${$REcdsnum}[1]} and $RECDSphases{$REmrnaid}{${$REcdsnum}[1]}=~/^[0-2]{1}$/) {
												${$REcds}{$REmrnaid}{'phase'}{${$REcdsnum}[0]}{${$REcdsnum}[1]}=$RECDSphases{$REmrnaid}{${$REcdsnum}[1]};
											}
											else {
												$REtest_partialphase=1;
											}
										}
										else {
											$REtest_partialphase=1;
										}
									}
									else {
										print STDERR $REsubinfo, "Warning: invalid CDS coordinates\n";
										print $REsubinfo, "Test: \$REcdsnum\n"; print Dumper $REcdsnum; print "\n";
										$REtestcdsseq=1;
										last;
									}
								}
								if ($REtestcdsseq==1 or $REcdsseq eq '') {
									print STDERR $REsubinfo, "Error: failed to get cdsseq\n";
									next;
								}
								elsif ($REtest_partialphase==1) {
									my $REcountend=0;
									my $REtestphase=0;
									my @REphases=();
									my $REprot={};
									if (${$REcds}{$REmrnaid}{'strand'} eq '+') {
										$REprot=Frame3Translation($REcdsseq);
										my @REbestframe=();
										foreach (my $REframe=0; $REframe<3; $REframe++) {
											if (exists ${$REprot}{$REframe} and (${$REprot}{$REframe} =~/\S+/)) {
												unless (${$REprot}{$REframe} =~ /\S+\*\S+/) {
													push (@REbestframe, $REframe);
												}
												else {
													if ($debug) {
														print $REsubinfo, "Test: Frame $REframe internal stop: "; print ${$REprot}{$REframe}, "\n";
													}
												}
											}
										}
										if (scalar(@REbestframe)>1) {
											my @REtempframe=();
											foreach my $REframe (@REbestframe) {
												if (${$REprot}{$REframe}=~/\S+\*$/) {
													push (@REtempframe, $REframe);
												}
												elsif (${$REprot}{$REframe}=~/^M/) {
													push (@REtempframe, $REframe);
												}
											}
											@REbestframe=@REtempframe;
											@REtempframe=();
										}
										if (scalar(@REbestframe)==1) {
											$REtestphase=1;
											my $z=shift @REbestframe;
											for(my $REx=0; $REx<scalar(@{$REcdsarr}); $REx++) {
												push (@REphases, $z);
												my $REy=${$REcdsarr}[$REx];
												$z=3 - ($REy->[1] - $REy->[0] +1 - $z)%3;
												$z=0 if ($z==3);
											}
										}
										else {
											print STDERR $REsubinfo, "Warnings: can not guess CDS phase ",join("/", @REbestframe) ,": mRNA $REmrnaid\n";
											print $REsubinfo, "Test: mRNA $REmrnaid STRAND ",${$REcds}{$REmrnaid}{'strand'} , " Protein seq: \n"; print Dumper $REprot; print "\n"; 
											print $REsubinfo, "Test: mRNA $REmrnaid CDS: \n"; print Dumper ${$REcds}{$REmrnaid}{'cds'}; print "\n";
										}
#										OLDBLOCK1: {### detect start/stop condon
#										if ($REcdsseq =~/(tag$)|(taa$)|(tga$)/i) {
#											$REcountend=3;
#											my $z=0;
#											for(my $REx=scalar(@{$REcdsarr})-1; $REx>=0; $REx--) {
#												my $REy=${$REcdsarr}[$REx];
#												$z=($REy->[1] - $REy->[0] +1 + $z)%3;
#												$z=0 if ($z==3);
#												unshift (@REphases, $z);
#											}
#										}
#										elsif ($REcdsseq =~/^atg/i) {
#											$REcountend=5;
#											my $z=0;
#											for(my $REx=0; $REx<scalar(@{$REcdsarr}); $REx++) {
#												push (@REphases, $z);
#												my $REy=${$REcdsarr}[$REx];
#												$z=3 - ($REy->[1] - $REy->[0] +1 + $z)%3;
#												$z=0 if ($z==3);
#											}
#										}}#OLDBLOCK1
									}
									elsif (${$REcds}{$REmrnaid}{'strand'} eq '-') {
										$REcdsseq=SeqRevComp($REcdsseq);
										$REprot=Frame3Translation($REcdsseq);
										my @REbestframe=();
										foreach (my $REframe=0; $REframe<3; $REframe++) {
											if (exists ${$REprot}{$REframe} and (${$REprot}{$REframe} =~/\S+/)) {
												unless (${$REprot}{$REframe} =~ /\S+\*\S+/) {
													push (@REbestframe, $REframe);
												}
												else {
													if ($debug) {
														print $REsubinfo, "Test: Frame $REframe internal stop: "; print ${$REprot}{$REframe}, "\n";
													}
												}
											}
										}
										if (scalar(@REbestframe)>1) {
											my @REtempframe=();
											foreach my $REframe (@REbestframe) {
												if (${$REprot}{$REframe}=~/\S+\*$/) {
													push (@REtempframe, $REframe);
												}
												elsif (${$REprot}{$REframe}=~/^M/) {
													push (@REtempframe, $REframe);
												}
											}
											@REbestframe=@REtempframe;
											@REtempframe=();
#											print "Test: \@REbestframe: @REbestframe\n";
										}
										if (scalar(@REbestframe)==1) {
											$REtestphase=1;
											my $z=shift @REbestframe;
											for(my $REx=0; $REx<scalar(@{$REcdsarr}); $REx++) {
												unshift (@REphases, $z);
												my $REy=${$REcdsarr}[$REx];
												$z=3 - ($REy->[1] - $REy->[0] +1 - $z)%3;
												$z=0 if ($z==3);
											}
										}
										else {
											print STDERR $REsubinfo, "Warnings: can not guess CDS phase ",join("/", @REbestframe) ,": mRNA $REmrnaid\n";
											print $REsubinfo, "Test: mRNA $REmrnaid STRAND ",${$REcds}{$REmrnaid}{'strand'} , " Protein seq: \n"; print Dumper $REprot; print "\n"; 
										}
#										OLDBLOCK2: {### detect start/stop condon
#										if ($REcdsseq =~/(^cta)|(^tta)|(^tca$)/i) {
#											$REcountend=3;
#											my $z=0;
#											for(my $REx=0; $REx<scalar(@{$REcdsarr}); $REx++) {
#												my $REy=${$REcdsarr}[$REx];
#												$z=($REy->[1] - $REy->[0] +1 + $z)%3;
#												$z=0 if ($z==3);
#												push (@REphases, $z);
#											}
#										}
#										elsif ($REcdsseq =~/cat$/i) {
#											$REcountend=5;
#											my $z=0;
#											for(my $REx=scalar(@{$REcdsarr})-1; $REx>=0; $REx--) {
#												unshift (@REphases, $z);
#												my $REy=${$REcdsarr}[$REx];
#												$z=3 - ($REy->[1] - $REy->[0] +1 + $z)%3;
#												$z=0 if ($z==3);
#											}
#										}}#OLDBLOCK2
									}
									else {
										print STDERR $REsubinfo, "Warnings: unknown CDS strand: $REmrnaid\n";
										next;
									}
#									if ($REcountend==0) {
#										print STDERR $REsubinfo, "Warnings: CDS phase unknown: $REmrnaid\n";
#										next;
#									}
									if ($REtestphase==1) {
										if (scalar(@REphases) eq scalar(@{$REcdsarr})) {
											for(my $REx=0; $REx<scalar(@{$REcdsarr}); $REx++) {
												${$REcds}{$REmrnaid}{'phase'}{${$REcdsarr}[$REx][0]}{${$REcdsarr}[$REx][1]}=$REphases[$REx];
											}
										}
										else {
											print $REsubinfo, "Test: \@REphases\n"; print Dumper \@REphases;print "\n";
											print $REsubinfo, "Test: \$REcdsarr\n"; print Dumper $REcdsarr;print "\n";
											print STDERR $REsubinfo, "Warnings: unequal phase: $REmrnaid\n";
											next;
										}
									}
								}
								else {
									###
								}
							}
						}
						else {
							print STDERR $REsubinfo, "Warnings: mRNA got no EXONs: $REmrnaid\n";
						}
					}
				}
				else {
					print STDERR $REsubinfo, "Warnings: mRNA got no CDS: $REmrnaid\n";
					next;
				}
				unless (exists ${$REmrnas}{$REmrnaid}{'start'} and 
					${$REmrnas}{$REmrnaid}{'start'} =~/^\d+$/ and 
					$REmrnaminstart !=0 and 
					${$REmrnas}{$REmrnaid}{'start'} == $REmrnaminstart
				) {### correcting mRNA start
					if ($debug) {
						print STDERR $REsubinfo, "Warnings: Correcting mRNA start : mRNA ",$REmrnaid , " ";
						if (exists ${$REmrnas}{$REmrnaid} and exists ${$REmrnas}{$REmrnaid}{'start'}) {
							print STDERR ${$REmrnas}{$REmrnaid}{'start'}, " => ";
						}
						print STDERR $REmrnaminstart, "\n";
					}
					${$REmrnas}{$REmrnaid}{'start'}=$REmrnaminstart;
				}
				unless (exists ${$REmrnas}{$REmrnaid}{'end'} and 
					${$REmrnas}{$REmrnaid}{'end'} =~/^\d+$/ and 
					$REmrnamaxend !=0 and
					${$REmrnas}{$REmrnaid}{'end'} == $REmrnamaxend
				) {### correcting mRNA ends
					if ($debug) {
						print STDERR $REsubinfo, "Warnings: Correcting mRNA end : mRNA ",$REmrnaid , " ";
						if (exists ${$REmrnas}{$REmrnaid} and exists ${$REmrnas}{$REmrnaid}{'end'}) {
							print STDERR ${$REmrnas}{$REmrnaid}{'end'}, " => ";
						}
						print STDERR $REmrnamaxend, "\n";
					}
					${$REmrnas}{$REmrnaid}{'end'}=$REmrnamaxend;
				}
			}
			$$REref2posgene{$RErefseq}{$REgeneminstart}{$REgeneid}++;
			unless (exists ${$REgene}{$REgeneid}{'start'} and 
				${$REgene}{$REgeneid}{'start'} =~/^\d+$/ and 
				$REgeneminstart !=0 and
				${$REgene}{$REgeneid}{'start'} == $REgeneminstart
			) {
				if ($debug) {
					print STDERR $REsubinfo, "Warnings: Correcting gene start : GENE ",$REgeneid , " ";
					if (exists ${$REgene}{$REgeneid} and exists ${$REgene}{$REgeneid}{'start'}) {
						print STDERR ${$REgene}{$REgeneid}{'start'}, " => ";
					}
					print STDERR $REgeneminstart, "\n";
				}
				${$REgene}{$REgeneid}{'start'} = $REgeneminstart;
			}
			unless (exists ${$REgene}{$REgeneid}{'end'} and 
				${$REgene}{$REgeneid}{'end'} =~/^\d+$/ and 
				$REgenemaxend !=0 and
				${$REgene}{$REgeneid}{'end'} == $REgenemaxend
			) {
				if ($debug) {
					print STDERR $REsubinfo, "Warnings: Correcting gene end : GENE ",$REgeneid , " ";
					if (exists  ${$REgene}{$REgeneid} and exists  ${$REgene}{$REgeneid}{'end'}) {
						print STDERR ${$REgene}{$REgeneid}{'end'}, " => ";
					}
					print STDERR $REgenemaxend, "\n";
				}
				${$REgene}{$REgeneid}{'end'} = $REgenemaxend;
			}
		}
	}
	unless ($REseqterminal==1) {
		die $REsubinfo, "Error: No EMBL ends //\n";
	}
	unlink $REtemp_fasta if (-e $REtemp_fasta);
	unlink $REtemp_fasta."index" if (-e $REtemp_fasta."index");
	return (1, $REref2posgene, $REgene2mrna, $REgene, $REmrnas, $REexons, $REcds);
}



sub CheckOverlap {
	my $COborderhash=shift;
	
	my $COsubinfo='SUB(CheckOverlap)';
	
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
				print STDERR $COsubinfo, "Error: feature overlap ($COarr[$COi][0] $COarr[$COi][1]) and ($COarr[$COj][0] $COarr[$COj][1])\n";
				return 0;
			}
		}
	}
	return (1, \@COarr);
}


sub GetGeneHash {
	my $GGHarr1=shift;
	
	my $GGHsubinfo='SUB(GetGeneHash)';
	my $GGHret_hash={};
	
	for (my $GGHi=2; $GGHi<scalar(@{$GGHarr1}); $GGHi++) {
		if (defined ${$GGHarr1}[$GGHi] and ${$GGHarr1}[$GGHi]=~/^\/(\w+)=(.*)$/) {
			my $GGHindex=$1;
			my $GGHvalue=$2;
			$GGHvalue=~s/^"//;$GGHvalue=~s/"$//;
			if (exists ${$GGHret_hash}{$GGHindex}) {
				print Dumper $GGHarr1;
				die $GGHsubinfo, "Error: repeated feature index\n";
			}
			${$GGHret_hash}{$GGHindex}=$GGHvalue;
		}
		else {
			print Dumper $GGHarr1;
			die $GGHsubinfo, "Error: unknown feature\n";
		}
	}
	return $GGHret_hash;
}
