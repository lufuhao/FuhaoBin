#!/usr/bin/env perl
use strict;
use warnings;
use List::Util qw (min max);
use Getopt::Long;
use constant USAGE=><<EOH;

SYNOPSIS:

perl $0 --input my.vcf [Options]
Version: LUFUHAO20150217

Requirements:
	Programs:
	Modules: Getopt::Long, constant

Descriptions:
	filter vcf by specified parameters

Options: -[cdefghinoqrstvx]
	--help|-h
		Print this help/usage;
	--input|-i	<File>
		[Msg] VCF file
	--output|-o	<File>
		[Msg] output vcf name
	--readdepth|-d	<'Int-Int'>
		[Opt] Total read depth, Default: '0-99999'
	--allelefrequency|-f	<'Float-Float'>
		[Opt] Allele (Ref and Var) frequency, Default '0-1'
	--allelecount|-c	<'Int-Int'>
		[Opt] Each alternative allele depth, Default: '0-99999'
	--refallelecount|-e	<'Int-Int'>
		[Opt] Refence allele depth, Default: '0-99999'
	--numallele|-n	<'Int-Int'>
		[Opt] Total number of alleles, Default: '0-99999'
	--meanmapq|-q	<'Int-Int'>
		[Opt] Mean mapping quality (6th Column)
		Default: '0-99999'
	--type|-t	<'Int/Int/Int/Int/Int'>
		[Opt] output variation filter: snp/mnp/insert/deletion/compleat
		Default: '1/1/1/1/1'; Change 1 to 0 for not output
	--numgenocol|-g	<INT>
		[Opt] Number of domains (':' delimiter) in GENO column, Default: 7
	--references|-s	<File>
		[Opt] Reference ID file, per ID perl line
	--gtfgff|-x <File>
		[Opt] Feature file; use Col1(ChrID), col3(Feature), Col3(Start) and Col4(End)
	--feature|-y	<feature list>
		[Opt] Feature list in Col3 of gtf/gff, delimited by comma
	--reverse|-r
		[Opt] Reverse selection of Reference ID above
	--verbose
		[Opt] Detailed output for trouble-shooting;
	--version|v!
		[Opt] Print current SCRIPT version;

	*float might be 5 digits after dot

Example:
	vcf-filter.pl -i 3DL.st3.calmd.vcf -o 3DL.st3.calmd.mnp.ref3 -d '3' -e '3' -t '1/1/0/0/0'

Author:
	Fu-Hao Lu
	Post-Doctoral Scientist in Micheal Bevan laboratory
	Cell and Developmental Department, John Innes Centre
	Norwich NR4 7UH, United Kingdom
	E-mail: Fu-Hao.Lu\@jic.ac.uk

EOH
###HELP ends#########################################################
die USAGE unless @ARGV;



###Receving parameter################################################
our ($help, $verbose, $ver, $debug);
our ($input, $output);
our ($read_depth, $allele_frequency, $allele_count, $ref_allele_count, $num_allele, $mean_mapq, $var_type, $num_geno_col);
our ($refid_list, $reverse_selection);
our ($gtfgff, $feature);
GetOptions(
	"help|h!" => \$help,
	"input|i:s" => \$input,
	"output|o:s" => \$output,
	"readdepth|d:s" => \$read_depth,
	"allelefrequency|f:s" => \$allele_frequency,
	"allelecount|c:s" => \$allele_count,
	"refallelecount|e:s" => \$ref_allele_count,
	"numallele|n:s" => \$num_allele,
	"meanmapq|q:s" => \$mean_mapq,
	"type|t:s" => \$var_type,
	"numgenocol|g:i" => \$num_geno_col,
	"references|s:s" => \$refid_list,
	"gtfgff|x:s" => \$gtfgff,
	"feature|y:s"	=> \$feature,
	"reverse|r!" => \$reverse_selection,
	"debug!" => \$debug,
	"verbose!" => \$verbose,
	"version|v!" => \$ver) or die USAGE;
($help or $ver) and die USAGE;



### Defaults ########################################################
$num_geno_col=6 unless (defined $num_geno_col);
$read_depth='0-99999' unless (defined $read_depth);
$allele_frequency='0-1' unless (defined $allele_frequency);
$allele_count='0-99999' unless (defined $allele_count);
$num_allele='0-99999' unless (defined $num_allele);
$ref_allele_count='0-99999' unless (defined $ref_allele_count);
$mean_mapq='0-99999' unless (defined $mean_mapq);
$var_type='1/1/1/1/1' unless (defined $var_type);

$reverse_selection=0 unless (defined $reverse_selection);
$verbose=0 unless (defined $verbose);

our @gtfgfffeatures=();
our %gtffeatures=();
if (defined $feature and $feature ne '') {
#	print "Input feature: ".$feature."\n";###for test###
	@gtfgfffeatures=split(/,/, $feature);
	if (scalar(@gtfgfffeatures) >=1) {
		die "Error: --feature|-y needs to specify --gtfgff|-x option\n" unless (defined $gtfgff and -s $gtfgff);
		map {$gtffeatures{$_}++} @gtfgfffeatures;
	}
	else {
		die "Error: invalid --feature|-y option\n";
	}
	print "Total number of ".scalar(@gtfgfffeatures)." detected: \n";
	map {print $_."\n"} @gtfgfffeatures;
}
### input and output ################################################
die "Error:VCF input not found\n" unless (defined $input and -s $input);
unlink $output if (defined $output and -e $output);

### Main ############################################################
&FilterVcf($input, $read_depth, $allele_frequency, $allele_count, $ref_allele_count, $num_allele, $var_type, $mean_mapq, $output);



#####################################################################
###                         sub functions                         ###
#####################################################################
###FilterVCF(vcf_file, readdepth)
###$FVvar_type=(snp/mnp/insertion/deletion/complex)
###Global: $refid_list, $reverse_selection, $gtfgff, %gtffeatures
sub FilterVcf {
	my ($FVfile_vcf, $FVread_depth, $FVallele_frequency, $FVallele_count, $FVref_allele_count, $FVnum_allele, $FVvar_type, $FVmean_mapq, $FVoutput)=@_;
###COMMIT: setup total read depth
	my ($FVmin_read_depth, $FVmax_read_depth)=(0, 99999);
	if (defined $FVread_depth) {
		if ($FVread_depth=~m/-/) {
			my @FVarr1=split(/-/,$FVread_depth);
			$FVmin_read_depth=$FVarr1[0] if (defined $FVarr1[0] and $FVarr1[0] ne '' and $FVarr1[0]>=0);
			$FVmax_read_depth=$FVarr1[1] if (defined $FVarr1[1] and $FVarr1[1] ne '' and $FVarr1[1]>=0);
		}
		else {
			$FVmin_read_depth=$FVread_depth;
		}
		if ($FVmin_read_depth > $FVmax_read_depth or $FVmin_read_depth !~ m/^\d+$/ or $FVmax_read_depth !~m/^\d+$/) {
			die "SUB(FilterVCF)Error: wrong read depth setting\n";
		}
	}
###COMMIT: setup Allele frequency
	my ($FVmin_allele_frequency, $FVmax_allele_frequency)=(0,1);
	if (defined $FVallele_frequency) {
		if ($FVallele_frequency=~m/-/) {
			my @FVarr2=split(/-/,$FVallele_frequency);
			$FVmin_allele_frequency=$FVarr2[0] if (defined $FVarr2[0] and $FVarr2[0] ne '' and $FVarr2[0]>=0);
			$FVmax_allele_frequency=$FVarr2[1] if (defined $FVarr2[1] and $FVarr2[1] ne '' and $FVarr2[1]>=0);
		}
		else {
			$FVmin_allele_frequency=$FVallele_frequency;
		}
		if ($FVmin_allele_frequency > $FVmax_allele_frequency or $FVmin_allele_frequency !~ m/^\d+$/ or $FVmax_allele_frequency !~m/^\d+$/) {
			die "SUB(FilterVCF)Error: wrong Allele frequency setting\n";
		}
	}
###COMMIT: setup alternative allele depth
	my ($FVmin_allele_count, $FVmax_allele_count)=(0,99999);
	if (defined $FVallele_count) {
		if ($FVallele_count=~m/-/) {
			my @FVarr4=split(/-/,$FVallele_count);
			$FVmin_allele_count=$FVarr4[0] if (defined $FVarr4[0] and $FVarr4[0] ne '' and $FVarr4[0]>=0);
			$FVmax_allele_count=$FVarr4[1] if (defined $FVarr4[1] and $FVarr4[1] ne '' and $FVarr4[1]>=0);
		}
		else {
			$FVmin_allele_count=$FVallele_count;
		}
		if ($FVmin_allele_count > $FVmax_allele_count or $FVmin_allele_count !~ m/^\d+$/ or $FVmax_allele_count !~m/^\d+$/) {
			die "SUB(FilterVCF)Error: wrong allele count setting\n";
		}
	}
###COMMIT: setup total number of alleles
	my ($FVmin_num_allele, $FVmax_num_allele)=(0,99999);
	if (defined $FVnum_allele) {
		if ($FVnum_allele=~m/-/) {
			my @FVarr3=split(/-/,$FVnum_allele);
#			print "@FVarr3\n";##For test###
			$FVmin_num_allele=$FVarr3[0] if (defined $FVarr3[0] and $FVarr3[0] ne '' and $FVarr3[0]>=0);
			$FVmax_num_allele=$FVarr3[1] if (defined $FVarr3[1] and $FVarr3[1] ne '' and $FVarr3[1]>=0);
		}
		else {
			$FVmin_num_allele=$FVnum_allele;
		}
		if ($FVmin_num_allele > $FVmax_num_allele or $FVmin_num_allele !~ m/^\d+$/ or $FVmax_num_allele !~m/^\d+$/) {
			die "SUB(FilterVCF)Error: wrong number of alleles setting\n";
		}
	}
#	print $FVnum_allele, "\t", $FVmin_num_allele,"\t", $FVmax_num_allele."\n"; ###For test###
###COMMIT: setup mean mapping quality
	my ($FVmin_mean_mapq, $FVmax_mean_mapq)=(0,99999);
	if (defined $FVmean_mapq) {
		if ($FVmean_mapq=~m/-/) {
			my @FVarr5=split(/-/,$FVmean_mapq);
			$FVmin_mean_mapq=$FVarr5[0] if (defined $FVarr5[0] and $FVarr5[0] ne '' and $FVarr5[0]>=0);
			$FVmax_mean_mapq=$FVarr5[1] if (defined $FVarr5[1] and $FVarr5[1] ne '' and $FVarr5[1]>=0);
		}
		else {
			$FVmin_mean_mapq=$FVmean_mapq;
		}
		if ($FVmin_mean_mapq > $FVmax_mean_mapq or $FVmin_mean_mapq !~ m/^\d+$/ or $FVmax_mean_mapq !~m/^\d+$/) {
			die "SUB(FilterVCF)Error: wrong number of alleles setting\n";
		}
	}
###COMMIT: setup reference allele depth
	my ($FVmin_ref_allele_count, $FVmax_ref_allele_count)=(0,99999);
	if (defined $FVref_allele_count) {
		if ($FVref_allele_count=~m/-/) {
			my @FVarr6=split(/-/,$FVref_allele_count);
			$FVmin_ref_allele_count=$FVarr6[0] if (defined $FVarr6[0] and $FVarr6[0] ne '' and $FVarr6[0]>=0);
			$FVmax_ref_allele_count=$FVarr6[1] if (defined $FVarr6[1] and $FVarr6[1] ne '' and $FVarr6[1]>=0);
		}
		else {
			$FVmin_ref_allele_count=$FVref_allele_count;
		}
		if ($FVmin_ref_allele_count > $FVmax_ref_allele_count or $FVmin_ref_allele_count !~ m/^\d+$/ or $FVmax_ref_allele_count !~m/^\d+$/) {
			die "SUB(FilterVCF)Error: wrong ref_allele_count setting\n";
		}
	}
###COMMIT: set output variation type
	my ($FVout_snp, $FVout_mnp, $FVout_insert, $FVout_deletion, $FVout_complex)=(1, 1, 1, 1, 1);
	if (defined $FVvar_type) {
		my @FVarr7=split(/\//, $FVvar_type);
		$FVout_snp=$FVarr7[0] if (defined $FVarr7[0] and $FVarr7[0] ne '' and $FVarr7[0] >=0);
		$FVout_mnp=$FVarr7[1] if (defined $FVarr7[1] and $FVarr7[1] ne '' and $FVarr7[1] >=0);
		$FVout_insert=$FVarr7[2] if (defined $FVarr7[2] and $FVarr7[2] ne '' and $FVarr7[2] >=0);
		$FVout_deletion=$FVarr7[3] if (defined $FVarr7[3] and $FVarr7[3] ne '' and $FVarr7[3] >=0);
		$FVout_complex=$FVarr7[4] if (defined $FVarr7[4] and $FVarr7[4] ne '' and $FVarr7[4] >=0);
	}

### Read ID list
	my %FVid_list=();
	my $FVfilterVCF_true=0;
	if (defined $refid_list and -s "$refid_list") {
		open (IDLIST, "<$refid_list") || die "SUB(FilterVCF)Error: ID list file open\n";
		while (my $FVline2=<IDLIST>) {
			chomp $FVline2;
#			print $FVline2."\n";###for Test###
			(my $FVseqid=$FVline2)=~m/^(\S+)\s*.*$/;
#			print $FVseqid."\n";###for Test###
			$FVid_list{"$FVseqid"}++ if (defined $FVseqid and $FVseqid ne '');
		}
		close IDLIST;
		print "SUB(FilterVCF)Info: Number of filter IDs: ". scalar(keys %FVid_list)."\n";
		$FVfilterVCF_true=1 if (scalar(keys %FVid_list)>0);
	}
	if ($FVfilterVCF_true and $debug) {
		print "SUB(FilterVCF)Info: filter IDs:\n";
		map {print $_."\n"} (keys %FVid_list);
	}
	
### gtfgff feature
	my %FVfeatures=();
	my $filter_featre=0;
	if ( defined $gtfgff and -s $gtfgff) {
		$filter_featre=1;
		open (GTFGFF, "<$gtfgff") || die "SUB(FilterVCF)Error: open gtfgff file\n";
		while (my $FVline3=<GTFGFF>) {
			chomp $FVline3;
			next if ($FVline3 =~/(^#)|(^$)/);
			my @FVarr8=split(/\t/, $FVline3);
			next if (scalar(@FVarr8)<5);
#			print $FVarr8[0]."\n";
			next unless (defined $FVarr8[0] and $FVarr8[0] ne '');
			if (scalar (keys %gtffeatures) >=1) {
				next unless (exists $gtffeatures{$FVarr8[2]});
			}
			my $FVfeat_min=min($FVarr8[3], $FVarr8[4]);
			my $FVfeat_max=max($FVarr8[3], $FVarr8[4]);
			${${$FVfeatures{$FVarr8[0]}}{$FVfeat_min}}{$FVfeat_max}++;
		}
		close GTFGFF;
		print "SUB(FilterVCF)Info: number of chrID by feature: ".scalar(keys %FVfeatures)."\n";
		die "SUB(FilterVCF)Error: number of chrID by feature: ".scalar(keys %FVfeatures)."\n" if (scalar(keys %FVfeatures)<1);
	}
###Summary of filters
	print STDERR "\n\n\n##### Filter summary #####\n";
	print STDERR "VCFfile:          $FVfile_vcf\n";
	print STDERR "TotalReadDepth:   $FVmin_read_depth ~ $FVmax_read_depth\n";
	print STDERR "AlleleFreuency:   $FVmin_allele_frequency ~ $FVmax_allele_frequency\n";
	print STDERR "RefAlleleDepth:   $FVmin_ref_allele_count ~ $FVmax_ref_allele_count\n";
	print STDERR "AltAlleleDepth:   $FVmin_allele_count ~ $FVmax_allele_count\n";
	print STDERR "NumAllele:        $FVmin_num_allele ~ $FVmax_num_allele\n";
	print STDERR "VariationType:    SNP/MNP/Ins/Del/Oth = $FVout_snp/$FVout_mnp/$FVout_insert/$FVout_deletion/$FVout_complex\n";
	print STDERR "MeanMAPQ:         $FVmin_mean_mapq ~ $FVmax_mean_mapq\n";
	print STDERR "\n\n\n";
### Parse VCF
	if ($FVfile_vcf=~/\.vcf$/i) {
		open (FVVCF, "cat $FVfile_vcf |") || die "SUB(FilterVCF)Error: vcf file open\n";
	}
	elsif ($FVfile_vcf=~/\.gz$/i) {
		open (FVVCF, "zcat $FVfile_vcf |") || die "SUB(FilterVCF)Error: vcf file open\n";
	}
	else {
		die "Can not guess VCF input format\n";
	}
	open (FVOUT, ">$FVoutput") || die "SUB(FilterVCF)Error: write vcf output\n";
	my $FVvcf_line_num=0;
	while (my $FVline=<FVVCF>) {
		$FVvcf_line_num++;
		chomp $FVline;
		if ($FVline=~/^#/) {#ignore comments
			print FVOUT $FVline."\n";
			next;
		}
		my @FVarr=split(/\t/, $FVline);
		if (scalar(@FVarr)<10) {#Ignore lines with column number less that 10
			print STDERR "SUB(FilterVCF)Error: $FVarr[0]\t$FVarr[1]\tLess10col\n";
			next;
		}
#Filter ID list
		if ($FVfilterVCF_true and defined $FVarr[0] and $FVarr[0] ne '') {##filter ids from id list
			if (exists $FVid_list{$FVarr[0]}) {
				if ($reverse_selection) {
					print STDERR "SUB(FilterVCF)Info: $FVarr[0]\t$FVarr[1]\tReverseSelection\n" if ($verbose);
					next;
				}
			}
			else {
				unless ($reverse_selection) {
					print STDERR "SUB(FilterVCF)Info: $FVarr[0]\t$FVarr[1]\tNotListed\n" if ($verbose);
					next;
				}
			}
		}
#Filter GTF.GFF feature
		if ($filter_featre) {
			my $FVkeep_this_record=0;
			next unless (exists $FVfeatures{$FVarr[0]});
			BLOCK3: {foreach my $FVstart (sort {$a<=> $b} keys %{$FVfeatures{$FVarr[0]}}) {
				if ($FVarr[1] >= $FVstart) {
					foreach my $FVend (sort {$b<=>$a} keys %{${$FVfeatures{$FVarr[0]}}{$FVstart}}) {
						if ($FVarr[1] <= $FVend) {
							$FVkeep_this_record=1;
							last BLOCK3;
						}
					}
				}
			}}
			next if ($FVkeep_this_record ==0);
		}
#Filter total number of alleles
		my $FVref_allele=$FVarr[3];
		my @FVvar_allele=split(/,/,$FVarr[4]);
		my $FVindv_total_num_alleles=scalar(@FVvar_allele)+1;
		if ($FVindv_total_num_alleles < $FVmin_num_allele or $FVindv_total_num_alleles > $FVmax_num_allele) {
			print STDERR "SUB(FilterVCF)Info: $FVarr[0]\t$FVarr[1]\ttotal number of alleles :$FVindv_total_num_alleles, $FVmin_num_allele, $FVmax_num_allele\n";
			next;
		}
#filter mean MAPQ
		if (defined $FVarr[5]) {
			if ($FVarr[5]<$FVmin_mean_mapq or $FVarr[5]>$FVmax_mean_mapq) {
				print STDERR "SUB(FilterVCF)Info: $FVarr[0]\t$FVarr[1]\tmapq\n";
				next;
			}
		}
		else {
			print STDERR "SUB(FilterVCF)Error: no QUAL at line $FVvcf_line_num; Ignored...\n";
			next;
		}
#filter variation type
		my ($FVthis_is_snp, $FVthis_is_mnp, $FVthis_is_ins, $FVthis_is_del, $FVthis_is_complex)=(0, 0, 0, 0, 0);
		foreach my $FVindv_var_alleles (@FVvar_allele) {
			BLOCK1: {
				if ($FVindv_var_alleles eq '-' or length($FVref_allele) > length($FVindv_var_alleles)) { $FVthis_is_del=1; last BLOCK1;}
				if (length($FVref_allele)==1 and length($FVindv_var_alleles)==1) {$FVthis_is_snp=1; last BLOCK1;}
				if (length($FVref_allele)>1 and length($FVref_allele) == length($FVindv_var_alleles)) {$FVthis_is_mnp=1; last BLOCK1;}
				if (length($FVref_allele) < length($FVindv_var_alleles)) {$FVthis_is_ins=1; last BLOCK1;}
			}
		}
		$FVthis_is_complex=$FVthis_is_snp+$FVthis_is_mnp+$FVthis_is_ins+$FVthis_is_del;
		if ($FVthis_is_complex != 1) {
			unless ($FVout_complex) {
				print STDERR "SUB(FilterVCF)Info: $FVarr[0]\t$FVarr[1]\tComplex\n" if ($verbose or $debug);
				next;
			}
		}
		elsif ($FVthis_is_snp){
			unless ($FVout_mnp){
				print STDERR "SUB(FilterVCF)Info: $FVarr[0]\t$FVarr[1]\tsnp\n" if ($verbose or $debug);
				next;
			}
		}
		elsif ($FVthis_is_mnp) {
			unless ($FVout_mnp) {
				print STDERR "SUB(FilterVCF)Info: $FVarr[0]\t$FVarr[1]\tmnp\n" if ($verbose or $debug);
				next;
			}
		}
		elsif ($FVthis_is_ins) {
			unless ($FVout_insert) {
				print STDERR "SUB(FilterVCF)Info: $FVarr[0]\t$FVarr[1]\tinsertion\n" if ($verbose or $debug);
				next;
			}
		}elsif ($FVthis_is_del) {
			unless ($FVout_deletion) {
				print STDERR "SUB(FilterVCF)Info: $FVarr[0]\t$FVarr[1]\tdeletion\n" if ($verbose or $debug);
				next;
			}
		}
#filter read_depth, AlleleFreuency
		my ($FVtest_readdepth, $FVtest_allelefreq, $FVtest_refalleledepth, $FVtest_allele_count)=(0,0,0,0);
		BLOCK2: {
			for (my $FVi=9; $FVi<scalar(@FVarr); $FVi++) {
				my @FVarr2=split(/:/, $FVarr[$FVi]);
				if (scalar(@FVarr2)==$num_geno_col) {
					unless ($FVarr2[1]<$FVmin_read_depth or $FVarr2[1]>$FVmax_read_depth) {
						$FVtest_readdepth=1;
					}
					unless ($FVarr2[2] < $FVmin_ref_allele_count or $FVarr2[2] > $FVmax_ref_allele_count) {
						$FVtest_refalleledepth=1;
					}
					my @FVarr3=split(/,/, $FVarr2[4]);
					unshift (@FVarr3, $FVarr2[2]);
					foreach my $FVao (@FVarr3) {
						unless (($FVao/$FVarr2[1])<$FVmin_allele_frequency or ($FVao/$FVarr2[1])>$FVmax_allele_frequency) {
							$FVtest_allelefreq=1;
						}
						unless ($FVao<$FVmin_allele_count or $FVao > $FVmax_allele_count) {
							$FVtest_allele_count=1;
						}
					}
					
				}
				else {
					print STDERR "SUB(FilterVCF)Error: colnum not equal to $num_geno_col at line $FVvcf_line_num; Ignored...\n";
					last BLOCK2;
				}
			}
		}
		unless ($FVtest_readdepth == 1) {
			print STDERR "SUB(FilterVCF)Info: $FVarr[0]\t$FVarr[1]\tread depth\n" if ($verbose or $debug);
			next;
		}
		unless ($FVtest_allelefreq == 1) {
			print STDERR "SUB(FilterVCF)Info: $FVarr[0]\t$FVarr[1]\tAllele frequency\n" if ($verbose or $debug);
			next;
		}
		unless ($FVtest_refalleledepth ==1) {
			print STDERR "SUB(FilterVCF)Info: $FVarr[0]\t$FVarr[1]\tRefAlleleFrequency\n" if ($verbose or $debug);
			next;
		}
		unless ($FVtest_allele_count ==1) {
			print STDERR "SUB(FilterVCF)Info: $FVarr[0]\t$FVarr[1]\tAlelle depth\n" if ($verbose or $debug);
			next;
		}
		print FVOUT $FVline."\n";
	}
	close FVOUT;
	close FVVCF;
}

#196_3dl 112     .       GTG     TTT,ATG 671.416 .       AB=0.333333,0.555556;ABP=9.52472,3.73412;AC=1,1;AF=0.5,0.5;AN=2;AO=9,15;
#CIGAR=1X1M1X,1X2M;DP=27;DPB=27;DPRA=0,0;EPP=5.18177,10.1038;EPPR=0;GTI=0;LEN=3,1;
#MEANALT=4,4;MQM=40,40;MQMR=0;NS=1;NUMALT=2;ODDS=49.0607;PAIRED=0.666667,0.8;
#PAIREDR=0;PAO=0,0;PQA=0,0;PQR=0;PRO=0;QA=310,553;QR=0;RO=0;RPL=3,8;RPP=5.18177,3.15506;
#RPPR=0;RPR=6,7;RUN=1,1;SAF=6,12;SAP=5.18177,14.7363;SAR=3,3;SRF=0;SRP=0;SRR=0;
#TYPE=complex,snp;technology.Illumina=1,1      GT:DP:RO:QR:AO:QA:GL    1/2:27:0:0:9,15:310,553:-76.5495,-51.3671,-48.6578,-31.2855,0,-26.7701
#	AB=Allele balance at heterozygous sites
#AC=Total number of alternate alleles in called genotypes
#AF=Estimated allele frequency in the range (0,1]
#AN=Total number of alleles in called genotypes
#	DP=Total read depth at the locus
#	GT=Genotype
#	GQ=Genotype Quality, the Phred-scaled marginal (or unconditional) probability of the called genotype
#	GL=Genotype Likelihood, log10-scaled likelihoods of the data given the called genotype for each possible genotype generated from the reference and alternate alleles given the sample ploidy
#	RO=Reference allele observation count
#	QR=Sum of quality of the reference observations
#	AO=Alternate allele observation count
#	QA=Sum of quality of the alternate observations
#MQM=Mean mapping quality of observed alternate alleles
#PAIRED=Proportion of observed alternate alleles which are supported by properly paired read fragments
#PAIREDR=Proportion of observed reference alleles which are supported by properly paired read fragments
