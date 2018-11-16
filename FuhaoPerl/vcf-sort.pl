#!/usr/bin/env perl
# modified from VCFlib (Author: petr.danecek@sanger)
#
use strict;
use warnings;
use Carp;
use constant USAGE =><<EOH;

usage:
  vcf-sort.pl file.vcf > out.vcf
  cat file.vcf | vcf-sort.pl > out.vcf

v20181115

Options:
  -c, --chromosomal-order       Use natural ordering (1,2,10,MT,X) rather 
                                than the default (1,10,2,MT,X). This requires 
                                new version of the unix "sort" command 
                                which supports the --version-sort option
  -p, --parallel <int>          Change the number of sorts run concurrently
                                to <int>
  -t, --temporary-directory     Use a directory other than /tmp as the 
                                temporary directory for sorting
  -h, -?, --help                This help message

Author:
    Fu-Hao Lu
    Post-Doctoral Scientist in Micheal Bevan laboratory
    Cell and Developmental Department, John Innes Centre
    Norwich NR4 7UH, United Kingdom
    E-mail: Fu-Hao.Lu\@jic.ac.uk
EOH

die USAGE if (scalar(@ARGV) !=3 or $ARGV[0] eq '-h' or $ARGV[0] eq '--help');

my $opts = parse_params();
sort_vcf($opts);

exit;

#--------------------------------

sub PrintError {
	my (@msg) = @_;
	if ( scalar @msg ) {
		croak @msg;
	}
}

sub parse_params {
	my $opts = {};
	while (my $arg=shift(@ARGV)) {
		if ( $arg eq '-p' || $arg eq '--parallel-sort' ) { $$opts{parallel_sort}=shift(@ARGV); next;}
		if ( $arg eq '-c' || $arg eq '--chromosomal-order' ) { $$opts{chromosomal_order}=1; next;}
		if ( $arg eq '-?' || $arg eq '-h' || $arg eq '--help' ) { die USAGE;}
		if ( $arg eq '-t' || $arg eq '--temporary-directory' ) { $$opts{temp_dir}=shift(@ARGV); next;}
		if ( -e $arg ) { $$opts{file}=$arg; next;}
		die "Error: Unknown parameter \"$arg\". Run -h for help.\n";
	}
	return $opts;
}

sub sort_vcf {
	my ($opts) = @_;

	my $fh;
	if (exists($$opts{file}) ) {
		if ( $$opts{file}=~/\.gz$/i) {
			open($fh,"gunzip -c $$opts{file} |") or PrintError("$$opts{file}: $!");
		}
		else {
			open($fh,'<',$$opts{file}) or PrintError("$$opts{file}: $!");
		}
	}
	else { $fh = *STDIN; }

	my $sort_opts = check_sort_options($opts);
	my $cmd;

	if ( exists($$opts{temp_dir}) ) {
		$cmd = "sort $sort_opts -T $$opts{temp_dir} -k2,2n";    
	}
	else {
		$cmd = "sort $sort_opts -k2,2n";
	}
	print STDERR "$cmd\n";
	open(my $sort_fh,"| $cmd") or PrintError("$cmd: $!");

	my $unflushed = select(STDOUT);
	$| = 1;
	while (my $line=<$fh>) {
		if ($line=~/^#/) {
			print $line;
			next;
		}
		print $sort_fh $line;
		last;
	}
	select($unflushed);
	while (my $line=<$fh>) {
		print $sort_fh $line;
	}
}

sub check_sort_options {
	my ($opts) = @_;

	my $sort_opts='-k1,1d';
	my $sort_help = join('',`sort --help`);

	if ($$opts{chromosomal_order}) {
		my $has_version_sort = ( $sort_help=~/\s+--version-sort\s+/ ) ? 1 : 0;
		if ($has_version_sort) {
			$sort_opts='-k1,1V';
		}
		else {
			print STDERR "Warnings: Old version of sort command installed, please run without the -c option.\n";
		}
	}
	
	if ($$opts{parallel_sort}) {
		my $has_parallel_sort = ( $sort_help=~/\s+--parallel=/ ) ? 1 : 0;
		if ($has_parallel_sort) {
			$sort_opts .= " --parallel $$opts{parallel_sort}";
		}
		else {
			print STDERR "Old version of sort command installed, please run without the -p option.\n";
		}
	}

	return $sort_opts;
}

