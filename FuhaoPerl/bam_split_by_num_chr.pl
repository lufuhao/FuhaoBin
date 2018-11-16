#!/usr/bin/env perl
use warnings;
use strict;
use constant USAGE =><<EOH;

Usage: perl $0 input num_chr output_prefix
min matchlength=36
Version: 20150409

EOH
die USAGE if (scalar(@ARGV)!=3);


###Default
my $path_samtools='samtools';
my $minmatch=36;



###input
my $bam_input=$ARGV[0];
chomp $bam_input;
die USAGE unless (-s $bam_input);
my $num_chr=$ARGV[1];
chomp $num_chr;
my $bam_outputprefix=$ARGV[2];
chomp $bam_outputprefix;



if (&SplitBam($bam_input, $num_chr, $bam_outputprefix)) {
	print STDERR "Bam splitting failed\n";
	exit 1;
}
else {
	print STDERR "Bam splitting succeeded\n";
	exit 0;
}


###split bams according to number of chromosomes
###&SplitBam(bam_input, num_chr, output_prefix)
###Global: 
###Dependency:
###Note: 
sub SplitBam {
	my ($SBban_input, $SBnum_chr, $SBbam_output_prefix)=@_;
	my $SBcmd='';
	my $SBbam_header=&RetrvNoExt($SBban_input);
	$SBbam_header.='.header';
	$SBcmd="$path_samtools view -H $SBban_input > $SBbam_header";
	if (&exec_cmd_return($SBcmd) or ! -s $SBbam_header) {
		print STDERR "SUB(SplitBam)Error: BAM header extract\n";
		return 1;
	}
	unless (open (BAMINPUT, "$path_samtools view $SBban_input |")) {
		print STDERR "SUB(SplitBam)Error: BAM open\n";
		return 1;
	}
	my $SBsuffix=0;
	my %SBseqids=();
	my %SBidlist=();
	my $SBtest_newout=1;
	while (my $SBline=<BAMINPUT>) {
		chomp $SBline;
		my @SBarr=split(/\t/, $SBline);
		$SBseqids{$SBarr[2]}++;
		if (scalar(keys %SBseqids)>$SBnum_chr) {
			$SBtest_newout=1;
			%SBseqids=();
			$SBseqids{$SBarr[2]}++;
			%SBidlist=();
			$SBsuffix++;
		}
		if ($SBtest_newout) {
			if ($SBsuffix != 0) {
			 	close BAMOUT;
				close SEQLIST;
			}
			my $SBindv_output=$SBbam_output_prefix.'_'.$SBsuffix.'.bam';
			my $SBseqID_output=$SBbam_output_prefix.'_'.$SBsuffix.'.id';
			unlink ($SBindv_output) if (-s $SBindv_output);
			unlink ($SBseqID_output) if (-s $SBseqID_output);
			unless (open (SEQLIST, ">$SBseqID_output")) {
				print STDERR "SUB(SplitBam)Error: ID output $SBindv_output\n";
				return 1;
			}
			unless (open (BAMOUT, "| $path_samtools view -bS - > $SBindv_output")) {
				print STDERR "SUB(SplitBam)Error: BAM output $SBindv_output\n";
				return 1;
			}
			unless (open (HEADER, "<$SBbam_header")) {
				print STDERR "SUB(SplitBam)Error: BAM header open for output $SBindv_output\n";
				return 1;
			}
			while (my $SBline2=<HEADER>) {
				print BAMOUT $SBline2;
			}
			close HEADER;
			$SBtest_newout=0;
		}
		unless ($SBarr[1] & 2048) {
			my $SBcount=0;
			while ($SBarr[5]=~/(\d+)M/g) {
				$SBcount+=$1;
			}
			next if ($SBcount<$minmatch);
		}
		print BAMOUT $SBline."\n";
		unless (exists $SBidlist{$SBarr[2]}) {
			print SEQLIST $SBarr[2]."\n";
			$SBidlist{$SBarr[2]}++;
		}
	}
	close BAMOUT;
	close SEQLIST;
	close BAMINPUT;
	unlink ($SBbam_header) if (-s $SBbam_header);
	return 0;
}



#&mytime()
sub mytime {
	my($sec,$min,$hour,$day,$mon,$year,$wday,$yday,$isdst)=localtime();
	$year += 1900;
	$mon  += 1;
	my $btime = sprintf("%04d%02d%02d %02d:%02d:%02d",$year,$mon,$day,$hour,$min,$sec);
	return $btime;
}
###Process command
#&exec_cmd(cmd)
sub exec_cmd {
	my ($cmd) = @_;
	print "#####\n".&mytime()."CMD: $cmd\n";
	my $start_time = time();
	my $return_code = system($cmd);
	my $end_time = time();
#	if ($return_code == -1) {
#		print “failed to execute: $!\n”;
#	}
#	elsif ($return_code & 127) {
#		printf “child died with signal %d, %s coredump\n”, ($return_code & 127),  ($return_code & 128) ? ‘with’ : ‘without’;
#	}
#	else {
#		printf “child exited with value %d\n”, $return_code >> 8;
#	}
	if ($return_code) {
#		print "Error, cmd: $cmd died with ReturnCode $return_code\n";
		die "SUB(exec_cmd)Error, cmd: $cmd died with ReturnCode $return_code\n";
		return $return_code;
	}
	else {
		print STDERR "Finished command: $cmd\nat ".&mytime()."\nRunning time:(".($end_time - $start_time)." seconds) with Returncode: $return_code\n";
		return $return_code;
	}
}
sub exec_cmd_return {
	my ($cmd) = @_;
	print "#####\n".&mytime()."CMD: $cmd\n";
	my $start_time = time();
	my $return_code = system($cmd);
	my $end_time = time();
#	if ($return_code == -1) {
#		print “failed to execute: $!\n”;
#	}
#	elsif ($return_code & 127) {
#		printf “child died with signal %d, %s coredump\n”, ($return_code & 127),  ($return_code & 128) ? ‘with’ : ‘without’;
#	}
#	else {
#		printf “child exited with value %d\n”, $return_code >> 8;
#	}
	if ($return_code) {
		print STDERR "SUB(exec_cmd_return)Error, cmd: $cmd died with ReturnCode $return_code\n";
#		die "Error, cmd: $cmd died with ReturnCode $return_code\n";
		return $return_code;
	}
	else {
		print STDERR "Finished command: $cmd\nat ".&mytime()."\nRunning time:(".($end_time - $start_time)." seconds) with Returncode: $return_code\n";
		return $return_code;
	}
}



###Get filename without extension
###&RetrvNoExt(file)
###Golobal: none
###Dependancy: none
sub RetrvNoExt {
	my $RNEfile=shift @_;
	chomp $RNEfile;
	my ($RNEreturn, $RNEfilename)=('', '');
	($RNEfilename=$RNEfile)=~s/.*\///s;
	($RNEreturn=$RNEfilename)=~s/^(\S+)\.\w+$/$1/;
	die "SUB(RetrvNoExt)Error: empty file name\n" if (!defined $RNEreturn or ($RNEreturn eq ''));
	return $RNEreturn;
}


=history
v20150409
	*Add min match control for each alignment except FLAG & 2048
v20150126
	*
=cut
