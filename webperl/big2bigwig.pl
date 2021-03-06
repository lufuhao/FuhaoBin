#!/usr/bin/perl -w
use Getopt::Std;
use strict;
use File::Basename;

# 1.0 changes .. 
# - results now in an output directory
# - remove dependency on samtools depth to calculate per-base depths due to hard-coded limit of 8000 by samtools depth

# 1.1 changes ..
# add a switch, -m,  to calculate normalized reads to reads per million before calculating depths

# 1.2 .. fixed option -t so that it defaults to 0, and no longer requires integers.

my $version = "1.2"; 
my $usage = "bam2wig
Version $version
USAGE:
bam2wig [options] [.bam]
DEPENDENCIES:
samtools [in PATH]
wigToBigWig [in PATH -- optional]
OPTIONS:
-t [float] : Threshold for reporting. Defaults to 0 (all depths reported)
-s [top/bottom/both] : Strand to report. Defaults to both. 'bottom' will quantify in negative numbers. Used in conjunction with option -g will subset the loci to those on the indicated strand.
-r [string] : Limit analysis to specified read group
-l [integer] : Minimum size of read to report. Default: no minimum.
-L [integer] : Maximum size of read to report. Default: no maximum.
-c [chr:start-stop] : Limit analysis to the indicated locus. 
-g [string] : Path to a GFF3 file. Loci in the GFF3 will be reported. Negates any setting of -c.
-m : Normalize reads to reads per million before calculating depths
-d : Treat as degradome\/PARE data .. tally only the depths of 5' ends. Option -s must be either \'bottom\' or \'top\'
-D : Name of output directory. Defaults to [/YourBamsPath/YourBamsBasename_bam2wig/].
-h : Get help \(print this message\)
-v : Print version number and quit
DOCUMENTATION:
type \'perldoc bam2wig\'
";

# if nothing in ARGV, quit with help message
unless($ARGV[0]) {
    die "$usage";
}

########### Options

our($opt_t,$opt_s,$opt_r,$opt_l,$opt_L,$opt_c,$opt_g,$opt_h,$opt_v,$opt_d,$opt_D,$opt_m);
$opt_t = 0; ## default
$opt_s = "both"; ## default

getopts('mdhvt:s:r:l:L:c:g:D:');

if($opt_h) {
    die "$usage";
}
if($opt_v) {
    die "bam2wig $version\n";
}

############ file checks
unless(-r $ARGV[-1]) {
    die "$usage";
}

### Validate options

unless($opt_t >= 0) {
    die "FATAL: Invalid value for option -t. It must be a number of 0 or more\n";
}

unless(($opt_s eq "both") or
       ($opt_s eq "top") or
       ($opt_s eq "bottom")) {
    die "FATAL: Option -s must be \'top\', \'bottom\', or \'both\'\n";
}

if($opt_l) {
    if($opt_l =~ /^\d+/) {
	if($opt_l >= 1) {
	} else {
	    die "FATAL: Option -l must be an integer >= 1\n";
	}
    } else {
	die "FATAL: Option -l must be an integer >= 1\n";
    }
}

if($opt_L) {
    if($opt_L =~ /^\d+/) {
	if($opt_L >= 1) {
	} else {
	    die "FATAL: Option -L must be an integer >= 1\n";
	}
    } else {
	die "FATAL: Option -L must be an integer >= 1\n";
    }
}

if(($opt_l) and ($opt_L)) {
    unless($opt_l <= $opt_L) {
	die "FATAL: Option l must be less than or equal to option L\n";
    }
}

if ($opt_c) {
    unless($opt_c =~ /^(\S+):(\d+)-(\d+)$/) {
	die "FATAL: Invalid format for option -c. Correct format is Chr:start-stop.\n";
    }
}

if($opt_g) {
    unless(-r $opt_g) {
	die "FATAL: Could not read GFF3 file given in option -g\n";
    }
    if($opt_c) {
	print STDERR "\nWARNING: use of option -g is over-riding use of option -c\n";
	$opt_c = '';
    }
}

if($opt_d) {
    if($opt_s eq "both") {
	die "FATAL: When degradome mode is activated by option -d, option -s must be either top or bottom\n";
    }
}

# do the output directory last
my $bamfilepath = pop @ARGV;
my ($bamfilename,$bamfilejustpath,$bamfilesuffix) = fileparse($bamfilepath,qr/\.[^.]*/);
if($opt_D) {
    # do nothing
} else {
    $opt_D = "$bamfilejustpath" . "$bamfilename" . "_bam2wig";
}
# create directory if needed
unless(-d $opt_D) {
    system "mkdir $opt_D";
}

########### Report on run initiation

print STDERR "\nbam2wig version $version\n";
print STDERR `date`;
print STDERR "Host: ";
print STDERR `hostname`;
print STDERR "Working Directory: ";
print STDERR `pwd`;
### Check dependencies
print STDERR "Checking dependencies\n";
#######
#######
# Check on samtools

open(SAMCHECK, "which samtools |");
my $sam_check = <SAMCHECK>;
close SAMCHECK;
if($sam_check) {
    print STDERR "\tsamtools: present at $sam_check";
} else {
    die "\tsamtools: NOT FOUND\! FATAL -- aborting\n";
}
# Check on wigToBigWig 

open(WIGTOBIGWIG, "which wigToBigWig |");
my $w2bw_check = <WIGTOBIGWIG>;
close WIGTOBIGWIG;
if($w2bw_check) {
    print STDERR "\twigToBigWig: present at $w2bw_check";
} else {
    print STDERR "\twigToBigWig: not found .. a bigwig file will NOT be created\n";
}



print STDERR "BAM file: $bamfilepath\n";
print STDERR "Output Directory: $opt_D\n";
my ($chrs,$chr_lens,$read_groups) = checkbam($bamfilepath);  ### references to arrays are returned from this sub-routine
print STDERR "Threshold: $opt_t read depth
Strand: $opt_s\n";
print STDERR "Read Group: ";
if($opt_r) {
    print STDERR "$opt_r\n";
} else {
    print STDERR "None -- all reads used\n";
}
print STDERR "Min read size: ";
if($opt_l) {
    print STDERR "$opt_l\n";
} else {
    print STDERR "No minimum\n";
}
print STDERR "Max read size: ";
if($opt_L) {
    print STDERR "$opt_L\n";
} else {
    print STDERR "No maximum\n";
}

print STDERR "Locus\/Loci: ";
if($opt_c) {
    print STDERR "$opt_c\n";
} elsif ($opt_g) {
    print STDERR "$opt_g\n";
} else {
    print STDERR "None \(whole genome\)\n";
}

if($opt_d) {
    print STDERR "Degradome mode activated .. tally depths of 5' ends only\n";
}

if($opt_m) {
    print STDERR "Reads normalized to reads per million mapped before calculating depths\n";
}

## get loci organized and checked, if warranted
my @loci = ();
if($opt_c) {
    push(@loci, $opt_c);
} elsif ($opt_g) {
    @loci = parse_gff3();
    unless($loci[0]) {
	die "FATAL: Failed to find any valid loci after initial parse of file $opt_g\n";
    }
}



my @loci_ok = ();
if($loci[0]) {
    @loci_ok = sort_loci(\@loci, \@$chrs);
    unless($loci_ok[0]) {
	die "FATAL: failed to find any valid loci after sorting\n";
    }
}



############ Write the wiggle file
print STDERR "\nPhase One: Generating depths\n";

# open the WIG output stream
my $wigfile = "$opt_D" . "\/" . "$bamfilename";
# Add notes to the file name if there were non-default options
if($opt_t != 0) {
    $wigfile .= "-t$opt_t";
}
if($opt_s ne "both") {
    $wigfile .= "-s$opt_s";
}
if($opt_r) {
    $wigfile .= "-r$opt_r";
}
if($opt_l) {
    $wigfile .= "-l$opt_l";
}
if($opt_L) {
    $wigfile .= "-L$opt_L";
}
if($opt_c) {
    my $safe_oc = $opt_c;
    $safe_oc =~ s/\:/_/g;  ## remove the colon, no colon in file names please
    $wigfile .= "-c$safe_oc";
}
if($opt_g) {
    my ($gffname,$gffpath,$gffsuffix) = fileparse($opt_g,qr/\.[^.]*/);
    $wigfile .= "-g$gffname";
}
if($opt_d) {
    $wigfile .= "-d";
}
if($opt_m) {
    $wigfile .= "-m";
}
$wigfile .= ".wig";
(open(WIG, "> $wigfile")) || die "Failed to open wiggle output file $wigfile for writing\n";

# open a BAM --> depth output stream
my $depthfile = "$opt_D" . "\/" . "$bamfilename" . "_depth.txt";

########## DEFUNCT
# degradome mode for bottom strand is a special case where the bam needs to be re-sorted...
#if(($opt_d) and ($opt_s eq "bottom")) {
#    (open(DEPTH, "| samtools view -S -b - 2> /dev/null | samtools sort -o - tempbam 2> /dev/null | samtools depth /dev/stdin 2> /dev/null > $depthfile")) || die "failed to open samtools depth stream\n";
#} else {
#    (open(DEPTH, "| samtools view -S -b - 2> /dev/null | samtools depth /dev/stdin 2> /dev/null > $depthfile")) || die "failed to open samtools depth stream\n";
#}
#
# feed in the header
#(open(HEADER, "samtools view -H $bamfilepath |")) || die "failed to open header of $bamfilepath\n";
#while (<HEADER>) {
#    print DEPTH "$_";
#}
#close HEADER;
#################

## If normalizing with option -m, get the number of aligned reads
my $n_reads;
if($opt_m) {
    (open(READS, "samtools view -F 0x4 $bamfilepath | wc -l |")) || die "Failed to count reads in file\n";
    $n_reads = <READS>;
    close READS;
    $n_reads =~ s/\s//g;
    $n_reads =~ s/\n//g;
    print STDERR "\nFound $n_reads mapped reads in file $bamfilepath\n";
}

(open(DEPTH, ">$depthfile")) || die "failed to open depth file $depthfile for writing\n";

my $view_call = "samtools view";
# strand
if($opt_s eq "top") {
    $view_call .= " -F 0x14";  ## filter non-mapped and minus-strand mapped reads
} elsif ($opt_s eq "bottom") {
    $view_call .= " -F 0x4 -f 0x10"; ## filter non-mapped and required minus-strand mapped reads
} else {
    $view_call .= " -F 0x4";  ## filter non-mapped reads
}
# read group
if($opt_r) {
    $view_call .= " -r $opt_r";

}

# add the bamfile
$view_call .= " $bamfilepath";


# initialize variables
my @sam_fields = ();
my $sam_out;
my $last_chr = "NULL";
my %hash = ();
my @sorted = ();
my $pos;

if($loci_ok[0]) {
    # multiple queries
    my $view_call_temp;
    foreach my $lo (@loci_ok) {
	$view_call_temp = $view_call;
	$view_call_temp .= " $lo";
	(open(SAM, "$view_call_temp |")) || die "Failed to open view_call $view_call_temp\n";
	while (<SAM>) {
	    $sam_out = check_sam($_);
	    if($sam_out) {
		# no header                                                                                                                                        
		if($_ =~ /^\@/) {
		    next;
		}
		chomp;
		@sam_fields = split ("\t", $_);
		# no unmapped reads                                                                                                                                
		if($sam_fields[1] & 4) {
		    next;
		}
		
		# record info for current read                                                                                                                     
		for (my $i = $sam_fields[3]; $i < ((length $sam_fields[9]) + $sam_fields[3]); ++$i) {
		    if($opt_m) {
			$hash{$i} += ((1/$n_reads)*1E6);
		    } else {
			++$hash{$i};
		    }
		}
		$last_chr = $sam_fields[2];
	    }
	}
	close SAM;
	@sorted = sort {$a <=> $b} keys %hash;
	foreach $pos (@sorted) {
	    print DEPTH "$last_chr\t";
	    print DEPTH "$pos\t";
	    print DEPTH "$hash{$pos}\n";
	}
	%hash = ();
    }
} else {
    # single query
    (open(SAM, "$view_call |")) || die "Failed to open view_call $view_call\n";
    while (<SAM>) {
	$sam_out = check_sam($_);
	if($sam_out) {
	    # no header                                                                                                                                        
	    if($_ =~ /^\@/) {
		next;
	    }
	    chomp;
	    @sam_fields = split ("\t", $_);
	    # no unmapped reads                                                                                                                                
	    if($sam_fields[1] & 4) {
		next;
	    }
	    
	    if($sam_fields[2] ne $last_chr) {
		unless($last_chr eq "NULL") {
		    # flush all entries from the last chr                                                                                                      
		    @sorted = sort {$a <=> $b} keys %hash;
		    foreach $pos (@sorted) {
			print DEPTH "$last_chr\t";
			print DEPTH "$pos\t";
			print DEPTH "$hash{$pos}\n";
		    }
		    %hash = ();
		}
	    }
	    # record info for current read                                                                                                                     
	    for (my $i = $sam_fields[3]; $i < ((length $sam_fields[9]) + $sam_fields[3]); ++$i) {
		if($opt_m) {
		    $hash{$i} += ((1/$n_reads)*1E6);
		} else {
		    ++$hash{$i};
		}
	    }
	    $last_chr = $sam_fields[2];
	}
    }
    close SAM;
    # final flush for last chr
    @sorted = sort {$a <=> $b} keys %hash;
    foreach $pos (@sorted) {
	print DEPTH "$last_chr\t";
	print DEPTH "$pos\t";
	print DEPTH "$hash{$pos}\n";
    }
    %hash = ();
}
close DEPTH;
print STDERR "\nPhase Two: Processing depths\n";

# (re)initialize variables
$last_chr = "NULL";
my $last_depth;
my $last_position;
my $last_initial_position;
my $chr;
my $position;
my $depth;
my $span;
my @dfields = ();

(open(DEPTH, "$depthfile")) || die "Failed to open depthfile $depthfile\n";
while (<DEPTH>) {
    chomp;
    @dfields = split ("\t", $_);
    $chr = $dfields[0];
    $position = $dfields[1];
    $depth = $dfields[2];
	
    ## fully ignore lines not meeting threshold
    if($depth < $opt_t) {
	next;
    }
    # special case first line crossing threshold, nothing ever written
    if($last_chr eq "NULL") {
	$last_chr = $chr;
	$last_position = $position;
	$last_initial_position = $position;
	$last_depth = $depth;
	$span = 1;
	next;
    }
    # new chromosome .. write info for last span on old chr and reinitialize
    if(($last_chr ne $chr) and ($last_chr ne "NULL")) {
	print WIG "variableStep chrom=$last_chr span=$span\n";
	if($opt_s eq "bottom") {
	    $last_depth = 0 - $last_depth;
	}
	print WIG "$last_initial_position $last_depth\n";
	$last_chr = $chr;
	$last_position = $position;
	$last_initial_position = $position;
	$last_depth = $depth;
	$span = 1;
	next;
    }
    # same chromosome, but discontonous in terms of position .. write last bin and reinitialize
    if(($last_chr eq $chr) and ($position > ($last_position + 1))) {
	print WIG "variableStep chrom=$last_chr span=$span\n";
	if($opt_s eq "bottom") {
	    $last_depth = 0 - $last_depth;
	}
	print WIG "$last_initial_position $last_depth\n";
	$last_chr = $chr;
	$last_position = $position;
	$last_initial_position = $position;
	$last_depth = $depth;
	$span = 1;
	next;
    }
    
    # same chromosome, continous in terms of position, but different depth ... write last bin and reinitialize
    if(($last_chr eq $chr) and ($position == ($last_position + 1)) and ($depth != $last_depth)) {
	print WIG "variableStep chrom=$last_chr span=$span\n";
	if($opt_s eq "bottom") {
	    $last_depth = 0 - $last_depth;
	}
	print WIG "$last_initial_position $last_depth\n";
	$last_chr = $chr;
	$last_position = $position;
	$last_initial_position = $position;
	$last_depth = $depth;
	$span = 1;
	next;
    }
    
    # same chromosome, continuous in position and same depth .. do not update last_initial_position, increase span
    if(($last_chr eq $chr) and ($position == ($last_position + 1)) and ($depth == $last_depth)) {
	++$span;
	$last_chr = $chr;
	$last_position = $position;
	$last_depth = $depth;
	next;
    }
    
    ## all cases covered above. If you are here, code is bad
    die "FATAL: Bad code in loop at line $_\n";
}
close DEPTH;
# report last interval
print WIG "variableStep chrom=$last_chr span=$span\n";
if($opt_s eq "bottom") {
    $last_depth = 0 - $last_depth;
}
print WIG "$last_initial_position $last_depth\n";
close WIG;

system "rm -f $depthfile";
print STDERR "\nWiggle file written to $wigfile\n";

if($w2bw_check) {
    print STDERR "\nPhase Three: Writing bigwig file\n";
    # write a chrom length file from the BAM information
    (open(LEN, ">chr_lens_temp.txt")) || die "FATAL: Failed to open a temp file for chrom lengths\n";
    foreach my $cname (@$chrs) {
	my $clen = shift @$chr_lens;
	print LEN "$cname\t$clen\n";
    }
    close LEN;
    my $bigwigfile = $wigfile;
    $bigwigfile =~ s/\.wig/\.bigwig/g;
    
    system "wigToBigWig $wigfile chr_lens_temp.txt $bigwigfile";
    system "rm -f chr_lens_temp.txt";
    print STDERR "bigwig file at $bigwigfile\n";
}

print STDERR "\nCompleted at ";
print STDERR `date`;
exit;



############
# SUB-ROUTINES
sub check_sam {
    my($samline) = @_;
    unless(($opt_l) or ($opt_L) or ($opt_d)) {
	return $samline;
    }
    my @fields = split ("\t", $samline);
    
    # test
    #print STDERR "in check_sam, length is ";
    #my $xxx = length $fields[9];
    #print STDERR "$xxx\n";
    
    my $newseq;
    if($opt_l) {
	if((length $fields[9]) < $opt_l) {
	    # test
	    #print STDERR "\tFailed less than opt_l $opt_l\n";
	    return 0;
	}
    }
    if($opt_L) {
	if((length $fields[9]) > $opt_L) {
	    # test
	    #print STDERR "\tFailed more than than opt_L $opt_L\n";
	    return 0;
	}
    }
    if($opt_d) {
	# process for degradome ... pretend it is a one-nt read
	if($fields[1] & 16) {
	    # minus strand process .. position is the END position ...
	    my $newpos = $fields[3] + (length $fields[9]) - 1;
	    $samline = "$fields[0]\t$fields[1]\t$fields[2]\t$newpos\t255\t1M\t\*\t0\t0\t";
	    # seq is the first letter of the revcomp
	    my $rc = reverse $fields[9];
	    $rc =~ tr/ACTG/TGAC/;
	    $newseq = substr($rc,0,1);
	    $samline .= "$newseq\t\*\n";
	    return $samline;
	} else {
	    # plus strand process
	    $samline = "$fields[0]\t$fields[1]\t$fields[2]\t$fields[3]\t255\t1M\t\*\t0\t0\t";
	    # seq is the first letter
	    $newseq = substr($fields[9],0,1);
	    $samline .= "$newseq\t\*\n";
	    return $samline;
	}
    } else {
	# if not degradome, must have passed the opt_l and/or opt_L filters already, so return it
	return $samline;
    }
    # why are you here?
    return 0;
}

sub sort_loci {
    my($in,$chrs) = @_; # by reference. All arrays
    
    # first, define the allowable names of chromosomes
    my %chr_names = ();
    foreach my $c (@$chrs) {
	$chr_names{$c} = 1;
    }
    
    # hash the clusters by chromosome name
    my %hash = ();
    my $no_suffix;
    foreach my $clus_name (@$in) {
	if($clus_name =~ /^(\S+):(\d+-\d+)$/) {
	    if(exists($chr_names{$1})) {
		push(@{$hash{$1}}, $2);
	    } else {
		die "FATAL: Chromosome $1 is not defined in your alignment\!\n";
	    }
	} else {
	    die "FATAL in sub-routine sort_clusters .. could not parse cluster name $clus_name\n";
	}
    }
    # define the output array
    my @output = ();
    # holder for any warnings .. for overlapping clusters
    my @warnings = ();
    # go chromosome by chromosome
    foreach my $o_c (@$chrs) {
	if(exists($hash{$o_c})) {
	    my %lefts = ();
	    foreach my $full (@{$hash{$o_c}}) {
		my @full_f = split ("-", $full);
		$lefts{$full_f[0]} = $full_f[1];
	    }
	    # sort ascending by left-most coordinate
	    my @lefts_sorted = sort {$a <=> $b} keys %lefts;
	    # file them away, flagging any cases where the left-most of the next cluster is <= the right-most of the previous
	    my $last_right = -2;
	    my $last_name = "NULL";
	    foreach my $left (@lefts_sorted) {
		my $name = "$o_c" . ":" . "$left" . "-" . "$lefts{$left}";
		push(@output, $name);
		if($left <= $last_right) {
		    push(@warnings,$last_name);
		    push(@warnings,$name);
		}
		$last_right = $lefts{$left};
		$last_name = $name;
	    }
	}
    }
    # FATAL if warnings. No overlaps allowed by this program
    if(@warnings) {
	print STDERR "FATAL: The loci defined in your GFF3 file overlap with each other. This is not allowed. Offending loci:\n";
	while (@warnings) {
	    my $name1 = shift @warnings;
	    my $name2 = shift @warnings;
	    print STDERR "$name1 and $name2\n";
	    exit;
	}
    }
    return @output;
}



sub parse_gff3 {
    (open(GFF3, "$opt_g")) || die "FATAL: in sub-routine parse_gff3. Could not open file $opt_g\n";
    my @fields = ();
    my @loci = ();
    my $locus;
    while (<GFF3>) {
	if ($_ =~ /^\#/) {
	    next;
	} else {
	    @fields = split ("\t", $_);
	    if((($opt_s eq "top") and ($fields[6] eq "+")) or
	       (($opt_s eq "bottom") and ($fields[6] eq "-")) or
	       ($opt_s eq "both")) {
		$locus = "$fields[0]" . ":" . "$fields[3]" . "-" . "$fields[4]";
		push(@loci,$locus);
	    }
	}
    }
    close GFF3;
    return @loci;
}

sub checkbam {
    my($bamfilepath) = @_;
    unless(-r $bamfilepath) {
	return "0";
    }
    my $baifilepath = "$bamfilepath" . ".bai";
    unless(-r $baifilepath) {
	print STDERR "BAM index $baifilepath was not found. Generating using samtools index\n";
	system "samtools index $bamfilepath";
	# check for success
	unless (-r $baifilepath) {
	    print STDERR "Failed to build BAM index\n";
	    return "0";
	}
    }
    
    # Validate header ... hash the sequences and lengths, and check for read groups
    my @chrs = ();
    my @chr_lens = ();
    my @read_groups = ();
    (open(HEADER, "samtools view -H $bamfilepath |")) || return 0;
    while (<HEADER>) {
	chomp;
	if($_ =~ /^\@SQ/) {
	    # a SQ line. Capture SN And LN
	    if($_ =~ /SN:(\S+)/) {
		push(@chrs, $1);
	    }
	    if($_ =~ /LN:(\d+)/) {
		push(@chr_lens, $1);
	    }
	} elsif ($_ =~ /^\@RG/) {
	    # a read-group line
	    if($_ =~ /ID:(\S+)/) {
		push(@read_groups, $1);
	    }
	}
    }
    close HEADER;
    # Report or fail
    unless($chrs[0]) {
	die "FATAL: No \@SQ information found in the BAM header. BAM header must include full header information about the chromosomes\!\n";
    }
    print STDERR "\tChromosomes in genome from $chrs[0] to $chrs[-1]\n";
    print STDERR "\tRead Groups in BAM file: ";
    if($read_groups[0]) {
	my $rg = join (",", @read_groups);
	print STDERR "$rg\n";
    } else {
	print STDERR "None\n";
    }
    
    return(\@chrs,\@chr_lens,\@read_groups);
    
}

__END__

=head1 SYNOPSIS
bam2wig
Conversion of a BAM alignment to wiggle and bigwig coverage files, with flexible reporting options.
=head1 LICENSE
Copyright (C) 2014 Michael J. Axtell                                                                                                                                                           
                                                                                                                                                                                                     
This program is free software: you can redistribute it and/or modify                                                                                                                                 
it under the terms of the GNU General Public License as published by                                                                                                                                 
the Free Software Foundation, either version 3 of the License, or                                                                                                                                    
(at your option) any later version.                                                                                                                                                                  
                                                                                                                                                                                                     
This program is distributed in the hope that it will be useful,                                                                                                                                      
    but WITHOUT ANY WARRANTY; without even the implied warranty of                                                                                                                                   
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the                                                                                                                                        
GNU General Public License for more details.                                                                                                                                                         
                                                                                                                                                                                                     
You should have received a copy of the GNU General Public License                                                                                                                                    
along with this program.  If not, see <http://www.gnu.org/licenses/>
=head1 CITATION
None yet..
=head1 AUTHOR
Michael J. Axtell, mja18@psu.edu
=head1 INSTALL
bam2wig is a perl program and so it requires perl to compile. It expects that perl is located at /usr/bin/perl on your system. If not, please modify line 1 (the hashbang) of the script accordingly. It was built and tested using perl 5.
Depending on your system, you may need super-user privileges to install bam2wig and/or its dependencies.
samtools <http://www.htslib.org/> is a required dependency. Install samtools to your PATH before running bam2wig.
wigToBigWig <http://genome.ucsc.edu/goldenPath/help/bigWig.html> is an optional dependency. If wigToBigWig is in your PATH, bam2wig will make a bigwig file in addition to a wiggle file.
Once the dependencies are installed to your PATH, install bam2wig to your PATH for convenience
=head1 USAGE
bam2wig [options] [alignments.bam]
=head1 OPTIONS
-t [float] : Threshold for reporting. Defaults to 0 (all depths reported)
-s [top/bottom/both] : Strand to report. Defaults to both. 'bottom' will quantify in negative numbers. Used in conjunction with option -g will subset the loci to those on the indicated strand.
-r [string] : Limit analysis to specified read group
-l [integer] : Minimum size of read to report. Default: no minimum.
-L [integer] : Maximum size of read to report. Default: no maximum.
-c [chr:start-stop] : Limit analysis to the indicated locus. 
-g [string] : Path to a GFF3 file. Loci in the GFF3 will be reported. Negates any setting of -c.
-D [string] : Name of output directory. Defaults to [/YourBamsPath/YourBamsBasename_bam2wig/]
-m : Count each as (1/(n mapped reads)) * 1E6. In other words, normalize to reads per million.
-d : Treat as degradome/PARE data .. tally only the depths of 5' ends. Option -s must be either 'bottom' or 'top'
-h : Get help (print this message)
-v : Print version number and quit
=head1 OUTPUT
A wiggle file is written to the working directory. The wiggle file LACKS a track line at the top. If wigToBigWig is installed, a bigwig file will also be written to the working directory. The file name will include information on non-default settings. For instance, file name 'mydata-t10-l20-L24.wig' indicates an analysis using option t set to 10, option l set to 20, and option L to to 24, with all other options left at their defaults.
=head1 BAM FILE REQUIREMENTS
The input BAM file must be sorted by chromosome, indexed (use samtools index), and possess a valid header with all @SQ lines present and the sequence names and lengths indicated by SN and LN, respectively (see SAM format specification). If option -r is used, the read groups must also be specified in the header using the @RG lines, per the SAM spec.
=head1 USING A GFF3 FILE
Using option -g allows the output to be restricted to the loci specified in gff3 file. When used under default setting of option -s (the strand set to 'both'), all loci from the GFF3 will be reported, and the values will represent the read coverage on BOTH strands, regardless of whether the gff3 locus is annotated as strand '+', '-', or '.'. When used with option -s set to "bottom", only GFF3 loci whose strand is annotated as "-" in the file will be used, and only the reads on the - strand of the reference analyzed. Conversely, when used with option -s set to "top", only GFF3 loci whose strand is annotated as "+" in the file will be used, and only the reads on the + strand of the reference analyzed.
=head1 DEGRADOME / PARE MODE
Option -d triggers activation of degradome / PARE mode. In this mode, the depth of coverage of 5' ends of reads is tallied, instead of total depth of coverage across the whole bodies of the reads. This mode must have option -s (strand) set to either "top" or "bottom".
=cut
