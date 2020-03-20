#!/usr/bin/env perl
use strict;
use warnings;
print '#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long;
use Cwd;
use FindBin qw($Bin);
use constant USAGE=><<EOH;

SYNOPSIS:

perl $0 --input my.fa [Options]
Version: v20200320

Requirements:
    Programs:
    Modiles: Scalar::Util, Cwd, Getopt::Long, FindBin

Descriptions:
    Determine the insert size given pairs of seqing data by
    mapping them to a reference.

Options:
    --help|-h
        Print this help/usage;
    --input|-i
        
    --verbose
        Detailed output for trouble-shooting;
    --version|v!
        Print current SCRIPT version;

Example:
    perl $0 

Author:
    Fu-Hao Lu
    Professor, PhD
    State Key Labortory of Crop Stress Adaptation and Improvement
    College of Life Science
    Jinming Campus, Henan University
    Kaifeng 475004, P.R.China
    E-mail: lufuhao@henu.edu.cn
EOH
###HELP ends#########################################################
die USAGE unless @ARGV;



###Receving parameter################################################
my ($help, $verbose, $debug, $version);
my ($input, $output);

GetOptions(
	"help|h!" => \$help,
	"input|i:s" => \$input,
#	"output|o:s" => \$output,
#	!:s:i
	"debug!" => \$debug,
	"verbose!" => \$verbose,
	"version|v!" => \$version) or die USAGE;
($help or $version) and die USAGE;



### Defaults ########################################################
$debug=0 unless (defined $debug);
$verbose=0 unless (defined $verbose);



### input and output ################################################



### Main ############################################################




#####################################################################
###                         sub functions                         ###
#####################################################################
### ReadSam
###&ReadSam(sam,ref, 1/2/3)
###Global:
###Dependency:
###Note:

';

#    Fu-Hao Lu
#    Post-Doctoral Scientist in Micheal Bevan laboratory
#    Cell and Developmental Department, John Innes Centre
#    Norwich NR4 7UH, United Kingdom
#    E-mail: Fu-Hao.Lu\@jic.ac.uk
