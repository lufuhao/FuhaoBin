#!/usr/bin/perl-w
use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;
use File::Basename;
# Date: 02-02-2010
# Modified on 01-02-1011 to parse r output as p-value = 0.005 or p-value < 0.005 or p-value > 0.005. priviously it was only looking for p-value = 0.005
# Modified on 11-10-2010 to print non-SNP position also in summary p-value file.
# Modified on 03-09-2010 to stop the program in case the maximum number of simulations allowed by R is reached before re-simulating (otherwise it re-simulates using 1,000,000,000 simulations which is beyond the limit of R)
# Modified on 29-07-2010 to write wiggle output and re-calculate p-value for last best 100 or 1000 SNP on higher number of simulation.
# Modified on 29-06-2010 to allow multi round of filtering of SNPs based on p-value.
# Modified on 08-03-2010 to fix bug in SNP calling function and discarding line with no coverage information.
# Modified on 09-03-2010 to allow snp wise or window wise p-value calulation.
# Modified on 06-04-2010 to implement uniform Window sliding module.

# Author: Ram Vinay Pandey 

# Define the variables
my $input;
my $output="";
my $help=0;
my $test=0;
my $verbose=1;

my $windowsize=1000;
my $step=1000;
my $windowunit="bp";
my $mincount=2;
my $mincoverage=4;
my $maxcoverage=1000;
my $numsimulations=100;
my $minCoverageFraction=0.6;
my $minpvalue=0.00001;
my $min_snp_fraction=0.8;
my $min_best_snp_loci=100;
my $magnitude=10;
my $trackname="p-value";
my $ucscfilter="";
my $ucscprepend="chr";
my $population;
my $nonsnp=0;
# This upper limit is defined due to R limitation: R does not support more than 100 million simulations.
my $max_simulation_limit = 100000000;
#my $max_simulation_limit = 100;

my $usage="perl $0 --min-count 2 --min-coverage 4 --max-coverage 1000 --min-snp-fraction 0.8 --min-covered-fraction 0.6 --number-simulations 100 --min-best-snp-loci 100 --window-size 100 --step-size 10 --window-unit snp --input synchronized-file.txt --output snp-probability.txt --trackname unknown --ucsc-filter 2L --ucsc-prepend dmel\n";

GetOptions(
    "input=s"       =>\$input,
    "output=s"      =>\$output,
    "min-count=i"   =>\$mincount,
    "min-coverage=i"=>\$mincoverage,
    "max-coverage=i"=>\$maxcoverage,
    "min-snp-fraction=f"=>\$min_snp_fraction,
    "number-simulations=i"=>\$numsimulations,
    "min-best-snp-loci=i"=>\$min_best_snp_loci,
    "window-unit=s"  =>\$windowunit,
    "window-size=i"  =>\$windowsize,
    "ucsc-filter=s"     =>\$ucscfilter,
    "ucsc-prepend=s"    =>\$ucscprepend,
    "trackname=s"       =>\$trackname,
    "step-size=i"   =>\$step,
    "min-covered-fraction=f"=>\$minCoverageFraction,
    "population=s"          =>\$population,
    "also-print-nonsnp"          =>\$nonsnp,
    "test"          =>\$test,
    "help"          =>\$help
) or pod2usage(-msg=>"Wrong options",-verbose=>1);

pod2usage(-verbose=>2) if $help;
Test::runTests() if $test;

my $paramfile=$output.".params";
open my $pfh, ">",$paramfile or die "Could not open $paramfile\n";
print $pfh "Using input\t$input\n";
print $pfh "Using output\t$output\n";
print $pfh "Using min-count\t$mincount\n";
print $pfh "Using min-coverage\t$mincoverage\n";
print $pfh "Using max-coverage\t$maxcoverage\n";
print $pfh "Using min-snp-fraction\t$min_snp_fraction\n";
print $pfh "Using number-simulations\t$numsimulations\n";
print $pfh "Using min-best-snp-loci\t$min_best_snp_loci\n";
print $pfh "Using window-unit\t$windowunit\n";
print $pfh "Using window-size\t$windowsize\n";
print $pfh "Using ucsc-filter\t$ucscfilter\n";
print $pfh "Using ucsc-preprend\t$ucscprepend\n";
print $pfh "Using trackname\t$trackname\n";
print $pfh "Using step-size\t$step\n";
print $pfh "Using min-covered-fraction\t$minCoverageFraction\n";
print $pfh "Using population\t$population\n";
print $pfh "Using also-print-nonsnp\t$nonsnp\n";
print $pfh "Using test\t$test\n";
print $pfh "Using help\t$help\n";
close $pfh;



$minpvalue = 2/$numsimulations;

my @simulations = ();
my $next_sim = $numsimulations;
for my $i (1..50) {
    
    $next_sim = $next_sim*$magnitude;
    push(@simulations,$next_sim);
    
}


#pod2usage(-msg=>"Minimum count for SNP calling must be 1 or larger than 1\n",-verbose=>1) if $mincount<1;
#pod2usage(-msg=>"Minimum coverage <1 not allowed",-verbose=>1) if $mincoverage<1;
#pod2usage(-msg=>"Number of simulation must be larger than 0\n",-verbose=>1) if $numsimulations<1;
    
pod2usage(-msg=>"A input file has to be provided\n",-verbose=>1) unless -e $input;
pod2usage(-msg=>"A output file has to be provided\n",-verbose=>1) unless $output;
pod2usage(-msg=>"Number of simulation must be atleast 100\n",-verbose=>1) if $numsimulations<100;
pod2usage(-msg=>"min_snp_fraction should be between 0 - 1\n",-verbose=>1) if (($min_snp_fraction<0) or ($min_snp_fraction>1));

my ( $name, $path, $extension ) = File::Basename::fileparse ( $output, '\..*' );
my ($to_delete_path,$to_delete_name) = ($path,$name);
my $output1 = $path.$name."1" . $extension;
my $output2 = $path.$name."-discarded-1" . $extension;
my $rbatchfile1 = $path.$name."1.r";
my $routputfile1 = $path.$name. "1-r-outputfile.txt";

my $summary_pvalue_output = $path.$name."-summary-pvalue" . $extension;
my $summary_pvalue_wiggle_output = $path.$name."-summary-pvalue.wig";
#=cut

open my $ofh1, ">$output1" or die "Could not open output file";
open my $ofh2, ">$output2" or die "Could not open output file";
open my $orfh1, ">$rbatchfile1" or die "Could not open output file";

#write_pvalue_in_wiggle_format($wiggle,$windowsize,$trackname,$ucscfilter,$ucscprepend);
#print "$input,$windowsize,$step,$mincount,$mincoverage,$numsimulations\n";

# population count
my $popcountObj = PopCount->new($input);
my $popCount=$popcountObj->count_samples();


####### Taking population defined by user

### Taking populations
my @population = ();
if ($population) {
    @population = split('\,',$population);
}
my $input_population = 0;
$input_population = scalar(@population);

if (($input_population<2) and ($popCount<2)) {
        
        pod2usage(-msg=>"invalid population count in syncronized file. Population count should be atleast 2 ",-verbose=>1);
        
}
elsif(($input_population<=1) and ($popCount>=2)) {

        @population = (1..$popCount);
}

elsif($input_population>$popCount) {
        
        @population = (1..$popCount);
}

#print "@population,,",scalar(@population),"\t$popCount\n";

#exit;
#######################################

my $reader;
if("bp" eq lc($windowunit))
{
    $reader=BpSlider->new($input,$windowsize,$step,$mincount,$mincoverage,$maxcoverage,$numsimulations);
}
elsif("snp"eq lc($windowunit))
{
    $windowsize = 1;
    $step = 1;
    $reader=BpSlider->new($input,$windowsize,$step,$mincount,$mincoverage,$maxcoverage,$numsimulations);
}
else
{
    pod2usage(-msg=>"invalid window unit; has to be bp or snp",-verbose=>1);
}

print "Started at: \t". localtime() .".....\n";
# population count
#my $popCount=$reader->count_samples();

my %snp_loci=();


while(my $window=$reader->nextWindow())
{
        my $chr=$window->{chr};
        my $pos=$window->{middle};
        my $win=$window->{window};
        my $above=$window->{count_covered};
        my $snpcount=$window->{countpuresnp};
        my $data=$window->{data};
        
        next unless @$data;
        #print "HI: $chr\t$pos\n";
        if ($nonsnp) {
            if ("bp" eq lc($windowunit)) {
                    $snp_loci{"$chr\t$pos"} = "1\t1\t$snpcount";
            }
            elsif("snp"eq lc($windowunit)) {
                    $snp_loci{"$chr\t$pos"} = "1\t$snpcount";
            }
        }
        my $snps= [grep {$_->{ispuresnp}} @$data];
        die "data not fitting" unless @$snps == $snpcount;

        my $coveredFrac=$above/$win;

        if ($coveredFrac>=$minCoverageFraction) {
            
            Utility::calculatePvalues($snps,$popCount,$pos,$chr,$orfh1);
        
        }
        else {
            if ($nonsnp) {
                if ("bp" eq lc($windowunit)) {
                    $snp_loci{"$chr\t$pos"} = "1\t1\t$snpcount";
                }
                elsif("snp"eq lc($windowunit)) {
                    $snp_loci{"$chr\t$pos"} = "1\t$snpcount";
                }
            }
        }
 
}

        

Utility::runR($rbatchfile1,$routputfile1);

%snp_loci = Utility::Parse_R_output($rbatchfile1,$routputfile1,$windowunit,$ofh1,$ofh2,$minpvalue,$min_snp_fraction,\%snp_loci);
close $orfh1;
close $ofh1;
close $ofh2;



##### If first output itself has less snp than select min-best-snp-loci then skip the loop of simulation and directly jump to simulation with maximum number of simulation.
my $snp_count_in_privious_file=0;
my $actual_simulation = 0;
my $ct=1;
my $ct1=0;
my %file_count = ();

my $yesno = "";
my $pvalue = 0;
my $nosimulation = 0;




my ($privious_output,$privious_rbatchfile,$privious_routputfile) = ("","","");
my ($next_output,$next_output_discarded,$next_rbatchfile,$next_routputfile) = ("","","","");

my $processed_simulation = 0;

# See if the privious fisher-exact-test.pl output file exists and contains something
        die "\n\n\t=> The privious output file: $output1 does not exists.\n\n" unless (-e $output1);
        die "\n\n\t=> The privious output file: $output1 appears to be empty!\n\n" if (-z $output1);
        
$snp_count_in_privious_file = Utility::snpCount_in_file($output1);

#print "Hi:$snp_count_in_privious_file\t$min_best_snp_loci\n";
#exit;
if ($snp_count_in_privious_file<=$min_best_snp_loci) {

    $privious_output=$output1;
    $actual_simulation = $numsimulations*$magnitude;
    $nosimulation = $max_simulation_limit;
    #$actual_simulation = $actual_simulation/$magnitude;

    $ct1=1;
    goto LEBEL1;
}

#exit;




foreach my $ns (@simulations) {
    
    $ct++;
    $ct1++;
    
    $file_count{$ct} = $ct;
    $file_count{$ct1} = $ct1;
    
    $privious_output = $path.$name.$ct1. $extension;
    $privious_rbatchfile = $path.$name.$ct1.".r";
    $privious_routputfile = $path.$name. $ct1."-r-outputfile.txt";
    
    $next_output = $path.$name.$ct. $extension;
    $next_output_discarded = $path.$name."-discarded-".$ct. $extension;
    $next_rbatchfile = $path.$name.$ct.".r";
    $next_routputfile = $path.$name. $ct."-r-outputfile.txt";
    
    $nosimulation = $ns;
    chomp($nosimulation);
    $pvalue = 2/$nosimulation;
    
    $snp_count_in_privious_file = Utility::snpCount_in_file($privious_output);

    if ($snp_count_in_privious_file==0) {
        
        Utility::writeSummary_pvalue();
        Utility::write_pvalue_in_wiggle_format();
        
    }
    elsif ($snp_count_in_privious_file<=$min_best_snp_loci) {
        goto LEBEL;
    }
    

        # See if the privious fisher-exact-test.pl output file exists and contains something
        die "\n\n\t=> The privious output file: $privious_output does not exists.\n\n" unless (-e $privious_output);
        die "\n\n\t=> The privious output file: $privious_output appears to be empty!\n\n" if (-z $privious_output);
        
        
        
        
        #print "$ct1\t$ct\t$pvalue\t$nosimulation\t$name\t$path\t$extension\n\n";
        
        #print "$privious_output\t$privious_rbatchfile\t$privious_routputfile\t$next_output\t$next_rbatchfile\t$next_routputfile\n\n";
        
        open my $ofh2, ">$next_output" or die "Could not open output file";
        open my $ofh3, ">$next_output_discarded" or die "Could not open output file";
        open my $orfh2, ">$next_rbatchfile" or die "Could not open output file";
    
        NextTire::calculatePvalues($privious_output,$privious_rbatchfile,$nosimulation,$orfh2);
        
        Utility::runR($next_rbatchfile,$next_routputfile);
        
        %snp_loci = Utility::Parse_R_output($next_rbatchfile,$next_routputfile,$windowunit,$ofh2,$ofh3,$pvalue,$min_snp_fraction,\%snp_loci);
        
        $processed_simulation = $nosimulation;
        
        if ($min_best_snp_loci<1) {
            
            if ($snp_count_in_privious_file==0) {
                Utility::writeSummary_pvalue();
                Utility::write_pvalue_in_wiggle_format();
                #last;
            }
            
        }
        elsif($nosimulation>=$max_simulation_limit) {
            Utility::writeSummary_pvalue();
            Utility::write_pvalue_in_wiggle_format();
            last;
        }
        elsif($nosimulation>=$max_simulation_limit) {
            Utility::writeSummary_pvalue();
            Utility::write_pvalue_in_wiggle_format();
            last;
        }
        else {
            if ($snp_count_in_privious_file<=$min_best_snp_loci) {
                goto LEBEL;
            }
        }

}



LEBEL:

#print "CT: $ct\t$ct1\n";
my $last_file_count = $ct;

while ($ct1>=1) {
    
    $ct--;
    $ct1--;
    #print "CT: $ct\t$ct1\n";
    
    $privious_output = $path.$name.$ct1. $extension;
    my $nosimulation1 = $nosimulation/$magnitude;
    $snp_count_in_privious_file = Utility::snpCount_in_file($privious_output);
    
    if ($snp_count_in_privious_file>=$min_best_snp_loci) {
        
        #print "$ct\t$ct1\t$nosimulation\t$snp_count_in_privious_file\t$privious_output\n";
        $actual_simulation = $nosimulation1;
        goto LEBEL1;
    }
    
    
}

#$nosimulation = $max_simulation_limit;
#$actual_simulation = $actual_simulation/$magnitude;

LEBEL1:

$nosimulation = $max_simulation_limit;
$actual_simulation = $actual_simulation/$magnitude;

#print "PRIV: $privious_output\t$ct\t$ct1\t$nosimulation\t$actual_simulation\n\n";

    $pvalue = 2/$actual_simulation;
    #print "PRIV: $privious_output\t$ct\t$ct1\t$nosimulation\t$actual_simulation\t$pvalue\n\n";

    $privious_output = $path.$name.$ct1. $extension;
    $privious_rbatchfile = $path.$name.$ct1.".r";
    $privious_routputfile = $path.$name. $ct1."-r-outputfile.txt";
    
    $next_output = $path.$name."-final". $extension;
    $next_output_discarded = $path.$name."-final-discarded-". $extension;
    $next_rbatchfile = $path.$name."-final.r";
    $next_routputfile = $path.$name."-final-r-outputfile.txt";
    
    # See if the privious fisher-exact-test.pl output file exists and contains something
    die "\n\n\t=> The privious output file: $privious_output does not exists.\n\n" unless (-e $privious_output);
    die "\n\n\t=> The privious output file: $privious_output appears to be empty!\n\n" if (-z $privious_output);
        
        
    open my $ofh22, ">$next_output" or die "Could not open output file";
    open my $ofh33, ">$next_output_discarded" or die "Could not open output file";
    open my $orfh22, ">$next_rbatchfile" or die "Could not open output file";
    
    NextTire::calculatePvalues($privious_output,$privious_rbatchfile,$nosimulation,$orfh22);
        
    Utility::runR($next_rbatchfile,$next_routputfile);
        
    %snp_loci = Utility::Parse_R_output($next_rbatchfile,$next_routputfile,$windowunit,$ofh22,$ofh33,$pvalue,$min_snp_fraction,\%snp_loci);
        
    $processed_simulation = $nosimulation;
        
    
    
    #print "$last_file_count\t$processed_simulation\n";
    
    Utility::writeSummary_pvalue();
    Utility::write_pvalue_in_wiggle_format();
    
    #unlink($next_rbatchfile);
    #unlink($next_rbatchfile);
            
#print "$ct\t$ct1\t$nosimulation\t$snp_count_in_privious_file\n";

foreach my $key (keys %snp_loci) {
    #print "P-VALUE:\t$key\t\t$snp_loci{$key}\n";
}

#print "size: ",scalar(keys %snp_loci),"\n";
=cut

do {
    $ct++;
    $ct1++;
    
    $file_count{$ct} = $ct;
    $file_count{$ct1} = $ct1;
    
    $privious_output = $path.$name.$ct1. $extension;
    $privious_rbatchfile = $path.$name.$ct1.".r";
    $privious_routputfile = $path.$name. $ct1."-r-outputfile.txt";
    
    $next_output = $path.$name.$ct. $extension;
    $next_output_discarded = $path.$name."-discarded-".$ct. $extension;
    $next_rbatchfile = $path.$name.$ct.".r";
    $next_routputfile = $path.$name. $ct."-r-outputfile.txt";
    
    print "You have already run till $processed_simulation similations, do you want to run next time? [yes/no] : ";
    $yesno = <STDIN>;
    chomp($yesno);

    if (($yesno eq "yes") or ($yesno eq "Yes") or ($yesno eq "YES") or ($yesno eq "Y") or ($yesno eq "y")) {
        
        
        # See if the privious fisher-exact-test.pl output file exists and contains something
        die "\n\n\t=> The privious output file: $output1 does not exists.\n\n" unless (-e $privious_output);
        die "\n\n\t=> The privious output file: $output1 appears to be empty!\n\n" if (-z $privious_output);
    
        #print "Enter p-value cutoff for next tier of fisher's exact test : ";
        #$pvalue = <STDIN>;
        #chomp($pvalue);
        print "Enter number of simulation for next tier of fisher's exact test : ";
        $nosimulation = <STDIN>;
        chomp($nosimulation);
        
        $pvalue = 2/$nosimulation;
        #print "$ct1\t$ct\t$pvalue\t$nosimulation\t$name\t$path\t$extension\n\n";
        
        #print "$privious_output\t$privious_rbatchfile\t$privious_routputfile\t$next_output\t$next_rbatchfile\t$next_routputfile\n\n";
        
        open my $ofh2, ">$next_output" or die "Could not open output file";
        open my $ofh3, ">$next_output_discarded" or die "Could not open output file";
        open my $orfh2, ">$next_rbatchfile" or die "Could not open output file";
    
        NextTire::calculatePvalues($privious_output,$privious_rbatchfile,$nosimulation,$orfh2);
        
        Utility::runR($next_rbatchfile,$next_routputfile);
        
        Utility::Parse_R_output($next_rbatchfile,$next_routputfile,$windowunit,$ofh2,$ofh3,$pvalue);
        
    }
    #else {
        #exit;
    #}

} until (($yesno eq "no") or ($yesno eq "No") or ($yesno eq "NO") or ($yesno eq "N") or ($yesno eq "n") or ($yesno eq ""));

=cut
    if (scalar(keys %file_count)>0) {
        foreach my $key (keys %file_count) {
            my $rbatchfile = $to_delete_path.$to_delete_name.$key.".r";
            my $routputfile = $to_delete_path.$to_delete_name. $key."-r-outputfile.txt";
            my $considered_output = $path.$name.$key. $extension;
            my $discarded_output = $path.$name."-discarded-".$key. $extension;
            #unlink($rbatchfile);
            #unlink($routputfile);
            
            if ($key>$ct) {
                #unlink($considered_output);
                #unlink($discarded_output);
            }

        }
    }

print "Complete at: \t". localtime() .".....\n";



exit;




{
    use warnings;
    use strict;
    package BpSlider;

    sub new
    {
        my $class=shift;
        my $file=shift;
        my $window=shift;
        my $step=shift;
        my $mincount=shift;
        my $mincov=shift;
        my $maxcov=shift;
        my $nosimulations=shift;
        
        open my $fh,"<$file" or die "Could not open file handle";
        
        return bless {
            lower=>0,
            upper=>$window,
            window=>$window,
            step=>$step,
            fh=>$fh,
            mincount=>$mincount,
            mincov=>$mincov,
            maxcov=>$maxcov,
            nosim=>$nosimulations,
            curwin=>[],
            buffer=>[]
        },__PACKAGE__;
    }
    
    sub count_samples
    {
        my $self=shift;
        my $l=$self->_nextline();
        my $p=Utility::_parseLightwight($l,$self->{mincount},$self->{mincov},$self->{maxcov});
        my $c=scalar(@{$p->{samples}});
        $self->_bufferline($l);
        return $c;
    }
    
    sub nextWindow
    {
        my $self=shift;
        my $mincount=$self->{mincount};
        my $mincov=$self->{mincov};
        my $maxcov=$self->{maxcov};
        my $nosimulations=$self->{nosim};
        
        #get the current window, and the current chromosome
        my $curwin=$self->{curwin};
        
        my $curChr="";
        $curChr=$curwin->[0]{chr} if @$curwin;
        
        my $resetchr=0;
        
        # empty unnecessary entries
        EMPTY: while(@$curwin)
        {
            my $e=shift @$curwin;
            if($e->{pos}>$self->{lower})
            {
                unshift @$curwin, $e;
                last EMPTY;
            }
            
        }
        
        # fill with novel entries
        my $line;
        FILL:while($line=$self->_nextline)
        {
            my $e=Utility::_parseLightwight($line,$mincount,$mincov,$maxcov,$nosimulations);
            $curChr=$e->{chr} unless $curChr;
            
            
            if($e->{chr} eq $curChr && $e->{pos} <= $self->{upper})
            {
                push @$curwin,$e;
            }
            else
            {
                $resetchr=1 if $e->{chr} ne $curChr;
                $self->_bufferline($line);
                last FILL;
            }
        }
        
        return undef unless $curChr;
        
        
        my $toret=Utility::_annotateWindow($curwin,$curChr,$self->{lower},$self->{upper},$self->{window});
        
        if($resetchr or not defined($line))
        {
            # we transgressed the boundaries to the next chromosome
            # reset the windows and the current buffer
            $self->{lower}=0;
            $self->{upper}=$self->{window};
            $self->{curwin}=[];
        }
        else
        {
            # next time we will still be in the same chromosome
            # increase the upper and lower boundaries by the stepsize and set the current buffer
            $self->{upper}+=$self->{step};
            $self->{lower}+=$self->{step};
            $self->{curwin}=$curwin;
        }

        return $toret;
    }
    
    
    
    sub _nextline
    {
        my $self=shift;
        my $fh=$self->{fh};
        my $buffer=$self->{buffer};
        
        return shift @$buffer if @$buffer;
        return <$fh>;
    }
    
    sub _bufferline
    {
        my $self=shift;
        my $line=shift;
        push @{$self->{buffer}},$line;
    }
    
    
}



{
    package Utility;
    use strict;
    use warnings;
    use List::Util qw[min max];
    
    sub writeSummary_pvalue {
        
        #my $ofh=shift;
        open my $summary, ">$summary_pvalue_output" or die "Could not open output file";
        my @snp_loci = ();
        my @vals = ();
        foreach my $key (keys %snp_loci) {
            my @key = split("\t",$key);
            my @value = split("\t",$snp_loci{$key});
            $vals[0] = $key[0];
            $vals[1] = $key[1];
            
            if("bp" eq lc($windowunit))
            {
                $vals[2] = $value[0];
                $vals[3] = $value[1];
                $vals[4] = $value[2];
            }
            elsif("snp"eq lc($windowunit))
            {
                $vals[2] = $value[0];
                $vals[4] = $value[1];
            }
            
            push @snp_loci,[@vals];
            @vals = ();
            
        }
        
        @snp_loci =  sort { $a->[0] cmp $b->[0] or $a->[1] <=> $b->[1] } @snp_loci;
        
        for (my $i=0;$i<scalar(@snp_loci);$i++) {
            
            print $summary "$snp_loci[$i][0]\t$snp_loci[$i][1]";
            if("bp" eq lc($windowunit))
            {
                print $summary "\t$snp_loci[$i][2]\t$snp_loci[$i][3]\t$snp_loci[$i][4]\n";
            }
            elsif("snp"eq lc($windowunit))
            {
                print $summary "\t$snp_loci[$i][2]\t$snp_loci[$i][4]\n";
            }
            
        }
    }
    
    sub write_pvalue_in_wiggle_format {
        
        open my $summary, "<$summary_pvalue_output" or die "Could not open output file";
        open my $wiggle, ">$summary_pvalue_wiggle_output" or die "Could not open output file";
        
        #calculate the offset; wiggle is 1-based
        my $offset=int($windowsize/2)-1;

        my $ucfilt;
        $ucfilt= {map{ $_,1} split /\s/,$ucscfilter} if $ucscfilter;
        
        #print track
        print $wiggle "track type=wiggle_0 name=\"$trackname\" visibility=full\n";
        
        my $activechr="";
        my $chrHash={};
        my ($chr,$pos,$average_pvalue,$product_pvalue,$snpcount) = ("",0,"","",0);
        my ($print_pvalue) = ("");
        while(my $l=<$summary>)
        {
            chomp $l;
            my($chr,$pos,undef,undef,$measure)=split /\t/,$l;
            
            if("bp" eq lc($windowunit))
            {
                ($chr,$pos,$average_pvalue,$product_pvalue,$snpcount)=split /\t/,$l;
                $print_pvalue = "$average_pvalue\t$product_pvalue";
            }
            elsif("snp"eq lc($windowunit))
            {
                ($chr,$pos,$average_pvalue,$snpcount)=split /\t/,$l;
                $print_pvalue = "$average_pvalue";
            }
            
            if($ucfilt)
            {
                next unless exists($ucfilt->{$chr});
            }
            if($chr ne $activechr)
            {
                my $printchr=$ucscprepend.$chr;
                print $wiggle "variableStep chrom=$printchr span=$windowsize\n";
                die "Chromosome $chr occured several times; Please make sure the file is properly sorted" if exists($chrHash->{$chr});
                $chrHash->{$chr}=1;
                $activechr=$chr;
            }
            
            next if $average_pvalue eq "na";
            $pos-=$offset;
            print $wiggle "$pos\t$print_pvalue\n";   
        }

    
    }
    
    sub snpCount_in_file {
        my $file=shift;
        my $snp_count_infile=0;
        
        open (IN,"<$file") or die "Could not open input file $file\n";
        while (<IN>)
        {
            my $line=$_;
            chomp($line);
            
            # discard blank line
            if ($line =~ m/^\s*$/g) {
                next;
            }
            else {
                $snp_count_infile++;
            }
        }
        close IN;
        
        return $snp_count_infile;
    }
    sub Parse_R_output {
        my $rbatchfile=shift;
        my $routputfile=shift;
        my $windowunit=shift;
        my $ofh1=shift;
        my $ofh2=shift;
        my $minpvalue=shift;
        my $min_snp_fraction=shift;
        my $snp_loci=shift;
        my %snp_loci = %{$snp_loci};
        
        open my $fh,"<$routputfile" or die "Could not open file handle";
        
        open (IN,"<$rbatchfile") or die "Could not open input file";

        while (my $line = <IN>) {
            chomp $line;
            
            if ($line =~ m/^\#\#/i) {
                
                $line =~ s/\#//gi;
                #my ($chr,$pos,$refnuc,$snpcount) = split("\t",$line);
                my @a = split("\t",$line);
                my $chr = shift @a;
                my $pos = shift @a;
                my $refnuc = shift @a;
                my $snpcount = shift @a;
                my $ns = shift @a;
                
                my $string = Utility::get_snp_frequency_string(\@a);
                #print "BEFORE: $chr\t$pos\t$refnuc\t\t$snpcount\n";
                my ($average_pvalue,$product_pvalue,$fh) = Utility::get_allpvalues($snpcount,$fh);
                if ("bp" eq lc($windowunit)) {
                    
                    my $min_snp_count = 0;
                    $min_snp_count = int($snpcount*$min_snp_fraction);
                    if ($min_snp_count==0) {
                        $min_snp_count=1;
                    }
                    else {
                        $min_snp_count = $min_snp_count;
                    }
                    
                    my $minpvalue1 = 1;
                    for my $i (1..$min_snp_count) {
                        $minpvalue1 = $minpvalue1*$minpvalue;
                    }
                    
                    if ($minpvalue1 > 0) {
                        $minpvalue1 = log($minpvalue1);
                    }
                    else {
                        $minpvalue1 = $minpvalue1;
                    }
            
                    if ($average_pvalue > 0) {
                        $average_pvalue = log($average_pvalue);
                    }
                    else {
                        $average_pvalue = $average_pvalue;
                    }
                    
                    if (($average_pvalue<=$minpvalue1) or ($product_pvalue<=$minpvalue1)) {
                        print $ofh1 "$chr\t$pos\t$refnuc\t$average_pvalue\t$product_pvalue\t$snpcount\t$minpvalue1\t$ns\t$min_snp_fraction\t$string\n";
                        print "$chr\t$pos\t$refnuc\t$average_pvalue\t$product_pvalue\t$snpcount\t$minpvalue1\t$ns\t$min_snp_fraction\t$string\n";
                        if(exists($snp_loci{"$chr\t$pos"})) {
                            $snp_loci{"$chr\t$pos"} = "$average_pvalue\t$product_pvalue\t$snpcount";
                        }
                    }
                    else {
                        print $ofh2 "$chr\t$pos\t$refnuc\t$average_pvalue\t$product_pvalue\t$snpcount\t$minpvalue1\t$ns\t$min_snp_fraction\t$string\n";
                        
                        if(exists($snp_loci{"$chr\t$pos"})) {
                            $snp_loci{"$chr\t$pos"} = "$average_pvalue\t$product_pvalue\t$snpcount";
                        }
                    }
                }
                elsif("snp"eq lc($windowunit)) {
                     if ($average_pvalue<=$minpvalue) {
                        print $ofh1 "$chr\t$pos\t$refnuc\t$average_pvalue\t$snpcount\t$minpvalue\t$ns\t$string\n";
                        print "$chr\t$pos\t$refnuc\t$average_pvalue\t$snpcount\t$minpvalue\t$ns\t$string\n";
                        
                        if(exists($snp_loci{"$chr\t$pos"})) {
                            $snp_loci{"$chr\t$pos"} = "$average_pvalue\t$snpcount";
                        }
                     }
                     else {
                        print $ofh2 "$chr\t$pos\t$refnuc\t$average_pvalue\t$snpcount\t$minpvalue\t$ns\t$string\n";
                        
                        if(exists($snp_loci{"$chr\t$pos"})) {
                            $snp_loci{"$chr\t$pos"} = "$average_pvalue\t$snpcount";
                        }
                        #else {
                            #$snp_loci{"$chr\t$pos"} = "1\t0";
                            #print "$chr\t$pos\t$refnuc\t1\t$snpcount\t$minpvalue\t0\t$string\n";
                        #}
                     }
                }
                #print "AFTER: $chr\t$pos\t$refnuc\t$average_pvalue\t$snpcount\n";
                #print "$chr\t$pos\t$refnuc\t$average_pvalue\t$snpcount\n";
        
            }
        }
        
        return %snp_loci;
    
    
    }
    
    sub get_snp_frequency_string {
        my $ref = shift;
        my @snp_frequency = @{$ref};
        
        my $popCount = shift @snp_frequency;
        my $string="";
        
        #print "SNP:\t@snp_frequency\n";
        for my $snp (@snp_frequency) {
            #print "SNP:\t$popCount\t$snp\n";
            my @freq = split('\,',$snp);
            $string .= $snp;
            $string .= "|";
        }
        
        $string =~ s/\|$//;
        return $string;
    }
    
    sub get_allpvalues {
    
            my $snpcount = shift;
            my $fh = shift;
            
            my $ct = 0;
            my $total_pvalue = 0;
            my $product_pvalue = 1;
            while(1) {
               my ($l,$fh) = Utility::read_rfile($fh);
               exit if $l eq "lastline";
               if ($l =~m/p-value\s+\=/i) {
                    $ct++;
                    chomp($l);
                    
                    my @val = split("=",$l);
                                        
                    my $pvalue=$val[1];
                    chomp($pvalue);
                    if($pvalue =~ m/^\s+/g) {
                        $pvalue =~ s/^\s+//g;
                    }
                                
                    if($pvalue =~ m/\s+$/g) {
                        $pvalue =~ s/\s+$//g;
                    }
                    $total_pvalue += $pvalue;
                    $product_pvalue *= $pvalue;
                    #print "$ct\t$l\t$pvalue\n";
               }
               elsif ($l =~m/p-value\s+\</i) {
                    $ct++;
                    chomp($l);
                    
                    my @val = split('\<',$l);
                                        
                    my $pvalue=$val[1];
                    chomp($pvalue);
                    if($pvalue =~ m/^\s+/g) {
                        $pvalue =~ s/^\s+//g;
                    }
                                
                    if($pvalue =~ m/\s+$/g) {
                        $pvalue =~ s/\s+$//g;
                    }
                    $total_pvalue += $pvalue;
                    $product_pvalue *= $pvalue;
                    #print "$ct\t$l\t$pvalue\n";
               }
               elsif ($l =~m/p-value\s+\>/i) {
                    $ct++;
                    chomp($l);
                    
                    my @val = split('\>',$l);
                                        
                    my $pvalue=$val[1];
                    chomp($pvalue);
                    if($pvalue =~ m/^\s+/g) {
                        $pvalue =~ s/^\s+//g;
                    }
                                
                    if($pvalue =~ m/\s+$/g) {
                        $pvalue =~ s/\s+$//g;
                    }
                    $total_pvalue += $pvalue;
                    $product_pvalue *= $pvalue;
                    #print "$ct\t$l\t$pvalue\n";
               }
               
               last if ($ct==$snpcount);
               last if ($l eq "");
            }
            
            my $average_pvalue = $total_pvalue/$snpcount;
            #$product_pvalue = log($product_pvalue);
            if ($product_pvalue > 0) {
                        $product_pvalue = log($product_pvalue);
            }
            else {
                $product_pvalue = $product_pvalue;
            }
                    
            #print "total pvalue: $total_pvalue\t$average_pvalue\n";
            return ($average_pvalue,$product_pvalue,$fh);
        }

        sub read_rfile {
            my $fh = shift;
            my $line = "";

            if (defined ($line = <$fh>)) {
                $line = $line;
            }
            elsif (eof){ $line = "lastline";}
            else {
                    $line = "";
            }
            return ($line,$fh);
        }
    
    
    sub calculatePvalues {
        
        my $data=shift;
        my $popcount=shift;
        my $pos=shift;
        my $chr=shift;
        my $orfh=shift;
        
        my $snpcount=@$data;
        my $refchar = "";
        my $total_pvalue = 0;
        my $average_pvalue = 0;
        my $result_string = "";
        #$refchar = $$data[0]->{refchar};
        
        $popcount = scalar(@population);
        #my %snp_loci = ();
        
        my $samp=[];
        
        my $ct=0;
        
        my @results = ();
        my @vals = ();
        my $ns = 0;
        foreach my $d (@$data)
        {
            next unless $d->{ispuresnp};
            $ct++;
            $refchar = $d->{refchar};
            $ns = $d->{nosim};
            if ($ct==1) {
                #print $orfh "##$chr\t$pos\t$refchar\t$snpcount\n";
                $vals[0] = $chr;
                my $ns = $d->{nosim};
                
                push (@results,"$chr\t$pos\t$refchar\t$snpcount\t$ns\t$popcount");
                #print "$chr\t$pos\t$refchar\t$snpcount\n";
            }
            
            #my $value = "";
            my $pvalue = 2;
            #my $flag = 1;
            my $cosa=$d->{samples}; #countsamples
            my @table = ();
            my @vals = ();
            #print "$d->{ispuresnp}\t$popcount\t$d->{nosim}\t$d->{file}\t$pos\t$d->{refchar}\t$snpcount\n";
            
            my $freqsa=[];
            for(my $i=0; $i<$popcount; $i++)
            #foreach my $p(@population)
            {
                #my $i=$p-1;
                my $t=$cosa->[$i];
                my $cov=$t->{eucov};
                #print "Cov:$t->{eucov}\tA:$t->{A}\tT:$t->{T}\tC:$t->{C}\tG:$t->{G}\n";
                #$value.="$t->{A},$t->{T},$t->{C},$t->{G},";
                $vals[0] = $t->{A};
                $vals[1] = $t->{T};
                $vals[2] = $t->{C};
                $vals[3] = $t->{G};
                push @table, [@vals];
                @vals=();
                
            }
            
            
            my $value = "";
             
            for(my $i=0;$i<4;$i++) {
                for(my $j=0;$j<scalar(@table);$j++) {
                    $value.="$table[$j][$i], ";
                }
            }
            
               
            $value =~ s/\,\s{1}$//;

            push (@results,$value);
            # Writing R-Batch file
            #open (RBATCH,">$rbatchfile") || die ("could not open $rbatchfile for writing $!\n\n");
            #print RBATCH "a<-(fisher.test(matrix(c($value),nrow=$no_population),alternative=\"two.sided\",simulate.p.value=T,B=$no_simulation))\na";
            #print $orfh "a<-(fisher.test(matrix(c($value), nrow = $popcount),alternative=\"two.sided\",simulate.p.value=T,B=$d->{nosim},conf.int=T,conf.level=0.95))\na\n";
            #print $orfh "fisher.test(matrix(c($value), nrow = $popcount),alternative=\"two.sided\",simulate.p.value=T,B=$d->{nosim},conf.int=T,conf.level=0.95)\n";

            @table = ();
        }
        $ct=0;
        

        if (scalar(@results)>1) {
            my $firstline = shift(@results);
            my ($chr1,$pos1,$refchar1,$snpcount1,$ns1,$popcount1) = split("\t",$firstline);
            $snp_loci{"$chr1\t$pos1"} = $pos1;
            my $value_string = join("\t", @results);
            print $orfh "##$firstline\t$value_string\n";
            print "$firstline\n";
            for (my $i=0;$i<scalar(@results);$i++) {
                
                print $orfh "fisher.test(matrix(c($results[$i]), nrow = $popcount1),alternative=\"two.sided\",simulate.p.value=T,B=$ns1,conf.int=T,conf.level=0.95)\n";
            }
        }
        
        @results = ();
       
       #return %snp_loci;

    }
    
    sub _annotateWindow
    {
        my $curwin=shift;
        my $chr=shift;
        my $start=shift;
        my $end=shift;
        my $window=shift;

        my $snps=0;
        my $aboveCoverage=0;
        foreach(@$curwin)
        {
            $snps++ if $_->{ispuresnp};
            $aboveCoverage++ if $_->{iscov};
        }

        return
        {
            chr=>$chr,
            start=>$start,
            end=>$end,
            middle=>int(($end+1+$start)/2),
            countpuresnp=>$snps,
            count_covered=>$aboveCoverage,
            window=>$window,
            data=>$curwin      
        };
    }
    
    
    
    sub _parseLightwight
    {
        my $line=shift;
        my $mincount=shift;
        my $mincov=shift;
        my $maxcov=shift;
        my $nosimulations=shift,
        
        chomp $line;
        my @a=split /\s+/,$line;
        my $chr=shift @a;
        my $pos=shift @a;
        my $rc=shift @a;
        
        my @samp;
        my $taintedsnp=0;
        #for(my $i=0; $i<@a; $i++)
        foreach my $p(@population)
        {
            my $i=$p-1;
            my $col=$a[$i];
            my $e;
            
            if($col=~/-/)
            {
                $e={index=>$i,eucov=>0,totcov=>0,A=>0,T=>0,C=>0,G=>0,N=>0,del=>0};
            }
            else
            {
                my @parts=split /:/,$col;
                
                die "failed parsing $col; does not have the correct number of entries" unless @parts ==6;
                
                $e = {
                     A=>$parts[0],
                     T=>$parts[1],
                     C=>$parts[2],
                     G=>$parts[3],
                     N=>$parts[4],
                     del=>$parts[5]
                    };
                
                #$taintedsnp=1 if $e->{del} > 0;
                $e->{eucov}  = ($e->{A}+$e->{T}+$e->{C}+$e->{G});
                $e->{totcov} = ($e->{eucov}+$e->{N}+$e->{del});
                $e->{index}  = $i;
            }
            push @samp,$e;
        }
=cut check if the position is a snp in any population (privious strategy)     changed on 16-09-2010   
        # check if the position is a snp
        my $issnp=0;
        my $chash={};
        my $is_suficient_covered=1;
        
        foreach(@samp)
        {
            $is_suficient_covered=0 if($_->{eucov}<$mincov || $_->{eucov}>=$maxcov);
            
            $chash->{A}++ if  $_->{A} >= $mincount;
            $chash->{T}++ if  $_->{T} >= $mincount;
            $chash->{C}++ if  $_->{C} >= $mincount;
            $chash->{G}++ if  $_->{G} >= $mincount;
        }
        
        $issnp=1 if values(%$chash)>1;          # a SNP is when there are differences in any data file.
        $issnp=0 unless $is_suficient_covered; # no SNP will ever be allowed at a position which is not sufficiently covered in all data files!!
=cut check if the position is a snp in any population (privious strategy)     changed on 16-09-2010
        
        # New SNP identification strategy (to check minor allele among all population)
        
        # check if the position is a snp
        my $issnp=0;
        my ($ca,$ct,$cc,$cg,$cn,$cdel)=(0,0,0,0,0,0);
        my $is_suficient_covered=1;
        foreach(@samp)
        {
            $is_suficient_covered=0 if($_->{eucov}<$mincov || $_->{eucov}>=$maxcov);
            $ca+= $_->{A};
            $ct+= $_->{T};
            $cc+= $_->{C};
            $cg+= $_->{G};
            $cn+= $_->{N};
            $cdel+= $_->{del};
        }
        
        my $allelecount=0;
        $allelecount++ if $ca>=$mincount;
        $allelecount++ if $ct>=$mincount;
        $allelecount++ if $cc>=$mincount;
        $allelecount++ if $cg>=$mincount;
        $issnp=1 if $allelecount>1;          # a SNP is when there are differences in any data file.
        
        # unset the snp if not sufficiently covered
        $issnp=0 unless $is_suficient_covered; # no SNP will ever be allowed at a position which is not sufficiently covered in all data files!!
        
        # the SNP is tainted if there are any deletions at the given position
        $taintedsnp=$cdel>0?1:0;
        
        my $en={
            chr=>$chr,
            pos=>$pos,
            refchar=>$rc,
            issnp=>$issnp,
            iscov=>$is_suficient_covered,
            nosim=>$nosimulations,
            samples=>\@samp
        };
        $en->{ispuresnp} = ($en->{issnp} and not $taintedsnp)?1:0;
        
        return $en;
    }
    
    sub runR {
        
        my $rbatchfile=shift;
        my $routputfile=shift;
        #my $rbatchfile1 = "/Volumes/Temp/popgen-tools/24-06-2010/test-16-07-2010bp2.r";
        #my $routputfile= "/Volumes/Temp/popgen-tools/24-06-2010/rrrrrrr.txt";
        print "Running R to calculate p-value Started at: \t". localtime() .".....\n";
            system("`R --vanilla --slave < $rbatchfile > $routputfile`");
        print "Running R to calculate p-value Completed at: \t". localtime() .".....\n";

    }
    
}





{
    package NextTire;
    use strict;
    use warnings;
    use List::Util qw[min max];

    
    sub calculatePvalues {
        my $output=shift;
        my $rbatchfile=shift;
        my $nosimulation=shift;
        my $orfh=shift;
        
        my @results = ();
        my $rline = "";
        
        open my $fh,"<$rbatchfile" or die "Could not open file handle";
        
        open (IN,"<$output") or die "Could not open input file";
        while (<IN>)
        {
            my $line=$_;
            chomp($line);
            
            # discard blank line
            if ($line =~ m/^\s*$/g) {
                next;
            }
            else {
                my @record = split("\t",$line);
                my $chr = shift @record;
                my $pos = shift @record;
                my $renuc = shift @record;
                
                ($rline,$fh) = NextTire::read_privious_r_batch_file($chr,$pos,$renuc,$fh);
                if ($rline ne "") {
                    ##2L        5500    N       38      10000   2
                    @results = split("\t",$rline);
                    my $chr1 = shift @results;
                    my $pos1 = shift @results;
                    my $refchar1 = shift @results;
                    my $snpcount1 = shift @results;
                    my $ns1 = shift @results;
                    my $popcount1 = shift @results;

                    if (scalar(@results)>0) {
                        
                        my $value_string = join("\t", @results);
                        
                        print $orfh "##$chr1\t$pos1\t$refchar1\t$snpcount1\t$nosimulation\t$popcount1\t$value_string\n";
                        print "$chr1\t$pos1\t$refchar1\t$snpcount1\t$nosimulation\t$popcount1\n";
                        for (my $i=0;$i<scalar(@results);$i++) {
                           
                            print $orfh "fisher.test(matrix(c($results[$i]), nrow = $popcount1),alternative=\"two.sided\",simulate.p.value=T,B=$nosimulation,conf.int=T,conf.level=0.95)\n";
                        }
                    }
                }
                
                @results = ();
                
                
                
            }
              
        }
        
        close IN;

    }
    
    
    sub read_privious_r_batch_file {
        my $chr=shift;
        my $pos=shift;
        my $renuc=shift;
        my $fh=shift;
        my $line = "";
        
        
        while (my $line = <$fh>) {
            chomp $line;
            
            if ($line =~ m/^\#\#/i) {
                
                $line =~ s/\#//gi;
                my @record = split("\t",$line);
                my $chr1 = shift @record;
                my $pos1 = shift @record;
                my $renuc1 = shift @record;
                if (($chr1 eq $chr) and ($pos1==$pos) and ($renuc1 eq $renuc)) {
                    return ($line,$fh);
                    exit;
                }       
            }
        }
    }
    
} # package NextTier finished


{
    use warnings;
    use strict;
    package PopCount;

    sub new
    {
        my $class=shift;
        my $file=shift;
        
        open my $fh1,"<$file" or die "Could not open file handle";
        
        return bless {
            fh=>$fh1,
            curwin=>[],
            buffer=>[]
        },__PACKAGE__;
        
        
    }
    
        sub count_samples
        {
            my $self=shift;
            my $l=$self->_nextline();
            my $p=_parseLightwight($l);
            my $c=scalar(@{$p->{samples}});
            $self->_bufferline($l);
           
            return $c;
        }
        
    sub _nextline
    {
        my $self=shift;
        my $fh=$self->{fh};
        my $buffer=$self->{buffer};
        
        return shift @$buffer if @$buffer;
        return <$fh>;
    }
    
    sub _bufferline
    {
        my $self=shift;
        my $line=shift;
        push @{$self->{buffer}},$line;
    }
    
    
    sub _parseLightwight
    {
        my $line=shift;

        chomp $line;
        my @a=split /\s+/,$line;
        my $chr=shift @a;
        my $pos=shift @a;
        my $rc=shift @a;
        
        my @samp;

        for(my $i=0; $i<@a; $i++)
        {
            my $col=$a[$i];
            my $e;
            
            if($col=~/-/)
            {
                $e={index=>$i,eucov=>0,totcov=>0,A=>0,T=>0,C=>0,G=>0,N=>0,del=>0};
            }
            else
            {
                my @parts=split /:/,$col;
                
                die "failed parsing $col; does not have the correct number of entries" unless @parts ==6;
                
                $e = {
                     A=>$parts[0],
                     T=>$parts[1],
                     C=>$parts[2],
                     G=>$parts[3],
                     N=>$parts[4],
                     del=>$parts[5]
                    };
                
                #$taintedsnp=1 if $e->{del} > 0;
                $e->{eucov}  = ($e->{A}+$e->{T}+$e->{C}+$e->{G});
                $e->{totcov} = ($e->{eucov}+$e->{N}+$e->{del});
                $e->{index}  = $i;
            }
            push @samp,$e;
        }
        
        my $en={
            chr=>$chr,
            pos=>$pos,
            refchar=>$rc,
            samples=>\@samp
        };
       
        return $en;
    }
    
}  # PopCount package finished

{
    package Test;
    use strict;
    use warnings;
    #use Test::More;
    #use Test::Simple;
    
    my $testcounter;
    
    sub _getBPSliderForString
    {
        my $str=shift;
        my $window=shift;
        my $step=shift;
        my $mincount=shift;
        my $mincov=shift;
        my $maxcov=shift;
        my $nosimulations=shift;
        open my $ofh,"<",\$str or die "could not open string filehandle";
        my $cr=bless {   
            lower=>0,
            upper=>$window,
            window=>$window,
            step=>$step,
            fh=>$ofh,
            mincount=>$mincount,
            mincov=>$mincov,
            maxcov=>$maxcov,
            nosim=>$nosimulations,
            curwin=>[],
            buffer=>[]
            
        },"BpSlider";
        return $cr;
    }
    
    sub runTests
    {
        $testcounter=1;
        testBpSlider();
        test_calculatePvalues();
        #test_callR();
        
        exit;
    }
    
    sub testArrangeSNPfrequenciesInColumns
    {

    }
    
    
    sub testBpSlider
    {
        my $teststr=
        "chr1\t1\tN\t0:0:6:0:0:0\t0:0:0:7:0:0\t0:0:8:0:0:0\n".
        "chr1\t2\tN\t0:4:8:0:0:0\t0:0:8:0:0:0\t0:0:8:0:0:0\n".
        "chr1\t3\tN\t0:0:3:0:0:0\t0:0:0:8:0:0\t0:0:13:0:0:0\n".
        "chr1\t4\tN\t0:0:5:0:0:0\t0:0:4:0:0:0\t0:0:6:0:0:0\n".
        "chr1\t5\tN\t0:0:8:0:0:0\t0:0:0:8:0:0\t0:0:8:0:0:0\n".
        "chr2\t3\tN\t0:0:11:0:0:0\t0:0:0:8:0:0\t0:0:20:0:0:0\n".
        "chr3\t1\tN\t0:0:8:0:0:0\t0:0:8:0:0:0\t0:0:8:0:0:0\n".
        "chr3\t2\tN\t0:0:1:0:0:0\t0:0:0:5:0:0\t0:0:8:0:0:0\n".
        "chr4\t4\tN\t0:0:8:0:0:0\t0:0:8:0:0:0\t0:0:8:0:0:0\n";
        my $bpsl=_getBPSliderForString($teststr,3,1,2,4,1000,100);
        
        my $w=$bpsl->nextWindow();
        
        is($w->{chr},"chr1","test BpSlider, correct chromosome");
        is($w->{count_covered},2,"test BpSlider, correct number of sufficiently covered regions");
        is($w->{countpuresnp},2,"test BpSlider, correct number of snps in region");
        is($w->{start},0,"test BpSlider, correct start position");
        is($w->{end},3,"test BpSlider, correct end position");
        is($w->{window},3,"test BpSlider, correct window length");
        is(scalar(@{$w->{data}}),3,"Correct number of data entries");
        is($w->{data}[0]{ispuresnp},1,"test BpSlider, correct identification of a pure SNP at the first position");
        is($w->{data}[1]{ispuresnp},1,"test BpSlider, correct identification of a pure SNP at the second position");
        is($w->{data}[2]{ispuresnp},0,"test BpSlider, correct, there is no pure SNP at position two");
        is($w->{data}[0]{issnp},1,"test BpSlider, correct identification of a SNP at the first position");
        is($w->{data}[1]{issnp},1,"test BpSlider, correct identification of a SNP at the second position");
        is($w->{data}[2]{issnp},0,"test BpSlider, correct, there is no SNP at position two");
        is($w->{data}[0]{iscov},1,"test BpSlider, first position is sufficiently covered");
        is($w->{data}[1]{iscov},1,"test BpSlider, second position is sufficiently covered for a SNP");
        is($w->{data}[2]{iscov},0,"test BpSlider, third position is not sufficiently covered for a SNP -> no SNP possible at this position");
        is($w->{data}[0]{pos},1,"test BpSlider, position is correct");
        is($w->{data}[1]{pos},2,"test BpSlider, position is correct");
        is($w->{data}[2]{pos},3,"test BpSlider, position is correct");
        
        is($w->{data}[0]{samples}[0]{eucov},6,"test BpSlider, coverage is correct");
        is($w->{data}[0]{samples}[1]{eucov},7,"test BpSlider, coverage is correct");
        is($w->{data}[0]{samples}[2]{eucov},8,"test BpSlider, coverage is correct");
        is($w->{data}[0]{samples}[0]{totcov},6,"test BpSlider, coverage is correct");
        is($w->{data}[0]{samples}[1]{totcov},7,"test BpSlider, coverage is correct");
        is($w->{data}[0]{samples}[2]{totcov},8,"test BpSlider, coverage is correct");
        is($w->{data}[0]{samples}[0]{index},0,"test BpSlider, index is correct");
        is($w->{data}[0]{samples}[1]{index},1,"test BpSlider, index is correct");
        is($w->{data}[0]{samples}[2]{index},2,"test BpSlider, index is correct");

        is($w->{data}[0]{samples}[0]{A},0,"test BpSlider, count of A is correct");
        is($w->{data}[0]{samples}[0]{T},0,"test BpSlider, count of T is correct");
        is($w->{data}[0]{samples}[0]{C},6,"test BpSlider, count of C is correct");
        is($w->{data}[0]{samples}[0]{G},0,"test BpSlider, count of G is correct");
        
        is($w->{data}[1]{samples}[0]{eucov},12,"test BpSlider, coverage is ok");
        is($w->{data}[1]{samples}[0]{A},0,"test BpSlider, count of A is correct");
        is($w->{data}[1]{samples}[0]{T},4,"test BpSlider, count of T is correct");
        is($w->{data}[1]{samples}[0]{C},8,"test BpSlider, count of C is correct");
        is($w->{data}[1]{samples}[0]{G},0,"test BpSlider, count of G is correct");
        
        is($w->{data}[2]{samples}[2]{eucov},13,"test BpSlider, coverage is ok");
        is($w->{data}[2]{samples}[2]{A},0,"test BpSlider, count of A is correct");
        is($w->{data}[2]{samples}[2]{T},0,"test BpSlider, count of T is correct");
        is($w->{data}[2]{samples}[2]{C},13,"test BpSlider, count of C is correct");
        is($w->{data}[2]{samples}[2]{G},0,"test BpSlider, count of G is correct");
        
        
        $w=$bpsl->nextWindow();
        is($w->{chr},"chr1","test BpSlider, correct chromosome");
        is($w->{count_covered},2,"test BpSlider, correct number of sufficiently covered regions");
        is($w->{countpuresnp},1,"test BpSlider, correct number of snps in region");
        is($w->{start},1,"test BpSlider, correct start position");
        is($w->{end},4,"test BpSlider, correct end position");
        is($w->{middle},3,"test BpSlider, correct end position");
        is($w->{window},3,"test BpSlider, correct window length");
        is(scalar(@{$w->{data}}),3,"Correct number of data entries");
        is($w->{data}[0]{pos},2,"test BpSlider, position is correct");
        is($w->{data}[1]{pos},3,"test BpSlider, position is correct");
        is($w->{data}[2]{pos},4,"test BpSlider, position is correct");
        is($w->{data}[0]{ispuresnp},1,"test BpSlider, correct identification of a SNP at the first position");
        is($w->{data}[1]{ispuresnp},0,"test BpSlider, correct there is no SNP at the second position");
        is($w->{data}[2]{ispuresnp},0,"test BpSlider, correct, there is no SNP at position three");
        is($w->{data}[0]{iscov},1,"test BpSlider, first position is sufficiently covered");
        is($w->{data}[1]{iscov},0,"test BpSlider, second position is not sufficiently covered for a SNP");
        is($w->{data}[2]{iscov},1,"test BpSlider, third position is sufficiently covered for a SNP");
        
        is($w->{data}[0]{samples}[0]{eucov},12,"test BpSlider, coverage is ok");
        is($w->{data}[0]{samples}[0]{A},0,"test BpSlider, count of A is correct");
        is($w->{data}[0]{samples}[0]{T},4,"test BpSlider, count of T is correct");
        is($w->{data}[0]{samples}[0]{C},8,"test BpSlider, count of C is correct");
        is($w->{data}[0]{samples}[0]{G},0,"test BpSlider, count of G is correct");
        
        is($w->{data}[1]{samples}[2]{eucov},13,"test BpSlider, coverage is ok");
        is($w->{data}[1]{samples}[2]{A},0,"test BpSlider, count of A is correct");
        is($w->{data}[1]{samples}[2]{T},0,"test BpSlider, count of T is correct");
        is($w->{data}[1]{samples}[2]{C},13,"test BpSlider, count of C is correct");
        is($w->{data}[1]{samples}[2]{G},0,"test BpSlider, count of G is correct");
        
        $w=$bpsl->nextWindow();
        is($w->{chr},"chr1","test BpSlider, correct chromosome");
        is($w->{count_covered},2,"test BpSlider, correct number of sufficiently covered regions");
        is($w->{countpuresnp},1,"test BpSlider, correct number of snps in region");
        is($w->{start},2,"test BpSlider, correct start position");
        is($w->{end},5,"test BpSlider, correct end position");
        is($w->{middle},4,"test BpSlider, correct end position");
        is($w->{window},3,"test BpSlider, correct window length");
        is(scalar(@{$w->{data}}),3,"Correct number of data entries");
        is($w->{data}[0]{pos},3,"test BpSlider, position is correct");
        is($w->{data}[1]{pos},4,"test BpSlider, position is correct");
        is($w->{data}[2]{pos},5,"test BpSlider, position is correct");
        
        # switch to new chromosome
        $w=$bpsl->nextWindow();
        is($w->{chr},"chr2","test BpSlider, correct chromosome");
        is($w->{count_covered},1,"test BpSlider, correct number of sufficiently covered regions");
        is($w->{countpuresnp},1,"test BpSlider, correct number of snps in region");
        is($w->{start},0,"test BpSlider, correct start position");
        is($w->{end},3,"test BpSlider, correct end position");
        is($w->{middle},2,"test BpSlider, correct end position");
        is($w->{window},3,"test BpSlider, correct window length");
        
        is(scalar(@{$w->{data}}),1,"Correct number of data entries");
        is($w->{data}[0]{pos},3,"test BpSlider, position is correct");
        is($w->{data}[0]{iscov},1,"test BpSlider, correct position is sufficiently covered");
        is($w->{data}[0]{ispuresnp},1,"test BpSlider, correct position is a snp");
        is($w->{data}[0]{samples}[0]{eucov},11,"test BpSlider, coverage is ok");
        is($w->{data}[0]{samples}[0]{A},0,"test BpSlider, count of A is correct");
        is($w->{data}[0]{samples}[0]{T},0,"test BpSlider, count of T is correct");
        is($w->{data}[0]{samples}[0]{C},11,"test BpSlider, count of C is correct");
        is($w->{data}[0]{samples}[0]{G},0,"test BpSlider, count of G is correct");
        is($w->{data}[0]{samples}[2]{eucov},20,"test BpSlider, coverage is ok");
        is($w->{data}[0]{samples}[2]{A},0,"test BpSlider, count of A is correct");
        is($w->{data}[0]{samples}[2]{T},0,"test BpSlider, count of T is correct");
        is($w->{data}[0]{samples}[2]{C},20,"test BpSlider, count of C is correct");
        is($w->{data}[0]{samples}[2]{G},0,"test BpSlider, count of G is correct");
        
        $w=$bpsl->nextWindow();
        is($w->{chr},"chr3","test BpSlider, correct chromosome");
        is($w->{count_covered},1,"test BpSlider, correct number of sufficiently covered regions");
        is($w->{countpuresnp},0,"test BpSlider, correct number of snps in region");
        is($w->{start},0,"test BpSlider, correct start position");
        is($w->{end},3,"test BpSlider, correct end position");
        is($w->{middle},2,"test BpSlider, correct end position");
        is($w->{window},3,"test BpSlider, correct window length");
        is(scalar(@{$w->{data}}),2,"Correct number of data entries");
        is($w->{data}[0]{pos},1,"test BpSlider, position is correct");
        is($w->{data}[1]{pos},2,"test BpSlider, position is correct");
        
        $w=$bpsl->nextWindow();
        is($w->{chr},"chr4","test BpSlider, correct chromosome");
        is($w->{count_covered},0,"test BpSlider, correct number of sufficiently covered regions");
        is($w->{countpuresnp},0,"test BpSlider, correct number of snps in region");
        is($w->{start},0,"test BpSlider, correct start position");
        is($w->{end},3,"test BpSlider, correct end position");
        is($w->{middle},2,"test BpSlider, correct end position");
        is($w->{window},3,"test BpSlider, correct window length");
        is(scalar(@{$w->{data}}),0,"Correct number of data entries");
        
        $w=$bpsl->nextWindow();
        is($w->{chr},"chr4","test BpSlider, correct chromosome");
        is($w->{count_covered},1,"test BpSlider, correct number of sufficiently covered regions");
        is($w->{countpuresnp},0,"test BpSlider, correct number of snps in region");
        is($w->{start},1,"test BpSlider, correct start position");
        is($w->{end},4,"test BpSlider, correct end position");
        is($w->{middle},3,"test BpSlider, correct end position");
        is($w->{window},3,"test BpSlider, correct window length");
        is(scalar(@{$w->{data}}),1,"Correct number of data entries");
        
        $w=$bpsl->nextWindow();
        not_exists($w,"correct end of file reached");
        $w=$bpsl->nextWindow();
        not_exists($w,"correct end of file reached");
        
        
        # weired characters
        $teststr=
        "chr1\t1\tN\t0:0:6:0:1:0\t0:0:0:7:3:0\t0:0:8:0:2:0\n".
        "chr1\t2\tN\t0:4:8:0:0:0\t0:0:8:0:0:0\t0:0:8:0:0:1\n".
        "chr1\t3\tN\t0:0:3:0:0:0\t0:0:0:8:0:0\t0:0:13:0:2:2\n".
        "chr1\t4\tN\t0:4:8:0:0:0\t-\t0:0:8:0:2:1\n";
        $bpsl=_getBPSliderForString($teststr,3,1,2,4,1000,100);
        
        $w=$bpsl->nextWindow();
        
        is($w->{chr},"chr1","test BpSlider, correct chromosome");
        is($w->{count_covered},2,"test BpSlider, correct number of sufficiently covered regions");
        is($w->{countpuresnp},1,"test BpSlider, correct number of snps in region");
        is($w->{start},0,"test BpSlider, correct start position");
        is($w->{end},3,"test BpSlider, correct end position");
        is($w->{window},3,"test BpSlider, correct window length");
        is($w->{data}[0]{issnp},1,"test BpSlider, snp correct");
        is($w->{data}[0]{ispuresnp},1,"test BpSlider, snp correct");
        is($w->{data}[0]{iscov},1,"test BpSlider, snp correct");
        is($w->{data}[1]{issnp},1,"test BpSlider, snp correct");
        is($w->{data}[1]{ispuresnp},0,"test BpSlider, snp correct");
        is($w->{data}[1]{iscov},1,"test BpSlider, snp correct");
        is($w->{data}[2]{issnp},0,"test BpSlider, snp correct");
        is($w->{data}[2]{ispuresnp},0,"test BpSlider, snp correct");
        is($w->{data}[2]{iscov},0,"test BpSlider, snp correct");
        
        is($w->{data}[0]{samples}[0]{eucov},6,"test BpSlider, coverage is ok");
        is($w->{data}[0]{samples}[0]{totcov},7,"test BpSlider, coverage is ok");
        is($w->{data}[0]{samples}[1]{eucov},7,"test BpSlider, coverage is ok");
        is($w->{data}[0]{samples}[1]{totcov},10,"test BpSlider, coverage is ok");
        
        is($w->{data}[1]{samples}[1]{eucov},8,"test BpSlider, coverage is ok");
        is($w->{data}[1]{samples}[1]{totcov},8,"test BpSlider, coverage is ok");
        is($w->{data}[1]{samples}[2]{eucov},8,"test BpSlider, coverage is ok");
        is($w->{data}[1]{samples}[2]{totcov},9,"test BpSlider, coverage is ok");
        is($w->{data}[1]{samples}[2]{del},1,"test BpSlider, deletion count is ok");
        is($w->{data}[1]{samples}[1]{del},0,"test BpSlider, deletion count is ok");
        is($w->{data}[1]{samples}[2]{A},0,"test BpSlider, A count is ok");
        is($w->{data}[1]{samples}[2]{T},0,"test BpSlider, T count is ok");
        is($w->{data}[1]{samples}[2]{C},8,"test BpSlider, C count is ok");
        is($w->{data}[1]{samples}[2]{G},0,"test BpSlider, G count is ok");
        is($w->{data}[1]{samples}[2]{N},0,"test BpSlider, N count is ok");
        is($w->{data}[0]{samples}[2]{N},2,"test BpSlider, N count is ok");
        
        
        $w=$bpsl->nextWindow();
        is($w->{data}[2]{pos},4,"test BpSlider, position correct");
        is($w->{data}[2]{issnp},0,"test BpSlider, snp correct");
        is($w->{data}[2]{ispuresnp},0,"test BpSlider, snp correct");
        is($w->{data}[2]{iscov},0,"test BpSlider, snp correct");
        is($w->{data}[2]{samples}[1]{del},0,"test BpSlider, deletion count is ok");
        is($w->{data}[2]{samples}[1]{A},0,"test BpSlider, A count is ok");
        is($w->{data}[2]{samples}[1]{T},0,"test BpSlider, T count is ok");
        is($w->{data}[2]{samples}[1]{C},0,"test BpSlider, C count is ok");
        is($w->{data}[2]{samples}[1]{G},0,"test BpSlider, G count is ok");
        is($w->{data}[2]{samples}[1]{N},0,"test BpSlider, N count is ok");
        is($w->{data}[2]{samples}[2]{del},1,"test BpSlider, deletion count is ok");
        is($w->{data}[2]{samples}[2]{A},0,"test BpSlider, A count is ok");
        is($w->{data}[2]{samples}[2]{T},0,"test BpSlider, T count is ok");
        is($w->{data}[2]{samples}[2]{C},8,"test BpSlider, C count is ok");
        is($w->{data}[2]{samples}[2]{G},0,"test BpSlider, G count is ok");
        is($w->{data}[2]{samples}[2]{N},2,"test BpSlider, N count is ok");
    }
    
    sub test_calculatePvalues
    {
        
        my $rbatchfile="test/fet-rbatch-file-19082010.r";
        my $routputfile="test/fet-routput-19082010.txt";
        my $considered = "test/fet-considered-19082010.txt";
        my $discarded = "test/fet-discarded-19082010.txt";
        
        open my $ofh1, ">$considered" or die "Could not open output file";
        open my $ofh2, ">$discarded" or die "Could not open output file";

        
        my %snp_loci=();
        
        Utility::runR($rbatchfile,$routputfile);
        %snp_loci = Utility::Parse_R_output($rbatchfile,$routputfile,"snp",$ofh1,$ofh2,0.8,0.8,\%snp_loci);
        
        close $ofh1;
        open(FILE, "<$considered") || die "Cannot open input file : '$considered'\n$!";
        
        while (<FILE>) {
              
              chomp;
              my $line = $_;
              # discard blank line
              if ($line =~ m/^\s*$/g) {
                    next;
              }
              # discard comment line
              elsif($line =~ m/^\s*#/g) {
                    next;
              }
              else {

                my @record = split("\t",$line);
                ok($record[3],"test calculatePvalues, p-value $record[3] is correct");
              }
        }
    
    }
    
    
    sub is
    {

        my $a=shift;
        my $b=shift;
        my $msg=shift;
        
        
        if($a eq $b)
        {
            print "$testcounter: OK $a = $b; $msg\n";
        }
        else
        {
            print "$testcounter: FAILED $a = $b; $msg\n";
        }
        $testcounter++;
    }
    
    sub not_exists
    {
        my $a=shift;
        my $msg=shift;
        if($a)
        {
            print "$testcounter FAILED; $msg\n";
        }
        else
        {
            print "$testcounter OK: $msg\n";
        }
        $testcounter++;
    }

    sub ok
    {
        my $a=shift;
        my $msg=shift;
        if($a)
        {
            print "$testcounter OK: $msg\n";
        }
        else
        {
            print "$testcounter FAILED; $msg\n";
        }
        $testcounter++;
    }
    
    
}



    #"input=s"      =>\$input,
    #"output=s"     =>\$output,
    #"min-count=i"   =>\$mincount,
    #"min-coverage=i"=>\$mincoverage,
    #"number-simulations=i"=>\$numsimulations,
    #"window-unit=s"  =>\$windowunit,
    #"window-size=i"  =>\$windowsize,
    #"step-size=i"   =>\$step,
    #"test"          =>\$test,
    #"help"         =>\$help
    
=head1 NAME

fisher-exact-test-multi-tier.pl - This script calculates for each SNP a fisher exact test based on the allele counts. Significant Fisher exact tests reflect SNPs for which the populations compared are significantly different from each other. This perl script requires R for the computation of the tests. Instead of computing exact p-values, the script uses R to simulate a number of contingency tables defined by the user from which an estimated p-value is calculated. 

=head1 SYNOPSIS

 perl fisher-exact-test-multi-tier.pl --min-count 2 --min-coverage 4 --max-coverage 1000 --min-snp-fraction 0.75 --number-simulations 100 --min-best-snp-loci 100 --window-size 100 --step-size 10 --window-unit snp --min-pvalue 0.00001 --population 1,2,3,4 --input synchronized-file.txt --output snp-probability.txt --trackname unknown --ucsc-filter 2L --ucsc-prepend dmel

=head1 OPTIONS

=over 4

=item B<--min-count>

the minimum count of the minor allele. used for SNP identification. SNPs will be identified considering all populations simultanously. default=2

=item B<--min-coverage>

the minimum coverage; used for SNP identification, the coverage in ALL populations has to be higher or equal to this threshold, otherwise no SNP will be called. default=4

=item B<--max-coverage>

the maximum coverage; used for SNP identification, the coverage in ALL populations has to be lower or equal to this threshold, otherwise no SNP will be called. default=1000

=item B<--min-snp-fraction>

the minimum snp fraction to be significant in window wise analysis for p-value cutoff. default=0.8

=item B<--number-simulations>

the number of simulation to randomly generate contingency table with fixed horizontal and vertical marginals. default=100

=item B<--min-best-snp-loci>

the minimum number of best snp to simulate with much higher number of simulation to get more significant p value. default=100

=item B<--window-size>

the size of the sliding window. Measured in C<--window-unit>; default=100

=item B<--step-size>

the size of the sliding window steps. Measured in C<--window-unit>; default=10

=item B<--window-unit>

the unit of the sliding window, may either be base pairs (bp) or SNPs (snp); [bp/snp]; default=bp

=item B<--min-covered-fraction>

the minimum fraction of a window being between min-coverage and max-coverage in ALL populations; float; default=0.6

=item B<--population>

the C<--population> parameter allows you to pick any populations among all for p-value caclulation [Optional parameter].
For example if a synchronized file contains 7 population and you want to calculate p-value for 1st, 2nd, 4th, and 6th population then give C<--population 1,2,4,6>.

This option is especially useful when user has a single synchronized file with 7 population then to get all pairwise p-value. user can use this parameter as following:

C<--population 1,2> to calculate p-value between population 1 and 2.
C<--population 1,3> to calculate p-value between population 1 and 3.
C<--population 1,4> to calculate p-value between population 1 and 4.
C<--population 1,5> to calculate p-value between population 1 and 5.
C<--population 1,6> to calculate p-value between population 1 and 6.
C<--population 1,7> to calculate p-value between population 1 and 7.
C<--population 2,3> to calculate p-value between population 2 and 3.
C<--population 2,4> to calculate p-value between population 2 and 4.
C<--population 2,5> to calculate p-value between population 2 and 5.
C<--population 2,6> to calculate p-value between population 2 and 6.
C<--population 2,7> to calculate p-value between population 2 and 7.
C<--population 3,4> to calculate p-value between population 3 and 4.
C<--population 3,5> to calculate p-value between population 3 and 5.
C<--population 3,6> to calculate p-value between population 3 and 6.
C<--population 3,7> to calculate p-value between population 3 and 7.
C<--population 4,5> to calculate p-value between population 4 and 5.
C<--population 4,6> to calculate p-value between population 4 and 6.
C<--population 4,7> to calculate p-value between population 4 and 7.
C<--population 5,6> to calculate p-value between population 5 and 6.
C<--population 5,7> to calculate p-value between population 5 and 7.
C<--population 6,7> to calculate p-value between population 6 and 7.
C<--population 1,2,5,7> to calculate p-value among populations 1,2,5 and 7.

=item B<--also-print-nonsnp>

This parameter should be used only if you want to print non-snp positions also in summary-pvalue file; default=not used

=item B<--input>

The input file has to be synchronized pileup file. Mandatory parameter

=item B<--output>

The output file. Mandatory parameter

=item B<--trackname>

This parameter is used to create a wiggle track format (WIG) output file for p-values. The name of the track. This information may for example be displayed in the IGV viewer. default=unknown

=item B<--ucsc-filter>

This parameter is used to create a wiggle track format (WIG) output file for p-values. The UCSC genome browser only accepts certain chromosome ids. Any attempts loading a wiggle file containing not recognised chromosome ids will result in an error; this option allows to filter for certain user provided
chromosome ids; for example C<--ucsc-filter "2L 2R 3L 3R 4 X"; default=""

=item B<--ucsc-prepend>

This parameter is used to create a wiggle track format (WIG) output file for p-values. For the UCSC genome browser, chromsomes IDs have to be formated according to the requirements of ucsc. For example chr2L has to be provided instead of 2L for D.mel.
This option allows to prepend the specified string to every chromosome id; This step is applied after the filtering using C<--ucsc-filter>; default=""

=item B<--test>

Run the unit tests for this script. 

=item B<--help>

Display help for this script

=back

=head1 Details

=head2 Input

Input is a single tab delimited file which contains a lightwight representation of every pileup file.
Every pileup file represents a population and will be parsed into a list of A-count:T-count:C-count:G-count:N-count:*-count

 2L     5943    C       0:0:10:0:0:0    0:0:10:0:0:0    0:0:10:0:0:0
 2L     5944    T       0:8:0:0:0:0     0:8:0:0:0:0     0:8:0:0:0:0
 2L     5945    G       0:0:0:8:0:0     0:0:0:8:0:0     0:0:0:8:0:0
 2L     5946    C       0:0:9:0:0:0     0:0:9:0:0:0     0:0:9:0:0:0
 2L     5947    T       0:7:0:0:0:0     0:7:0:0:0:0     0:7:0:0:0:0
 
 col1: reference contig (chromosome)
 col2: position in the reference contig
 col3: reference character
 col4: population 1
 col5: population 2
 coln: population n
 
 population data are in the form
 A:T:C:G:N:*
 A: count of character A
 T: count of character T
 C: count of character C
 G: count of character G
 N: count of character N
 *: deletion, count of deletion
 
=head2 Output

This program creates 2 outputs

Output1: An output of this program looks like in the given example (when window-unit snp):

 2L     10452   T       0.009901        1       0.02    100     0, 0, 79, 177, 0, 0, 8, 0
 2L     45207   A       0.009901        1       0.02    100     20, 6, 0, 0, 0, 0, 84, 177
 2L     51931   A       0.009901        1       0.02    100     53, 118, 0, 0, 0, 0, 8, 0
 2L     69287   A       0.009901        1       0.02    100     53, 64, 0, 0, 0, 0, 46, 110
 2L     101241  A       0.009901        1       0.02    100     74, 148, 0, 0, 0, 0, 14, 3
 2L     111842  T       0.009901        1       0.02    100     0, 0, 110, 246, 5, 0, 0, 0
 2L     111844  G       0.009901        1       0.02    100     6, 0, 0, 0, 0, 0, 110, 248
 2L     119806  A       0.009901        1       0.02    100     17, 12, 0, 0, 0, 0, 98, 202
 2L     160545  T       0.009901        1       0.02    100     0, 0, 61, 182, 8, 0, 0, 0
 2L     164590  T       0.009901        1       0.02    100     0, 0, 31, 17, 0, 0, 138, 194

 col 1: reference chromosome
 col 2: position in the reference chromosome
 col 3: reference genome base
 col 4: p-value for each SNP
 col 5: number of SNP
 col 6: p-value cutoff (p-value cutoff was calculated as 2/number of simulation)
 col 7: number of simulation
 col 8: SNP frequency table (A1,A2,A3...T1,T2,T3...C1,C2,C3....G1,G2,G3.....)

Output1: An output of this program looks like in the given example (when window-unit bp):

 2L     111400  C       -2.76938851503629       -10.3798246577252       3       -7.82404601085629       100     0.75    0, 0, 0, 6, 66, 133, 0, 0|0, 0, 110, 246, 5, 0, 0, 0|6, 0, 0, 0, 0, 0, 110, 248
 2L     111500  C       -2.6227857657937        -10.2079744007986       3       -7.82404601085629       100     0.75    0, 0, 0, 6, 66, 133, 0, 0|0, 0, 110, 246, 5, 0, 0, 0|6, 0, 0, 0, 0, 0, 110, 248
 2L     111600  C       -3.00576914148472       -11.2805101980636       3       -7.82404601085629       100     0.75    0, 0, 0, 6, 66, 133, 0, 0|0, 0, 110, 246, 5, 0, 0, 0|6, 0, 0, 0, 0, 0, 110, 248
 2L     111700  C       -2.45572850251064       -10.6674057251772       3       -7.82404601085629       100     0.75    0, 0, 0, 6, 66, 133, 0, 0|0, 0, 110, 246, 5, 0, 0, 0|6, 0, 0, 0, 0, 0, 110, 248
 2L     111800  C       -2.82344982933196       -11.0728708332853       3       -7.82404601085629       100     0.75    0, 0, 0, 6, 66, 133, 0, 0|0, 0, 110, 246, 5, 0, 0, 0|6, 0, 0, 0, 0, 0, 110, 248
 2L     111900  C       -2.71809043599663       -10.9550877976289       3       -7.82404601085629       100     0.75    0, 0, 0, 6, 66, 133, 0, 0|0, 0, 110, 246, 5, 0, 0, 0|6, 0, 0, 0, 0, 0, 110, 248
 2L     112000  C       -2.34653349742236       -9.89431684194352       3       -7.82404601085629       100     0.75    0, 0, 0, 6, 66, 133, 0, 0|0, 0, 110, 246, 5, 0, 0, 0|6, 0, 0, 0, 0, 0, 110, 248
 2L     112100  T       -4.20972174426658       -8.53719285812341       2       -3.91202300542815       100     0.75    0, 0, 110, 246, 5, 0, 0, 0|6, 0, 0, 0, 0, 0, 110, 248
 2L     112200  T       -4.20972174426658       -8.53719285812341       2       -3.91202300542815       100     0.75    0, 0, 110, 246, 5, 0, 0, 0|6, 0, 0, 0, 0, 0, 110, 248


 col 1: reference chromosome
 col 2: middle position of sliding window in the reference chromosome
 col 3: reference genome base
 col 4: log(average of p-values in a given window)
 col 5: log(product of p-values in a given window)
 col 6: number of SNP
 col 7: log(p-value cutoff)
 col 8: number of simulation
 col 9: C<--min-snp-fraction> the minimum snp fraction to be significant in window wise analysis for p-value cutoff.
 col 10: SNP frequency table (A1,A2,A3...T1,T2,T3...C1,C2,C3....G1,G2,G3.....)

Output2: Two summary output files are generated by this script ie. output-summary-pvalue.txt (tab delemited file) and  output-summary-pvalue.wig (wiggle file to display in IGV or UCSC genome browser). Both summary files contains minimum p-value for all SNP.

 2L     9229    0.9307  1
 2L     9385    0.6832  1
 2L     9418    0.5644  1
 2L     9507    0.4356  1
 2L     9533    1       1
 2L     9644    0.5842  1
 
 
 col 1: reference chromosome
 col 2: position in the reference chromosome
 col 3: p-value for each SNP
 col 4: number of SNP

=head1 DESCRIPTION

=item B<--min-pvalue>

the minimum p-value cut off  to filter all snp with > min-pvalue cutoff.
If user defines C<--min-pvalue> then this program will create two output file 1) output-considered.txt and 2) output-discarded.txt

=head1 AUTHORS

Ram vinay pandey

Robert Kofler

Pablo Orozco terWengel

Christian Schloetterer

=cut
