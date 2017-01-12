#!/usr/bin/perl

# Correlation
# Didier Gonze
# Updated: 28/4/2004

##########################################################################################


&ReadArguments;

&ReadData;

&Correlation;


##########################################################################################
### Read arguments from the command line

sub ReadArguments {

    $verbo=0;
    $infile = "";
    $col1=1; # 1st column
    $col2=2; # 2nd column
    $title=0;
    
    foreach my $a (0..$#ARGV) {

    ### help
    if ($ARGV[0] eq "-h") {
    	die "Syntax: correlation.pl -i filename -col # #\n";
    }
    elsif ($ARGV[0] eq "-help") {
    	&PrintHelp;
    }
            
    ### input file
    elsif ($ARGV[$a] eq "-i") {
	$ok=1;
    	$infile = $ARGV[$a+1];
    }

    ### column with the data
    elsif ($ARGV[$a] eq "-col") {
    	$col1 = $ARGV[$a+1];
    	$col2 = $ARGV[$a+2];
    }
	
    ### verbosity 
    elsif ($ARGV[$a] eq "-v") {
	$verbo=1;
    }
	
    }
	
    if ($infile eq "") {
    	die "STOP! You have to give the name of the input file!\n";
    }
	
}  #  End of ReadArguments

##########################################################################################
### Print help


sub PrintHelp {
  open HELP, "| more";
  print <<EndHelp;
NAME
        correlation.pl

DESCRIPTION
	Calculate the (Pearson) correlation coefficient between two columns
	of data. The formula are as described in Wolfram website:
	http://mathworld.wolfram.com/CorrelationCoefficient.html

AUTHOR
	Didier Gonze (dgonze\@ulb.ac.be)  

OPTIONS
	-i file_name
		Specify the input file containing the data. 
		This argument is obligatory (except if using option -h).
	       
	-col # #
		Specify the columns containing the data (default: 1 2).
		
	-v 
		Verbosity: print detailed informations during the 
		process.
	       
	-h 
		Give syntax. This argument must be the first.

	-help 
		Give detailed help (print this message). This argument 
		must be the first.

EXAMPLE
        perl correlation.pl -i datafile -col 2 3

EndHelp

close HELP;
die "\n";

}  #  End of PrintHelp


##########################################################################################
### Read the data and fill the data vector

sub ReadData {

my $i=0;
my $badlines=0;

open inf, $infile or die "STOP! File $infile not found.\n";
if ($verbo==1) {print "Open input file: $infile\n";}
if ($verbo==1) {print "Select columns: $col1 and $col2\n";}

foreach $line (<inf>){
    chomp $line;
    @line=split /\t/,$line;
    
    if ($line[$col1-1] =~ /\d+/ and $line[$col2-1] =~ /\d+/){
       $i++;
       $x[1][$i]=$line[$col1-1];
       $x[2][$i]=$line[$col2-1];
    }
    else{
       $badlines++;
    }
}

$nbdata=$i;

if ($nbdata==0) {die "STOP! No numeric data...\n";}

if ($verbo==1) {print "Total number of data = $nbdata\n";}
if ($verbo==1) {print "Number of rejected lines (no numeric data) = $badlines\n";}

close inf;

}  # End of ReadData


##########################################################################################
### Correlation

sub Correlation {

$mean[1]=&Mean(1);
$mean[2]=&Mean(2);

if ($verbo==1) {
  $xmean1=sprintf("%.3f",$mean[1]);
  $xmean2=sprintf("%.3f",$mean[2]);
  print "Mean($col1) = $xmean1\n";
  print "Mean($col2) = $xmean2\n";
}

$ssxx=&SS(1,1);
$ssyy=&SS(2,2);
$ssxy=&SS(1,2);

if ($verbo==1) {
  $xssxx=sprintf("%.3f",$ssxx);
  $xssyy=sprintf("%.3f",$ssyy);
  $xssxy=sprintf("%.3f",$ssxy);
  print "SS($col1)= $xssxx\n";
  print "SS($col2) = $xssyy\n";
  print "SS($col1,$col2) = $xssxy\n";
}

$correl=&Correl($ssxx,$ssyy,$ssxy);

$xcorrel=sprintf("%.4f",$correl);

print "Correlation = $xcorrel\n";

}  # End of Correlation


##########################################################################################
### Mean

sub Mean {

my ($a)=@_;
my ($i,$sum)=(0,0);

for ($i=1;$i<=$nbdata;$i++){
  $sum=$sum+$x[$a][$i];
}
$mu=$sum/$nbdata;

return $mu;

}

##########################################################################################
### SS = sum of squared deviations to the mean

sub SS {

my ($a,$b)=@_;
my ($i,$sum)=(0,0);

for ($i=1;$i<=$nbdata;$i++){
  $sum=$sum+($x[$a][$i]-$mean[$a])*($x[$b][$i]-$mean[$b]);
}

return $sum;

}

##########################################################################################
### Correlation

sub Correl {

my($ssxx,$ssyy,$ssxy)=@_;

$sign=$ssxy/abs($ssxy);

$correl=$sign*sqrt($ssxy*$ssxy/($ssxx*$ssyy));

return $correl;

}
