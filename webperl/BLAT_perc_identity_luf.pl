#!usr/bin/env perl
#http://happyjlj1988.blog.163.com/blog/static/182475182201452704634716
#代码只能使用于mRNA-DNA的比对
#Blat的"percent identity"的Perl代码
use warnings;
use strict;
my ($file,$line,$pid,$query,$subject);
open $file,"alignment_up_gar.psl";
<$file>; <$file>;  <$file>; <$file>; <$file>;
print join('	',qw(query subject p_identify(%))),"\n";   
while($line=<$file>)
{
$line=trim($line);
($query,$subject,$pid) = pslCalcMilliBad($line);
print join('	',($query,$subject,$pid)),"\n";   
} 

 
sub trim
{
    my $value=shift;
chomp($value);
$value=~s/^\s*//;
$value=~s/\s*$//;
return $value;	
}
 

sub pslCalcMilliBad 
#/* Calculate badness in parts per thousand. */
{
    my $line=shift;   
my $sizeMul =1; #not protein
my $isMrna=1;  # mRNA
my @cols=split(/\t/,$line);
my ($qAliSize, $tAliSize, $aliSize);
my $milliBad = 0;
my $sizeDif;
my $insertFactor;
my $total;

# cols[0]  matches 
# cols[1]  misMatches 
# cols[2]  repMaches 
# cols[4]  qNumInsert 
# cols[6]  tNumInsert 
# cols[11] qStart 
# cols[12] qEnd 
# cols[15] tStart 
# cols[16] tEnd  
# cols[13] subject
# cols[9] query
my ($matches, $misMatches,$repMaches,$qNumInsert,$tNumInsert,$qStart,$qEnd,$tStart,$tEnd,$subject,$query)=($cols[0],$cols[1],$cols[2],$cols[4],$cols[6],$cols[11],$cols[12],$cols[15],$cols[16],$cols[13],$cols[9]);
$qAliSize = $sizeMul * ($qEnd - $qStart);
$tAliSize = $tEnd - $tStart;
my $aliSize; 
if($qAliSize < $tAliSize)
{
  $aliSize = $qAliSize;
}
else
{
$aliSize = $tAliSize; 
}

if ($aliSize <= 0)
{
return ($query,$subject,100);
}
$sizeDif = $qAliSize - $tAliSize;
if ($sizeDif < 0)
{
if ($isMrna)
{
$sizeDif = 0;
}
else
{
$sizeDif = -($sizeDif);
}
}

$insertFactor = $qNumInsert;
if (!$isMrna)
{
$insertFactor += $tNumInsert;
}

$total = ($sizeMul * ($matches + $repMaches + $misMatches));
if ($total != 0)
{
$milliBad = (1000 * ($misMatches*$sizeMul + $insertFactor + 
round(3*log(1+$sizeDif)))) / $total;
}	 
return ($query,$subject,(100.0 - ($milliBad * 0.1 ) ));
}
 
sub round {
    my $number = shift;
    return int( $number + .5 );
}
