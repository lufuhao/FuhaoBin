#!/usr/bin/env perl
use warnings;
use strict;
 
my $file = $ARGV;
my $offset = &test_qualities($file);
exit;
 
sub iterator{
	my $handle = shift || die $!;
	my %return;
	return sub{ ## actual iterator
		my %return;
		$return{'head'} = readline($handle) || return; #if the next line exists , get it otherwise return null
		$return{'seq'} = readline($handle);
		$return{'head2'} = readline($handle);
		$return{'quals'} = readline($handle);
		map {chomp $return{$_}} keys %return;  
		return \%return;
	};
}
 
sub test_qualities{
	my $file = shift;
	open(my $TEST, '<', $file) || die $!;
	my $iterator = &iterator($TEST);
	while (my $record = $iterator->()){
		map {if ($_ > 74) {return 64;} elsif ($_ < 59) {return 33}} unpack("W*" , $record->{'quals'});
	}
	close $TEST;
	die "Was unable to determine quality offset , you'll need to do this yourself\n";
}

=reference
  SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS.....................................................
  ..........................XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX......................
  ...............................IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII......................
  .................................JJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJ......................
  LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL....................................................
  !"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_`abcdefghijklmnopqrstuvwxyz{|}~
  |                         |    |        |                              |                     |
 33                        59   64       73                            104                   126
  0........................26...31.......40                                
                           -5....0........9.............................40 
                                 0........9.............................40 
                                    3.....9.............................40 
  0........................26...31........41                               

 S - Sanger        Phred+33,  raw reads typically (0, 40)
 X - Solexa        Solexa+64, raw reads typically (-5, 40)
 I - Illumina 1.3+ Phred+64,  raw reads typically (0, 40)
 J - Illumina 1.5+ Phred+64,  raw reads typically (3, 40)
     with 0=unused, 1=unused, 2=Read Segment Quality Control Indicator (bold) 
     (Note: See discussion above).
 L - Illumina 1.8+ Phred+33,  raw reads typically (0, 41)
 =cut
