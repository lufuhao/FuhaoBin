#!/usr/bin/env perl
# This script can be used to find the Reverse Complement in a DNA Sequence.

# While executing this script it asks for the file name of the DNA sequence. 
# If the DNA sequence file is not in the same directory of this script, enter the file name with its full path.
# Example: 
# In windows:  c:\rnafile.txt 
# In Linux  : /home/user/sequence/rnafile.txt

# Send your comments and suggessions to techcuriosity @gmail.com
# Website : http://www.techcuriosity.com

use File::Path;
print "\n\t#################### GET REVERSE COMPLEMENT OF DNA ####################\n\n";
print "PLEASE TYPE THE FILENAME OF THE DNA SEQUENCE := ";
$dnafilename = <STDIN>;
chomp $dnafilename;
unless ( open(DNAFILE, $dnafilename) ) {
	print "Error: can not open file \"$dnafilename\"\n\n";
	exit;
}
@DNA = <DNAFILE>;
close DNAFILE;
$DNA = join( '', @DNA);
$DNA =~ s/\s//g;
print "\nTHE ORIGINAL DNA SEQUENCE :=\n$DNA\n\n";
@DNA = split( '', $DNA );
print"REVERSE COMPLEMENT OF THE DNA SEQUENCE :=\n";

foreach $nucleotide(reverse(@DNA)) {
	if ($nucleotide =~ /a/i) {
		print "T";
		print WRITE "T";
	}
	elsif ($nucleotide =~ /t/i) {
		print "A";
		print WRITE "A";
	}
	elsif ($nucleotide =~ /g/i) {
		print "C";
		print WRITE "C";
	}
	elsif ($nucleotide =~ /c/i) {
		print "G";
		print WRITE "G";
	}
	else {
		die "$0:  Bad nucleotide!  [$nucleotide]\n";
	}
}

