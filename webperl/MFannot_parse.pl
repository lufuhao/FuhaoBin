#!/usr/bin/perl -w
use strict;
use warnings;
use Cwd;
use Getopt::Long;
use FindBin qw($Bin);
#http://blog.sina.com.cn/s/blog_83f77c940102wjqe.html
#用来解析在线线粒体注释工具MFannot的输出工具，主要的输入文件是tbl和输入的线粒体的fasta
my($tbl,$fasta,$outdir);
GetOptions(
     "tbl:s"=>\$tbl,
     "fa:s"=>\$fasta,
     "o:s"=>\$outdir,
           );
$outdir||=getcwd;

sub usage{
	print qq {
This script will parse mfannot(http://megasun.bch.umontreal.ca/cgi-bin/mfannot/mfannotInterface.pl) result.

usage:
perl $0 -fa genome.fasta -tbl mfannot.tbl -o outputdirectory
-fa                 the genome (fasta file)
-p                  the prefix of output
-o                  /path/to/out/directoy
-tbl                tbl file from mfannot output
Email:fanyucai1\@126.com
2016.12.7
};
	exit;
}



if(!$tbl ||!$fasta) {
	&usage();
}



my (%gene,%tRNA);
open(FA,"$fasta");
my $seq;
while() {
	chomp;
	if($_!~"\>") {
		$seq.=$_;
	}
}
close FA;



open(TBL,"$tbl");
my (%start,%end,$num,%anno);
$num=-1;
my %hash1;
my %hash2;
my %hash3;#"+"or "-"
my %name;#genename
while() {
	chomp;
	my @array=split(/\t/,$_);
	if($#array==2 && $array[2] eq "gene") {#get the position
		$num++;
		$start{$num}=$array[0];
		$end{$num}=$array[1];
	}
	if( $#array>3 && ($array[3] eq "gene") ) {#get the gene name
		$name{$num}=$array[4];
		if($array[4]=~ /^rnl/ || $array[4] =~ /^rns/) {
			$hash1{$num}="rRNA";
		}
		elsif($array[4]=~ /^trn/ ) {
			$hash1{$num}="tRNA";
		}
		else {
			$hash1{$num}="gene";
		}
	}
	if($#array==1) {#get the sequence when face the exon
		if($array[1]>$array[0]) {
			$hash2{$num}.=substr($seq,$array[0]-1,$array[1]-$array[0]+1);
			$hash3{$num}="+";
		}
		else {
			$hash2{$num}.=substr($seq,$array[1]-1,$array[0]-$array[1]+1);
			$hash3{$num}="-";
		}
		if($name{$num}=~/cox/i) {
			print "cox\t$array[0]\t$array[1]\t$name{$num}\n";
		}
		elsif($name{$num}=~/orf/i) {
			print "orf\t$array[0]\t$array[1]\t$name{$num}\n";            
		}
		elsif($name{$num}=~/cob/i) {
			print "cob\t$array[0]\t$array[1]\t$name{$num}\n"; 
		}
		elsif($name{$num}=~/nad/i) {
			print "nad\t$array[0]\t$array[1]\t$name{$num}\n"; 
		}
		else {
			print "$hash1{$num}\t$array[0]\t$array[1]\t$name{$num}\n";
		}
    }
    if($#array==2 && ($_=~"CDS" ||$_=~"RNA")) {#get the sequence
        if($array[1]>$array[0]) {
            $hash2{$num}.=substr($seq,$array[0]-1,$array[1]-$array[0]+1);
            $hash3{$num}="+";
        }
        else {
            $hash2{$num}.=substr($seq,$array[1]-1,$array[0]-$array[1]+1);
            $hash3{$num}="-";
        }
        if($name{$num}=~/cox/i) {
            print "cox\t$array[0]\t$array[1]\t$name{$num}\n";
        }
        elsif($name{$num}=~/orf/i) {
            print "orf\t$array[0]\t$array[1]\t$name{$num}\n";            
        }
        elsif($name{$num}=~/cob/i) {
            print "cob\t$array[0]\t$array[1]\t$name{$num}\n"; 
        }
        elsif($name{$num}=~/nad/i) {
            print "nad\t$array[0]\t$array[1]\t$name{$num}\n"; 
        }
        else {
            print "$hash1{$num}\t$array[0]\t$array[1]\t$name{$num}\n";
        }
    }
    if($#array>3 && $array[3] ne "note" &&  $array[3] ne "number" && $array[3] ne "gene") {
        for (my $i=3;$i<=$#array;$i++) {
            $anno{$num}.=$array[$i];
        }
    }
}
close TBL;



open(TT,">$outdir/tRNA.fasta");
open(RR,">$outdir/rRNA.fasta");
open(GN,">$outdir/gene.fasta");
for(my $k=0;$k<=$num;$k++) {
    if($hash3{$k} eq "-") {
        $hash2{$k}=~ tr/ACGTacgt/TGCAtgca/;
        $hash2{$k}=reverse($hash2{$k});
    }
    if($hash1{$k} eq "gene") {
       print GN ">",$name{$k},"\t",$hash3{$k},"\t",$start{$k},"\t",$end{$k},"\t",$anno{$k},"\n",$hash2{$k},"\n";
    }
    if($hash1{$k} eq "tRNA") {
       print TT ">","$name{$k}","\t\t",$hash3{$k},"\t",$start{$k},"\t",$end{$k},"\t",$anno{$k},"\n",$hash2{$k},"\n";
    }
    if($hash1{$k} eq "rRNA") {
       print RR  ">",$name{$k},"\t",$hash3{$k},"\t",$start{$k},"\t",$end{$k},"\t",$anno{$k},"\n",$hash2{$k},"\n";
    }
}
close TT;
close RR;
close GN

