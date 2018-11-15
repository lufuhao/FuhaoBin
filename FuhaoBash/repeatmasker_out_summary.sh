#!/bin/bash
### https://www.biostars.org/p/86669/

### Exit if command fails
set -o errexit
### Set readonly variable
#readonly passwd_file=”/etc/passwd”
### exit when variable undefined
set -o nounset
### Script Root
RootDir=$(cd `dirname $(readlink -f $0)`; pwd)
### MachType
if [ ! -z $(uname -m) ]; then
	machtype=$(uname -m)
elif [ ! -z "$MACHTYPE" ]; then
	machtype=$MACHTYPE
else
	echo "Warnings: unknown MACHTYPE" >&2
fi

#export NUM_THREADS=`grep -c '^processor' /proc/cpuinfo 2>/dev/null || sysctl -n hw.ncpu 2>/dev/null || echo 1`;
ProgramName=${0##*/}
echo "MachType: $machtype"
echo "RootPath: $RootDir"
echo "ProgName: $ProgramName"
RunPath=$PWD
echo "RunDir: $RunPath"

###echo color
#Black        0;30     Dark Gray     1;30
#Red          0;31     Light Red     1;31
#Green        0;32     Light Green   1;32
#Brown/Orange 0;33     Yellow        1;33
#Blue         0;34     Light Blue    1;34
#Purple       0;35     Light Purple  1;35
#Cyan         0;36     Light Cyan    1;36
#Light Gray   0;37     White         1;37
#RED='\033[0;31m'
#NC='\033[0m' # No Color
#printf "I ${RED}love${NC} Stack Overflow\n"

################# help message ######################################
help() {
cat<<HELP

$0 --- Brief Introduction

Version: 20181018

Requirements:
	Linux

Descriptions:
	Summarize RepeatMasker output

Example:
  $0 repeatmasker.out output_prefix

Author:
  Fu-Hao Lu
  Post-Doctoral Scientist in Micheal Bevan laboratory
  Cell and Developmental Department, John Innes Centre
  Norwich NR4 7UH, United Kingdom
  E-mail: Fu-Hao.Lu@jic.ac.uk
HELP
exit 0
}
[[ -z "$@" ]] && help
[ -z "$1" ] && help
[ "$1" = "-h" ] || [ "$1" = "--help" ] && help
#################### Environments ###################################

BedIn=$1
OutPrefix=$2

cat $BedIn | tr -s ' ' | sed 's/^ *//g' | tr ' ' '\t' > $OutPrefix.tabout

BedTabIn="$OutPrefix.tabout"

#grep -v '^SW   perc perc perc\|score   div. del. ins\|^$' $BedTabIn | awk '{print $7-$6+1,$11}' | awk '{group[$2]}; {count[$2]+=$1}; END {for (i in group) print i"\t"count[i]}' | sort -k1,1 > $OutPrefix.subclass.bed
grep -P ^"\d+\t" $BedTabIn | awk '{print $7-$6+1,$11}' | awk '{group[$2]}; {count[$2]+=$1}; END {for (i in group) print i"\t"count[i]}' | sort -k1,1 > $OutPrefix.subclass.bed

OutPfx="SUM"
echo -n -e "$OutPfx\t";
grep -P ^"\d+\t" $BedTabIn | awk '{print $5,"\t",$6-1,"\t",$7,"\t",$11}' | perl -lne 's/ //g;print;' | sort -k1,1 -k2,2n -k3,3n | cut -f 1,2,3 > $OutPrefix.$OutPfx.bed
bedtools merge -i $OutPrefix.$OutPfx.bed | perl -lane '$sum=$sum+$F[2]-$F[1];END {print $sum;}'

OutPfx="DNA"
echo -n -e "$OutPfx\t";
grep -P ^"\d+\t" $BedTabIn | awk '{print $5,"\t",$6-1,"\t",$7,"\t",$11}' | perl -lne 's/ //g;@F=split(/\t/);print if($F[3]=~/^DNA/);' | sort -k1,1 -k2,2n -k3,3n| cut -f 1,2,3  | perl -lne 's/ //g;print;'> $OutPrefix.$OutPfx.bed
bedtools merge -i $OutPrefix.$OutPfx.bed | perl -lane '$sum=$sum+$F[2]-$F[1];END {print $sum;}'

OutPfx="DNA-En-Spm"
echo -n -e "$OutPfx\t";
grep -P ^"\d+\t" $BedTabIn | awk '{print $5,"\t",$6-1,"\t",$7,"\t",$11}' | perl -lne 's/ //g;@F=split(/\t/);print if($F[3]=~/^DNA\/En-Spm/ or $F[3]=~/^DNA\/CMC-EnSpm/);' | sort -k1,1 -k2,2n -k3,3n| cut -f 1,2,3  | perl -lne 's/ //g;print;'> $OutPrefix.$OutPfx.bed
bedtools merge -i $OutPrefix.$OutPfx.bed | perl -lane '$sum=$sum+$F[2]-$F[1];END {print $sum;}'
#OutPfx="DNA-CMC-EnSpm"
#echo -n -e "$OutPfx\t";
#grep -P ^"\d+\t" $BedTabIn | awk '{print $5,"\t",$6-1,"\t",$7,"\t",$11}' | perl -lne 's/ //g;@F=split(/\t/);print if($F[3]=~/^DNA\/CMC-EnSpm/);' | sort -k1,1 -k2,2n -k3,3n| cut -f 1,2,3  | perl -lne 's/ //g;print;'> $OutPrefix.$OutPfx.bed
#bedtools merge -i $OutPrefix.$OutPfx.bed | perl -lane '$sum=$sum+$F[2]-$F[1];END {print $sum;}'
OutPfx="LINE"
echo -n -e "$OutPfx\t";
grep -P ^"\d+\t" $BedTabIn | awk '{print $5,"\t",$6-1,"\t",$7,"\t",$11}' | perl -lne 's/ //g;@F=split(/\t/);print if($F[3]=~/^LINE/);' | sort -k1,1 -k2,2n -k3,3n| cut -f 1,2,3  | perl -lne 's/ //g;print;'> $OutPrefix.$OutPfx.bed
bedtools merge -i $OutPrefix.$OutPfx.bed | perl -lane '$sum=$sum+$F[2]-$F[1];END {print $sum;}'
OutPfx="Low_complexity"
echo -n -e "$OutPfx\t";
grep -P ^"\d+\t" $BedTabIn | awk '{print $5,"\t",$6-1,"\t",$7,"\t",$11}' | perl -lne 's/ //g;@F=split(/\t/);print if($F[3]=~/^Low_complexity/);' | sort -k1,1 -k2,2n -k3,3n| cut -f 1,2,3  | perl -lne 's/ //g;print;'> $OutPrefix.$OutPfx.bed
bedtools merge -i $OutPrefix.$OutPfx.bed | perl -lane '$sum=$sum+$F[2]-$F[1];END {print $sum;}'

OutPfx="LTR"
echo -n -e "$OutPfx\t";
grep -P ^"\d+\t" $BedTabIn |    awk '{print $5,"\t",$6-1,"\t",$7,"\t",$11}' | perl -lne 's/ //g;@F=split(/\t/);print if($F[3]=~/^LTR/);' | sort -k1,1 -k2,2n -k3,3n| cut -f 1,2,3  | perl -lne 's/ //g;print;'> $OutPrefix.$OutPfx.bed
bedtools merge -i $OutPrefix.$OutPfx.bed | perl -lane '$sum=$sum+$F[2]-$F[1];END {print $sum;}'
OutPfx="LTR-Copia"
echo -n -e "$OutPfx\t";
grep -P ^"\d+\t" $BedTabIn |    awk '{print $5,"\t",$6-1,"\t",$7,"\t",$11}' | perl -lne 's/ //g;@F=split(/\t/);print if($F[3]=~/^LTR\/Copia/);' | sort -k1,1 -k2,2n -k3,3n| cut -f 1,2,3  | perl -lne 's/ //g;print;'> $OutPrefix.$OutPfx.bed
bedtools merge -i $OutPrefix.$OutPfx.bed | perl -lane '$sum=$sum+$F[2]-$F[1];END {print $sum;}'
OutPfx="LTR-Gypsy"
echo -n -e "$OutPfx\t";
grep -P ^"\d+\t" $BedTabIn |    awk '{print $5,"\t",$6-1,"\t",$7,"\t",$11}' | perl -lne 's/ //g;@F=split(/\t/);print if($F[3]=~/^LTR\/Gypsy/);' | sort -k1,1 -k2,2n -k3,3n| cut -f 1,2,3  | perl -lne 's/ //g;print;'> $OutPrefix.$OutPfx.bed
bedtools merge -i $OutPrefix.$OutPfx.bed | perl -lane '$sum=$sum+$F[2]-$F[1];END {print $sum;}'
OutPfx="MobileElement"
echo -n -e "$OutPfx\t";
grep -P ^"\d+\t" $BedTabIn |    awk '{print $5,"\t",$6-1,"\t",$7,"\t",$11}' | perl -lne 's/ //g;@F=split(/\t/);print if($F[3]=~/^MobileElement/);' | sort -k1,1 -k2,2n -k3,3n| cut -f 1,2,3  | perl -lne 's/ //g;print;'> $OutPrefix.$OutPfx.bed
bedtools merge -i $OutPrefix.$OutPfx.bed | perl -lane '$sum=$sum+$F[2]-$F[1];END {print $sum;}'
OutPfx="Other"
echo -n -e "$OutPfx\t";
grep -P ^"\d+\t" $BedTabIn |    awk '{print $5,"\t",$6-1,"\t",$7,"\t",$11}' | perl -lne 's/ //g;@F=split(/\t/);print if($F[3]=~/^Other/);' | sort -k1,1 -k2,2n -k3,3n| cut -f 1,2,3  | perl -lne 's/ //g;print;'> $OutPrefix.$OutPfx.bed
bedtools merge -i $OutPrefix.$OutPfx.bed | perl -lane '$sum=$sum+$F[2]-$F[1];END {print $sum;}'
OutPfx="RC"
echo -n -e "$OutPfx\t";
grep -P ^"\d+\t" $BedTabIn |    awk '{print $5,"\t",$6-1,"\t",$7,"\t",$11}' | perl -lne 's/ //g;@F=split(/\t/);print if($F[3]=~/^RC/);' | sort -k1,1 -k2,2n -k3,3n| cut -f 1,2,3  | perl -lne 's/ //g;print;'> $OutPrefix.$OutPfx.bed
bedtools merge -i $OutPrefix.$OutPfx.bed | perl -lane '$sum=$sum+$F[2]-$F[1];END {print $sum;}'
OutPfx="rRNA"
echo -n -e "$OutPfx\t";
grep -P ^"\d+\t" $BedTabIn |    awk '{print $5,"\t",$6-1,"\t",$7,"\t",$11}' | perl -lne 's/ //g;@F=split(/\t/);print if($F[3]=~/^rRNA/);' | sort -k1,1 -k2,2n -k3,3n| cut -f 1,2,3  | perl -lne 's/ //g;print;'> $OutPrefix.$OutPfx.bed
bedtools merge -i $OutPrefix.$OutPfx.bed | perl -lane '$sum=$sum+$F[2]-$F[1];END {print $sum;}'
OutPfx="Simple_repeat"
echo -n -e "$OutPfx\t";
grep -P ^"\d+\t" $BedTabIn |    awk '{print $5,"\t",$6-1,"\t",$7,"\t",$11}' | perl -lne 's/ //g;@F=split(/\t/);print if($F[3]=~/^Simple_repeat/);' | sort -k1,1 -k2,2n -k3,3n| cut -f 1,2,3  | perl -lne 's/ //g;print;'> $OutPrefix.$OutPfx.bed
bedtools merge -i $OutPrefix.$OutPfx.bed | perl -lane '$sum=$sum+$F[2]-$F[1];END {print $sum;}'
OutPfx="SINE"
echo -n -e "$OutPfx\t";
grep -P ^"\d+\t" $BedTabIn |    awk '{print $5,"\t",$6-1,"\t",$7,"\t",$11}' | perl -lne 's/ //g;@F=split(/\t/);print if($F[3]=~/^SINE/);' | sort -k1,1 -k2,2n -k3,3n| cut -f 1,2,3  | perl -lne 's/ //g;print;'> $OutPrefix.$OutPfx.bed
bedtools merge -i $OutPrefix.$OutPfx.bed | perl -lane '$sum=$sum+$F[2]-$F[1];END {print $sum;}'

