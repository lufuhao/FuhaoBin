#!/bin/bash
RunDir=$(cd `dirname $(readlink -f $0)`; pwd)
MachType=$(uname -m)

################# help message ######################################
help () {
cat<<HELP

$0 blast.outfmt6 ...

Version: 20170918

Requirement:
	blast.filter.pl
	blast92gff3.pl

Descriptions:
  Filter blast hsp length >=250
  Convert blast tab to gff3

Author:
  Fu-Hao Lu
  Post-Doctoral Scientist in Micheal Bevan laboratory
  Cell and Developmental Department, John Innes Centre
  Norwich NR4 7UH, United Kingdom
  E-mail: Fu-Hao.Lu@jic.ac.uk
HELP
exit 0
}
[ -z "$1" ] && help
[ "$1" = "-h" ] && help
[ "$1" = "--help" ] && help



#################### Subfuctions ####################################
###Detect command existence
CmdExists () {
  if command -v $1 >/dev/null 2>&1; then
    return 0
  else
    return 1
  fi
}



####################################################################
if [ $(CmdExists blast.filter.pl) -ne 0 ]; then
	echo "Error: script 'blast.filter.pl'  is required but not found.  Aborting..." >&2
	exit 100
fi
if [ $(CmdExists blast92gff3.pl) -ne 0 ]; then
	echo "Error: script 'blast92gff3.pl'  is required but not found.  Aborting..." >&2
	exit 100
fi



###################################################################
for i in $@; do
	echo $i
	perl /usr/users/celldev/luf/local/perl/blast.filter.pl $i 250 $i.fil250
	perl /usr/users/celldev/luf/local/perl/blast92gff3.pl < $i.fil250 > $i.fil250.gff3
done 
