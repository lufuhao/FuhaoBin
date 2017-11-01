#!/bin/bash
################# help message ######################################
help() {
cat<<HELP

$0 --- Brief Introduction

Version: 20170619

Requirements:
	Linux: wget, md5sum, tar

Descriptions:
	download and setup NCBI NR database

Options:
  -h    Print this help message

Example:
  $0 download_dir

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
[ "$1" = "-h" ] || [ "$1" = "--help" ] && help
#################### Environments ###################################
echo -e "\n######################\nProgram $ProgramName initializing ...\n######################\n"
#echo "Adding $RunDir/bin into PATH"
#export PATH=$RunDir/bin:$RunDir/utils/bin:$PATH

#################### Initializing ###################################
#opt_s=0
#opt_a=0
#opt_t=1
#################### Parameters #####################################
#while [ -n "$1" ]; do
#  case "$1" in
#    -h) help;shift 1;;
#    -i) opt_i=$2;shift 2;;
#    -t) opt_t=$2;shift 2;;
#    -s) opt_s=1;shift 1;;
#    -a) opt_a=1;shift 1;;
#    --) shift;break;;
#    -*) echo "error: no such option $1. -h for help" > /dev/stderr;exit 1;;
#    *) break;;
#  esac
#done


#################### Subfuctions ####################################
###Detect command existence
CmdExists () {
  if command -v $1 >/dev/null 2>&1; then
    return 0
  else
    return 1
  fi
}





#################### Command test ###################################

if [[ $(CmdExists 'wget') -ne 0 ]]; then
	echo "Error: CMD 'wget' is required but not found.  Aborting..." >&2 
	exit 127
fi
if [[ $(CmdExists 'tar') -ne 0 ]]; then
	echo "Error: CMD 'tar' is required but not found.  Aborting..." >&2 
	exit 127
fi
if [[ $(CmdExists 'md5sum') -ne 0 ]]; then
	echo "Error: CMD 'md5sum' is required but not found.  Aborting..." >&2 
	exit 127
fi



#################### Defaults #######################################




#################### Input and Output ###############################




#################### Main ###########################################

### NR fasta: ftp://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nr.gz
### NR db: ftp://ftp.ncbi.nlm.nih.gov/blast/db/





echo "Downloading NCBI nr protein database"

if [ -z "$1" ]; then
	installdir=$(pwd);
else
	installdir=$1;
	cd $installdir
fi

downloadtest=0;
dbpath='ftp://ftp.ncbi.nlm.nih.gov/blast/db/'
numdigit=2
firstnumber=0



while [[ $downloadtest -eq 0 ]]; do
	curnum=$(printf "%.${numdigit}d" $firstnumber)
	if [ ! -s nr.$curnum.tar.gz ]; then
		wget ftp://ftp.ncbi.nlm.nih.gov/blast/db/nr.$curnum.tar.gz > /dev/null 2>&1
		if [[ $? -ne 0 ]] || [ ! -s nr.$curnum.tar.gz ]; then
			echo "Error: Failed to download nr.$curnum.tar.gz" >&2
			downloadtest=1
			exit 0;
		fi
		wget ftp://ftp.ncbi.nlm.nih.gov/blast/db/nr.$curnum.tar.gz.md5 > /dev/null 2>&1
		if [[ $? -ne 0 ]] || [ ! -s nr.$curnum.tar.gz.md5 ]; then
			echo "Error: Failed to download nr.$curnum.tar.gz.md5" >&2
			downloadtest=1
			exit 0;
		fi
		md5sum -c nr.$curnum.tar.gz.md5
		if [[ $? -ne 0 ]]; then
			echo "Error: md5 not match nr.$curnum.tar.gz.md5" >&2
			exit 0;
		fi
		rm nr.$curnum.tar.gz.md5
	fi
	tar xvf nr.$curnum.tar.gz >/dev/null 2>&1
	rm nr.$curnum.tar.gz
	((firstnumber++))
done
