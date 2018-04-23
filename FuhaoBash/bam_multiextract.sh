#!/bin/bash
RunDir=$(cd `dirname $(readlink -f $0)`; pwd)
MachType=$(uname -m)

################# help message ######################################
help()
{
cat<<HELP

$0 --- Extract multiple sequences from multiple bams

Version: 20150728

Descriptions:
Extract subset of sequences from BAMs

Options:
  -h    Print this help message
  -i    BAM list file, 1 file/line
  -b	comma delimited BAM files
  -n    number of BAM files
  -f    Sequence ID file, 1 ID/line
  -s    Comma-delimited sequence IDs
  -o    Output
  -d    Debug mode, keep temporary files.

Example:
  $0  [Options]

  BAMlist      1 BAM file per line
  sequenceids  comma-delimited seqids or ID file
               1 ID perl line
  BAMout       Final BAM output name

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
#################### Defaults ########################################
debug=0
bamlist=''
bamline=''
idfile=''
idlist=''
bamout=''
numbam=0

while [ -n "$1" ]; do
  case "$1" in
    -h) help;shift 1;;
    -i) bamlist=$2;shift 2;;
    -b) bamline=$2;shift 2;;
    -n) numbam=$2;shift 2;;
    -f) idfile=$2;shift 2;;
    -s) idlist=$2;shift 2;;
    -o) bamout=$2;shift 2;;
    -d) debug=1;shift 1;;
    --) shift;break;;
    -*) echo "error: no such option $1. -h for help" > /dev/stderr;exit 1;;
    *) break;;
  esac
done



#################### Subfuctions #####################################
###Detect command existence
CmdExists () {
	if command -v $1 >/dev/null 2>&1; then
		echo 0
	else
		echo 1
	fi
}

# Read the file in parameter and fill the array named "array"
#getArray() {
#    local i=0
#    local -a array=()
#    while read line # Read a line
#	for line in `cat $1`; do
#		echo $line
#        array[i]=$line # Put it into the array
#        i=$(($i + 1))
#    done
#    echo ${#array[@]} >&2  ###For test###
#    echo ${array[*]}
#}

# split string into Array
splitString () {
	local -a array=()
	local str2split=$1
	local str2delimiter=$2
	array=(`echo $str2split | sed -e "s/$str2delimiter/\n/g"`)
	echo ${array[@]}
}



######################check cmd#######################################
if [ $(CmdExists samtools) -eq 1 ]; then
	echo "Error: CMD 'samtools' in PROGRAM 'samtools' is required but not found.  Aborting..." >&2
	exit 127
fi
if [ $(CmdExists bamverify) -eq 1 ]; then
	echo "Error: script 'bamverify' is required but not found.  Aborting..." >&2
	exit 127
fi



#################### Input and Output ################################
###read bam files
test_bam=0;
if [ -z "$bamlist" ] || [ ! -s $bamlist ] || [[ "$bamlist" =~ ^\s+$ ]]; then
	echo "Warning: not use bamlist" >&2
else
	test_bam=1
	bamfile1=($(cat $bamlist))
fi
if [ -z "$bamline" ] || [[ "$bamline" =~ ^\s+$ ]]; then
	echo "Warning: not use delimited BAMs" >&2
else
	test_bam=1
	bamfile2=($(echo "$bamline" | perl -ne '@arr=split(/,/); map {print $_."\n"} @arr;'))
fi
if [ $test_bam -eq 0 ]; then
	echo "Error: invalid bamline or bamlist" >&2
	exit 1
fi
declare -a bamcombined=("${bamfile1[@]}" "${bamfile2[@]}")
declare -a bamfiles=($(echo "${bamcombined[@]}" | perl -ne '$line=$_;@arr=split(/\s+/,$line);map {$file{$_}++ if (-s $_)} @arr; foreach (keys %file){print $_."\n";}'))
if [ ${#bamfiles[@]} -lt 1 ]; then
	echo "Error: Empty bam lists" >&2
	exit 1
elif [ ${#bamfiles[@]} -ne $numbam ]; then
	echo "Error: numbam $numbam != number of bam files ${#bamfiles[@]}" >&2
	exit 1
fi
if [ $debug -ne 0 ]; then
    echo "######BAMfiles (Total: ${#bamfiles[@]}) detected:"
    for indbam in ${bamfiles[@]}; do
		if [ -s $indbam ]; then
			echo "Info: BAM detected: $indbam"
		else
		    echo "Error: Invalid BAM file $indbam" >&2
		    exit 1
		fi
	done
	echo -e "\n\n\n"
fi



###read ID names to be extracted
declare -a arrayid1=()
if [ ! -z "$idfile" ] && [ -s $idfile ]; then
	if [ $debug -ne 0 ]; then
			echo "Info: ID list in a file"
	fi
	arrayid1=($(cat $idfile))
fi
declare -a arrayid2=()
if [ ! -z "$idlist" ]; then
	if [ $debug -ne 0 ]; then
			echo "Info: ID list in a list"
	fi
#	if [ -z "${idlist##*,*}" ]; then	
	if [[ "$idlist" =~ ',' ]]; then
		if [ $debug -ne 0 ]; then
			echo "Info: ID list is multiple seqIDs"
		fi
		arrayid2=(`echo $idlist | sed -e "s/,/\n/g"`)
	else
		if [ $debug -ne 0 ]; then
			echo "Info: ID list is a single seqID"
		fi
		arrayid2[0]=$idlist
	fi
fi
declare -a combined=()
declare -a idnames=()
combined=("${arrayid1[@]}" "${arrayid2[@]}")
idnames=$(echo "${combined[@]}" | tr ' ' '\n' | sort -u | tr '\n' ' ')
if [ ${#idnames[@]} -lt 1 ]; then
	echo "Error: empty seq ID list" >&2
	exit 1
fi
declare mergeid=''
if [ $debug -ne 0 ]; then
	echo "###Read IDs to be extracted"
fi
for indseqid in ${idnames[@]}; do
	if [ $debug -ne 0 ]; then
		echo "$indseqid"
	fi
	mergeid="$mergeid $indseqid"
done
if [ $debug -ne 0 ]; then
	echo "Merged IDs: $mergeid"
	echo -e "\n\n\n"
fi


### detect BAM output
if [ -e $bamout ]; then
	echo "Warnings: BAM output exists" >&2
	rm $bamout
fi
if [ $debug -ne 0 ]; then
	echo -e "###Bamoutput\n$bamout\n\n\n"
fi



#################### Main ###########################################
if [ ${#bamfiles[@]} -eq 1 ]; then
	echo ${bamfiles[0]}
	if [ ! -s ${bamfiles[0]}.bai ]; then
		samtools index ${bamfiles[0]}
		if [ $? -ne 0 ] || [ ! -s ${bamfiles[0]}.bai ]; then
			echo "Error: Samtools index: ${bamfiles[0]}" >&2
			exit 1
		fi
	fi
	samtools view -bh ${bamfiles[0]} $mergeid > $bamout 2> /dev/null
	if [ $? -ne 0 ] || [ ! -s $bamout ]; then
		echo "Error: Samtools extract: ${bamfiles[0]} $mergeid > $bamout" >&2
		exit 1
	fi
elif [ ${#bamfiles[@]} -ge 2 ]; then
	tempfiles=''
	declare -i j=0
	for indvbam in ${bamfiles[@]}; do
		if [ ! -s $indvbam.bai ]; then
			samtools index $indvbam > /dev/null 2>&1
			if [ $? -ne 0 ] || [ ! -s $indvbam.bai ]; then
				echo "Error: samtools index failed: $indvbam" >&2
				exit 1
			fi
		fi
		if [[ "$indvbam" =~ \.[sS][aA][mM]$  ]]; then
			samtools view -bSh $indvbam $mergeid > temp_extractbam.$j.bam 2>/dev/null
		elif [[ "$indvbam" =~ \.[bB][aA][mM]$  ]]; then
			samtools view -bh $indvbam $mergeid > temp_extractbam.$j.bam 2>/dev/null
		else
			echo "Error: can not guess SAM/BAM format: $indvbam" >&2
			exit 1
		fi
		tempfiles="$tempfiles temp_extractbam.$j.bam"
		if [ $? -ne 0 ] || [ ! -s temp_extractbam.$j.bam ]; then
			echo "Error: samtools view extract: $indvbam" >&2
			rm $tempfiles > /dev/null 2>&1
			exit 1
		fi
		j=$(($j + 1))
	done
	echo "Info: Merging $j BAM files"
	samtools merge $bamout $tempfiles > /dev/null 2>&1
	if [ $? -ne 0 ] || [ ! -s $bamout ]; then
		echo "Error: samtools merge $tempfiles > $bamout"
		exit 1
	fi
	if [ $debug -eq 0 ]; then
		echo "Deleting temporary files: $tempfiles"
		rm $tempfiles > /dev/null 2>&1
	fi
fi
samtools sort $bamout $bamout.st
if [ $? -ne 0 ] || [ ! -s "$bamout.st.bam" ]; then
	echo "Error: samtools sort $bamout $bamout.st" >&2
	exit 1
fi
rm $bamout
mv $bamout.st.bam $bamout
samtools index $bamout
if [ $? -ne 0 ] || [ ! -s $bamout ]; then
	echo "Error: samtools sort $tempfiles > $bamout" >&2
	exit 1
fi
exit 0
