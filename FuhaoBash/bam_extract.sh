#!/bin/bash
RootDir=$(cd `dirname $(readlink -f $0)`; pwd)
if [ ! -z $(uname -m) ]; then
	machtype=$(uname -m)
elif [ ! -z "$MACHTYPE" ]; then
	machtype=$MACHTYPE
else
	echo "Warnings: unknown MACHTYPE" >&2
fi
abs2rel () { perl -MFile::Spec -e 'print(File::Spec->abs2rel($ARGV[1], $ARGV[0]), "\n")' "$@"; }
export NUM_THREADS=`grep -c '^processor' /proc/cpuinfo 2>/dev/null || sysctl -n hw.ncpu 2>/dev/null || echo 1`;
ProgramName=${0##*/}
echo "MachType: $machtype"
echo "RootPath: $RootDir"
echo "ProgName: $ProgramName"

################# help message ######################################
help() {
cat<<HELP

$0 --- Brief Introduction

Version: 20150628

Requirements:
	samtools

Descriptions:
	Extract BAM using region file OR bed file
	remove RG info and merge 4 bams

Options:
  -h    help
  -1    bamfile1
  -2    bamfile2
  -3    bamfile3
  -4    bamfile4
  -r    region2extract
  -b    bedfile file
  -e    reheader file
  -o1   outbamfile NOT reheadered
  -o2   outbam2 reheadered
  -m    outmerge
  -d    debug
  -g    runmerge

Example:
  $0 -i ./chr1.fa -t 10

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
echo -e "\n######################\nProgram initializing ...\n######################\n"
#echo "Adding $RunDir/bin into PATH"
#export PATH=$RunDir/bin:$RunDir/utils/bin:$PATH

#################### Initializing ###################################
region2extract=''
runrename=0
### AABBDD
bamfile1=''
bamfile2=''
bamfile3=''
bamfile4=''
outbamfile=''
outreheader=''
outmerge=''
bedfile=
packagename="$ProgramName"
debug=0
runmerge=0
#################### Parameters #####################################
while [ -n "$1" ]; do
  case "$1" in
    -h) help;shift 1;;
    -1) bamfile1=$2;shift 2;;
    -2) bamfile2=$2;shift 2;;
    -3) bamfile3=$2;shift 2;;
    -4) bamfile4=$2;shift 2;;
    -r) region2extract=$2;shift 2;;
    -b) bedfile=$2;shift 2;;
    -e) outreheader=$2;shift 2;;
    -o1) outbamfile=$2;shift 2;;
    -o2) outbam2==$2;shift 2;;
    -m) outmerge=$2;shift 2;;
    -d) debug=1;shift 1;;
    -g) runmerge=1;shift 1;;
    --) shift;break;;
    -*) echo "error: no such option $1. -h for help" > /dev/stderr;exit 1;;
    *) break;;
  esac
done



#################### Defaults #######################################
declare -a tempfiles=()
declare -a finalbam=()

if [ ! -z "$outmerge" ]; then
	echo "${packagename}Info: -m outmerge specified, auto enable -g runmerge"
	runmerge=1
fi
if [ $runmerge -eq 1 ]; then
	if [ -z "$bamfile1" ] || [ -z "$bamfile2" ] || [ -z "$bamfile3" ] || [ -z "$bamfile4" ] || [ -z "$outreheader" ] || [ -z "$outmerge" ]; then
		echo "${packagename}Error: runmerge needs --bamfile1/2/3/4 BAMs and -e header"
		exit 1
	fi
fi
if [ ! -z $outbam2 ]; then
	if [ -z "$bamfile1" ] || [ -z "$outreheader" ] || [ -z "$outbamfile" ]; then
		echo "${packagename}Error: -o2 reRG.bam needs --bamfile1 BAMs, -e header and -o1 ourbam1file"
		exit 1
	fi
fi
if [ -z "$region2extract" ]; then
	if [ -z "$bedfile" ]; then
		echo "${packagename}Error: Neither region or bed file not specified" >&2
		exit 1
	elif [ -s "$region2extract" ]; then
		cmd_samtools="samtools view -L $region2extract"
	fi
else
	if [ ! -z "$bedfile" ]; then
		echo "${packagename}Error: both region and bed file specified" >&2
		exit 1
	fi
fi


#################### Input and Output ###############################


#################### Subfuctions ####################################
###Detect command existence
CmdExists () {
  if command -v $1 >/dev/null 2>&1; then
    echo 0
  else
#    echo "I require $1 but it's not installed.  Aborting." >&2
    echo 1
  fi
#  local cmd=$1
#  if command -v $cmd >/dev/null 2>&1;then
#    echo >&2 $cmd "  :  "`command -v $cmd`
#    exit 0
#  else
#    echo >&2 "Error: require $cmd but it's not installed.  Exiting..."
#    exit 1
#  fi
}
###Usage: array=(`split delimiter string`)
split () {
	local separator=$1
	local mystring=$2
	echo $mystring | sed -e "s/$separator/\n/g"
}
#str='ni,hai,a'
#a=(`SplitString ',' $str`)
#echo ${#a[@]} ${a[0]} ${a[1]} ${a[2]}
#Usage: string=$(join delimiter array)
join () {
        local separator=$1
        shift 1
        local -a array=(`echo $@`)
        local returnstr=$(printf "$separator%s" "${array[@]}")
        returnstr=${returnstr:1}
        echo $returnstr
}
DeleteTemp () {
	local -a file2delete=("$@")
	for ind_file in ${file2delete[@]}; do
		if [ -e $ind_file ]; then
			rm $ind_file
		fi
	done
}
CleanCache () {
	local -a file2clean=("$@")
	for ind_file2 in ${file2clean[@]}; do
		if [ -e $ind_file2 ]; then
			rm $ind_file2
		fi
	done
	`DeleteTemp ${tempfile[@]}`
}
### ExtractBam (BAMin, region, header, BAMout)
### Return: 0=success; 1=fail; 2=empty;
ExtractBam () {
	local EBbam=$1
	local EBregion=$2
	local EBheader=$3
	local EBbamout=$4

	local EBshinfo="${packagename}::SH(ExtractBam)"
	
	local EBtestbam=$($cmd_samtools -F 4 $EBbam $EBregion | perl -e '$line=<>; if (defined $line and $line =~/\S+/) {print "0\n";exit 0;}else {print "1\n";exit 0}')
	if [ $EBtestbam -eq 1 ]; then
		echo 2 && exit 2
	else
		if [ -z "$EBheader" ]; then
			$cmd_samtools -b -h -F 4 $EBbam $EBregion > $EBbamout
		elif [ -s $EBheader ]; then
			$cmd_samtools -F 4 $EBbam $EBregion | cat $EBheader - | samtools view -b -S -h - > $EBbamout
		else
			echo 1 && exit 1
		fi
		
		if [ $? -ne 0 ] || [ ! -s $EBbamout ]; then
			echo "${EBshinfo}Error: extract error: BAM ($EBbam) Region($EBregion)" >&2
			if [ -e "$EBbamout" ]; then
				rm "$EBbamout"
			fi
			echo 1 && exit 1
		else
			echo 0 && exit 0
		fi
	fi
}



### MergeBam (BAMout bam1 bam2 bam3 ...)
### Return: 0=success; 1 = fail
MergeBam () {
	local MBbamout=$1
	shift 1
	local -a MBbams2merge=("$@")

	local MBshinfo="${packagename}::SH(MergeBam)"
	local MBind_file=''
	
	if [ -z "${MBbams2merge[@]}" ] || [ -z "$MBbamout" ] || [ ${#MBbams2merge[@]} -lt 1 ]; then
		echo 1 && exit 1
	elif [ ${#MBbams2merge[@]} -eq 1 ]; then
		mv ${MBbams2merge[0]} $MBbamout
		echo 0 && exit 0
	elif [ ${#MBbams2merge[@]} -gt 1 ]; then
		samtools merge $MBbamout "${MBbams2merge[@]}"
		if [ $? -ne 0 ] || [ ! -s $MBbamout ]; then
			echo "${MBshinfo}Error: samtools merge error: $MBbamout" >&2
			if [ -e $MBbamout ]; then
				rm $$MBbamout
			fi
			echo 1 && exit 1
		else
			for MBind_file in ${MBbams2merge[@]}; do
				if [ -e $MBind_file ]; then
					rm $MBind_file
				fi
			done
			echo 0 && exit 0
		fi
	else
		echo 1 && exit 1
	fi
}



#################### Command test ###################################
#echo "Error: CMD/script 'samtools' in PROGRAM 'SAMtools' is required but not found.  Aborting..." >&2 
if [ $(CmdExists 'samtools') -ne 0 ]; then
	echo "Error: CMD 'samtools' in package 'SAMtools' is required but not found.  Aborting..." >&2 
	exit 127
fi



#################### Main ###########################################
if [ ! -s $outbamfile ]; then
	###extract AABBDD bam
	declare -a bamarr1=($(echo "$bamfile1" | perl -e '$line=<>;@arr=split(/,/, $line); map {print $_."\n"} @arr;'))
	declare -a bamext=()
	tempi=0
	for ind1bam in ${bamarr1[@]}; do
		let tempi+=1
		if [ $debug -ne 0 ]; then
			echo "${packagename}Info: Extracting AABBDD BAM: $ind1bam ..."
		fi
		test_code=$(ExtractBam $ind1bam $region2extract ' ' $outbamfile.AABBDD.$tempi.bam)
		if [ $test_code -eq 0 ]; then
			bamext+=("$outbamfile.AABBDD.$tempi.bam")
			tempfiles+=("$outbamfile.AABBDD.$tempi.bam")
		elif [ $test_code -eq 1 ]; then
			echo "${packagename}Error: failed to extract BAM: $ind1bam ..." >&2
			`DeleteTemp "${tempfiles[@]}"`
			exit 1
		elif [ $test_code -eq 2 ]; then
			continue
		fi
	done
	if [ ${#bamext[@]} -lt 1 ]; then
		echo "${packagename}Error: AABBDD no BAM to merge" >&2
		`DeleteTemp "${tempfiles[@]}"`
		exit 1
	else
		test_code2=$(MergeBam "$outbamfile" "${bamext[@]}")
		if [ test_code2 -ne 0 ]; then
			echo "${packagename}Error: samtools AABBDD merge error" >&2
			`DeleteTemp "${tempfiles[@]}"`
			exit 1
		fi
	fi
	if [ $runrename -eq 1 ]; then
		mv $outbamfile $outbamfile.temp.bam
		samtools view -h $outbamfile.temp.bam | perl -ne 'if (/^@/){ print; next;}else{chomp; if (/\tRG:Z:(\S+)/) {$rg=$1;$_=s/^[a-zA-Z0-9-]+:[a-zA-Z0-9-]+:[a-zA-Z0-9-]+/$rg/;}print;}' | samtools view -b -h -S - > $outbamfile
		if [ $? -ne 0 ] || [ ! -s $outbamfile ]; then
			echo "${packagename}Error: samtools AABBDD rename error" >&2
			`DeleteTemp "${tempfiles[@]}"`
			exit 1
		fi
		tempfiles+=("$outbamfile.temp.bam")
	fi
	samtools index $outbamfile
	if [ $? -ne 0 ] || [ ! -s $outbamfile.bai ]; then
		samtools sort $outbamfile $outbamfile.st
		if [ $? -ne 0 ] || [ ! -s $outbamfile.st.bam ]; then
			echo "${packagename}Error: samtools AABBDD sort error" >&2
			`DeleteTemp "${tempfiles[@]}"`
			exit 1
		fi
		rm $outbamfile
		mv $outbamfile.st.bam $outbamfile
		samtools index $outbamfile
		if [ $? -ne 0 ] || [ ! -s $outbamfile.bai ]; then
			echo "${packagename}Error: samtools AABBDD index error" >&2
			`DeleteTemp "${tempfiles[@]}"`
			exit 1
		fi
		if [ $runmerge -eq 0 ]; then
			echo "${packagename}Info: Ignoring Final merge"
			exit 0
		fi
	fi
fi

if [ ! -z "$outbam2" ] || [ ! -z "$outmerge" ]; then
	samtools view $outbamfile | perl -ne 'if (/\tRG:Z:\S+/) {$_=~s/RG:Z:\S+/RG:Z:AABBDD/;print;}else {chomp; $_.="\tRG:Z:AABBDD\n";print;}' | cat $outreheader - | samtools view -b -h -S - > $outbamfile.rehead.bam
	if [ $? -ne 0 ] || [ ! -s $outmerge.AABBDD.bam ]; then
		echo "${packagename}Error: AABBDD replace RG error2" >&2
		`DeleteTemp "${tempfiles[@]}"`
		exit 1
	fi
	finalbam=("$outbamfile.rehead.bam")
	tempfiles+=("$outbamfile.rehead.bam")
	if [ -z "$outmerge" ]; then
		mv $outbamfile.rehead.bam $outbam2
		`DeleteTemp "${tempfiles[@]}"`
		exit 0
	fi
fi

### create AABB bam
declare -a bamarr2=($(echo "$bamfile2" | perl -e '$line=<>;@arr=split(/,/, $line); map {print $_."\n"} @arr;'))
bamext=()
tempi=0
test_AABB=0
for ind2bam in ${bamarr2[@]}; do
	let tempi+=1
	if [ $debug -ne 0 ]; then
		echo "${packagename}Info: Extracting AABB BAM: $ind2bam ..."
	fi
	test_code=$(ExtractBam $ind2bam $region2extract ' ' "$outmerge.AABB.$tempi.bam")
	if [ $test_code -eq 0 ]; then
		bamext+=("$outmerge.AABB.$tempi.bam")
		tempfiles+=("$outmerge.AABB.$tempi.bam")
	elif [ $test_code -eq 1 ]; then
		echo "${packagename}Error: failed to extract BAM: $ind2bam ..." >&2
		`DeleteTemp "${tempfiles[@]}"`
		exit 1
	elif [ $test_code -eq 2 ]; then
		continue
	fi
done
if [ ${#bamext[@]} -lt 1 ]; then
	echo "${packagename}Error: AABB no BAM to merge" >&2
	`DeleteTemp "${tempfiles[@]}"`
else
	test_code2=$(MergeBam "$outmerge.AABB.merge.bam" "${bamext[@]}")
	if [ test_code2 -ne 0 ]; then
		echo "${packagename}Error: samtools AABB merge error" >&2
		`DeleteTemp "${tempfiles[@]}"`
		exit 1
	else
		tempfiles+=("$outmerge.AABB.merge.bam")
		samtools view $outmerge.AABB.merge.bam | perl -ne 'if (/\tRG:Z:\S+/) {$_=~s/RG:Z:\S+/RG:Z:AABB/;print;}else {chomp; $_.="\tRG:Z:AABB\n";print;}' | cat $outreheader - | samtools view -b -h -S - > $outmerge.AABB.bam
		if [ $? -ne 0 ] || [ ! -s $outmerge.AABB.bam ]; then
			echo "${packagename}Error: AABB replace RG error" >&2
			exit 1
		else
			tempfiles+=("$outmerge.AABB.bam")
			finalbam+=("$outmerge.AABB.bam")
			test_AABB=1
		fi
	fi
fi




### create AA bam
declare -a bamarr3=($(echo "$bamfile3" | perl -e '$line=<>;@arr=split(/,/, $line); map {print $_."\n"} @arr;'))
bamext=()
tempi=0
test_AA=0
for ind3bam in ${bamarr3[@]}; do
	let tempi+=1
	if [ $debug -ne 0 ]; then
		echo "${packagename}Info: Extracting AA BAM: $ind3bam ..."
	fi
	test_code=$(ExtractBam $ind3bam $region2extract ' ' "$outmerge.AA.$tempi.bam")
	if [ $test_code -eq 0 ]; then
		bamext+=("$outmerge.AA.$tempi.bam")
		tempfiles+=("$outmerge.AA.$tempi.bam")
	elif [ $test_code -eq 1 ]; then
		echo "${packagename}Error: failed to extract BAM: $ind3bam ..." >&2
		`DeleteTemp "${tempfiles[@]}"`
		exit 1
	elif [ $test_code -eq 2 ]; then
		continue
	fi
done
if [ ${#bamext[@]} -lt 1 ]; then
	echo "${packagename}Error: AA no BAM to merge" >&2
	`DeleteTemp "${tempfiles[@]}"`
	exit 1
else
	test_code2=$(MergeBam "$outmerge.AA.merge.bam" "${bamext[@]}")
	if [ test_code2 -ne 0 ]; then
		echo "${packagename}Error: samtools AA merge error" >&2
		`DeleteTemp "${tempfiles[@]}"`
		exit 1
	else
		tempfiles+=("$outmerge.AA.merge.bam")
		samtools view $outmerge.AA.merge.bam | perl -ne 'if (/\tRG:Z:\S+/) {$_=~s/RG:Z:\S+/RG:Z:AA/;print;}else {chomp; $_.="\tRG:Z:AA\n";print;}' | cat $outreheader - | samtools view -b -h -S - > $outmerge.AA.bam
		if [ $? -ne 0 ] || [ ! -s $outmerge.AA.bam ]; then
			echo "${packagename}Error: AA replace RG error" >&2
			exit 1
		else
			tempfiles+=("$outmerge.AA.bam")
			finalbam+=("$outmerge.AA.bam")
			test_AA=1
		fi
	fi
fi



### create DD bam
declare -a bamarr4=($(echo "$bamfile4" | perl -e '$line=<>;@arr=split(/,/, $line); map {print $_."\n"} @arr;'))
bamext=()
tempi=0
test_DD=0
for ind4bam in ${bamarr4[@]}; do
	let tempi+=1
	if [ $debug -ne 0 ]; then
		echo "${packagename}Info: Extracting DD BAM: $ind4bam ..."
	fi
	test_code=$(ExtractBam $ind4bam $region2extract ' ' "$outmerge.DD.$tempi.bam")
	if [ $test_code -eq 0 ]; then
		bamext+=("$outmerge.DD.$tempi.bam")
		tempfiles+=("$outmerge.DD.$tempi.bam")
	elif [ $test_code -eq 1 ]; then
		echo "${packagename}Error: failed to extract BAM: $ind4bam ..." >&2
		`DeleteTemp "${tempfiles[@]}"`
		exit 1
	elif [ $test_code -eq 2 ]; then
		continue
	fi
done
if [ ${#bamext[@]} -lt 1 ]; then
	echo "${packagename}Error: DD no BAM to merge" >&2
	`DeleteTemp "${tempfiles[@]}"`
	exit 1
else
	test_code2=$(MergeBam "$outmerge.DD.merge.bam" "${bamext[@]}")
	if [ test_code2 -ne 0 ]; then
		echo "${packagename}Error: samtools DD merge error" >&2
		`DeleteTemp "${tempfiles[@]}"`
		exit 1
	else
		tempfiles+=("$outmerge.DD.merge.bam")
		samtools view $outmerge.DD.merge.bam | perl -ne 'if (/\tRG:Z:\S+/) {$_=~s/RG:Z:\S+/RG:Z:DD/;print;}else {chomp; $_.="\tRG:Z:DD\n";print;}' | cat $outreheader - | samtools view -b -h -S - > $outmerge.DD.bam
		if [ $? -ne 0 ] || [ ! -s $outmerge.DD.bam ]; then
			echo "${packagename}Error: DD replace RG error" >&2
			exit 1
		else
			tempfiles+=("$outmerge.DD.bam")
			finalbam+=("$outmerge.DD.bam")
			test_DD=1
		fi
	fi
fi



###create AABBDD bam
if [ $test_AABB -eq 0 ] || [ $test_AA -eq 0 ] || [ $test_DD -eq 0 ]; then
	if [ -z "$outbam2" ]; then
		mv ${finalbam[0]} $outmerge
	else
		cp ${finalbam[0]} $outmerge
	fi
	samtools index $outmerge
	if [ $? -ne 0 ] || [ ! -s $outmerge.bai ]; then
		echo "${packagename}Error: AABBDD noRG index error" >&2
		`DeleteTemp "${tempfiles[@]}"`
		exit 1
	fi
	exit 0
else
	test_code=$(MergeBAM $outmerge "${finalbam[@]}")
	if [ $test_code -ne 0 ]; then
		echo "${packagename}Error: Final Merge error" >&2
		`DeleteTemp "${tempfiles[@]}"`
		exit 1
	fi
	tempfiles+=("$outmerge")
fi
samtools index $outmerge
if [ $? -ne 0 ] || [ ! -s $outmerge.bai ]; then
	echo "${packagename}Error: Final merge noRG index error" >&2
	`DeleteTemp "${tempfiles[@]}"`
	exit 1
fi
exit 0
