#!/bin/bash
### Source: http://data.biostarhandbook.com/scripts/wonderdump.sh
set -ue
### Exit if command fails
#set -o errexit
### Set readonly variable
#readonly passwd_file=”/etc/passwd”
### exit when variable undefined
#set -o nounset
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

$0 --- download SRR using prefetch in SRA-toolkit

Version: v20200723

Requirements:
	SRA-toolkit

Descriptions:
  Download SRA files directly using SRA-toolkit
    should be better than wget or curl

Options: 
  -h    Print this help message
  -f    SRR list file
  -i    SRR names, comma delimited
  -d    Output Path
  -n    Options for prefetch
          Default:  --max-size 50G

Example:

  $0 -f SRR_list.file -i SRRxxxxxx,SRRyyyyyyy -d $HOME/puclic/SRA

Author:
  Fu-Hao Lu
  Professor, PhD
  State Key Labortory of Crop Stress Adaptation and Improvement
  College of Life Science
  Jinming Campus, Henan University
  Kaifeng 475004, P.R.China
  E-mail: lufuhao@henu.edu.cn
HELP
exit 2
}
[ $# -lt 1 ] && help
[ "$1" = "-h" ] || [ "$1" = "--help" ] && help
#################### Environments ###################################
echo -e "\n######################\nProgram $ProgramName initializing ...\n######################\n"
#echo "Adding $RunDir/bin into PATH"
#export PATH=$RunDir/bin:$RunDir/utils/bin:$PATH

#################### Initializing ###################################
declare -a SRRlist=()
opt_f=""
opt_d=$PWD
opt_n=" --max-size 50G "
#################### Parameters #####################################
while [ -n "$1" ]; do
  case "$1" in
    -h) help;shift 1;;
    -i) SRRlist=($(echo $2 | tr "," "\n"));shift 2;;
    -f) opt_f=$(readlink -f $2);shift 2;;
    -d) opt_d=$2;shift 2;;
    -n) opt_n=$2;shift 2;;
    --) shift;break;;
    -*) echo "Error: no such option $1. -h for help" > /dev/stderr;exit 1;;
    *) break;;
  esac
done


#################### Subfuctions ####################################
###Detect command existence
CmdExists () {
  if command -v $1 >/dev/null 2>&1; then
    return 0
  else
    return 1
  fi
}




SrrDownloadUsingAscp() {
	local $SDUCsrr=$1
	
	local SDUCsub="(SrrDownloadUsingCurl)"
	# Create the full path to the file.
	local Sra_File="$opt_d/$SDUCsrr.sra"
	local Tmp_File="$opt_d/$SDUCsrr.tmp"
	
	echo "($SDUCsub)Info: Getting SRR run: $SDUCsrr"

	cd $opt_d
	# Download only if it does not exist.
	if [ ! -f $Sra_File ]; then
#		PATH1=${SDUCsrr:0:6}
#		PATH2=${SDUCsrr:0:10}
#		URL="anonftp@ftp-private.ncbi.nlm.nih.gov:/sra/sra-instant/reads/ByRun/sra/SRR/${PATH1}/${PATH2}/${SDUCsrr}.sra"
#		echo "($SDUCsub)Info: Downloading: $URL"
		echo "($SDUCsub)Info: Saving to: $Sra_File"
		prefetch $opt_n -O $opt_d --output-file $Tmp_File $SDUCsrr
		if [ $? -ne 0 ]; then
			echo "($SDUCsub)Error: ascp failed to download SRR: $SDUCsrr" >&2
			echo "($SDUCsub)Error: CMD used: ascp $opt_n $URL $opt_d" >&2
			return 100;
		else # Move to local file only if successful.
			mv $Tmp_File $Sra_File
			if [ -e $Sra_File.curl.err ]; then
				rm $Sra_File.curl.err >/dev/null 2>&1
			fi
		fi
	else
		echo "($SDUCsub)Warnings: existing SRA file found: $Sra_File"
	fi

	return 0
}
###ascp的用法：ascp [参数] 目标文件 目标地址，在线文档
#-v verbose mode 唠叨模式，能让你实时知道程序在干啥，方便查错。有些作者的程序缺乏人性化，运行之后，只见光标闪，压根不知道运行到哪了
#-T 取消加密，否则有时候数据下载不了
#-i 提供私钥文件的地址，地址一般是$HOME/.aspera/connect/etc/asperaweb_id_dsa.openssh文件
#-l 设置最大传输速度，一般200m到500m，如果不设置，反而速度会比较低，可能有个较低的默认值
#-k 断点续传，一般设置为值1
#-Q 不懂，一般加上它
#-P 提供SSH port，一般是33001，反正我不懂

### SRA数据库下载
#首先记住，数据的存放地址是ftp-private.ncbi.nlm.nih.gov，SRA在Aspera的用户名是anonftp，下载举例：
#如果我想下载SRR949627.sra文件，首先我需要找到地址，去ncbi ftp-private或者ncbi faspftp，一层层寻找，直至找到，然后记下链接地址，就可以开始下载了：
#ascp -v -i ~/.aspera/connect/etc/asperaweb_id_dsa.openssh -k 1 -T -l200m anonftp@ftp-private.ncbi.nlm.nih.gov:/sra/sra-instant/reads/ByRun/sra/SRR/SRR949/SRR949627/SRR949627.sra ~/biostar/aspera/
#    注意：anonftp@ftp-private.ncbi.nlm.nih.gov后面是:号，不是路径/！
#    一般来说，NCBI的sra文件前面的地址都是一样的/sra/sra-instant/reads/ByRun/sra/SRR/...，那么写脚本批量下载也就不难了！



#################### Command test ###################################
CmdExists 'ascp'
if [ $? -ne 0 ]; then
	echo "Error: CMD 'ascp' in PROGRAM 'Aspera Connect' is required but not found.  Aborting..." >&2 
	exit 127
fi



#################### Defaults #######################################



#################### Input and Output ###############################
# Make the directory if it does not exist.
if [ ! -d $opt_d ]; then
	mkdir -p $opt_d
fi



#################### Main ###########################################

cd $opt_d
if [ ! -z "$opt_f" ] && [ -s "$opt_f" ]; then
	echo ""
	echo "### Using SRR list file: $opt_f"
	echo ""
	while read SrrID; do
		echo "###     SRR: $SrrID"
		if SrrDownloadUsingAscp $SrrID; then
			echo "###     SRR: $SrrID    OK"
		else
			echo "###     SRR: $SrrID    failed" >&2
		fi
	done < $opt_f
fi

if [[ ${SRRlist[@]} -gt 0 ]]; then
	echo ""
	echo "### Using SRR comma list"
	echo ""
	for SrrID in ${SRRlist[@]}; do
		echo "###     SRR: $SrrID"
		if SrrDownloadUsingAscp $SrrID; then
			echo "###     SRR: $SrrID    OK"
		else
			echo "###     SRR: $SrrID    failed" >&2
		fi
	done
fi

exit 0
