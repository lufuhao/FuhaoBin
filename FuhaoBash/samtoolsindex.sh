#!/bin/sh
RunDir=$(cd `dirname $(readlink -f $0)`; pwd)
MachType=$(uname -m)

################# help message ######################################
help()
{
cat<<HELP

Synposis
  samtoolsindex path|bam_files

Version: 20181101

Descriptions:
  samtools index all path(/*.sorted.bam) files
  Accept both path and files

    NOTE: BAM files must be sorted

Options:
  -h    Print this help message
  -n    Ignore indexed files

Example:
  ## index all BAMs (.bam) in folders
  samtoolsindex.sh BAM_folder1 BAM_folder2
  
  ### index multiple bam files
  samtoolsindex.sh in1.bam in2.bam 
  
  ### Mixed
  samtoolsindex.sh BAM_folder1 in2.bam

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
#################### Defaults #######################################

opt_g=0
#################### parameters #####################################
while [ -n "$1" ]; do
  case "$1" in
    -h) help;shift 1;;
    -n) opt_g=1;shift 1;;
    --) shift;break;;
    -*) echo "Error: no such option $1. -h for help";exit 100;;
    *) break;;
  esac
done

subsamtoolsindex() {
	local file_to_index=$1
	
	if [ -e "$file_to_index.bai" ]; then
		if [ $opt_g -eq 1 ]; then
			echo "    Info: $file_to_index already indexed" >&2
			return 0;
		else
			echo "    Info: $file_to_index already indexed, try to re-generate ..." >&2
			rm "$file_to_index.bai" >/dev/null 2>&1
		fi
	fi
	echo "    Indexing $file_to_index"
	samtools index $file_to_index >/dev/null 2>&1
	if [ $? -ne 0 ] || [ ! -s "$file_to_index.bai" ]; then
		echo "    Error: Samtools index error: $file_to_index" >&2
		return 1;
	fi
	
	return 0;
}

for indv_bam in $@; do
  if [ -d $indv_bam ]; then
    echo "$file is a directory. Indexing..."
    indv_path=`echo $indv_bam | sed "s/\/$//"`
    for indv2_bam in `ls $indv_path/*.bam`; do
      subsamtoolsindex $indv2_bam
    done
  elif [ -f $indv_bam ]; then
    echo "$indv_bam is a file. "
    subsamtoolsindex $indv_bam
  else
    echo "$indv_bam is not file or folder"
    exit 100;
  fi
done


exit 0
