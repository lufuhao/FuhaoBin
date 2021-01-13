#!/bin/bash
#!/bin/bash
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

$0 --- Brief Introduction

Version: v20210105

Requirements:
	Aspera
	perl && File::Spec
	taxonkit
	Linux: gzip, tar
	csvtk

Descriptions:
	Extract subset frm NCBI NR database

Options:
  -h    Print this help message
  -i    NR.gz file
  -d    Output path
  -td   NCBI /pub/taxonomy/taxdump.tar.gz file
  -pat  NCBI /pub/taxonomy/accession2taxid/prot.accession2taxid.gz file
  -aid  default: $HOME/.aspera/connect/etc/asperaweb_id_dsa.openssh


Example:
  $0 -i ./chr1.fa -t 10

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
opt_d=$PWD
opt_i=""
opt_aid="$HOME/.aspera/connect/etc/asperaweb_id_dsa.openssh"
opt_taxdump=""
opt_t=1
#################### Parameters #####################################
while [ -n "$1" ]; do
  case "$1" in
    -h) help;shift 1;;
#    -i) FastQR1Arr=($(echo $2 | tr  "\n"));shift 2;;
    -i) opt_i=$2;shift 2;;
    -d) opt_d=$2;shift 2;;
    -td) opt_taxdump=$2;shift 2;;
    -pat) opt_pat=$2; shift 2;;
    -aid) opt_aid=$2;shift 2;;
    -t) opt_t=$2;shift 2;;
#    -1) seq_rfn=(${seq_rfn[@]} "$2");shift 2;;
#    -s) opt_s=1;shift 1;;
#    -a) opt_a=1;shift 1;;
    --) shift;break;;
    -*) echo "error: no such option $1. -h for help" > /dev/stderr;exit 1;;
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





#################### Command test ###################################
CmdExists 'samtools'
if [ $? -ne 0 ]; then
	echo "Error: CMD/script 'samtools' in PROGRAM 'SAMtools' is required but not found.  Aborting..." >&2 
	exit 127
fi
if [[ $(CmdExists 'mum.stat') -eq 1 ]]; then
	echo "Error: script 'mum.stat' is required but not found.  Aborting..." >&2 
	exit 127
fi


#################### Defaults #######################################




#################### Input and Output ###############################
if [ ! -d $opt_d ]; then
	mkdir -p $opt_d
fi
if [ ! -s "$HOME/.aspera/connect/etc/asperaweb_id_dsa.openssh" ]; then
	echo "Error: not found: $HOME/.aspera/connect/etc/asperaweb_id_dsa.openssh" >&2
	exit 100
fi



#################### Main ###########################################
if [ $? -ne 0 ] || [ ! -s $gffout ]; then
	echo "GFFSORT_Error: sort error" >&2
	echo "CMD used: bedGraphToBigWig $opt_bg $opt_fai $opt_o" >&2
	exit 100
fi


### Step1: 从NCBI下载Nr数据库和分类数据库文件
opt1database="$opt_d/1.database"
if [ "$opt_i" == "" ]; then
	cd $opt_d
# 下载Nr数据库（FASTA文件）
	if [ ! -d $opt1database ]; then
		mkdir -p $opt1database
	fi
		ascp -T -l 200M -i $opt_aid --host=ftp.ncbi.nih.gov --user=anonftp --mode=recv /blast/db/FASTA/nr.gz $opt1database
		if [ $? -ne 0 ] || [ ! -s "$opt1database/nr.gz" ]; then
			echo "Error: ascp downloading failed: nr.gz" >&2
			echo "CMD used: ascp -T -l 200M -i $HOME/.aspera/connect/etc/asperaweb_id_dsa.openssh --host=ftp.ncbi.nih.gov --user=anonftp --mode=recv /blast/db/FASTA/nr.gz $opt1database" >&2
			exit 100
		fi
		opt_i="$opt1database/nr.gz"
	fi
elif [ -s "$opt_i" ]; then
	echo "Info: using existing NR.gz: $opt_i"
else
	echo "Error: invalid NR.gz: $opt_i" >&2
	exit 100
fi
# 下载NCBI的分类数据库文件
if [ "$opt_taxdump" == "" ]; then
	cd $opt_d
# 下载taxdump数据库
	if [ ! -d $opt1database ]; then
		mkdir -p $opt1database
	fi
		ascp -T -l 200M -i $opt_aid --host=ftp.ncbi.nih.gov --user=anonftp --mode=recv /pub/taxonomy/taxdump.tar.gz $opt1database
		if [ $? -ne 0 ] || [ ! -s "$opt1database/taxdump.tar.gz" ]; then
			echo "Error: ascp downloading failed: nr.gz" >&2
			echo "CMD used: ascp -T -l 200M -i $opt_aid --host=ftp.ncbi.nih.gov --user=anonftp --mode=recv /pub/taxonomy/taxdump.tar.gz $opt1database" >&2
			exit 100
		fi
		opt_taxdump="$opt1database/taxdump.tar.gz"
	fi
elif [ -s "$opt_taxdump" ]; then
	echo "Info: using existing taxdump.tar.gz: $opt_taxdump"
else
	echo "Error: invalid taxdump.tar.gz: $opt_taxdump" >&2
	exit 100
fi
if [ "$opt_pat" == "" ]; then
	cd $opt_d
# 下载taxdump数据库
	if [ ! -d $opt1database ]; then
		mkdir -p $opt1database
	fi
	ascp -T -l 200M -i $opt_aid --host=ftp.ncbi.nih.gov --user=anonftp --mode=recv /pub/taxonomy/accession2taxid/prot.accession2taxid.gz $opt1database
	if [ $? -ne 0 ] || [ ! -s "$opt1database/prot.accession2taxid.gz" ]; then
		echo "Error: ascp downloading failed: nr.gz" >&2
		echo "CMD used: ascp -T -l 200M -i $opt_aid --host=ftp.ncbi.nih.gov --user=anonftp --mode=recv /pub/taxonomy/accession2taxid/prot.accession2taxid.gz $opt1database" >&2
		exit 100
	fi
	opt_pat="$opt1database/prot.accession2taxid.gz"
elif [ -s "$opt_pat" ]; then
	echo "Info: using existing prot.accession2taxid.gz: $opt_pat"
else
	echo "Error: invalid prot.accession2taxid.gz: $opt_pat" >&2
	exit 100
fi



### Step2: decomprass
opt2decompress="$opt_d/2.decompress"
opt2protaccession2taxid="$opt2decompress/prot.accession2taxid"
if [ ! -d $opt2decompress ]; then
	mkdir -p $opt2decompress
fi

# 解压缩两个NCBI的分类数据库文件
gzip -dc $opt_pat > $opt2protaccession2taxid
if [ ! -d $HOME/.taxonkit ]; then
	mkdir -p $HOME/.taxonkit
fi
tar zxf taxdump.tar.gz -C $HOME/.taxonkit
### 其主要有效文件有两个：
### names.dmp 记录物种名及其分类编号
### nodes.dmp 记录分类编号的节点信息
### 查看~/.taxonkit/names.dmp文件，使用关键词检索得到目标类的分类编号，例如：
### fungi 4751             # grep -P "\|\s+[fF]ungi\w*\s*\|" ~/.taxonkit/names.dmp
### plants 3193            # grep -P "\|\s+[pP]lant\w*\s*\|" ~/.taxonkit/names.dmp
### animals 33208          # grep -P "\|\s+[aA]nimal\w*\s*\|" ~/.taxonkit/names.dmp



### Step3: 下载并安装NCBI分类数据库解析软件TaxonKitTaxonKit，并解析nodes.dmp文件的物种节点信息，得到指定类的所有物种列表信息。
opt2taxid="$opt_d/3.taxid"
if [ ! -d $opt2taxid ]; then
	mkdir -p $opt2taxid
fi
# 下载并安装NCBI分类数据库解析软件TaxonKit
#wget https://github.com/shenwei356/taxonkit/releases/download/v0.2.4/taxonkit_linux_amd64.tar.gz
#tar zxvf taxonkit_linux_amd64.tar.gz
# 提取古菌(2157)、细菌(2)和病毒(10239) 4751 fungi这几个大类下的所有物种编号。
taxonkit list --ids 2,2157,10239 --threads $opt_t -o $opt2taxid/sub.meta.list

# 再编写程序extract_sub_data_from_Nr.pl获得列表中物种在Nr数据库中的序列信息。
gzip -dc nr.gz | extract_sub_data_from_Nr.pl --sub_taxon $opt2taxid/sub.meta.list --acc2taxid $opt2decompress/prot.accession2taxid - > $opt2taxid nr_meta.fasta



#Listing all taxids below $id using taxonkit.
id=6656
# 6656 is the phylum Arthropoda
# echo 6656 | taxonkit lineage | taxonkit reformat
# 6656    cellular organisms;Eukaryota;Opisthokonta;Metazoa;Eumetazoa;Bilateria;Protostomia;Ecdysozoa;Panarthropoda;Arthropoda    Eukaryota;Arthropoda;;;;;
# 2    bacteria
# 2157 archaea
# 4751 fungi
# 10239 virus
# time: 3s
taxonkit list --ids $id --indent "" > $id.taxid.txt
# taxonkit list --ids 2,4751,10239 --indent "" > microbe.taxid.txt
wc -l $id.taxid.txt
# 518373 6656.taxid.txt
# time: 4min
#https://bioinf.shenwei.me/taxonkit/tutorial/
#Retrieving target accessions. There are two options:
pigz -dc prot.accession2taxid.gz \
    | csvtk grep -t -f taxid -P $id.taxid.txt \
    | csvtk cut -t -f accession.version,taxid \
    | sed 1d \
    > $id.acc2taxid.txt

cut -f 1 $id.acc2taxid.txt > $id.acc.txt

wc -l $id.acc.txt
# 8174609 6656.acc.txt


### Step3: 提取fungi/plants/animals子集
taxonkit list --ids 4751 --threads $opt_t > sub.fungi.list
taxonkit list --ids 3193 --threads $opt_t > sub.plants.list
taxonkit list --ids 33208 --threads $opt_t > sub.animals.list
taxonkit list --ids 10239 --threads $opt_t > sub.virus.list
gzip -dc nr.gz | perl extract_sub_data_from_Nr.pl --sub_taxon sub.fungi.list --acc2taxid prot.accession2taxid - > nr_fungi.fasta
gzip -dc nr.gz | perl extract_sub_data_from_Nr.pl --sub_taxon sub.plants.list --acc2taxid prot.accession2taxid - > nr_plants.fasta
gzip -dc nr.gz | perl extract_sub_data_from_Nr.pl --sub_taxon sub.animals.list --acc2taxid prot.accession2taxid - > nr_animals.fasta
gzip -dc nr.gz | perl extract_sub_data_from_Nr.pl --sub_taxon sub.virus.list --acc2taxid prot.accession2taxid - > nr_virus.fasta



### Step5: 使用makeblastdb创建blast本地数据库
makeblastdb -in nr_fungi.fasta -dbtype prot -title nr_fungi -parse_seqids -out nr_fungi_`date +%Y%m%d` -logfile nr_fungi_`date +%Y%m%d`.log
makeblastdb -in nr_plants.fasta -dbtype prot -title nr_plants -parse_seqids -out nr_plants_`date +%Y%m%d` -logfile nr_plants_`date +%Y%m%d`.log
makeblastdb -in nr_animals.fasta -dbtype prot -title nr_animals -parse_seqids -out nr_animals_`date +%Y%m%d` -logfile nr_animals_`date +%Y%m%d`.log
makeblastdb -in nr_virus.fasta -dbtype prot -title nr_virus -parse_seqids -out nr_virus_`date +%Y%m%d` -logfile nr_virus_`date +%Y%m%d`.log
