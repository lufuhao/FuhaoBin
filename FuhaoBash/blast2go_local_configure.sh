#!/bin/bash
RootDir=$(cd `dirname $(readlink -f $0)`; pwd)
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

################# help message ######################################
help() {
cat<<HELP

$0 --- Brief Introduction

Version: 20160730

Requirements:
	[Internet]
	Linux: wget

Descriptions:
	xxx

Options:
  -h    Print this help message
  -o    Database path to setup, default: ./
  -t    Number of threads, default: 1
  -s    Not run simulation
  -a    Not run assembly

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
databasepath=$PWD
opt_s=0
opt_a=0
opt_t=1
#################### Parameters #####################################
while [ -n "$1" ]; do
  case "$1" in
    -h) help;shift 1;;
    -i) opt_i=$2;shift 2;;
    -t) opt_t=$2;shift 2;;
    -s) opt_s=1;shift 1;;
    -a) opt_a=1;shift 1;;
    --) shift;break;;
    -*) echo "error: no such option $1. -h for help" > /dev/stderr;exit 1;;
    *) break;;
  esac
done


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

#Usage: string=$(join delimiter array)
join () {
        local separator=$1
        shift 1
        local -a array=(`echo $@`)
        local returnstr=$(printf "$separator%s" "${array[@]}")
        returnstr=${returnstr:1}
        echo $returnstr
}
abs2rel () { perl -MFile::Spec -e 'print(File::Spec->abs2rel($ARGV[1], $ARGV[0]), "\n")' "$@"; }


#################### Command test ###################################
if [ $(CmdExists 'santools') -eq 0 ]; then
	exit 0
else
	echo "Error: CMD/script 'samtools' in PROGRAM 'SAMtools' is required but not found.  Aborting..." >&2 
	exit 127
fi
if [ $(CmdExists 'santools') -eq 1 ]; then
	echo "Error: CMD/script 'samtools' in PROGRAM 'SAMtools' is required but not found.  Aborting..." >&2 
	exit 127
fi



#################### Defaults #######################################




#################### Input and Output ###############################

#Create path for database setup
if [ ! -d "$databasepath" ]; then
	mkdir -p "$databasepath"
fi



#################### Main ###########################################



#Step1. 下载uniprotKB（http://www.uniprot.org/）的swiss-prot和trEMBL
# swiss-prot:有报道过的蛋白质序列数据库，
# trEMBL: 翻译自EMBL核酸序列的蛋白质数据库，是未经报道过的。
# 在这三个FTP服务器的current_release目录下总是放置最新的uniprot版本，可从任一FTP下载
# Delete old
if [ -e "$databasepath/uniprot_sprot.fasta.gz" ]; then
	rm "$databasepath/uniprot_sprot.fasta.gz"
fi
if [ -e "$databasepath/uniprot_sprot.fasta" ]; then
	rm "$databasepath/uniprot_sprot.fasta"
fi
if [ -e "$databasepath/uniprot_trembl.fasta.gz" ]; then
	rm "$databasepath/uniprot_trembl.fasta.gz"
fi
if [ -e "$databasepath/uniprot_trembl.fasta" ]; then
	rm "$databasepath/uniprot_trembl.fasta"
fi
#1.1 swiss-prot
#1.1.1 uniprot FTP
wget -P "$databasepath" ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz
#1.1.2 EBI FTP
#wget -P "$databasepath" ftp://ftp.ebi.ac.uk/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz
#1.1.3 Expasy FTP
#wget -P "$databasepath" ftp://ftp.expasy.org/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz
if [ $? -ne 0 ] || [ ! -s "$databasepath/uniprot_sprot.fasta.gz" ]; then
	echo "Error: Step1 downloading swiss-prot error: uniprot_sprot.fasta.gz" >&2
	exit 1
fi
#1.2 trEMBL
#1.2.1 uniprot FTP
wget -P "$databasepath" ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_trembl.fasta.gz
#1.2.2 EBI FTP
#wget -P "$databasepath" ftp://ftp.ebi.ac.uk/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_trembl.fasta.gz
#1.2.3 Expasy FTP
#wget -P "$databasepath" ftp://ftp.expasy.org/databases/uniprot/current_release/knowledgebase/complete/uniprot_trembl.fasta.gz
if [ $? -ne 0 ] || [ ! -s "$databasepath/uniprot_trembl.fasta.gz" ]; then
	echo "Error: Step1 downloading trEMBL error: uniprot_trembl.fasta.gz" >&2
	exit 1
fi


#Step2. 下载uniref数据库
# 这是一个对uniprot蛋白序列进行聚类之后的数据库，去除了冗余数据，也常用于基因注释
# uniref100、uniref90和uniref50分别以100%，90%以及50%的相似度进行聚类。
#2.1 uniref100
if [ -e "$databasepath/uniref100.fasta.gz" ]; then
	rm "$databasepath/uniref100.fasta.gz"
fi
if [ -e "$databasepath/uniref100.fasta" ]; then
	rm "$databasepath/uniref100.fasta"
fi
wget -P "$databasepath" ftp://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref100/uniref100.fasta.gz
if [ $? -ne 0 ] || [ ! -s "$databasepath/uniref100.fasta.gz" ]; then
	echo "Error: Step2.1 downloading uniref100 error: uniref100.fasta.gz" >&2
	exit 1
fi
#2.2 uniref90
if [ -e "$databasepath/uniref90.fasta.gz" ]; then
	rm "$databasepath/uniref90.fasta.gz"
fi
if [ -e "$databasepath/uniref90.fasta" ]; then
	rm "$databasepath/uniref90.fasta"
fi
wget -P "$databasepath" ftp://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref90/uniref90.fasta.gz
if [ $? -ne 0 ] || [ ! -s "$databasepath/uniref90.fasta.gz" ]; then
	echo "Error: Step2.2 downloading uniref90 error: uniref90.fasta.gz" >&2
	exit 1
fi
#2.3 uniref50
if [ -e "$databasepath/uniref50.fasta.gz" ]; then
	rm "$databasepath/uniref50.fasta.gz"
fi
if [ -e "$databasepath/uniref50.fasta" ]; then
	rm "$databasepath/uniref50.fasta"
fi
wget -P "$databasepath" ftp://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref50/uniref50.fasta.gz
if [ $? -ne 0 ] || [ ! -s "$databasepath/uniref50.fasta.gz" ]; then
	echo "Error: Step2.3 downloading uniref50 error: uniref50.fasta.gz" >&2
	exit 1
fi



#Step3. uniprot的GOA（GO Gene Association）
# GOA是将uniprot蛋白质序列和GO ID关联的一个文件。GOA Slim则是前者的简要版本。通过该文件可以寻找到每个uniprot序列的GO注释。
#最新的GOA版本一般总会放在这两个地址下，是相对不变的。
if [ -e "$databasepath/gene_association.goa_uniprot.gz" ]; then
	rm "$databasepath/gene_association.goa_uniprot.gz"
fi
if [ -e "$databasepath/gene_association.goa_uniprot" ]; then
	rm "$databasepath/gene_association.goa_uniprot"
fi
#3.1 
wget -P "$databasepath" ftp://ftp.geneontology.org/pub/go/gene-associations/gene_association.goa_uniprot.gz
#3.2
#wget -P "$databasepath" ftp://ftp.ebi.ac.uk/pub/databases/GO/goa/UNIPROT/gene_association.goa_uniprot.gz
if [ $? -ne 0 ] || [ ! -s "$databasepath/gene_association.goa_uniprot.gz" ]; then
	echo "Error: Step3 downloading GOA error: gene_association.goa_uniprot.gz" >&2
	exit 1
fi


#Step4. uniprot GOA（GO Gene Association） GOA Slim版本
# goslim.map是关联GOA和GOA slim的文件
#4.1 GOA Slim
if [ -e "$databasepath/gene_association.goa_uniprot_slimmed.gz" ]; then
	rm "$databasepath/gene_association.goa_uniprot_slimmed.gz"
fi
if [ -e "$databasepath/gene_association.goa_uniprot_slimmed" ]; then
	rm "$databasepath/gene_association.goa_uniprot_slimmed"
fi
wget -P "$databasepath" ftp://ftp.ebi.ac.uk/pub/databases/GO/goa/goslim/gene_association.goa_uniprot_slimmed.gz
if [ $? -ne 0 ] || [ ! -s "$databasepath/gene_association.goa_uniprot.gz" ]; then
	echo "Error: Step4 downloading GOA-Sim error: gene_association.goa_uniprot.gz" >&2
	exit 1
fi
#4.2 GOA-Sim map
if [ -e "$databasepath/goaslim.map" ]; then
	rm "$databasepath/goaslim.map"
fi
wget -P "$databasepath" ftp://ftp.ebi.ac.uk/pub/databases/GO/goa/goslim/goaslim.map
if [ $? -ne 0 ] || [ ! -s "$databasepath/goaslim.map" ]; then
	echo "Error: Step4 downloading GOA-Sim-map error: goaslim.map" >&2
	exit 1
fi



#Step5 
# GO的OBO文件，可从GO的FTP上进行下载：
if [ -e "$databasepath/gene_ontology.obo" ]; then
	rm "$databasepath/gene_ontology.obo"
fi
wget -P "$databasepath" ftp://ftp.geneontology.org/pub/go/ontology/gene_ontology.obo
if [ $? -ne 0 ] || [ ! -s "$databasepath/gene_ontology.obo" ]; then
	echo "Error: Step5 downloading GO OBO error: gene_ontology.obo" >&2
	exit 1
fi



#首先，将预测基因在本地和uniprot数据库进行BLASTP比对，获得最相似的uniprot同源基因；然后，在uniprot GOA或者GOA Slim文件中提取相关uniprot基因的GO ID，一个Uniprot基因通常对应多个GO ID；之后，将关联的GO ID赋给预测基因，由此得到预测基因的GO功能注释；最后，可将预测基因的GO功能注释转换为wego所读取的格式，由wego进行GO功能分类，及作图
#将青枯雷尔氏菌RS91和RS98的GeneMark预测基因的蛋白质序列和本地的swiss-Prot数据库进行blastp比对，期望值e值设为1e-10，排除一些不相干的匹配项，取最佳匹配项，即只需输出一个最相似的匹配结果。
#因为之后需要对blastp的比对结果进行过滤，需要比较匹配片段的比例，原blastp的tab格式，没有query的长度，因此，在之前进行比对的时候，参数-outfmt采用默认值，即输出pair格式的blastp结果（可称之为bout文件）。
#通过perl编程处理bout文件，得到一个简要的列表文件（简而记之曰blst文件，这是笔者自己所定义的格式，信息更全面，类似于blast -fmtout 6所产生的列表文件）。针对blst文件进行blastp结果的过滤，过滤设定如下：相似性在30%以上，匹配片段占比对双方较短的序列的60%以上。这里分别用到了blastp-lst00.pl和filter-blst.pl两个perl脚本。
#RS91和RS98匹配到的uniprot基因在GOA文件中会有对应的若干个GO ID，这些GO ID是对该uniprot基因的注释，通过GO ID又可以从GO的OBO文件中提取具体的GO Term术语，即功能注释。因此，从Uniprot的GOA文件中提取GO ID，转而赋给匹配的RS91或RS98的预测基因，通过Perl编程read_GOA.pl脚本实现这一步骤。
#接着，可以对read_GOA.pl的结果文件进行统计。

# 配置以下7行
godbname=go_201307-assocdb-data   # 根据http://archive.geneontology.org/latest-full/下assocdb-data.gz文件更改
dbname=b2gdb                      # 数据库 名称，不用改
dbuser=root                       # 数据库 用户名
dbpass=passwordofroot             # 数据库 用户密码
dbhost=localhost                  # 数据库 所在ip
dbport=3306                       # 数据库 端口，3306是默认的，如果是无root权限安装的MySQL，一定要改为设置的端口，比如我的33060
path=/home/shenwei/Public/Data/local_b2g # 数据文件目录，注意路径末尾不要有“/”

# 如果已经下载数据文件，下列部分保持注释
### Download the GO database the NCBI mapping files and the PIR mapping
# wget http://archive.geneontology.org/latest-full/$godbname.gz
# wget ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/gene_info.gz
# wget ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2accession.gz
# wget ftp://ftp.pir.georgetown.edu/databases/idmapping/idmapping.tb.gz

# 如果已经下载并解压数据文件，下列部分保持注释
###unzip files
# gzip -dv $godbname.gz
# gzip -dv gene_info.gz
# gzip -dv gene2accession.gz
# gzip -dv idmapping.tb.gz


echo 1. Create the DB Tables and user
mysql -h$dbhost -P$dbport -u$dbuser -p$dbpass < b2gdb.sql

### Import data to the GO Database
echo 2. Import $godbname
mysql -h$dbhost -P$dbport -u$dbuser -p$dbpass $dbname < $godbname

echo 3. Import gene2accession
mysql -h$dbhost -P$dbport -u$dbuser -p$dbpass $dbname -e"LOAD DATA LOCAL INFILE '$path"/gene2accession"' INTO TABLE gene2accession FIELDS TERMINATED BY '\t' LINES TERMINATED BY '\n';"

echo 4. Import gene_info
mysql -h$dbhost -P$dbport -u$dbuser -p$dbpass $dbname -e"LOAD DATA LOCAL INFILE '$path"/gene_info"' INTO TABLE gene_info FIELDS TERMINATED BY '\t' LINES TERMINATED BY '\n';"

echo 5. Import idmapping.tb
java -cp .:mysql-connector-java-5.0.8-bin.jar: ImportIdMapping $path/idmapping.tb $dbhost:$dbport $dbname blast2go blast4it
echo All data imported.
