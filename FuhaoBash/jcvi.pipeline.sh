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

Version: 20200316

Requirements:
	LAST (http://last.cbrc.jp/)
	get_the_longest_transcripts.py (https://github.com/xuzhougeng/myscripts)
	seqkit (https://github.com/shenwei356/seqkit)
#	perl && File::Spec

Descriptions:
	

Options:
  -h    Print this help message
  -i1   Query fasta file in fa.gz or fa format
  -i2   Subject fasta file in fa.gz or fa format
  -g1   Query genome annotation in GFF format
  -g2   Subject genome annotation in GFF format
  -p1   Query Prefix
  -p2   Subject Prefix
  -d    Output directory, default: .
  -fmt   Molecular type: cds/pep, default:pep
  -type GFF feature to be extracted, default: mRNA
  -key  GFF feature to seqID: ID/Name, default: ID
  -cs   csscore for 'jcvi.compara.catalog ortholog'
  
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
[ $# -lt 1 ] && help
[ "$1" = "-h" ] || [ "$1" = "--help" ] && help
#################### Environments ###################################
echo -e "\n######################\nProgram $ProgramName initializing ...\n######################\n"
#echo "Adding $RunDir/bin into PATH"
#export PATH=$RunDir/bin:$RunDir/utils/bin:$PATH

#################### Initializing ###################################
opt_d=$PWD;
opt_type="mRNA";
opt_key="ID";
opt_fmt="pep"
opt_cscore=0.7


opt_s=0
opt_a=0
opt_t=1
#################### Parameters #####################################
while [ -n "$1" ]; do
  case "$1" in
    -h) help;shift 1;;
    -i1) opt_i1=$2;shift 2;;
    -i2) opt_i2=$2;shift 2;;
    -g1) opt_g1=$2;shift 2;;
    -g2) opt_g2=$2;shift 2;;
    -p1) opt_p1=$2;shift 2;;
    -p2) opt_p2=$2;shift 2;;
    -d) opt_d=$2;shift 2;;
    -type) opt_type=$2;shift 2;;
    -key)  opt_key=$2;shift 2;;
    -mt)   opt_mt=$2;shift 2;;
    -cc)   opt_cscore=$2;shift 2;;
    
    
    
#FastQR1Arr=($(echo $2 | tr  "\n"));shift 2;;
    -t) opt_t=$2;shift 2;;
    -1) seq_rfn=(${seq_rfn[@]} "$2");shift 2;;
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
    return 0
  else
    return 1
  fi
}





#################### Command test ###################################
CmdExists 'lastal'
if [ $? -ne 0 ]; then
	echo "Error: CMD 'lastal' in PROGRAM 'LAST' is required but not found.  Aborting..." >&2 
	exit 127
fi
CmdExists 'lastdb'
if [ $? -ne 0 ]; then
	echo "Error: CMD 'lastdb' in PROGRAM 'LAST' is required but not found.  Aborting..." >&2 
	exit 127
fi




if [ $? -ne 0 ]; then
	echo "Error: CMD/script 'samtools' in PROGRAM 'SAMtools' is required but not found.  Aborting..." >&2 
	exit 127
fi
if [[ $(CmdExists 'mum.stat') -eq 1 ]]; then
	echo "Error: script 'mum.stat' is required but not found.  Aborting..." >&2 
	exit 127
fi


#################### Defaults #######################################
step=0;
path_data="$opt_d/0.data"
path_ortholog="$opt_d/1.ortholog"
path_synteny="$opt_d/2.synteny"

mt_type=''
if [[ $opt_fmt =~ ^[pP][eE][pP]$ ]];then
	opt_fmt='pep'
	mt_type='prot'
elif [[ $opt_fmt =~ ^[cC][dD][sS]$ ]]
	opt_fmt='cds'
	mt_type='nucl'
else
	echo "Error: unknown molecular type" >&2
	exit 100
fi


#################### Input and Output ###############################
opt_i1
opt_i2
opt_g1
opt_g2
opt_p1
opt_p2
opt_d=$PWD;
opt_type="mRNA";
opt_key="ID";



#################### Main ###########################################

### Step1: prepare data
((step++));
if [ ! -d $path_data ]; then
	mkdir -p $path_data
else
	echo "Warnings: Step$step: existing path $path_data; Delete this folder if you have new data" >&2
fi
cd $path_data

# prepare fasta
if [ ! -s $path_data/$opt_p1.$opt_fmt ]; then
	if [[ $opt_i1 =~ ^.*\.[gG][zZ]$ ]]; then
		echo "#Step${step} Info: Query in fasta.gz format, decompressing"
		gunzip -c $opt_i1 > $path_data/$opt_p1.$opt_fmt
	else
		echo "#Step${step} Info: Query in fasta format, Linking"
		ln -sf $opt_i1 $path_data/$opt_p1.$opt_fmt
	fi
else
	echo "Warnings: Step${step}: use existing path $path_data; Delete this folder if you have new data" >&2
fi
if [ ! -s $path_data/$opt_p2.fa ]; then
	if [[ $opt_i2 =~ ^.*\.[gG][zZ]$ ]]; then
		echo "#Step${step} Info: Subject in fasta.gz format, decompressing"
		gunzip -c $opt_i2 > $path_data/$opt_p2.$opt_fmt
	else
		echo "#Step${step} Info: Subject in fasta format, Linking"
		ln -sf $opt_i2 $path_data/$opt_p2.$opt_fmt
	fi
else
	
fi

# prepare bed
if [ ! -s $path_data/$opt_p1.bed ]; then
	python2 -m jcvi.formats.gff bed --type=$opt_type --key=$opt_key $opt_g1 > $path_data/$opt_p1.bed
	if [ $? -ne 0 ] || [ ! -s $path_data/$opt_p1.bed ]; then
		echo "Error: Step${step}: Failed to convert query GFF to BED" >&2
		echo "CMD used: python2 -m jcvi.formats.gff bed --type=$opt_type --key=$opt_key $opt_g1 > $path_data/$opt_p1.bed"
		exit 100;
	fi
else
	echo "Warnings: Step${step}: using existing Query BED: $path_data/$opt_p1.bed" >&2
fi
echo "#Step${step} Info: top 5 lines pf query BED"
head -n 5 $path_data/$opt_p1.bed
echo -e "\n\n\n"

if [ ! -s $path_data/$opt_p2.bed ]; then
	python2 -m jcvi.formats.gff bed --type=$opt_type --key=$opt_key $opt_g2 > $path_data/$opt_p2.bed
	if [ $? -ne 0 ] || [ ! -s $path_data/$opt_p2.bed ]; then
		echo "Error: Step${step}: Failed to convert subject GFF to BED" >&2
		echo "CMD used: python2 -m jcvi.formats.gff bed --type=$opt_type --key=$opt_key $opt_g2 > $path_data/$opt_p2.bed"
		exit 100;
	fi
else
	echo "Warnings: Step${step}: using existing Subject BED: $path_data/$opt_p2.bed" >&2
fi
echo "#Step${step} Info: top 5 lines pf subject BED"
head -n 5 $path_data/$opt_p2.bed
echo -e "\n\n\n"


#zcat Alyrata_384_v2.1.gene.gff3.gz | python ~/scripts/python/get_the_longest_transcripts.py > aly_lst_gene.txt
#zcat Athaliana_167_TAIR10.gene.gff3.gz | python ~/scripts/python/get_the_longest_transcripts.py  > ath_lst_gene.txt
#sed -i 's/\.v2\.1//g' aly_lst_gene.txt
#sed -i 's/\.TAIR10//g' ath_lst_gene.txt
#seqkit grep -f  <(cut -f 2 ath_lst_gene.txt ) Athaliana_167_TAIR10.cds.fa.gz > ath.cds
#seqkit grep -f  <(cut -f 2 aly_lst_gene.txt ) Alyrata_384_v2.1.cds.fa.gz > aly.cds




### Step2: run JCVI
((step++));
if [ ! -d $path_ortholog ]; then
	mkdir -p $path_ortholog
else
	echo "Warnings: Step${step}: existing path $path_ortholog; Delete this folder if you have new data" >&2
fi
cd $path_ortholog


ln -sf $path_data/$opt_p1.bed $path_ortholog/$opt_p1.bed
ln -sf $path_data/$opt_p2.bed $path_ortholog/$opt_p2.bed
ln -sf $path_data/$opt_p1.$opt_fmt $path_ortholog/$opt_p1.$opt_fmt
ln -sf $path_data/$opt_p2.$opt_fmt $path_ortholog/$opt_p2.$opt_fmt

python -m jcvi.compara.catalog ortholog $opt_p1 $opt_p2 --no_strip_names --dbtype $mt_type --cscore=$opt_cscore
if [ $? -ne 0 ] || [ ! -s $path_ortholog/$opt_p1.$opt_p2.anchors ]; then
	echo "Error: Step${step}: jcvi.compara.catalog ortholog running error" >&2
	echo "CMD used: python -m jcvi.compara.catalog ortholog $opt_p1 $opt_p2 --no_strip_names --dbtype $mt_type --cscore=$opt_cscore" >&2
	exit 100
fi
rm *.bck *.des *.prj *.sds *.ssp *.suf *.tis


#Pairwise synteny visualization using dot plot.
python -m jcvi.graphics.dotplot $opt_p1.$opt_p2.anchors -o $opt_p1.$opt_p2.anchors.dotplot.pdf --dpi=600 --format=pdf --font=Arial


#We could also quick test if the synteny pattern is indeed 1:1, by running:
python -m jcvi.compara.synteny depth --histogram $opt_p1.$opt_p2.anchors



### Step3. synteny
((step++));
if [ ! -d $path_synteny ]; then
	mkdir -p $path_synteny
else
	echo "Warnings: Step${step}: existing path $path_synteny; Delete this folder if you have new data" >&2
fi
cd $path_synteny

ln -sf $path_ortholog/$opt_p1.$opt_p2.anchors $path_synteny/$opt_p1.$opt_p2.anchors
python -m jcvi.compara.synteny screen --minspan=30 --simple --minsize=10 $opt_p1.$opt_p2.anchors $opt_p1.$opt_p2.anchors.new




### Step4. plot

#seqids
chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19
Pp01,Pp02,Pp03,Pp04,Pp05,Pp06,Pp07,Pp08


#layout file, which tells the plotter where to draw what. The whole canvas is 0-1 on x-axis and 0-1 on y-axis. First, three columns specify the position of the track. Then rotation, color, label, vertical alignment (va), and then the genome BED file. Track 0 is now grape, track 1 is now peach. The next stanza specifies what edges to draw between the tracks. e, 0, 1 asks to draw edges between track 0 and 1, using information from the .simple file.
# y, xstart, xend, rotation, color, label, va,  bed
 .6,     .1,    .8,       0,      , Grape, top, grape.bed
 .4,     .1,    .8,       0,      , Peach, top, peach.bed
# edges
e, 0, 1, grape.peach.anchors.simple


python -m jcvi.graphics.karyotype --dpi=600 --format=png --font Arial seqids layout



### Step5. 1:1





if [ $? -ne 0 ] || [ ! -s $gffout ]; then
	echo "GFFSORT_Error: sort error" >&2
	echo "CMD used: bedGraphToBigWig $opt_bg $opt_fai $opt_o" >&2
	exit 100
fi
