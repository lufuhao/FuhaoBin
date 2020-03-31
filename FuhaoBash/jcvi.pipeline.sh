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
      python3
    seqkit (https://github.com/shenwei356/seqkit)
    jcvi (https://github.com/tanghaibao/jcvi)
      Biopython (http://www.biopython.org/)
      numpy (http://numpy.scipy.org/)
      matplotlib (http://matplotlib.org/)

Descriptions:
  This script calculate gene synteny and draw the synteny plot
    Source file needed
      1. Protein/CDS/cDNA multifasta file for both query and subject
      2. GFF3 file
    Steps
      1. Convert GFF3 to BED format
        4column: chr[tab]start[tab]end[tab]ID
          ID should be the same ith multifasta seqID
      2. run JCVI to get the raw synteny/anchors
        JCVI use LAST for similarity search
      3. Synteny screen to get high quality anchors
      4. Make synteny plot
      5. Filter 1:1 match for furthr analyais

Options:
  -h    Print this help message
  -i1   Query fasta file in fa.gz or fa format
  -i2   Subject fasta file in fa.gz or fa format
  -g1   Query genome annotation in GFF format
  -g2   Subject genome annotation in GFF format
  -p1   Query Prefix
  -p2   Subject Prefix
  -ul   Use longest transcript for each geneIDs both query and subject
  -ul1  Use longest transcript for each query geneIDs
  -ul2  Use longest transcript for each subject geneIDs 
  -d    Output directory, default: .
  -mt   Molecular type: cds/pep, default:pep
  -type GFF feature to be extracted, default: mRNA
  -key  GFF feature to seqID: ID/Name, default: ID
  -cc   csscore cutoff for jcvi.compara.catalog ortholog, default: 0.7
  -nsn  Do not strip alternative splicing (e.g. At5g06540.1 ->
          At5g06540) [default: False]
          --no_strip_names jcvi.compara.catalog ortholog
  -mp   minspan for jcvi.compara.synteny screen, default: 0
  -mz   minsize for jcvi.compara.synteny screen, default: 0
  -clean Clean Temporary files:
            * lastdb files



Example:
  $0 \\
      -i1 ./sp1.fa.gz -i2 ./sp2.fa.gz -g1 ./sp1.gff.gz -g2 ./sp2.gff.gz \\
      -p1 sp1 -p2 sp2 -fmt pep -type mRNA -key ID -mp 30 -mz 10

Author:
  Fu-Hao Lu
  Professor, PhD
  State Key Labortory of Crop Stress Adaptation and Improvement
  College of Life Science
  Jinming Campus, Henan University
  Kaifeng 475004, P.R.China
  E-mail: lufuhao@henu.edu.cn
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
opt_cc=0.7
opt_nsn=0
opt_mp=0
opt_mz=0
opt_clean=0
opt_useLongest=0;
opt_useQueryL=0;
opt_useSubjectL=0;

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
    -ul) opt_useLongest=1;shift;;
    -ul1) opt_useQueryL=1;shift;;
    -ul2) opt_useSubjectL=1;shift;;
    -d) opt_d=$2;shift 2;;
    -type) opt_type=$2;shift 2;;
    -key)  opt_key=$2;shift 2;;
    -mt)   opt_mt=$2;shift 2;;
    -cc)   opt_cc=$2;shift 2;;
    -nsn)  opt_nsn=1;shift;;
    -mp)   opt_mp=$2;shift 2;;
    -mz)   opt_mz=$2;shift 2;;
    -clean) opt_clean=1; shift;;
    
    
#FastQR1Arr=($(echo $2 | tr  "\n"));shift 2;;
#    -t) opt_t=$2;shift 2;;
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
CmdExists 'lastal'
if [ $? -ne 0 ]; then
	echo "Error: CMD 'lastal' in PROGRAM 'LAST' (http://last.cbrc.jp/) is required but not found.  Aborting..." >&2 
	exit 127
fi
CmdExists 'lastdb'
if [ $? -ne 0 ]; then
	echo "Error: CMD 'lastdb' in PROGRAM 'LAST' (http://last.cbrc.jp/) is required but not found.  Aborting..." >&2 
	exit 127
fi
CmdExists 'get_the_longest_transcripts.py'
if [ $? -ne 0 ]; then
	echo "Error: script 'get_the_longest_transcripts.py' on https://github.com/xuzhougeng/myscripts is required but not found.  Aborting..." >&2 
	exit 127
fi
CmdExists 'seqkit'
if [ $? -ne 0 ]; then
	echo "Error: CMD 'seqkit' in PROGRAM 'seqkit' (https://github.com/shenwei356/seqkit) is required but not found.  Aborting..." >&2 
	exit 127
fi



#################### Defaults #######################################
step=0;
path_data="$opt_d/1.data"
path_ortholog="$opt_d/2.ortholog"
path_synteny="$opt_d/3.synteny"
path_plot="$opt_d/4.plot"
path_1to1="$opt_d/5.1to1"

mt_type=''
if [[ $opt_fmt =~ ^[pP][eE][pP]$ ]];then
	opt_fmt='pep'
	mt_type='prot'
elif [[ $opt_fmt =~ ^[cC][dD][sS]$ ]]; then
	opt_fmt='cds'
	mt_type='nucl'
else
	echo "Error: unknown molecular type" >&2
	exit 100
fi
if [ $opt_useLongest -eq 1 ]; then
	echo "Info: -ul was set; both query and subject longest genes would b used"
	opt_useQueryL=1;
	opt_useSubjectL=1;
fi



#################### Input and Output ###############################



#################### Main ###########################################
echo "################# Parameter check ####################"
echo "    -i1    $opt_i1"
echo "    -i2    $opt_i2"
echo "    -g1    $opt_g1"
echo "    -g2    $opt_g2"
echo "    -p1    $opt_p1"
echo "    -p2    $opt_p2"
echo "    -ul    $opt_useLongest"
echo "    -ul1   $opt_useQueryL"
echo "    -ul2   $opt_useSubjectL"
echo "    -d     $opt_d"
echo "    -type  $opt_type"
echo "    -key   $opt_key"
echo "    -mt    $opt_mt"
echo "    -cc    $opt_cc"
echo "    -nsn   $opt_nsn"
echo "    -mp    $opt_mp"
echo "    -mz    $opt_mz"
echo "    -clean $opt_clean"
echo -e "\n\n\n"





### Step1: prepare data
((step++));
if [ ! -d $path_data ]; then
	mkdir -p $path_data
else
	echo "Warnings: Step$step: existing path $path_data; Delete this folder if you have new data" >&2
fi
cd $path_data

#
#zcat Athaliana_167_TAIR10.gene.gff3.gz | python ~/scripts/python/get_the_longest_transcripts.py  > ath_lst_gene.txt
#sed -i 's/\.v2\.1//g' aly_lst_gene.txt
#sed -i 's/\.TAIR10//g' ath_lst_gene.txt
#
#
# prepare bed
echo "Step${step}: Preparing Data"; echo "Step${step}: Preparing Data" >&2;
if [ ! -s $path_data/$opt_p1.bed ]; then
	if [ $opt_useLongest -eq 0 ]; then
		python2 -m jcvi.formats.gff bed --type=$opt_type --key=$opt_key $opt_g1 > $path_data/$opt_p1.bed
		if [ $? -ne 0 ] || [ ! -s $path_data/$opt_p1.bed ]; then
			echo "Error: Step${step}: Failed to convert query GFF to BED" >&2
			echo "CMD used: python2 -m jcvi.formats.gff bed --type=$opt_type --key=$opt_key $opt_g1 > $path_data/$opt_p1.bed"
			exit 100;
		fi
	else
		zcat  $opt_g1 | get_the_longest_transcripts.py | cut -f 2 > $path_data/$opt_p1.longest.genes
		python2 -m jcvi.formats.gff bed --type=$opt_type --key=$opt_key $opt_g1 > $path_data/$opt_p1.all.bed
		grep -f $path_data/$opt_p1.longest.genes $path_data/$opt_p1.all.bed > $path_data/$opt_p1.bed
	fi
else
	echo "Warnings: Step${step}: using existing Query BED: $path_data/$opt_p1.bed" >&2
fi
echo "#Step${step} Info: top 5 lines pf query BED"
head -n 5 $path_data/$opt_p1.bed

if [ ! -s $path_data/$opt_p2.bed ]; then
	if [ $opt_useLongest -eq 0 ]; then
		python2 -m jcvi.formats.gff bed --type=$opt_type --key=$opt_key $opt_g2 > $path_data/$opt_p2.bed
		if [ $? -ne 0 ] || [ ! -s $path_data/$opt_p2.bed ]; then
			echo "Error: Step${step}: Failed to convert subject GFF to BED" >&2
			echo "CMD used: python2 -m jcvi.formats.gff bed --type=$opt_type --key=$opt_key $opt_g2 > $path_data/$opt_p2.bed"
			exit 100;
		fi
	else
		zcat  $opt_g2 | get_the_longest_transcripts.py > $path_data/$opt_p2.longest.genes
		python2 -m jcvi.formats.gff bed --type=$opt_type --key=$opt_key $opt_g2 > $path_data/$opt_p2.all.bed
		grep -f $path_data/$opt_p2.longest.genes $path_data/$opt_p2.all.bed > $path_data/$opt_p2.bed
	fi
else
	echo "Warnings: Step${step}: using existing Subject BED: $path_data/$opt_p2.bed" >&2
fi
echo "#Step${step} Info: top 5 lines pf subject BED"
head -n 5 $path_data/$opt_p2.bed
echo -e "\n\n\n"

# prepare fasta
if [ ! -s $path_data/$opt_p1.$opt_fmt ]; then
	if [ $opt_useQueryL -eq 0 ]; then
		if [[ $opt_i1 =~ ^.*\.[gG][zZ]$ ]]; then
			echo "#Step${step} Info: Query in fasta.gz format, decompressing"
			gunzip -c $opt_i1 > $path_data/$opt_p1.$opt_fmt
		else
			echo "#Step${step} Info: Query in fasta format, Linking"
			ln -sf $opt_i1 $path_data/$opt_p1.$opt_fmt
		fi
	else
		echo "#Step${step} Info: Extracting longest Query transcript IDs"
		if [ ! -s $path_data/$opt_p1.longest.genes ]; then
			zcat  $opt_g1 | get_the_longest_transcripts.py | cut -f 2 > $path_data/$opt_p1.longest.genes
		fi
		seqkit grep -f $path_data/$opt_p1.longest.genes $opt_i1 > $path_data/$opt_p1.$opt_fmt
	fi
else
	echo "Warnings: Step${step}: use existing file: $path_data/$opt_p1.$opt_fmt; Delete this if you have new data" >&2
fi
if [ ! -s $path_data/$opt_p2.$opt_fmt ]; then
	if [ $opt_useSubjectL -eq 0 ]; then
		if [[ $opt_i2 =~ ^.*\.[gG][zZ]$ ]]; then
			echo "#Step${step} Info: Subject in fasta.gz format, decompressing"
			gunzip -c $opt_i2 > $path_data/$opt_p2.$opt_fmt
		else
			echo "#Step${step} Info: Subject in fasta format, Linking"
			ln -sf $opt_i2 $path_data/$opt_p2.$opt_fmt
		fi
	else
		echo "#Step${step} Info: Extracting longest Subject transcript IDs"
		if [ ! -s $path_data/$opt_p2.longest.genes ]; then
			zcat  $opt_g2 | get_the_longest_transcripts.py | cut -f 2 > $path_data/$opt_p2.longest.genes
		fi
		seqkit grep -f $path_data/$opt_p2.longest.genes $opt_i2 > $path_data/$opt_p2.$opt_fmt
	fi
else
	echo "Warnings: Step${step}: use existing file: $path_data/$opt_p2.$opt_fmt; Delete this if you have new data" >&2
fi
echo -e "\n\n\n"



### Step2: run JCVI
((step++));
if [ ! -d $path_ortholog ]; then
	mkdir -p $path_ortholog
else
	echo "Warnings: Step${step}: existing path $path_ortholog; Delete this folder if you have new data" >&2
fi
cd $path_ortholog
echo "Step${step}: run jcvi.compara.catalog ortholog"; echo "Step${step}: run jcvi.compara.catalog ortholog" >&2;
if [ ! -s $path_ortholog/$opt_p1.$opt_p2.anchors ]; then
	ln -sf $path_data/$opt_p1.bed $path_ortholog/$opt_p1.bed
	ln -sf $path_data/$opt_p2.bed $path_ortholog/$opt_p2.bed
	ln -sf $path_data/$opt_p1.$opt_fmt $path_ortholog/$opt_p1.$opt_fmt
	ln -sf $path_data/$opt_p2.$opt_fmt $path_ortholog/$opt_p2.$opt_fmt
	
	CatalogOptions=""
	if [ $opt_nsn -eq 1 ]; then
		CatalogOptions="--no_strip_names"
	fi
	python -m jcvi.compara.catalog ortholog $opt_p1 $opt_p2 $CatalogOptions --dbtype $mt_type --cscore=$opt_cc > $opt_p1.$opt_p2.compara.catalog.ortholog.log 2>&1
	if [ $? -ne 0 ] || [ ! -s $path_ortholog/$opt_p1.$opt_p2.anchors ]; then
		echo "Error: Step${step}: jcvi.compara.catalog ortholog running error" >&2
		echo "CMD used: python -m jcvi.compara.catalog ortholog $opt_p1 $opt_p2 $CatalogOptions --dbtype $mt_type --cscore=$opt_cc" >&2
		exit 100
	fi
fi
if [ $opt_clean -eq 1 ]; then
	rm *.bck *.des *.prj *.sds *.ssp *.suf *.tis
fi

#Pairwise synteny visualization using dot plot.
if [ ! -s $opt_p1.$opt_p2.anchors.dotplot.pdf ]; then
	python -m jcvi.graphics.dotplot $opt_p1.$opt_p2.anchors -o $opt_p1.$opt_p2.anchors.dotplot.pdf --dpi=600 --format=pdf --font=Arial > $opt_p1.$opt_p2.graphics.dotplot.log 2>&1
fi

#We could also quick test if the synteny pattern is indeed 1:1, by running:
if [ ! -s $opt_p1.$opt_p2.depth.pdf ]; then
	python -m jcvi.compara.synteny depth --histogram $opt_p1.$opt_p2.anchors > $opt_p1.$opt_p2.compara.synteny.depth.log 2>&1
fi
echo -e "\n\n\n"



### Step3. synteny
((step++));
if [ ! -d $path_synteny ]; then
	mkdir -p $path_synteny
else
	echo "Warnings: Step${step}: existing path $path_synteny; Delete this folder if you have new data" >&2
fi
cd $path_synteny
echo "Step${step}: run jcvi.compara.synteny screen"; echo "Step${step}: run jcvi.compara.synteny screen" >&2;
ln -sf $path_ortholog/$opt_p1.$opt_p2.anchors $path_synteny/$opt_p1.$opt_p2.anchors
ln -sf $path_data/$opt_p1.bed $path_synteny/$opt_p1.bed
ln -sf $path_data/$opt_p2.bed $path_synteny/$opt_p2.bed
if [ ! -s $opt_p1.$opt_p2.anchors.new ]; then
	python -m jcvi.compara.synteny screen --minspan=$opt_mp --simple --minsize=$opt_mz $opt_p1.$opt_p2.anchors $opt_p1.$opt_p2.anchors.new > $opt_p1.$opt_p2.compara.synteny.screen.log 2>&1
	if [ $? -ne 0 ] || [ ! -s $path_synteny/$opt_p1.$opt_p2.anchors.new ] || [ ! -s $path_synteny/$opt_p1.$opt_p2.anchors.simple ]; then
		echo "Error: Step${step}: jcvi.compara.synteny screen running error" >&2
		echo "CMD used: python -m jcvi.compara.synteny screen --minspan=$opt_mp --simple --minsize=$opt_mz $opt_p1.$opt_p2.anchors $opt_p1.$opt_p2.anchors.new" >&2
		exit 100
	fi
fi
echo -e "\n\n\n"



### Step4. plot
((step++));
if [ ! -d $path_plot ]; then
	mkdir -p $path_plot
else
	echo "Warnings: Step${step}: existing path $path_plot; Delete this folder if you have new data" >&2
fi
cd $path_plot
echo "Step${step}: run jcvi.graphics.karyotype"; echo "Step${step}: run jcvi.graphics.karyotype" >&2;


ln -sf $path_data/$opt_p1.bed $path_plot/$opt_p1.bed
ln -sf $path_data/$opt_p2.bed $path_plot/$opt_p2.bed
ln -sf $path_synteny/$opt_p1.$opt_p2.anchors.simple $path_plot/

#seqids
if [ ! -s $path_plot/$opt_p1.$opt_p2.seqids ]; then
	cut -f 1 $opt_p1.bed | sort -u | tr "\n" "," | sed 's/,$/\n/' > $path_plot/$opt_p1.$opt_p2.seqids
	if [ $? -ne 0 ] || [ ! -s $path_plot/$opt_p1.$opt_p2.seqids ]; then
		echo "Error: Step${step}: Failed to collect query seqIDs for $opt_p1" >&2
		exit 100
	fi
	cut -f 1 $opt_p2.bed | sort -u | tr "\n" "," | sed 's/,$/\n/' >> $path_plot/$opt_p1.$opt_p2.seqids
	if [ $? -ne 0 ] || [ ! -s $path_plot/$opt_p1.$opt_p2.seqids ]; then
		echo "Error: Step${step}: Failed to collect subject seqIDs for $opt_p2" >&2
		exit 100
	fi
fi

#layout
#The whole canvas is 0-1 on x-axis and 0-1 on y-axis. First, three columns specify the position of the track. Then rotation, color, label, vertical alignment (va), and then the genome BED file. Track 0 is now grape, track 1 is now peach. The next stanza specifies what edges to draw between the tracks. e, 0, 1 asks to draw edges between track 0 and 1, using information from the .simple file.
if [ ! -s $path_plot/$opt_p1.$opt_p2.layout ]; then
	echo "# y, xstart, xend, rotation, color, label, va, bed, label_va" > $path_plot/$opt_p1.$opt_p2.layout
	echo ".6,     .1,    .8,       0,      m, $opt_p1, top, $opt_p1.bed, center" >> $path_plot/$opt_p1.$opt_p2.layout
	echo ".4,     .1,    .8,       0,      k, $opt_p2, bottom, $opt_p2.bed, center" >> $path_plot/$opt_p1.$opt_p2.layout
	echo "# edges" >> $path_plot/$opt_p1.$opt_p2.layout
	echo "e, 0, 1, $opt_p1.$opt_p2.anchors.simple" >> $path_plot/$opt_p1.$opt_p2.layout
fi


#utils/webcolors.py
#perl -i -lane '$i="";$j=""; $line=$_; $F[0]=~s/^.*\*//; if ($F[0]=~/$TraesCS(\d+)[ABD]\d+G\d+$/) {$i=$1;}else{print STDERR "Error: no match1";} if ($F[2]=~/$TraesCS(\d+)[ABD]\d+G\d+$/) {$j=$1;}else{print STDERR "Error: no match2";} print STDERR "Info: i: $i; j : $j"; if ($i eq $j) {$line="red*".$line;}else{$line="black*".$line;} print $line;' aa.bb.anchors.simple

if [ ! -s $path_plot/$opt_p1.$opt_p2.pdf ]; then
	python -m jcvi.graphics.karyotype --dpi=600 --format=pdf --font=Arial $path_plot/$opt_p1.$opt_p2.seqids $path_plot/$opt_p1.$opt_p2.layout > $opt_p1.$opt_p2.graphics.karyotype.log 2>&1
	if [ $? -ne 0 ] || [ ! -s karyotype.pdf ]; then
		echo "Error: Step${step}: Failed to run jcvi.graphics.karyotype for seqids: $path_plot/$opt_p1.$opt_p2.seqids layout:$path_plot/$opt_p1.$opt_p2.layout" >&2
		echo "CMD used: python -m jcvi.graphics.karyotype --dpi=600 --format=pdf --font=Arial $path_plot/$opt_p1.$opt_p2.seqids $path_plot/$opt_p1.$opt_p2.layout"
		exit 100
	fi
	mv karyotype.pdf $path_plot/$opt_p1.$opt_p2.pdf
fi
echo -e "\n\n\n"



### Step5. 1:1
((step++));
if [ ! -d $path_1to1 ]; then
	mkdir -p $path_1to1
else
	echo "Warnings: Step${step}: existing path $path_1to1; Delete this folder if you have new data" >&2
fi
cd $path_1to1
echo "Step${step}: run  1:1 synteny"; echo "Step${step}: run 1:1 synteny" >&2;

#ln -sf $path_synteny/$opt_p1.$opt_p2.anchors.new $path_1to1/
SimpleFile=$path_synteny/$opt_p1.$opt_p2.anchors.new

if [ ! -s $opt_p1.$opt_p2.duplicated.list ]; then
	grep -v ^'#' $SimpleFile | cut -f 1 | sort | uniq -d > $opt_p1.$opt_p2.duplicated.list
	grep -v ^'#' $SimpleFile | cut -f 2 | sort | uniq -d >> $opt_p1.$opt_p2.duplicated.list
fi
if [ -s $opt_p1.$opt_p2.duplicated.list ]; then
	if [ ! -s $path_1to1/$opt_p1.$opt_p2.anchors.new.uniq ]; then
		grep -v -f $opt_p1.$opt_p2.duplicated.list $SimpleFile > $path_1to1/$opt_p1.$opt_p2.anchors.new.uniq
	fi
	LineNum1=$(grep -v ^'#' $SimpleFile | wc -l)
	LineNum2=$(cat $opt_p1.$opt_p2.duplicated.list | wc -l)
	LineNum3=$(grep -v ^'#' $path_1to1/$opt_p1.$opt_p2.anchors.new.uniq | wc -l)
	echo "Info: geneIDs with duplication: $LineNum1"
	echo "Info: duplication IDs         : $LineNum2"
	echo "Info: geneIDs uniq            : $LineNum3"
	echo "Info: no duplicated geneIDs detected, final uniq synteny: $path_1to1/$opt_p1.$opt_p2.anchors.new.uniq"
else
	grep -v ^'#' $SimpleFile | wc -l
	cat $opt_p1.$opt_p2.duplicated.list | wc -l
	echo "Info: no duplicated geneIDs detected, final uniq synteny: $SimpleFile"
fi

#list.merger.pl Final.AABBDD.uniq "undef" final.uniq dd.aa.anchors.new.uniq dd.bb.anchors.new.uniq

exit 0
