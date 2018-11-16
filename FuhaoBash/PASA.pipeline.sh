#!/bin/bash
RunDir=$(cd `dirname $(readlink -f $0)`; pwd)
MachType=$(uname -m)

################# help message ######################################
help()
{
cat<<HELP

$0 --- Setup Config and MySQL, and then run PASA pipeline

Version: 20150528

Descriptions:
	1. make PASA config in PWD
	2. check mySQL exists or not
	3. run PASA pipeline

Requirements:
    Linux: cp, perl, MySQL
    PASA
        \$PASAHOME/scripts/Launch_PASA_pipeline.pl
        PASAHOME
        

Options:
  -h    Print this help message
  -t    Transcript file in fasta
  -g    Genome file in fasta
  -f    Full-length cDNA file in fasta
  -db   MySQL database name
  -ru   MySQL readonly username
  -rp   MwSQL readonly password
  -wu   MySQL readwrite username
  -wp   MySQL readwrite password
  -m    PASAHOME
  -c    Number of CPUs
  --runpasa 
  -a    Aligners: gmap|blat|gmap,blat|blat,gmap
  
  

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
[ "$1" = "-h" ] && help
[ "$1" = "--help" ] && help
#################### Defaults #######################################
echo -e "\n######################\nProgram initializing ...\n######################\n"
#echo "Adding $RunDir/bin into PATH"
#export PATH=$RunDir/bin:$RunDir/utils/bin:$PATH

###Defaults
transcript=''
genome=''
fulllength=''
mysqldatabase=''
readonly_user='readonly'
readonly_pawd='123456'
readwrite_user='luf'
readwrite_pawd='123456'
runpasa=0
pasacmd=' -C -R '
numcpus=1
aligners='gmap,blat'
#################### parameters #####################################
while [ -n "$1" ]; do
  case "$1" in
    -h) help;shift 1;;
    -t) transcript=$2;shift 2;;
    -g) genome=$2;shift 2;;
    -f) fulllength=$2;shift 2;;
    -db) mysqldatabase=$2;shift 2;;
    -ru) readonly_user=$2;shift 2;;
    -rp) readonly_pawd=$2;shift 2;;
    -wu) readwrite_user=$2;shift 2;;
    -wp) readwrite_pawd=$2;shift 2;;
    -m)	PASAHOME=$2; shift 2;;
    -c) numcpus=$2; shift 2;;
    --runpasa) runpasa=1; shift 1;;
    -a) aligners=$2; shift 2;;
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



#################### Command test ###################################
#requirements: perl mysql
#perl
if [ $(CmdExists 'perl') -eq 1 ]; then
	echo "Error: CMD 'perl'  is required but not found.  Aborting..." >&2 
#	echo "Error: CMD/script \'samtools\' in PROGRAM \'samtools\'  is required but not found.  Aborting..." >&2 
	exit 127
fi
#mysql
if [ $(CmdExists 'mysql') -eq 1 ]; then
	echo "Error: CMD 'mysql'  is required but not found.  Aborting..." >&2 
	exit 127
fi
#PASA
if [ -z $PASAHOME ] || [ ! -d $PASAHOME ]; then
	echo "Error: invalid PASAHOME" >&2
	exit 127
else
	echo "PASAHOME=$PASAHOME"
fi
if [ $runpasa -eq 1 ]; then
	if [ ! -s $PASAHOME/scripts/Launch_PASA_pipeline.pl ]; then
		echo "Error: Script Launch_PASA_pipeline.pl not found: $PASAHOME/scripts/Launch_PASA_pipeline.pl" >&2
		exit 127
	fi
fi



#################### Defaults #######################################
finalpasacmd=" $pasacmd --CPU $numcpus --ALIGNERS $aligners "



#################### Input and Output ###############################
if [ -z $transcript ] || [ ! -s $transcript ]; then
	echo "Error: invalid transcript -t $transcript" >&2
	exit 1
else
	finalpasacmd="$finalpasacmd -t $transcript "
fi
if [ -z $genome ] || [ ! -s $genome ]; then
	echo "Error: invalid genome -g $genome" >&2
	exit 1
else
	finalpasacmd="$finalpasacmd -g $genome "
fi
if [ ! -z $fulllength ]; then
	if [ ! -s $fulllength ]; then
		echo "Error: invalid full length accessions -f $fulllength" >&2
		exit 1
	else
		finalpasacmd="$finalpasacmd -f $fulllength "
	fi
fi
if [ -z $mysqldatabase ]; then
	echo "Error: invalid mysql database name" >&2
	exit 1
fi
if [ -z $readonly_user ]; then
	echo "Error: invalid mysql readonly username" >&2
	exit 1
fi
if [ -z $readonly_pawd ]; then
	echo "Error: invalid mysql readonly password" >&2
	exit 1
fi
if [ -z $readwrite_user ]; then
	echo "Error: invalid mysql read&write username" >&2
	exit 1
fi
if [ -z $readwrite_pawd ]; then
	echo "Error: invalid mysql read&write password" >&2
	exit 1
fi
echo "PASA parameters: $finalpasacmd"


#################### Main ###########################################
###1.Configure PASA_conf
echo -e "\n\n\n###Check PASA config"
if [ ! -s $PASAHOME/pasa_conf/conf.txt ]; then
	if [ -s $PASAHOME/pasa_conf/pasa.CONFIG.template ]; then
		echo "Preparing $PASAHOME/pasa_conf/conf.txt"
		cp $PASAHOME/pasa_conf/pasa.CONFIG.template $PASAHOME/pasa_conf/conf.txt
		perl -p -i -e "s/MYSQL_RW_USER=.*/MYSQL_RW_USER=$readwrite_user/; s/MYSQL_RW_PASSWORD=.*/MYSQL_RW_PASSWORD=$readwrite_pawd/; s/MYSQL_RO_USER=.*/MYSQL_RO_USER=$readonly_user/; s/MYSQL_RO_PASSWORD=.*/MYSQL_RO_PASSWORD=$readonly_pawd/; " $PASAHOME/pasa_conf/conf.txt
	else
		echo "Error: Can not find PASA $PASAHOME/pasa_conf/pasa.CONFIG.template file" >&2
		exit 1
	fi
else
	echo "PASA config detected: $PASAHOME/pasa_conf/conf.txt"
fi

###2.Configure PASA assembly
#Test if mysql database exists
echo -e "\n\n\n###Check MySQL databases"
mysql -u $readwrite_user -p$readwrite_pawd --skip-column-names -e "use $mysqldatabase" >/dev/null 2>&1
if [ $? -ne 0 ]; then
	echo -e "\n\n\n###Preparing ./alignAssembly.config"
	if [ ! -s ./alignAssembly.config ]; then
		cp $PASAHOME/pasa_conf/pasa.alignAssembly.Template.txt alignAssembly.config
	fi
	perl -p -i -e "s/MYSQLDB=.*/MYSQLDB=$mysqldatabase/" alignAssembly.config
else
	echo "Error: MySQL database: $mysqldatabase already exists, exit..." >&2
	exit 1
fi

###3. run PASA
if [ $runpasa -eq 1 ]; then
	if [ $(CmdExists 'data') -eq 0 ]; then
		mytime=$(data +%Y%m%d-%H%M%S)
	else
		mytime='temp'
	fi
	echo "$PASAHOME/scripts/Launch_PASA_pipeline.pl -c ./alignAssembly.config $finalpasacmd > $mytime.log 2> $mytime.err"
	$PASAHOME/scripts/Launch_PASA_pipeline.pl -c ./alignAssembly.config $finalpasacmd  > $mytime.log 2> $mytime.err
else
	echo "PASA not running as --runpasa not spcified"
fi
exit 0
