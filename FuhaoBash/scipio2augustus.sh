#!/bin/sh
RunDir=$(cd `dirname $(readlink -f $0)`; pwd)
MachType=$(uname -m)

################# help message ######################################
help()
{
cat<<HELP

$0 --- NGS sIMulation PipeLinE

Version: 20150320

Requirements:
	PATH: new_species.pl, etraining, optimize_augustus.pl, filterGenes.pl (AUGUSTUS)
	PATH: perl

Descriptions:
From scipio gff to augustus training

Options:
  -h    Print this help message
  -i    CONFIG file
  -n    Species name ([a-zA-Z0-9_] and no space)
  -d    Keep temporary files

Example:
  $0 -i raw.gb -n triticum_aestivum

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
#################### Defaults #######################################
echo -e "\n######################\nNGSimple initializing ...\n######################\n"
echo "Adding $RunDir/bin into PATH"
export PATH=$RunDir/bin:$RunDir/utils/bin:$PATH

###Defaults
debug=0
#################### parameters #####################################
while [ -n "$1" ]; do
  case "$1" in
    -h) help;shift 1;;
    -i) scipiogb=$2;shift 2;;
    -n) speciesname=$2;shift 2;;
    -d) debug=1;shift 1;;
    --) shift;break;;
    -*) echo "error: no such option $1. -h for help" > /dev/stderr;exit 1;;
    *) break;;
  esac
done

###Default###
if [ -z "$AUGUSTUS_CONFIG_PATH" ]; then
  echo "Error: please setup AUGUSTUS_CONFIG_PATH"
  exit 1
fi
if [ ! -d "$AUGUSTUS_CONFIG_PATH/species" ]; then
  echo "Info: AUGUSTUS_CONFIG_PATH/species not exists, creating..."
  mkdir -p "$AUGUSTUS_CONFIG_PATH/species"
fi
if [ -z "$speciesname" ]; then
  echo "Error: invalid speciesname"
  exit 1
fi
if [ -d "$AUGUSTUS_CONFIG_PATH/species/$speciesname" ]; then
  echo "Warning: $AUGUSTUS_CONFIG_PATH/species/$speciesname exists. Exit..."
  exit 1
fi



###Input the output###
if [ ! -s $scipiogb ]; then
  echo "Error: invald scipio gb file"
  exit 1
fi
scipiogb_name=${scipiogb##*/}
scipiogb_base=${scipiogb_name%.*}



###Main###
###1. create parameters for new training set
echo -e "\n\n\nStep1: create training set"
new_species.pl --species=$speciesname
if [ $? -eq 0 ]; then
  echo "Step1: Training set created: $AUGUSTUS_CONFIG_PATH/species/$speciesname"
else
  echo "Step1: Error to create training set"
  exit 1
fi
###2. e-training to get problematic seq IDs
echo -e "\n\n\nStep2: remove problematic seqIDs"
etraining --species=$speciesname --stopCodonExcludedFromCDS=true $scipiogb 2> $scipiogb_base.train_err
cat $scipiogb_base.train_err | perl -pe 's/.*n\s+sequence (\S+):.*/$1/' | sort | uniq > $scipiogb_base.badgenes_lst
filterGenes.pl $scipiogb_base.badgenes_lst $scipiogb > $scipiogb_base.filter.gb
if [ $debug -eq 0 ]; then
  echo "Step2: delete temporary files"
  rm $scipiogb_base.train_err $scipiogb_base.badgenes_lst
fi
if [ ! -s $scipiogb_base.filter.gb ]; then
  echo "Step2: Error to filter raw genbank files"
  exit 1
fi
###3. e-training
echo -e "\n\n\nStep3: e-training"
etraining --species=$speciesname --stopCodonExcludedFromCDS=true $scipiogb_base.filter.gb 2> $scipiogb_base.filter.train_error2
if [ $? -eq 0 ]; then
  echo "Step3: successful etraining"
else
  echo "Step3: failed etraining"
  exit 1
fi
echo "IMPORTANT: please check $scipiogb_base.filter.train_error2 to make sure all the problematic sequences are excluded"
###4. optimize augustus
echo -e "\n\n\nStep4: optimize augustus"
optimize_augustus.pl --species=$speciesname --metapars=$AUGUSTUS_CONFIG_PATH/species/$speciesname/${speciesname}_metapars.cfg $scipiogb_base.filter.gb 1> $scipiogb_base.filter.optimize.log 2> $scipiogb_base.filter.optimize.err
if [ $? -eq 0 ]; then
  echo "Step4: optimize_augustus.pl"
else
  echo "Step4: optimize_augustus.pl"
  exit 1
fi
exit 0
