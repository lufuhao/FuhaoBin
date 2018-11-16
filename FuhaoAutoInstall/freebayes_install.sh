#!/bin/sh
RunDir=$(cd `dirname $(readlink -f $0)`; pwd)
MachType=$(uname -m)

################# help message ######################################
help()
{
cat<<HELP

install_freebayes.sh

Version: 20150309

Descriptions:
Install freebayes

Options:
  -h    Print this help message
  -d    [Msg] Installation ROOT
  		Should install ROOT/freebayes/version/$MachType/
  -p	[Opt] Program name (Default: freebayes)
  -m	[Opt] MachType (Default: uname -m)
  -v	[Opt] Version (Default: auto-detect)

Example:
  install_freebayes.sh -d install_path

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
echo -e "\n######################\nFreebayes installing ...\n######################\n"
#echo "Adding $RunDir/bin into PATH"
#export PATH=$RunDir/bin:$RunDir/utils/bin:$PATH

###Defaults
rootpath="$HOME/local"
program='freebayes'
if [ -z $MachType ]; then
  MachType='x86_64'
fi
version='ver'
#################### parameters #####################################
while [ -n "$1" ]; do
  case "$1" in
    -h) help;shift 1;;
    -d) rootpath=$2;shift 2;;
    -p) program=$2;shift 2;;
    -m) MachType=$2;shift 2;;
    -v) version=$2;shift2;;
    --) shift;break;;
    -*) echo "error: no such option $1. -h for help"  > /dev/stderr;exit 1;;
    *) break;;
  esac
done
if [ -z "$rootpath" ]; then
  echo "Error: missing -d option" > /dev/stderr
  exit 1
fi
if [ -z "$program" ]; then
  echo "Error: missing -p option" > /dev/stderr
  exit 1
fi
if [ -z "$MachType" ]; then
  echo "Error: missing -m option" > /dev/stderr
  exit 1
fi
if [ -z "$version" ]; then
  echo "Error: missing -v version" > /dev/stderr
  exit 1
fi



if [ ! -d $rootpath/$program ]; then
  mkdir -p $rootpath/$program
fi
cd $rootpath/$program
if [ -d $rootpath/$program/version ]; then
  echo "Path: $rootpath/$program/version exists" > /dev/stderr
  exit 1
else
  mkdir -p $rootpath/$program/version
fi
cd $rootpath/$program/version
echo -e "\n\n\n##### Downloading #####"
git clone --recursive git://github.com/ekg/freebayes.git
if [ ! -d $rootpath/$program/version/freebayes ]; then
  echo "git clone freebayes fails" > /dev/stderr
  exit 1
fi
mv freebayes $MachType
cd $rootpath/$program/version/$MachType
echo -e "\n\n\n##### Compiling #####"
make

if [ -s $rootpath/$program/version/$MachType/bin/freebayes ]; then
  version=$($rootpath/$program/version/$MachType/bin/freebayes --version | cut -f2 -d':' |sed 's/ //g')
else
  echo "Compilation fails" > /dev/stderr
  exit 1
fi
if [ -z "$version" ]; then
  echo "Version detection failed" > /dev/stderr
  exit 1
fi
cd $rootpath/$program
mv version $version
echo -e "\n\n\n"
echo "### Freebayes wrapper"
echo "PATH=$rootpath/$program/$version/$MachType/bin"':$PATH'
exit 0
