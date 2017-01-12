#!/bin/sh
#Require python (2.6<=version<3), c compiler, python-dev, Cython
###APT
#apt-get install python-dev cython
###YUM
#yum install python-devel cython
#package_x86=ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.2.29/ncbi-blast-2.2.29+-ia32-linux.tar.gz
#package_x64=ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.2.29/ncbi-blast-2.2.29+-x64-linux.tar.gz
package=ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.2.29/ncbi-blast-2.2.29+-x64-linux.tar.gz
RunDir=$(cd `dirname $0`; pwd)
MachType=$(uname -m)

if which makeblastdb 2>/dev/null; then
  echo "Blast++ already installed"
  exit 1
fi

if [ ! -d $RunDir/blastpp ]; then
  mkdir -p $RunDir/blastpp
fi
cd $RunDir/blastpp
if [ ! -s ncbi-blast-*-linux.tar.gz ]; then
  echo "Did not find ncbi-blast-*-linux.tar.gz in $RunDir/blastpp."
  echo "Trying to download from NCBI ..."
  wget $package
fi
if [ -s ncbi-blast-*-linux.tar.gz ]; then
  blastpp_package=`ls ncbi-blast-*-linux.tar.gz`
  blastpp_basename=${blastpp_package%+-*-linux.tar.gz}
  blastpp_version=${blastpp_basename#ncbi-blast-}
  blastpp_basename2=${blastpp_package%-*}
  blastpp_machtype=${blastpp_basename2##*-}
  if [ "v$blastpp_version" = "v" ]; then
    echo "Warning: Can not get BLAST++ version, try to continue..."
  else
    echo "BLAST++ version $blastpp_version detected"
  fi
  echo "BLAST++ MachType: $blastpp_machtype"
  if [ "$blastpp_machtype" = "ia32" ]; then
    blastpp_machtype='x86'
  elif [ "$blastpp_machtype" = "x64" ]; then
    blastpp_machtype='x64'
  else
    echo "Warning: Can not detect BLAST++ MachType ... Use System MachType instead"
    blastpp_machtype=$MachType
  fi
else
  echo "Downloading BLAST++ from NCBI failed"
  exit 0
fi
if [ -d $RunDir/blastpp/v$blastpp_version ]; then
 rm -rf $RunDir/blastpp/v$blastpp_version
fi
mkdir -p $RunDir/blastpp/v$blastpp_version/src
mv $blastpp_package $RunDir/blastpp/v$blastpp_version/src/
cd $RunDir/blastpp/v$blastpp_version/
tar xvf $RunDir/blastpp/v$blastpp_version/src/$blastpp_package > /dev/null
mv ncbi-blast-* $blastpp_machtype
if [ ! -s $RunDir/blastpp/v$blastpp_version/$blastpp_machtype/bin/makeblastdb ]; then
  echo "Installing BLAST++ failed"
  cd $RunDir/blastpp/
  mv $RunDir/blastpp/v$blastpp_version/src/$blastpp_package $RunDir/blastpp/
  exit 0
fi
echo "*******************************************"
echo "************ Important ********************"
echo "*******************************************"
echo "Add to /etc/profile or ~/.bashrc\n\n"
echo "###BLAST++"
echo "export PATH=$RunDir/blastpp/v$blastpp_version/$blastpp_machtype/bin"':$PATH'
echo "\n\n\n"
