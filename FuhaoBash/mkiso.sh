#!/bin/sh

 

#****************************************************

# Author: Muddyboot - toobyddum@gmail.com

# Last modified: 2007-09-11 21:21

# Filename: mkiso.sh

# Description: script for easy creating ISO image

#****************************************************

 

if [ $# -lt 3 ]; then

echo -e "\nUsage: `basename $0` source_dir output_iso cd_label \n"

exit 1

fi

 

source=$1

output=$2

label=$3

 

### extra mkiso argument

shift 3

for i in $@; do

extra_args="$extra_args $1 "

shift

done

 

if [ ! -e "$source" ]; then

echo -e "\nERR: Source file or directory does not exist ! \n"

exit 1

fi

 

## remove exists TRANS.TBL files

if [ -d "$source" ]; then

find $source -name TRANS.TBL | xargs rm -f

fi

 

### 制作ISO

mkisofs -J -T -R $extra_args \

-V $label -o $output $source

 

### 加入 MD5 校验信息

MD5_CHECKSUM=`whereis implantisomd5|awk -F': ' '{print $2}'`

 

if [ -z "$MD5_CHECKSUM" ]; then

echo -e "\n** WARNING: implantisomd5 not found, no md5sum added.\n"

else

echo -e "\n** Good, implantisomd5 program found."

echo "Adding md5sum information for ISO image ..."

implantisomd5 --force $output

fi

 

echo

echo "** ISO image $output created successfully ! "

echo
