#!/bin/bash

declare -a QuotaDirArr=("/home/hpcusers/asistant" "/home/hpcusers/doctor" "/home/hpcusers/master" "/home/hpcusers/postdoc/")

for inddir in ${QuotaDirArr[@]}; do
        for inddir2 in `find $inddir/ -maxdepth 1 -type d`; do
                size=$(du -sh $inddir2 2> /dev/null | perl -lne 's/\s+.*$//;if (/^(\d+\.*\d*)T$/){print $1 * 1000000000000;}elsif (/^(\d+\.*\d*)G$/){print $1 * 1000000000;}elsif (/^(\d+\.*\d*)M$/){print $1 * 1000000;}elsif (/^(\d+\.*\d*)K$/){print $1 * 1000;}elsif (/^\d+$/){print $1 * 1000;}else {print 0;}')
                if [ $size -ge 2000000000000 ]; then
                        du -sh $inddir2 2> /dev/null | mail -s " /home heavy users " -r 381717988@qq.com -aFrom:noreply@LUFHPC.com lufuhao@henu.edu.cn
                	du -sh $inddir2 2> /dev/null > /home/hpcsoft/ProdSoft/hpcutils/quotacheck.log
		fi
        done    
done
