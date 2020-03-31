#!/bin/bash

#GO annotation using eggnog-mapper

### Requirements
###		python2.7
###		40Gb disk
###

git clone https://github.com/eggnogdb/eggnog-mapper.git

#eggNOG数据库
python ./download_eggnog_data.py
#OR
wget http://eggnogdb.embl.de/download/emapperdb-5.0.0/eggnog.db.gz
wget http://eggnogdb.embl.de/download/emapperdb-5.0.0/eggnog_proteins.dmnd.gz
#下载好后移至eggnog-mapper安装目录的data文件夹下并解压。

python emapper.py -i pep.fa --output out -m diamond --cpu 12

#-i：输入蛋白序列。
#--output：输出文件前缀。
#-m diamond：使用DIAMOND进行序列比对。
#--cpu：使用的线程数。
#使用DIAMOND进行比对的速度非常快。30万条序列用12个线程注释花了5个多小时。
#注释完成后会输出两个文件，emapper.annotations为后缀的文件记录了注释结果。
#文件一共有22列：
#1. query_name 输入的ID
#2. seed eggNOG ortholog 在eggNOG中比对到的最佳结果
#3. seed ortholog evalue
#4. seed ortholog score
#5. Predicted taxonomic group
#6. Predicted protein name 预测得到的蛋白名
#7. Gene Ontology terms 注释到的GO terms
#8. EC number
#9. KEGG_ko 注释到的ko
#10. KEGG_Pathway 注释到的通路
#11. KEGG_Module
#12. KEGG_Reaction
#13. KEGG_rclass
#14. BRITE
#15. KEGG_TC
#16. CAZy
#17. BiGG Reaction
#18. tax_scope: eggNOG taxonomic level used for annotation
#19. eggNOG OGs
#20. bestOG (deprecated, use smallest from eggnog OGs)
#21. COG Functional Category
#22. eggNOG free text description

###大家可以根据自己的需求提取对应的信息。

#参考资料：
#https://github.com/eggnogdb/eggnog-mapper/wiki/eggNOG-mapper-v2
#https://www.jianshu.com/p/e646c0fa6443



