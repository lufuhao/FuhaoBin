#!/usr/bin/python
#coding:utf-8
'''使用方法
python get_species_placement_in_NCBI.py Organism_name.txt 2 placement.txt
输出结果
root,cellular organisms,Eukaryota,Viridiplantae,Streptophyta,Streptophytina,Embryophyta,Tracheophyta,Euphyllophyta,Spermatophyta,Magnoliophyta,Mesangiospermae,eudicotyledons,Gunneridae,Pentapetalae,rosids,malvids,Myrtales,Combretaceae,Lumnitzera,Lumnitzera littorea
root,cellular organisms,Eukaryota,Viridiplantae,Streptophyta,Streptophytina,Embryophyta,Tracheophyta,Euphyllophyta,Spermatophyta,Magnoliophyta,Mesangiospermae,eudicotyledons,Gunneridae,Pentapetalae,rosids,malvids,Myrtales,Lythraceae,Punica,Punica granatum
root,cellular organisms,Eukaryota,Viridiplantae,Streptophyta,Streptophytina,Embryophyta,Tracheophyta,Euphyllophyta,Spermatophyta,Magnoliophyta,Mesangiospermae,eudicotyledons,Gunneridae,Pentapetalae,rosids,malvids,Myrtales,Lythraceae,Heimia,Heimia myrtifolia
root,cellular organisms,Eukaryota,Viridiplantae,Streptophyta,Streptophytina,Embryophyta,Tracheophyta,Euphyllophyta,Spermatophyta,Magnoliophyta,Mesangiospermae,eudicotyledons,Gunneridae,Pentapetalae,rosids,malvids,Myrtales,Lythraceae,Sonneratia,Sonneratia alba
root,cellular organisms,Eukaryota,Viridiplantae,Streptophyta,Streptophytina,Embryophyta,Tracheophyta,Euphyllophyta,Spermatophyta,Magnoliophyta,Mesangiospermae,eudicotyledons,Gunneridae,Pentapetalae,rosids,malvids,Myrtales,Onagraceae,Onagroideae,Epilobieae,Epilobium,Epilobium ulleungensis'''
import sys
from ete3 import NCBITaxa
input_file = sys.argv[1]
colNum = sys.argv[2]
output_file = sys.argv[3]
ncbi = NCBITaxa()
fw = open(output_file,"w")
with open(input_file,"r") as fr:
    for line in fr:
        species_line = line.strip()
        species_arr=species_line.split("\t")
        species_name = species_arr[ int(colNum)-1 ]
        fw.write(species_line+"\t")
        name2taxid = ncbi.get_name_translator([species_name])
        for a,b in name2taxid.items():
            lineage = ncbi.get_lineage(b[0])
            names = ncbi.get_taxid_translator(lineage)
            i = 1
            for taxid in lineage:
                if i < len(lineage):
                    fw.write(names[taxid]+"\t")
                    i = i + 1
                else:
                    fw.write(names[taxid]+"\n")
        print(species_name + ":","OK")
fw.close()
