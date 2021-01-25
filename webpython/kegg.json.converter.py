#!/usr/bin/env python3
### https://zhuhaizhen.github.io/2020/06/30/KEGG-all-kos/
### download: https://www.kegg.jp/kegg-bin/get_htext?ko00001
### convert ko00001.json to ko.txt


import os
import json
import re

class KOParser(object):
	def __init__(self, map):
		self.map = map

	def json_parser(self):
		# check the existance of json file
		if os.path.exists(self.map):
			with open(self.map, 'r') as f:
				ko_list = []
				# convert json to dict
				map_dict = json.load(f)
				# print(type(map_dict)) dict
				maps = map_dict['children']
				# print(type(maps)) list
				for map in maps:
					map_name = map['name'][0:5] + '\t' + map['name'][6:]
					# print(map_name)
					# print(type(map['children']))
					for map_l_1 in map['children']:
						map_l_1_name = map_l_1['name'][0:5] + '\t' + map_l_1['name'][6:]
						# print(map_l_1_name)
						for pathway in map_l_1['children']:
							try:
								for genes in pathway['children']:
									pathway_name = pathway['name'][0:5] + '\t' + pathway['name'][6:]
									# print(genes['name'])
									k_num = genes['name'].split(sep='  ')[0]
									gene_name = genes['name'].split(sep='  ')[1].split(sep=';')[0]
									anno = genes['name'].split(sep='  ')[1].split(sep=';')[-1]
									try:
										pattern = re.compile('(.*)(\[EC:.*\])')
										product = re.search(pattern, anno).group(1)
										ec = re.search(pattern, anno).group(2)
										ko = k_num + '\t' + gene_name + '\t' + product + '\t' + ec
									except:
										ko = k_num + '\t' + gene_name + '\t' + anno
									# print(ko)
									info = map_name + '\t' + map_l_1_name + '\t' + pathway_name + '\t' + ko
									ko_list.append(info)
							except:
								continue
				return ko_list
		else:
			print('Error: Json file does not exist.')

	def save_data(self, result, file_name):
		# convert list to string for writing as file
		for info in result:
			kos = ''.join(info) + '\n'
			with open(file_name, 'a', encoding='utf-8') as f:
				f.write(kos)


if __name__ == '__main__':
	map = 'ko00001.json'
	file_name = 'ko.txt'
	ko_json = KOParser(map)
	ko = ko_json.json_parser()
	ko_json.save_data(ko, file_name)
