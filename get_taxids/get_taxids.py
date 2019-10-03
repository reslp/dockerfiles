#!/usr/bin/env python

import datetime
import sys

if len(sys.argv) != 3:
	print("get_taxids.py /absolute/path/to/mapping/prot.accessions2taxid /absolute/path/to/tabular_blast_results_file")
	sys.exit(1)


file1 = open(sys.argv[1], "r")
all_ids_dict={}
i = 0

sys.stderr.write(str(datetime.datetime.now())+" Reading all taxids\n")

for line in file1:
	elements = line.split("\t")
	all_ids_dict[elements[1]]=elements[2]
file1.close()

#print(all_ids_dict)
sys.stderr.write(str(datetime.datetime.now())+" Mapping taxids\n")
file2 = open(sys.argv[2], "r")

for line in file2:
	elements = line.strip("\n").split("\t")
	try:	
		taxid = all_ids_dict[elements[1]]
	except:
		continue
	print(elements[0],"\t",taxid,"\t", elements[11])

file2.close()
sys.stderr.write(str(datetime.datetime.now())+" done\n")
