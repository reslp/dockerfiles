#!/usr/bin/env python3
# written by Philipp Resl

import argparse
import sys
import random
import os


if sys.version_info[0] < 3:
    raise Exception("Must be using Python 3")

pars = argparse.ArgumentParser(prog="parse_dbcan.py", description = """This script parses the overview.txt outpfile from dbcan. It will only keep CAZyme family names and remove EC codes and other annotations.""", epilog = """written by Philipp Resl""")
pars.add_argument('-i', dest="inputfile", required=True, help="Path to overview.txt files, comma seperated. The name of the file is used as species name. e.g.: -f file1.txt,file2.txt")
pars.add_argument('--mergeby=', dest="mergeby", default="none", help="How to merge lower-level CAZyme annotations. Possible values: families, subfamilies, enzyme_class, none. Default: none")
pars.add_argument('--predictor=', dest="predictor", required=False, default="diamond", help="The results of which predictor should be used (Currently only diamond is implemented)? possible values: hmmer, dbcansub, diamond. default: diamond")
pars.add_argument('--for-itol', dest="itol", action="store_true", default=False, help="Create itol Datafile. default: false")
pars.add_argument('--seed', dest="seed", action="store_true", default="42", help="Random seed. Change to get different colors in itol output. default: 42")
pars.add_argument('--transpose', dest="transpose", action="store_true", default=False, help="Transpose output table.")
pars.add_argument('--include=', dest="include", default="all", help="List of CAZyme groups to be included in the output. Comma seperated. e.g. --include=GH5,GH43,GH16. Default: all")
args=pars.parse_args()
geneids = []
cazymes = []

piechart_template = """
DATASET_PIECHART
SEPARATOR COMMA

DATASET_LABEL,example multi bar chart

COLOR,#ff0000

FIELD_COLORS


FIELD_LABELS

DATA
# sp,position,radius,values,seperated,by,comma
"""
def process_cazy(cazyme_string, predictor):
	if predictor == "diamond":
		cazyme_string = cazyme_string.split("\t")[4]
	if cazyme_string.startswith("-"):
		return []
	if args.mergeby == "none":
		return [cazyme_string]
	if args.mergeby == "families":
		if "+" in cazyme_string:
			cazy = []
			for i in cazyme_string.split("+"):
				if "_" in i:
					cazy.append(i.split("_")[0])
				elif "." in i:
					#print("CAZyme annotation contains non-standard term:",i, ". Will skip this", file=sys.stderr)
					continue
				else:
					cazy.append(i)
			return cazy
		else:
			if "_" in cazyme_string:
				return [cazyme_string.split("_")[0]]
			else:
				return [cazyme_string]
	if args.mergeby == "enzyme_class":
		cazy = []
		for el in cazyme_string.split("+"):
			if el.startswith("GH"):
				cazy.append("GH")
			if el.startswith("CE"):
				cazy.append("CE")
			if el.startswith("CBM"):
				cazy.append("CBM")
			if el.startswith("GT"):
				cazy.append("GT")
			if el.startswith("PL"):
				cazy.append("PL")
			if el.startswith("AA"):
				cazy.append("AA")
		return cazy
	if args.mergeby == "subfamilies":
		if "+" in cazyme_string:
			cazy = []
			for i in cazyme_string.split("+"):
				if "." in i:
					#print("CAZyme annotation contains non-standard term:",i, ". Will skip this", file=sys.stderr)
					continue
				else:
					cazy.append(i)
			return cazy
		else:
			return [cazyme_string]

# main loop
outdict = {}
all_terms = []
for filename in args.inputfile.split(","):
	with open(filename, "r") as file:
		sp = os.path.basename(filename).split(".")[0]
		outdict[sp] = []
		for line in file:
			if line.startswith("Gene"): # skip first line
				continue
			line = line.strip()
			outdict[sp] += process_cazy(line, args.predictor)
		all_terms += list(set(outdict[sp])) 

# check if only specific items should be displayed
if args.include != "all":
	all_terms = args.include.split(",")
else:
	all_terms = sorted(list(set(all_terms)))

# output
if args.itol:
	random.seed(int(args.seed))
	colors = ["#"+''.join([random.choice('0123456789ABCDEF') for j in range(6)]) for i in range(len(all_terms))]

	field_labels = "FIELD_LABELS," + ",".join(all_terms)
	piechart_template = piechart_template.replace("FIELD_LABELS", field_labels)
	field_colors = "FIELD_COLORS," +  ",".join(colors)
	piechart_template = piechart_template.replace("FIELD_COLORS", field_colors)
	for sp in outdict.keys():
		sp_str = sp + ",-1,50"
		for term in all_terms:
			if term not in outdict[sp]:
				sp_str +=  ",0"
			else:
				sp_str += "," + str(outdict[sp].count(term))
		piechart_template = piechart_template + "\n" + sp_str
	print(piechart_template)
else:
	if args.transpose:
		header="cazyme"
		for key in outdict.keys():
			header += "\t" + key
		print(header)
		for term in all_terms:
			sp_str = term
			for key in outdict.keys():
				if term not in outdict[key]:
					sp_str += "\t0"
				else:
					sp_str += "\t" + str(outdict[key].count(term))
			print(sp_str)
	else:
		header="species"
		for term in all_terms:
			header += "\t" + term
		print(header)
		for key in outdict.keys():
			sp_str = key
			for term in all_terms:
				if term not in outdict[key]:
					sp_str +=  "\t0"
				else:
					sp_str += "\t" + str(outdict[key].count(term))
			print(sp_str)
