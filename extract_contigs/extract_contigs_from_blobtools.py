#!/usr/bin/env python

import datetime
import sys
from Bio import SeqIO
from collections import defaultdict

if len(sys.argv) != 3:
	print("extract_contigs_from_blobtools.py /absolute/path/to/assembly.fasta /absolute/path/to/blobtools_table.txt")
	sys.exit(1)

assembly_name = sys.argv[1].split("/")[-1]
blobfile = open(sys.argv[2], "r")

contig_taxon={}
for line in blobfile:
	if line.startswith("#"):
		continue
	elements = line.split("\t")
	contig_taxon[elements[0]] = elements[5]	

blobfile.close()


seqfile = open(sys.argv[1], "r")
out_dict = defaultdict(list)
i = 0
for seq_record in SeqIO.parse(seqfile, "fasta"):
	i += 1	
	print("Extracting sequence no.: %s" % str(i), end ='\r')
	if seq_record.id in contig_taxon.keys():
		out_dict[contig_taxon[seq_record.id]].append(seq_record.format("fasta"))
seqfile.close()
print()

for taxon, seqs in out_dict.items():
	if len(seqs) < 10: #skipping taxa with less than 10 sequences
		continue
	print("Writing:", assembly_name.replace(".fasta","")+"_blobtools_"+taxon.replace(" ", "_"), end='\r')
	outfile = open(assembly_name.replace(".fasta","")+"_blobtools_"+taxon.replace(" ", "_")+".fa", "w")
	for seq in seqs:
		outfile.write(seq)
	outfile.close()
print()


