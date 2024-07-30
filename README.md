# dockerfiles

This repository contains Dockerfiles I accumulated for different purposes.

Some of the images are bare-bone and should maybe be changed when incorporated into analysis pipelines.

The **working directory** in these containers is generally `/data`. This can be used as a bindpoint. Downloaded software is generally kept in `/software`. Some older containers here may however not follow these rules.

The **tags** given for each container correspond the the available versions of different programs. The Dockerfiles refer to the latest tags.

**Currently available containers:**

[abyss](#abyss)
[alitv](#alitv)
[angsd](#angsd)
[astral](#astral)
[augustus](#augustus)
[bamutil](#bamutil)
[bbmap](#bbmap)
[beast](#beast)
[bedops](#bedops)
[bedtools](#bedtools)
[bibcure](#bibcure)
[biopython_plus](#biopython_plus)
[blobtools](#blobtools)
[bmge](#bmge)
[bonito](#bonito)
[bowtie2](#bowtie2)
[busco](#busco)
[cafe](#cafe)
[cafe-viz](#cafe-viz)
[clustalo](#clustalo)
[colord](#colord)
[concoct](#concoct)
[consent](#consent)
[cupp](#cupp)
[deeploc](#deeploc)
[dentist](#dentist)
[diamond](#diamond)
[edta](#edta)
[eggnog-mapper](#eggnog-mapper)
[emboss](#emboss)
[ete3](#ete3)
[fasta3](#fasta3)
[fastcodeml](#fastcodeml)
[fastqc](#fastqc)
[fastuniq](#fastuniq)
[ffsend](#ffsend)
[flappie](#flappie)
[flash](#flash)
[fmlrc](#fmlrc)
[gatk](#gatk)
[genome-viz](#genome-viz)
[get_organelle](#get_organelle)
[gnuplot](#gnuplot)
[guidance2](#guidance2)
[hmmer](#hmmer)
[ima2](#ima2)
[imagemagick](#imagemagick)
[interproscan-wrapper](#interproscan-wrapper)
[ipyrad](#ipyrad)
[iqtree](#iqtree)
[julia](#julia)
[jupyter-lab](#jupyter-lab)
[keras](#keras)
[kmergenie](#kmergenie)
[macse](#macse)
[mafft](#mafft)
[magic-wormhole](#magic-wormhole)
[mamba](#mamba)
[maxbin](#maxbin)
[mcl](#mcl)
[metabat](#metabat)
[metaxa](#metaxa)
[miniBarcoder](#miniBarcoder)
[mitos](#mitos)
[mosdepth](#mosdepth)
[mrbayes](#mrbayes)
[multiqc](#multiqc)
[mummer](#mummer)
[muscle](#muscle)
[mutt](#mutt)
[ncbi-blast](#ncbi-blast)
[ncbi-datasets](#ncbi-datasets)
[necat](#necat)
[nim](#nim)
[nimbio_plus](#nimbio-plus)
[norgal](#norgal)
[novoplasty](#novoplasty)
[oma-standalone](#oma-standalone)
[ont-fast5-api](#ont-fast5-api)
[orthofinder](#orthofinder)
[paml](#paml)
[pal2nal](#pal2nal)
[pandoc](#pandoc)
[pass](#pass)
[pear](#pear)
[phykit](#phykit)
[phylobayes](#phylobayes)
[phylobayes-mpi](#phylobayes-mpi)
[phylopy](#phylopy)
[phylo-scripts](#phylo-scripts)
[picard](#picard)
[pigz](#pigz)
[purge_haplotigs](#purge_haplotigs)
[pyinstaller](#pyinstaller)
[pyscaf](#pyscaf)
[raxml-ng](#raxml-ng)
[rename](#rename)
[rmarkdown](#rmarkdown-based-on-r)
[repeatmodeler-repeatmasker](#repeatmodeler-repeatmasker)
[rphylogenetics](#rphylogenetics)
[samtools](#samtools)
[scrape_cazy](#scrape_cazy)
[seqtk](#seqtk)
[singularity](#singularity)
[snakemake](#snakemake)
[spades](#spades)
[sra-tools](#sra-tools)
[stringtie](#stringtie)
[table2asn](#table2asn)
[tbl2asn](#tbl2asn)
[treeshrink](treeshrink)
[trimal](#trimal)
[trimmomatic](#trimmomatic)
[typst](#typst)
[ubuntu](#ubuntu)
[viralmsa](#viralmsa)
[quast](#quast)
[quicktree](#quicktree)






## iqtree
tags: 2.0rc2, 2.0.7

```
docker run -v $(pwd):/data reslp/iqtree:2.0.7
```
 

## mrbayes-MPI

tags: 3.2.6

```
docker run -v $(pwd):/data -it reslp/mrbayes-mpi:3.2.6 /sbin/my_init -- mpirun --allow-run-as-root -np 8 mb /data/alignment.nex
```

## phylo-scripts
tags: latest

```
docker run reslp/phylo-scripts:latest
```


## Spades
tags: 3.12.0, 3.13.0, 3.15.3

WORKDIR /data

```
docker pull reslp/spades:3.13.0
```


## trimmomatic
tags: 0.38

```
docker run -t -v $(pwd):/data reslp/trimmomatic:0.38 trimmomatic PE -phred33 -threads 8 /data/pair1.fq.gz /data/pair2.fq.gz /data/trimmed_pair1.fq.gz /data/unpaired_pair1.fq.gz /data/trimmed_pair2.fq.gz /data/unpaired_pair2.fq.gz HEADCROP:20 ILLUMINACLIP:/data/adapters.fa:2:30:10 LEADING:30 TRAILING:30 SLIDINGWINDOW:4:15 MINLEN:80
```
or

```
docker run --user $(id -u):$(id -g) --rm -v $(pwd):/data reslp/trimmomatic trimmomatic
```

Important: The adapter file is not part of the docker image. It needs to passed together with the mounted directory.


## fastqc
tags: 0.11.7

```
docker run -v $(pwd):/data --rm reslp/fastqc:0.11.7 fastqc /data/illumina_reads_5_R1_trimmed.fq -o /data

docker run --user $(id -u):$(id -g) --rm -v $(pwd):/data reslp/fastqc:0.11.7 fastqc -h
```

## maxbin
tags: 2.2.6

`docker pull reslp/maxbin:2.2.6`

`docker run --rm reslp/maxbin:2.2.6`

## concoct
tags: 1.1

`docker pull reslp/concoct:1.1`

`docker run --rm reslp/concoct:1.1`

A complete example filtering a metagenome. The `metagenome.fasta` file should be in the current working directory.

```
docker run -t --user $(id -u):$(id -g) -v $(pwd):/data/ --rm reslp/concoct:1.1 cut_up_fasta.py /data/metagenome.fasta -c 10000 -o 0 --merge_last -b /data/metagenome.fasta_contigs_10K.bed > metagenome.fasta_contigs_10K.fa
docker run -t --user $(id -u):$(id -g) -v $(pwd):/data/ --rm reslp/concoct:1.1 concoct_coverage_table.py /data/metagenome.fasta_contigs_10K.bed /data/metagenome.fasta.bam > concoct_coverage_table.tsv
docker run -t --user $(id -u):$(id -g) -v $(pwd):/data/ --rm reslp/concoct:1.1 concoct --composition_file /data/metagenome.fasta_contigs_10K.fa --coverage_file /data/concoct_coverage_table.tsv -b /data/concoct/metagenome.fasta_concoct --threads 24
docker run --user $(id -u):$(id -g) -v $(pwd):/data/ --rm reslp/concoct:1.1 merge_cutup_clustering.py /data/concoct/metagenome.fasta_concoct_clustering_gt1000.csv > concoct/metagenome.fasta_concoct_clustering_merged.csv
docker run -t --user $(id -u):$(id -g) -v $(pwd):/data/ --rm reslp/concoct:1.1 extract_fasta_bins.py /data/metagenome.fasta /data/concoct/metagenome.fasta_concoct_clustering_merged.csv --output_path /data/concoct/bins
```

## metabat
tags: 2.13

```
docker run --rm reslp/metabat:2.13
```

## busco
tags: 3.0.2, 4.1.4, 5.2.1

The command for busco >4 is busco instead of run_busco on busco <4.

```
docker run --rm reslp/busco:3.0.2
docker run --rm reslp/busco:4.1.4
```

## blobtools
tags: 1.1.1

These commands expect that blast results are alread present:

```
docker run -t --user $(id -u):$(id -g) -v $(pwd):/data/ --rm reslp/blobtools:1.1.1 create -i /data/metagenome.fasta -b /data/metagenome.fasta.bam -t /data/blobtools/metagenome.fasta_diamond_matches_formatted -o /data/blobtools/metagenome.fasta_blobtools
docker run -t --user $(id -u):$(id -g) -v $(pwd):/data/ --rm reslp/blobtools:1.1.1 view -i /data/blobtools/metagenome.fasta_blobtools.blobDB.json -o /data/blobtools/
docker run -t --user $(id -u):$(id -g) -v $(pwd):/data/ --rm reslp/blobtools:1.1.1 plot -i /data/blobtools/metagenome.fasta_blobtools.blobDB.json -o /data/blobtools/
```


## QUAST
tags: 5.0.2

Run QUAST on all fasta files in current directory:

`docker run --rm -t -v $(pwd):/data/ reslp/quast:5.0.2 quast.py $(ls -d -1 "*.fasta" | sed 's/^/\/data\//') -o /data/quast --silent
`

## ncbi-blast
tags: 2.9.0, 2.15.0

To create a blast database and run a blastp search:
`my_sequences.fasta`and `sequences_to_search.fa`need to be in the current working directory. The blastdb will also be created there.

```
docker run --rm -v $(pwd):/files/ reslp/ncbi-blast:2.9.0 makeblastdb -in /files/my_sequences.fasta -dbtype 'prot' -hash_index -out /files/my_blastdb

docker run --rm -v $(pwd):/files/ reslp/ncbi-blast:2.9.0 blastp -db /files/my_blastdb -query /files/sequences_to_search.fa -outfmt 6 > blast_results.txt

```

## singularity
tags: 3.4.1, 3.6.3

To run the conatiner and to be able to run singularity containers inside it docker needs to be run with privileged mode:

```
docker run --privileged -it --rm -v $(pwd):/data --entrypoint /bin/bash reslp/singularity:3.4.1
```

Inside the container it is now possible to run singularity images. For example the singularity image for funannotate (which was downloaded and converted to .sif from Docker Hub).

```
singularity run funannotate.sif mask -i Trapelia_coarctata_sorted.fas -o Trapelia_coarctata_mask.fas -m repeatmasker
```


## eggnog-mapper

tags: 1.0.3, 2.0.1

```
docker pull reslp/eggnog-mapper:2.0.1
docker run reslp/eggnog-mapper:2.0.1
```

## ASTRAL
tags: 5.7.1

```
docker pull reslp/astral:5.7.1
docker run -it --rm reslp/astral:5.7.1 java -jar /ASTRAL-5.7.1/Astral/astral.5.7.1.jar
```

## pal2nal
tags: 14

```
docker pull reslp/pal2nal:14
docker run -it reslp/pal2nal:14
```

## mafft
tags: 7.464

```
docker run -it --rm -v $(pwd):/data reslp/mafft:7.464 mafft --auto --quiet /data/OG0000331.fa > OG0000331_aligned.fa
```

## scrape_cazy
tags: 1


Containerized version of a python script to download information from CAZY.org:

```
docker run -it --rm reslp/scrape_cazy:1
docker run --rm -it -v $(pwd):/data reslp/scrape_cazy:1 scrape_cazy.py -f GH3
```

## deeploc
tags: 1.0

This container needs some things set up properly to work:

1. Download DeepLoc: https://services.healthtech.dtu.dk/software.php
2. Change the first line in `bin/deeploc`to `#!/usr/bin/env python`
3. Mind that you correctly use the bindpoints of the container on execution:

```
docker run --rm -it -v $(pwd)/deeploc-1.0/bin:/external -v $(pwd)/deeploc-1.0/DeepLoc:/usr/lib/python3/dist-packages/DeepLoc -v $(pwd):/data reslp/deeploc:1.0 deeploc -f /data/deeploc-1.0/test.fasta
```

## CONSENT
Correct long reads: https://github.com/morispi/CONSENT

tags: v2.1

```
docker run --rm -v $(pwd):/corrected reslp/consent:v2.1 CONSENT-correct --in /corrected/Lmin_combined.fasta --out /corrected/Lmin_combined_corrected.fasta --type ONT
```

## fmlrc
Correct long reads: https://github.com/holtjma/fmlrc

tags: v1.0.0

```
docker run reslp/fmlrc:v1.0.0
```

## biopython_plus
This image should serve as a base environment to execute python scripts which use biopython and scipy. Tag version corresponds to biopython.

tags: 1.77

```
docker run reslp:biopython_plus:1.77
```


## miniBarcoder
This is a dockerfile for: https://github.com/asrivathsan/miniBarcoder

tags: 5e1dc3b

```
docker run  --rm -it --entrypoint /bin/bash reslp/minibarcoder:5e1dc3b
```

## metaxa
tags: 2.2

WORKDIR is /data
```
docker run --rm -it -v $(pwd):/data reslp/metaxa:2.2
```

## phylobayes
tags: 4.1c

WORDIR is /data
```
docker run --rm -it -v $(pwd):/data reslp/phylobayes:4.1c
```

## phylobayes-mpi
tags: 1.8b, git_dca7bdf, c12c4d7, 3c77656 (commit from Feb 2023)

WORDIR is /data
```
docker run --rm -it -v $(pwd):/data reslp/phylobayes-mpi:c12c4d7
```


## GUIDANCE2
tags: 2.02

WORDIR is /data

```
docker run --rm -it -v $(pwd):/data reslp/guidance2:2.02
```

## raxml-ng
tags: 1.0.0, 1.1.0

WORKDIR is /data

```
docker run --rm -it -v $(pwd):/data reslp/raxml-ng:1.1.0
```

## rmarkdown based on R
The tags versions correspond to the R version.

tags: 3.6.3, 4.0.3 

```
docker run --rm -it -v $(pwd):/data reslp/rmarkdown:3.6.3 Rscript -e "rmarkdown::render('./report.rmd')"
```
## mamba

The Fast Cross-Platform Package Manager (https://github.com/mamba-org/mamba)

This can serve as a base image for things installed using conda and mamba.

tags: 0.5.3, 0.6.0, 0.7.0, 0.15.3, 0.21.1, 0.23.1, 0.25.0, 1.4.1, 1.5.2(broken!)

```
docker run --rm -it reslp/mamba:1.4.1
```

## abyss
tags: 2.2.5

WORKDIR is /data

```
docker run --rm -it reslp/abyss:2.2.5
```

## kmergenie
tags: 1.7051

```
docker run --rm -it reslp/kmergenie:1.7051
```

## ubuntu
This is based on the official ubuntu image and containes build-essential, vim and wget and gdb.

tags: 20.04

```
docker run --rm -it reslp/ubuntu:20.04
```

## cafe
Software for Computational Analysis of gene Family Evolution 

tags: 4.2.1, 5.0.0b2, 5.1.0

WORKDIR /data (5.0.0b2, 5.1.0)
```
docker run --rm -it reslp/cafe:4.2.1
docker run --rm -it -v $(pwd):/data reslp/cafe:5.0.0b2 cafe
```

## hmmer
tags: 3.1b2

```
docker run --rm -it -v $(pwd):/data reslp/hmmer:3.1b2 /bin/bash
```

## emboss
tags: 6.6.0

```
docker run --rm -it -v $(pwd):/data reslp/emboss:6.6.0 /bin/bash
```

## augustus
tags: 3.3.3

```
docker run --rm -it $(pwd):/data reslp/augustus:3.3.3
```

## clustalo
tags: 1.2.4

```
docker run --rm -it reslp/clustalo:1.2.4
```

## bedtools
tags: 2.29.2

```
docker run --rm -it reslp/bedtools:2.29.2
```

## get_organelle
tags: 1.7.1

```
docker run --rm -it reslp/getorganelle:1.7.1 get_organelle_from_reads.py
```

## samtools
tags: 1.9, 1.11

```
docker run --rm -it reslp/samtools:1.11
```

## trimal
tags: 1.4.1

```
docker run --rm -it -v $(pwd):/data reslp/trimal:1.4.1 trimal -in /data/aligned_all_PF00083.fas -out /data/trimmed_all_PF00083.fas -gappyout
```

## table2asn
tags: 20201119, 1.26.678

Program had no good versioning while in beta, till it was released tags are build dates.

```
docker run --rm -it -v $(pwd):/data reslp/table2asn:1.26.678 table2asn -h
```

## tbl2asn
tags: 25.8

This is the legacy version of the NCBI tool.

```
docker run --rm -it -v $(pwd):/data reslp/tbl2asn:25.8
```




## flappie
tags: 4de542f

WORKDIR /data

For convenience the flappy container now also includes ont-fast5-api v. 3.1.6. 
```
# see also ont-fast5-api container on how to convert multi-read fast5 (which is the current read type) to single-read fast5.
docker run --rm -it -v $(pwd):/data reslp/flappie:4de542f flappie --help
docker run --rm -it -v $(pwd):/data reslp/flappie:4de542f flappie /data/fast5_single/0/ > 0.fastq
# using ont-fast5-api:
docker run --rm -it -v $(pwd):/data reslp/flappie:4de542f multi_to_single_fast5 -i /data/fast5_pass -s /data/fast5_single -t 4
```

## ont-fast5-api
tags: 3.1.6

```
# this will convert multi-read fast5 (newer Minion read type) to single-read fast5. See here: https://github.com/nanoporetech/flappie/issues/32
docker run --rm -it -v $(pwd):/data reslp/ont-fast5-api:3.1.6 multi_to_single_fast5 -i /data/fast5_pass -s /data/fast5_single -t 4
```

## mummer
tags: 4.0.0beta2

WORKDIR /data

```
docker run --rm -it -v $(pwd):/data reslp/mummer:4.0.0beta2
```

## bonito
tags: 0.3.2

WORKDIR /data

This container works only with CPUs not GPUs. To basecall Minion reads run this:

```
docker run --rm -i -v $(pwd):/data reslp/bonito:0.3.2 bonito basecaller dna_r9.4.1@v2 --use_openvino --device=cpu /data/fast5_pass 1> basecalls.fasta
docker run --rm -i -v $(pwd):/data reslp/bonito:0.3.2 bonito basecaller dna_r9.4.1@v3 --use_openvino --device=cpu /data/fast5_pass 1> basecalls.fasta
docker run --rm -i -v $(pwd):/data reslp/bonito:0.3.2 bonito basecaller dna_r9.4.1@v2 --use_openvino --device=cpu /data/fast5_pass --fastq 1> basecalls.fastq
```

## interproscan-wrapper
tags: 5.39-77.0, 5.48-83.0

This is a wrapper container for use with interproscan. Interproscan is huge so it can't be properly packed into a container. The tags refer to the tested versions of interproscan which work with this container.
The dockerfile provided may not be the same for different tags. For older versions, Ubuntu 18.04 is used as base due to MobiDBs incompatibility with python >3.8.2 (see here: https://github.com/ebi-pf-team/interproscan/issues/147). 
Newer versions of the container are based on Ubuntu 20.04 which has python 3.8.5. Older versions of interproscan (before 5.39-77.0) may not be compatible because they may have again different requirements.

```
docker pull reslp/interproscan-wrapper:5.39-77.0
```

## pyscaf
https://github.com/lpryszcz/pyScaf

tags: 0.12a4

WORKDIR /data

```
docker run -v $(pwd):/data -it --rm reslp/pyscaf:0.12a
```

## edta
Extensive de-novo TE annotator - https://github.com/oushujun/EDTA

tags: 1.9.6, 1.9.9, 2.0.1

WORKDIR /data

This container uses rmblast 2.9.0 to make the conda installed repeatmasker work, version 2.10.0 won't work. Also, for edta to work in singularity, the configure perl script of RepeatMasker needs to be executed again.

```
docker run --rm -it reslp/edta:1.9.6
```

## bedops
high performance genomic feature operations - https://github.com/bedops/bedops

tags: 2.4.39

WORKDIR /data

```
docker run -v $(pwd):/data --rm -it reslp/bedops:2.4.39
```

## norgal
Mitochondrial genome assembler - https://bitbucket.org/kosaidtu/norgal/src/master/

tags: 1.0

WORDIR /data

```
docker run -v $(pwd):/data --rm -it reslp/norgal:1.0 norgal
docker run -v $(pwd):/data --rm -it reslp/norgal:1.0 norgal.py
```

## seqtk
FASTA and FASTQ processing tool - https://github.com/lh3/seqtk

tags: 1.3

WORKDIR /data

```
docker run -v $(pwd):/data --rm -it reslp/seqtk:1.3 seqtk
```

## sra-tools
The SRA Toolkit and SDK from NCBI is a collection of tools and libraries for using data in the INSDC Sequence Read Archives. - https://github.com/ncbi/sra-tools

tags: 2.20.1, 2.10.9

WORKDIR /data
```
docker run --rm -it reslp/sra-tools:2.10.9 fastq-dump
```

## quicktree
Fast implementation of the neighbour-joining phylogenetic inference method - https://github.com/khowe/quicktree

tags: 2.5

WORKDIR /data

```
docker run --rm -it reslp/quicktree:2.5
```

## bbmap
Splice-aware global aligner for DNA and RNA sequencing reads - https://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/bbmap-guide/

tags: 38.90

WORKDIR /data

```
docker run --rm -it -v $(pwd):/data reslp/bbmap:38.90 bbmap.sh
```

## novoplasty
The organelle assembler and heteroplasmy caller - https://github.com/ndierckx/NOVOPlasty

tags: 4.2

WORKDIR /data

```
docker run --rm -it -v $(pwd):/data reslp/novoplasty:4.2
docker run --rm -it -v $(pwd):/data reslp/novoplasty:4.2 NOVOPlasty.pl
```

## mitos
MITOS is a tool for the annotation of metazoan mitochondrial genomes. - https://gitlab.com/Bernt/MITOS

tags: 2.0.8, 1.0.5

WORKDIR /data

```
docker run --rm -it -v $(pwd):/data reslp/mitos:2.0.8 runmitos.py
docker run --rm -it -v $(pwd):/data reslp/mitos:2.0.8
```

## rphylogenetics
This is an R container with some phylogenetic packages preinstalled. It is based on my rmarkdown container.
Packages: ape, phytools, wesanderson, RColorBrewer

tags: 4.0.3 (refers to used R version), 4.0.3-strict-versions (with more strict versioning of R packages)

WORKDIR /data

```
docker run --rm -it -v $(pwd):/data reslp/rphylogenetics:4.0.3
```

## Orthofinder
OrthoFinder: phylogenetic orthology inference for comparative genomics

tags: 2.3.3, 2.5.2, 2.5.4, 2.5.5

**Note**: The version 2.5.4 of orthofinder seems to have a problem with formatting trees. This seems to be related to [this](https://github.com/davidemms/OrthoFinder/issues/782). Although in my container the python version is 3.9. Up to stride orthofinder seems to run correctly.

WORKDIR /data

```
docker run --rm -it -v $(pwd):/data reslp/orthofinder:2.5.5
```

## Jupyter-Lab
This Jupyter Lab Notebook also has Tensorflow, Scikitlearn, pandas and numpy installed.
The tags refer to the tensorflow version.

tags: 2.4.1

To start a server (either jupyter lab or jupyter notebook):

```
docker run -d --rm -p 8888:8888 reslp/jupyter-tf:2.4.1 jupyter-lab
docker run -d --rm -p 8888:8888 reslp/jupyter-tf:2.4.1
```

This will return a hash. Run docker logs on that hash to see the access token. eg:

```
docker logs 406f5662e0a0fb544ed80d95054642e55439f6eb18fcae30487c8488159e36e4
```



## CUPP
CUPP is a new approach for sequence analysis based on the hypothesis that conserved peptide patterns can be used to identify similarities between proteins. -
https://www.bioengineering.dtu.dk/english/ResearchNy/Research-Sections/Section-for-Protein-Chemistry-and-Enzyme-Technology/Enzyme-Technology/CUPP

tags: 1.0.14

```
docker run --rm -it reslp/cupp:1.0.14 python /software/CUPP_v1.0.14/CUPPprediction.py --help

```

## ima2
IM is a program, written with Rasmus Nielsen, for the fitting of an isolation model with migration to haplotype data drawn from two closely related species or populations. - https://bio.cst.temple.edu/~tuf29449/software

tags: 8.27.12

WORKDIR /data


```
docker run --rm -it reslp/ima2:8.27.12
```


## alitv
AliTV is a program to visualize whole genome alignments: https://alitv.readthedocs.io/en/latest/manual.html

tags: d06bd45

WORKDIR /data

```
docker run --rm -it -v $(pwd):/data reslp/alitv:d06bd45
```

## speedtest-cli
This a command line program which performs network speedtests: https://www.speedtest.net/de/apps/cli

tags: 2.1.2, 2.1.3

```
docker run --rm reslp/speedtest-cli:2.1.3 speedtest-cli --csv >> /home/pi/speedtest/speedtest_data.txt
```

## bowtie2
Bowtie 2 is an ultrafast and memory-efficient tool for aligning sequencing reads to long reference sequences. - https://github.com/BenLangmead/bowtie2

tags: 2.3.5

WORKDIR /data

```
docker run reslp/bowtie2:2.3.5 bowtie2
```

## fastuniq
an ultrafast de novo duplicates removal tool for paired short DNA sequences - https://sourceforge.net/projects/fastuniq/

tags: 1.1

WORKDIR /data

```
docker run --rm -it reslp/fastuniq:1.1 fastuniq
```

## flash
Fast Length Adjustment of SHort reads - https://ccb.jhu.edu/software/FLASH/ 

tags: 1.2.11

WORKDIR /data

```
docker run --rm -it reslp/flash:1.2.11 flash
```

## bamutil
Programs that perform operations on SAM/BAM files, all built into a single executable, bam. - https://genome.sph.umich.edu/wiki/BamUtil

tags: 1.0.15

WORKDIR /data

```
docker run --rm -it reslp/bamutil:1.0.15 bam
```

## picard
Java tools for working with NGS data in the BAM format - http://broadinstitute.github.io/picard/

tags: 2.18.7

WORKDIR /data

```
docker run --rm -it reslp/picard:2.18.7 picard
```

## gatk
Genome analysis toolkit - https://gatk.broadinstitute.org/hc/en-us

tags: 3.8

WORKDIR /data

Hint: Executable name is gatk3 NOT gatk!

```
docker run --rm -it reslp/gatk:3.8 gatk3
```

## angsd
ANGSD is a software for analyzing next generation sequencing data - http://www.popgen.dk/angsd/index.php/ANGSD

tags: 0.935

```
docker run --rm -it reslp/angsd:0.935 angsd
```

## purge_haplotigs

Pipeline to help with curating heterozygous diploid genome assemblies - https://bitbucket.org/mroachawri/purge_haplotigs/src/master/

tags: 1.1.1

WORKDIR /data

```
docker run --rm -it $(pwd):/data reslp/purge_haplotigs:1.1.1
```

## fasta3

The FASTA package - protein and DNA sequence similarity searching and alignment programs

tags: 36.3.8

WORKDIR /data

*Hint:* All the programs of the FASTA3 package have 36 appended to the end of their name. So eg. glsearch is glsearch36 in this container.

```
docker run --rm -it -v $(pwd):/data reslp/fasta3:36.3.8 glsearch36 /data/query.fas /data/subset.fasta
```

# macse

MACSE: Multiple Alignment of Coding SEquences Accounting for Frameshifts and Stop Codons.

WORKDIR /data

tags: 1.2, 2.05

```
docker run --rm -it reslp/macse:1.2 macse
docker run --rm -it reslp/macse:2.05 macse -h
```

# julia

julia - The Julia Programming Language

This container includes julia and a number of packages from BioJulia and packages for data handling, plotting and report generation.

tags: 1.6.3

WORDKIR /data

```
docker run --rm -it reslp/julia:1.6.3 julia
```

# nim

Nim programming lanuage

tags: 1.4.8, 1.6.6, 1.6.12

```
docker run --rm -it reslp/nim:1.6.6
```

# mosdepth

fast BAM/CRAM depth calculation for WGS, exome, or targeted sequencing (https://github.com/brentp/mosdepth)

tags: 0.3.3

```
docker run --rm -it reslp/mosdepth:0.3.3
```

# viralmsa

Reference-guided multiple sequence alignment of viral genomes  (https://github.com/niemasd/ViralMSA)

Note: There is an offical container, but it lacks several dependencies and does not work (tested March 21, 2022) hence the container here which largely builds on the official version.

tags: 1.1.20

WORKDUT /data

```
docker run --rm -it reslp/viralmsa:1.1.20 ViralMSA.py -h
```

# dentist

Close assembly gaps using long-reads at high accuracy. (https://github.com/a-ludi/dentist)

Note: This container also includes snakemake 6.0.2 so it is possible to run the whole dentist pipeline within the container.

tags: 3.0.0

```
docker run --rm -it reslp/dentist:3.0.0 dentist
```

# muscle

Multiple sequence alignment with top benchmark scores scalable to thousands of sequences. (https://github.com/rcedgar/muscle)

tags: 5.1

WORKDIR /data

```
docker run --rm -it reslp/muscle:5.1 muscle
```

# necat

NECAT is an error correction and de-novo assembly tool for Nanopore long noisy reads. (https://github.com/xiaochuanle/NECAT)

tags: 0.0.1_update20200803

WORKDIR /data

The script in the container is called necat not necat.pl.

```
docker run --rm -it reslp/necat:0.0.1_update20200803 necat
```

# colord

A versatile compressor of third generation sequencing reads. (https://github.com/refresh-bio/colord)

tags: 1.1.0

WORKDIR /data

```
docker run --rm -it reslp/colord:1.1.0 colord
```

# bmge

BMGE (Block Mapping and Gathering with Entropy) is a program that selects regions in a multiple sequence alignment that are suited for phylogenetic inference. (https://bioweb.pasteur.fr/packages/pack@BMGE@1.12/)

tags: 1.12

WORKDIR /data

```
docker run --rm -it reslp/bmge:1.12 bmge
```

# ete3

A Python framework for the analysis and visualization of trees. - http://etetoolkit.org/

tags: 3.1.2

Hint: External programs don't work due to the omnipresent bug in ete3: https://github.com/etetoolkit/ete/issues/444

WORKDIR /data

```
docker run --rm -it reslp/ete3:3.1.2
```

# dendropy

DendroPy Phylogenetic Computing Library - https://dendropy.org/

tags: 4.5.2

WORKDIR /data

```
docker run --rm -it reslp/dendropy:4.5.2
```

# phylopy

A container with dendropy, ete3 and pandas based on python 3.7.

tags: 1

WORKDIR /data

```
docker run --rm -it reslp/phylopy:1
```

# ffsend

Easily and securely share files from the command line. A fully featured Firefox Send client. https://github.com/timvisee/ffsend

tags: 0.2.76

WORKDIR /data

```
docker run --rm -it -v $(pwd):/data reslp/ffsend:0.2.76 ffsend upload myfile.txt
```

# magic-wormhole

Get things from one computer to another, safely. https://github.com/magic-wormhole/magic-wormhole

tags: 0.12.0

WORKDIR /data

```
docker run --rm -it reslp/magic-wormhole:0.12.0
```

# pandoc

A universal document converter. https://pandoc.org/

tags: 2.18, 3.1.3

WORKDIR /data

```
docker run --rm -it -v $(pwd):/data reslp/pandoc:2.18 pandoc test.md -o test.pdf
```

# beast

Bayesian Evolutionary Analysis by Sampling Trees - https://github.com/CompEvol/beast2

tags: 2.6.7

WORKDIR /data

```
docker run --rm -it -v $(pwd):/data reslp/beast:2.6.7
```

# pigz

Parallel implementation of gzip

tags: 2.6

WORKDIR /data

```
docker run --rm -it -v $(pwd):/data reslp/pigz:2.6
```

# rename

Rename files

tags: 1.10-1

WORKDIR /data

```
docker run --rm -it -v $(pwd):/data reslp/rename:1.10-1 
```

# pyinstaller

PyInstaller bundles a Python application and all its dependencies into a single package.

tags: 5.1, 6.1

WORKDIR /data

```
docker run --rm -it -v $(pwd):/data reslp/pyinstaller:5.1
```

# pear

PEAR is an ultrafast, memory-efficient and highly accurate pair-end read merger

tags: 0.9.6

WORKDIR /data

```
docker run --rm -it reslp/pear:0.9.6
```

# taxonkit

A Practical and Efficient NCBI Taxonomy Toolkit - https://github.com/shenwei356/taxonkit

tags: 0.12.0

WORKDIR /data

```
docker run --rm -it reslp/taxonkit:0.12.0
```

# exonerate

Exonerate - A generic tool for pairwise sequence comparison / alignment

tags: 2.4.0

```
docker run --rm -it reslp/exonerate:2.4.0
```

# oma-standalone

Standalone version of the OMA Orthology Matrix - https://omabrowser.org/oma/home/

tags: 2.5.0

WORKDIR /data

```
docker run --rm -it reslp/oma-standalone:2.5.0 OMA
```

# repeatmodeler-repeatmasker

This container contains RepeatModeler and Repeatmasker including all dependencies (hopefully).

The tag consists of the repeatmodeler version and the repeatmasker version:

tags: 2.0.3-4.1.3

WORKDIR /data

```
docker run --rm -it reslp/repeatmodeler-repeatmasker:2.0.3-4.1.3
```

# mutt

Open source CLI email client

tags: 3.2.6

```
docker run --rm -it --user $(id -u):$(id -g) -v /etc/passwd:/etc/passwd:ro -v /etc/group:/etc/group:ro -v $HOME/:/home/$USER/ -e HOME=/home/$USER reslp/mutt:3.2.6 mutt
```

# PAML

Phylogenetic Analysis by Maximum Likelihood - https://github.com/abacus-gene/paml

tags: 4.9 (conda version, comes with biopython1.81 scipy), 4.10.6

WORKDIR /data

```
docker run --rm -it reslp/paml:4.10.6
```

# pass

the standard unix password manager - https://www.passwordstore.org/

HOME /data

GnuPG and pass needs to be configured correctly for this to work.

```
# linux
docker run --rm -it -v $HOME:/data --user $(id -u):$(id -g) reslp/pass:1.7.4

# mac
docker run --rm -it -v $HOME:/data reslp/pass:1.7.4
```

# snakemake

The Snakemake workflow management system - https://snakemake.readthedocs.io/en/stable/

tags: 6.0.2

```
docker run --rm -it reslp/snakemake:6.0.2
```

# nimbio_plus

My own nimbio library together with nim, glob, ggplotnim, clapfn, nimx 

tags: 0.0.1, 0.0.2, 0.0.3

```
docker run --rm -it reslp/nimbio_plus:0.0.2
```

# fastcodeml

FastCodeML is a software to infer positive selection along positions of a protein coding gene using the Branch-Site model of evolution - https://bioweb.pasteur.fr/packages/pack@FastCodeML@1.1.0

tags: 1.1.0

```
docker run --rm -it reslp/fastcodeml:1.1.0 fast
```

# diamond

Accelerated BLAST compatible local sequence aligner. - https://github.com/bbuchfink/diamond

tags: 0.9.22, 2.0.7, 2.1.7

WORKDIR /data

```
docker run --rm -it reslp/diamond:2.0.7
```

# MCL


MCL - a cluster algorithm for graphs - https://github.com/micans/mcl

tags: 14.137, 22.282

WORKDIR /data

```
docker run --rm -it reslp/mcl:14.137
```

# cafe-viz

Tools for visualizing cafe output. Currently includes:  
cafeplotter and CAFE_Fig  
https://github.com/moshi4/CafePlotter  
https://github.com/LKremer/CAFE_fig

tags: 1.0

WORKDIR /data

```
docker run --rm -it reslp/cafe-viz:1.0 cafeplotter
docker run --rm -it reslp/cafe-viz:1.0 python3 /usr/bin/CAFE_fig.py
```


# phykit

a UNIX shell toolkit for processing and analyzing multiple sequence alignments and phylogenies - https://github.com/jlsteenwyk/phykit

tags: 1.11.15

WORKDIR /data

```
docker run --rm -it reslp/phykit:1.11.15 phykit -h
```

# genome-viz

a python script to visualize genome assemblies - https://github.com/reslp/genomics

tags: 1.0

```
docker run --rm -it reslp/genome-viz:1.0 genome_viz.py
```

# bibcure

Bibcure helps in boring tasks by keeping your bibfile up to date and normalized...also allows you to easily download all papers inside your bibtex. - https://github.com/bibcure/bibcure

tags: 0.3.0

WORKDIR /data

```
docker run --rm -it -v $(pwd):/data reslp/bibcure:0.3.0 doi2bib -i literature_only_dois.txt -o bib.txt
```

# amplicon_sorter

Sorts amplicons from Nanopore sequencing data based on similarity - https://github.com/avierstr/amplicon_sorter

Tags are based on git commits, since there is no versioning. 
tags: 33f0d38

```
docker run --rm -it -v $(pwd):/data reslp/amplicon_sorter:33f0d38 amplicon_sorter.py
```

# typst

A new markup-based typesetting system that is powerful and easy to learn. - https://github.com/typst/typst

tags: 0.7.0, 0.9.0

WORKDIR /data

```
docker run --rm -it -v $(pwd):/data reslp/typst:0.9.0 typst watch /data/document.typ
```

# gnuplot

Gnuplot is a portable command-line driven graphing utility for Linux, OS/2, MS Windows, OSX, VMS, and many other platforms. - http://www.gnuplot.info/

tags: 5.4.9

```
docker run --rm -it reslp/gnuplot:5.4.9
```

# ncbi-acc-download

Download files from NCBI Entrez by accession - https://github.com/kblin/ncbi-acc-download/

tags: 0.2.8

WORKDIR /data

```
docker run --rm -it -v $(pwd):/data reslp/ncbi-acc-download:0.2.8 ncbi-acc-download
```

# stringtie

Transcript assembly and quantification for RNA-Seq - https://github.com/gpertea/stringtie

tags: 2.2.1

WORKDIR /data

```
docker run --rm -it -v $(pwd):/data reslp/stringtie:2.2.1 stringtie
```

# imagemagick

ImageMagick is a free, open-source software suite, used for editing and manipulating digital images. - https://imagemagick.org/index.php

tags: 7.1.0_46

WORKDIR /data

```
docker run --rm -it -v $(pwd):/data reslp/imagemagick:7.1.0_46 convert
```

# tectonic

Tectonic is a modernized, complete, self-contained TeX/LaTeX engine, powered by XeTeX and TeXLive. - https://tectonic-typesetting.github.io/en-US/

tags: 0.14.1

WORKDIR /data

```
docker run --rm -it -v $(pwd):/data reslp/tectonic:0.14.1 tectonic
```

# multiqc

Aggregate results from bioinformatics analyses across many samples into a single report  - https://multiqc.info/

tags: 1.6, 1.17

WORKDIR /data

```
docker run --rm -it reslp/multiqc:1.17 multiqc
```

# keras

These containers come with the deep-learing library KERAS together with several additional libraries such as pandas, numpy, sci-kit and scikit-learn, matplotlib etc.

tags: 3.0.0-torch, 3.0.0-tf

```
docker run --rm -it reslp/keras:3.0.0-tf
```

# treeshrink

TreeShrink is an algorithm for detecting abnormally long branches in one or more phylogenetic trees. - https://github.com/uym2/TreeShrink

tags: 1.3.9

WORKDIR /data

```
docker run --rm -it reslp/treeshrink:1.3.9 run_treeshrink -h
```

# ipyrad

Interactive assembly and analysis of RAD-seq data sets - https://github.com/dereneaton/ipyrad

tags: 0.9.93

WORKDIR /data

```
docker run --rm -it reslp/ipyrad:0.9.93 ipyrad
```

# ncbi-datasets

NCBI Datasets command-line clients - https://www.ncbi.nlm.nih.gov/datasets/docs/v2/download-and-install/

tags: 16.25.0

```
docker run --rm -it reslp/datasets:16.25.0 datasets
```


