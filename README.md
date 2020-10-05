# dockerfiles

This repository contains Dockerfiles I accumulated for different purposes.

Some of the images are bare-bone and should maybe be changed when incorporated into analysis pipelines.

The **working directory** in these containers is generally `/data`. This can be used as a bindpoint. Downloaded software is generally kept in `/software`. Some older containers here may however not follow this rules.

The **tags** given for each container correspond the the available versions of different programs.

**Currently available containers:**

[astral](#astral)
[biopython_plus](#biopython_plus)
[blobtools](#blobtools)
[busco](#busco)
[concoct](#concoct)
[consent](#consent)
[deeploc](#deeploc)
[eggnog-mapper](#eggnog-mapper)
[fastqc](#fastqc)
[fmlrc](#fmlrc)
[guidance2](#guidance2)
[iqtree](#iqtree)
[mafft](#mafft)
[mamba](#mamba)
[maxbin](#maxbin)
[metabat](#metabat)
[metaxa](#metaxa)
[miniBarcoder](#miniBarcoder)
[mrbayes](#mrbayes)
[pal2nal](#pal2nal)
[phylobayes](#phylobayes)
[phylobayes-mpi](#phylobayes-mpi)
[phylo-scripts](#phylo-scripts)
[raxml-ng](#raxml-ng)
[rmarkdown](#rmarkdown-based-on-r)
[scrape_cazy](#scrape_cazy)
[singularity](#singularity)
[spades](#spades)
[trimmomatic](#trimmomatic)
[quast](#quast)
[ncbi-blast](#ncbi-blast)






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
tags: 3.12.0, 3.13.0

`docker pull reslp/spades:3.13.0`



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
tags: 3.0.2

```
docker run --rm reslp/busco:3.0.2
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
tags: 2.9.0

To create a blast database and run a blastp search:
`my_sequences.fasta`and `sequences_to_search.fa`need to be in the current working directory. The blastdb will also be created there.

```
docker run --rm -v $(pwd):/files/ reslp/ncbi-blast:2.9.0 makeblastdb -in /files/my_sequences.fasta -dbtype 'prot' -hash_index -out my_blastdb

docker run --rm -v $(pwd):/files/ reslp/ncbi-blast:2.9.0 blastp -db /files/my_blastdb -query /files/sequences_to_search.fa -outfmt 6 > blast_results.txt

```

## singularity
tags: 3.4.1

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
docker run --rm reslp/scrape_vim cazy:1
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
tags: 1.8b

WORDIR is /data
```
docker run --rm -it -v $(pwd):/data reslp/phylobayes:1.8b
```


## GUIDANCE2
tags: 2.02

WORDIR is /data

```
docker run --rm -it -v $(pwd):/data reslp/guidance2:2.02
```

## raxml-ng
tags: 1.0.0

WORKDIR is /data

```
docker run --rm -it -v $(pwd):/data reslp/raxml-ng:1.0.0
```

## rmarkdown based on R
The tags versions correspond to the R version.

tags: 3.6.3 

WORKDIR is /data

```
docker run --rm -it -v $(pwd):/data reslp/rmarkdown:3.6.3 Rscript -e "rmarkdown::render('./report.rmd')"
```
## mamba
This can serve as a base image for things installed using conda and mamba.

tags: 0.5.3

```
docker run --rm -it reslp/mamba:0.5.3
```

