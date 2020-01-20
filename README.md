dockerfiles
=========

This repository contains Docker files I accumulated for different purposes. For these to work, first install [Docker](http://www.docker.com) or if you are on Windows install [Docker Toolbox](https://www.docker.com/products/docker-toolbox).

Some of the images are bare-bone and should maybe be changed when incorporated into analysis pipelines.



IQ-tree
===================



1. Install [Docker](http://www.docker.com) or if you are on Windows use [Docker Toolbox](https://www.docker.com/products/docker-toolbox).
2. Download the docker file: `curl https://raw.githubusercontent.com/reslp/dockerfiles/master/iqtree/Dockerfile > Dockerfile`. If you are using Docker Toolbox you will need to start the Virtual Machine first by clicking on the start.sh file in the Docker Toolbox Folder.
3. In the same directory you can now build the container: `docker build --tag iqtree .`
4. To execute `docker run -v <full_path_to_input_files>:/input_files/ iqtree`

Important: <full_path_to_input_files> has to be the full path to your files starting from the root of your file system. This folder should contain alignment and or partitiuon files for use with iqtree.

An example command for IQ-tree may look like this:
`docker run -v /user/home/philipp/alignments/:/input/ iqtree -nt 1 -s /input/alignment.fas`
The most basic command which will show you an overview of IQ-trees command line paramters is:
`docker run iqtree`

MrBayes
===================


1. Install [Docker](http://www.docker.com) or if you are on Windows use [Docker Toolbox](https://www.docker.com/products/docker-toolbox).
2. Download the docker file: `curl https://raw.githubusercontent.com/reslp/dockerfiles/master/mrbayes/Dockerfile > Dockerfile`. If you are using Docker Toolbox you will need to start the Virtual Machine first by clicking on the start.sh file in the Docker Toolbox Folder.
3. In the same directory you can now build the container: `docker build --tag mrbayes .`
4. To execute `docker run -v <path_to_local_folder>:/input/ -it mb_ubuntu /sbin/my_init -- mpirun --allow-run-as-root -np 8 mb /input/<nexus_file.nex>`

Important: <path_to_local_folder> has to be the full path to your files starting from the root of your file system. This folder should contain the NEXUS file used to run mrbayes.

Phylo-Scripts
===================

1. Install [Docker](http://www.docker.com) or if you are on Windows use [Docker Toolbox](https://www.docker.com/products/docker-toolbox).
2. Install the phylo-scripts container: `docker build github.com/reslp/phylo-scripts --tag phylo-scripts`. If you are using Docker Toolbox you will need to start the Virtual Machine first by clicking on the start.sh file in the Docker Toolbox Folder.
3. To create concatenated alignments execute `docker run -v <full_path_to_input_files>:/input_files/ phylo-scripts`

Important: <full_path_to_input_files> has to be the full path to your files starting from the root of your file system. This folder should contain all the single locus files which should be used in the multi gene alignment and the file with Sequence IDs (`IDs_used_for_tree.txt`). 

Spades 3.12.0
================================

`docker pull reslp/spades`



trimmomatic 0.38
===============================
`docke3r pull reslp/trimmomatic`

Example command:

```
docker run -t -v $(pwd):/data reslp/trimmomatic trimmomatic PE -phred33 -threads 8 /data/pair1.fq.gz /data/pair2.fq.gz /data/trimmed_pair1.fq.gz /data/unpaired_pair1.fq.gz /data/trimmed_pair2.fq.gz /data/unpaired_pair2.fq.gz HEADCROP:20 ILLUMINACLIP:/data/adapters.fa:2:30:10 LEADING:30 TRAILING:30 SLIDINGWINDOW:4:15 MINLEN:80
```
or

```
docker run --user $(id -u):$(id -g) --rm -v $(pwd):/data reslp/trimmomatic trimmomatic
```

Important: The adapter file is not part of the docker image. It needs to passed together with the mounted directory.


fastqc 0.11.7
===============================

Get the docker image:

`docker pull reslp/fastqc`

Example command for running the container:
```
docker run -v $(pwd):/data --rm reslp/fastqc fastqc /data/illumina_reads_5_R1_trimmed.fq -o /data
```
or 

```
docker run --user $(id -u):$(id -g) --rm -v $(pwd):/data reslp/fastqc fastqc -h
```
funannotate 1.7.2
======


### External dependencies:

Things not included in the Container, which need to be installed manually due to licensing incompatibility:

1. Signal-P 4.1: [www.cbs.dtu.dk/services/SignalP/](http://www.cbs.dtu.dk/services/SignalP/) 
2. GeneMark-ES: [exon.gatech.edu/GeneMark/](http://exon.gatech.edu/GeneMark/)
3. Repeatmasker libraries from RepBase

The way to get these programs into the container is to place them into a folder and then mount this folder to a specific point in the container by adding a certain flag to the docker run command:

	-v /local/location_of_programs:/root/

The docker container is set up in such a way as that it searches for specific folders in the root directory and adds them to the path. This way, funannotate running inside the container finds the desired programs. Currently the container is etup to add the following folders o the PATH hence the version names:

	/root/signalp-4.1
	/root/gm_et_linux_64

Also it expects the GeneMark license key file in /root/.




### Setting up the funannotate database:

Funannotate relies on a larger dataset of different kinds which it uses to add functional annotations, train gene finders etc. This database needs to be created with 

	funannotate setup

The docker container has the location folder for the database hardcoded (by setting the FUNANNOTATE_DB environment variable) to /root/database/. This folder needs to be overwritten with the local location of the database when funannotate is run. This is again done when docker run is invoked:

	-v /local/location/of/database:/root/database



### Where is my data?

Data will be stored in /data which can be mounted from an external folder as well like so:

	-v $(pwd):/data




### An example command for the funannotate docker container:

This command mount external dependencies and a database folder:

	docker run --rm -it -v $(pwd)/external:/root -v $(pwd)/database:/root/database -v $(pwd):/data reslp/funannotate check







maxbin 2.2.6
===============================

`docker pull reslp/maxbin`

`docker run --rm reslp/maxbin`

concoct 1.1
===============================

`docker pull reslp/concoct`

`docker run --rm reslp/concoct`

A complete example filtering a metagenome. The `metagenome.fasta` file should be in the current working directory.

```
docker run -t --user $(id -u):$(id -g) -v $(pwd):/data/ --rm reslp/concoct cut_up_fasta.py /data/metagenome.fasta -c 10000 -o 0 --merge_last -b /data/metagenome.fasta_contigs_10K.bed > metagenome.fasta_contigs_10K.fa
docker run -t --user $(id -u):$(id -g) -v $(pwd):/data/ --rm reslp/concoct concoct_coverage_table.py /data/metagenome.fasta_contigs_10K.bed /data/metagenome.fasta.bam > concoct_coverage_table.tsv
docker run -t --user $(id -u):$(id -g) -v $(pwd):/data/ --rm reslp/concoct concoct --composition_file /data/metagenome.fasta_contigs_10K.fa --coverage_file /data/concoct_coverage_table.tsv -b /data/concoct/metagenome.fasta_concoct --threads 24
docker run --user $(id -u):$(id -g) -v $(pwd):/data/ --rm reslp/concoct merge_cutup_clustering.py /data/concoct/metagenome.fasta_concoct_clustering_gt1000.csv > concoct/metagenome.fasta_concoct_clustering_merged.csv
docker run -t --user $(id -u):$(id -g) -v $(pwd):/data/ --rm reslp/concoct extract_fasta_bins.py /data/metagenome.fasta /data/concoct/metagenome.fasta_concoct_clustering_merged.csv --output_path /data/concoct/bins
```

metabat 2.13
===============================

`docker pull reslp/metabat`

`docker run --rm reslp/metabat`

busco 3.0.2
===============================

`docker pull reslp/busco`

`docker run --rm reslp/busco`

blobtools
=======
`docker pull reslp/blobtools`

These commands expect that blast results are alread present:

```
docker run -t --user $(id -u):$(id -g) -v $(pwd):/data/ --rm reslp/blobtools create -i /data/metagenome.fasta -b /data/metagenome.fasta.bam -t /data/blobtools/metagenome.fasta_diamond_matches_formatted -o /data/blobtools/metagenome.fasta_blobtools
docker run -t --user $(id -u):$(id -g) -v $(pwd):/data/ --rm reslp/blobtools view -i /data/blobtools/metagenome.fasta_blobtools.blobDB.json -o /data/blobtools/
docker run -t --user $(id -u):$(id -g) -v $(pwd):/data/ --rm reslp/blobtools plot -i /data/blobtools/metagenome.fasta_blobtools.blobDB.json -o /data/blobtools/
```


QUAST 
=========
`docker pull reslp/quast`

Run QUAST on all fasta files in current directory:

`docker run --rm -t -v $(pwd):/data/ reslp/quast quast.py $(ls -d -1 "*.fasta" | sed 's/^/\/data\//') -o /data/quast --silent
`

ncbi-blast 2.9.0
==========
`docker pull reslp/ncbi-blast`

To create a blast database and run a blastp search:
`my_sequences.fasta`and `sequences_to_search.fa`need to be in the current working directory. The blastdb will also be created there.

```
docker run --rm -v $(pwd):/files/ reslp/ncbi-blast makeblastdb -in /files/my_sequences.fasta -dbtype 'prot' -hash_index -out my_blastdb

docker run --rm -v $(pwd):/files/ reslp/ncbi-blast blastp -db /files/my_blastdb -query /files/sequences_to_search.fa -outfmt 6 > blast_results.txt

```

