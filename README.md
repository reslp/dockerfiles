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

Important: The adapter file is not part of the docker image. It needs to passed together with the mounted directory.


fastqc 0.11.7
===============================

`docker pull reslp/fastqc`

maxbin 2.2.6
===============================

`docker pull reslp/maxbin`
