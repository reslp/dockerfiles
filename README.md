dockerfiles
=========

This repository contains docker files I accumulated for different purposes.



IQ-tree
===================

Installation:
1. Install [Docker](http://www.docker.com) or if you are on Windows use [Docker Toolbox](https://www.docker.com/products/docker-toolbox).
2. Install the phylo-scripts container: `docker build github.com/reslp/dockerfiles/iqtree --tag iqtree`. If you are using Docker Toolbox you will need to start the Virtual Machine first by clicking on the start.sh file in the Docker Toolbox Folder.
3. To execute `docker run -v <full_path_to_input_files>:/input_files/ iqtree`

Important: <full_path_to_input_files> has to be the full path to your files starting from the root of your file system. This folder should contain alignment and or partitiuon files for use with iqtree.

An example command for IQ-tree may look like this:
`docker run -v /user/home/philipp/alignments/:/input/ iqtree -nt 1 -s /input/alignment.fas`




