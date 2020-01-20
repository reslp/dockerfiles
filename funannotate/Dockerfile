FROM ubuntu:18.04

MAINTAINER <resl.philipp@bio.lmu.de>

WORKDIR /software

#these two layers should take care of all the python and perl dependencies:
RUN apt-get update && apt-get -y upgrade && apt-get -y install wget && apt-get -y install python2.7 python-pip && pip install funannotate
RUN apt-get install -y bioperl && apt-get install -y cpanminus && cpanm File::Which Hash::Merge JSON Logger::Simple Parallel::ForkManager Scalar::Util::Numeric Text::Soundex

#Software dependencies:
#CodingQuarry
RUN wget https://sourceforge.net/projects/codingquarry/files/CodingQuarry_v2.0.tar.gz && tar xvfz CodingQuarry_v2.0.tar.gz && cd CodingQuarry_v2.0 && make 
ENV QUARRY_PATH="/software/CodingQuarry_v2.0/QuarryFiles"
ENV PATH="/software/CodingQuarry_v2.0:$PATH"

#Trinity
RUN apt-get -y install cmake && apt-get -y install zlib1g-dev && wget https://github.com/trinityrnaseq/trinityrnaseq/releases/download/v2.8.6/trinityrnaseq-v2.8.6.FULL.tar.gz && tar xvfz trinityrnaseq-v2.8.6.FULL.tar.gz && cd trinityrnaseq-v2.8.6 && make
ENV TRINITY_HOME="/software/trinityrnaseq-v2.8.6"
ENV PATH="/software/trinityrnaseq-v2.8.6:$PATH"

#Augustus
#RUN apt-get install -y augustus
#ENV AUGUSTUS_CONFIG_PATH="/usr/share/augustus/config/"
RUN apt-get install -y libboost-iostreams-dev zlib1g-dev libgsl-dev libboost-graph-dev libboost-all-dev libsuitesparse-dev liblpsolve55-dev libsqlite3-dev libgsl-dev libboost-graph-dev libboost-all-dev libsuitesparse-dev liblpsolve55-dev libmysql++-dev libbamtools-dev libboost-all-dev && wget https://github.com/Gaius-Augustus/Augustus/archive/3.3.2.tar.gz && tar xvfz 3.3.2.tar.gz && cd Augustus-3.3.2/auxprogs && sed -i 's/cd bam2wig; make;/#cd bam2wig; make;/g' Makefile && sed -i 's/cd bam2wig; make clean;/#cd bam2wig; make clean;/g' Makefile && cd .. && make 
ENV AUGUSTUS_CONFIG_PATH="/software/Augustus-3.3.2/config"
ENV PATH="/software/Augustus-3.3.2/bin:$PATH"

#bamtools
RUN apt-get install -y bamtools

#BLAT
RUN mkdir blat && cd blat && wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/blat/blat && chmod 755 ./blat
ENV PATH="/software/blat/:$PATH"

#FASTA36
RUN wget http://faculty.virginia.edu/wrpearson/fasta/fasta36/fasta-36.3.8g.tar.gz && tar xvfz fasta-36.3.8g.tar.gz && cd fasta-36.3.8g/src && make -f ../make/Makefile.linux_sse2 all && cp ../bin/fasta36 /usr/local/bin/fasta

#diamond
RUN wget http://github.com/bbuchfink/diamond/releases/download/v0.9.29/diamond-linux64.tar.gz && tar xvfz diamond-linux64.tar.gz && mv diamond /usr/local/bin

#glimmerhmm
# this is currently commented out, it seems to halt the predict step of funannotate for some strange reason
#RUN wget https://ccb.jhu.edu/software/glimmerhmm/dl/GlimmerHMM-3.0.4.tar.gz && tar xvfz GlimmerHMM-3.0.4.tar.gz && cd GlimmerHMM/bin && cp glimmerhmm_linux_x86_64 glimmerhmm
#ENV PATH="/software/GlimmerHMM/bin:$PATH"
#ENV PATH="/software/GlimmerHMM/train:$PATH"

#GMAP
RUN wget http://research-pub.gene.com/gmap/src/gmap-gsnap-2019-09-12.tar.gz && tar xvfz gmap-gsnap-2019-09-12.tar.gz && cd gmap-2019-09-12/ && ./configure && make
ENV PATH="/software/gmap-2019-09-12/src:$PATH"

#Java
RUN apt-get install -y default-jre

#HISAT2
RUN apt-get install -y hisat2

#minimap2
RUN wget https://github.com/lh3/minimap2/releases/download/v2.17/minimap2-2.17_x64-linux.tar.bz2 && tar -jxvf minimap2-2.17_x64-linux.tar.bz2
ENV PATH="/software/minimap2-2.17_x64-linux:$PATH"

#kallisto
RUN wget https://github.com/pachterlab/kallisto/releases/download/v0.46.1/kallisto_linux-v0.46.1.tar.gz && tar xvfz kallisto_linux-v0.46.1.tar.gz
ENV PATH="/software/kallisto:$PATH"

#Proteinortho
RUN wget https://gitlab.com/paulklemm_PHD/proteinortho/-/archive/v6.0.12/proteinortho-v6.0.12.tar.gz && tar xvfz proteinortho-v6.0.12.tar.gz
ENV PATH="/software/proteinortho-v6.0.12:$PATH"

#pslCDnaFilter
# this does currently not select a specific version!! It uses always the latest version
RUN wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/pslCDnaFilter && chmod 755 pslCDnaFilter && mv pslCDnaFilter /usr/local/bin

#salmon
RUN wget https://github.com/COMBINE-lab/salmon/releases/download/v1.1.0/salmon-1.1.0_linux_x86_64.tar.gz && tar xvfz salmon-1.1.0_linux_x86_64.tar.gz
ENV PATH="/software/salmon-latest_linux_x86_64/bin:$PATH"

#snap
RUN apt-get install -y git && git clone https://github.com/KorfLab/SNAP.git && cd SNAP && make 
ENV ZOE="/software/SNAP/Zoe"
ENV PATH="/software/SNAP:$PATH"

#stringtie
RUN wget http://ccb.jhu.edu/software/stringtie/dl/stringtie-2.0.6.Linux_x86_64.tar.gz && tar xvfz stringtie-2.0.6.Linux_x86_64.tar.gz
ENV PATH="/software/stringtie-2.0.6.Linux_x86_64:$PATH"

#tRNA-Scan
RUN wget http://trna.ucsc.edu/software/trnascan-se-2.0.5.tar.gz && tar xvfz trnascan-se-2.0.5.tar.gz && cd tRNAscan-SE-2.0 && ./configure && make && make install 
RUN wget http://eddylab.org/infernal/infernal-1.1.3-linux-intel-gcc.tar.gz && tar xvfz infernal-1.1.3-linux-intel-gcc.tar.gz && cd infernal-1.1.3-linux-intel-gcc/binaries && cp * /usr/local/bin

#tantan
RUN apt-get install -y tantan

#trimal
RUN wget http://trimal.cgenomics.org/_media/trimal.v1.2rev59.tar.gz && tar xvfz trimal.v1.2rev59.tar.gz && cd trimAl/source && make
ENV PATH="/software/trimAl/source:$PATH"

#trimmomatic
RUN apt-get install -y trimmomatic && touch trimmomatic && echo '#!/bin/bash' >> trimmomatic && echo 'java -jar /usr/share/java/trimmomatic-0.36.jar "$@"' >> trimmomatic && chmod 755 trimmomatic && mv trimmomatic /usr/bin

#PASA Pipeline
#I am not sure of this contains everything. this will have to be tested when funannotate is run
RUN cpanm install DBI && apt-get install -y mysql-server mysql-client libdbd-mysql-perl && wget https://github.com/PASApipeline/PASApipeline/releases/download/pasa-v2.4.1/PASApipeline.v2.4.1.FULL.tar.gz && tar xvfz PASApipeline.v2.4.1.FULL.tar.gz && cd PASApipeline.v2.4.1 && make
ENV PASAHOME="/software/PASApipeline.v2.4.1"

#Evidence Modeler
RUN wget https://github.com/EVidenceModeler/EVidenceModeler/archive/v1.1.1.tar.gz && tar xvfz v1.1.1.tar.gz
ENV EVM_HOME="/software/EVidenceModeler-1.1.1"

#ete3
# this installs a slightly older version, but the pyqt issue does not occur with it.
RUN apt-get install -y python-numpy python-qt4 python-lxml python-six && pip install --upgrade ete3==3.0.0b35

#Repeatmasker
RUN wget http://www.repeatmasker.org/RepeatMasker-4.1.0.tar.gz && tar xvfz RepeatMasker-4.1.0.tar.gz
ENV PATH="/software/RepeatMasker:$PATH"

#tbl2asn
RUN wget ftp://ftp.ncbi.nih.gov/toolbox/ncbi_tools/converters/by_program/tbl2asn/linux64.tbl2asn.gz && gunzip linux64.tbl2asn.gz && mv linux64.tbl2asn /usr/local/bin/tbl2asn && chmod 755 /usr/local/bin/tbl2asn

RUN rm *.gz

#export variables:
ENV PATH="/root/gm_et_linux_64:$PATH"
ENV GENEMARK_PATH="/root/gm_et_linux_64"
ENV PATH="/root/signalp-4.1:$PATH"
ENV FUNANNOTATE_DB="/root/database"
# to correct signalps path
# sed -i -e 's#$ENV{SIGNALP} = *+#$ENV{SIGNALP} = \/root\/signalp-4.1#g' /root/signalp-4.1/signalp

ENTRYPOINT ["funannotate"]
CMD ["-v"]


