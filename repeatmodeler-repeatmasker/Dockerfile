FROM ubuntu:20.04

ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get update && apt-get install -y wget build-essential zlib1g-dev libgomp1 \
    perl \
    python3-h5py \
    libfile-which-perl \
    libtext-soundex-perl \
    libjson-perl liburi-perl libwww-perl

WORKDIR /software

# RECON
RUN wget http://www.repeatmasker.org/RepeatModeler/RECON-1.08.tar.gz && \
	tar xvfz RECON-1.08.tar.gz && \
	cd RECON-1.08/src && \
	make && make install && \
	sed -i -e 's#\$path = \"\";#\$path = \"/software/RECON-1.08/bin\";#' /software/RECON-1.08/scripts/recon.pl
# RepeatScout
RUN wget http://www.repeatmasker.org/RepeatScout-1.0.6.tar.gz && \
	tar xvfz RepeatScout-1.0.6.tar.gz && \
	cd RepeatScout-1.0.6 && \
	make && \
	cp RepeatScout build_lmer_table /usr/local/bin
# TRF
RUN wget https://github.com/Benson-Genomics-Lab/TRF/releases/download/v4.09.1/trf409.linux64 && \
	chmod +x trf409.linux64 && \
	mv trf409.linux64 /usr/local/bin/trf
# rmblast
RUN wget http://www.repeatmasker.org/rmblast-2.9.0+-p2-x64-linux.tar.gz && \
	tar xvfz rmblast-2.9.0+-p2-x64-linux.tar.gz

RUN wget http://www.repeatmasker.org/RepeatModeler/RepeatModeler-2.0.3.tar.gz && \
 	tar xvfz RepeatModeler-2.0.3.tar.gz && \
 	cd RepeatModeler-2.0.3 && \
 	perl -i -0pe 's/^#\!.*/#\!\/usr\/bin\/env perl/g' \
 	configure \
 	BuildDatabase \
 	Refiner \
 	RepeatClassifier \
 	RepeatModeler \
 	TRFMask \
 	util/dfamConsensusTool.pl \
 	util/renameIds.pl \
 	util/viewMSA.pl && \
 	sed -i -e 's#/usr/local/RECON#/software/RECON-1.08/bin#g' RepModelConfig.pm && \
 	sed -i -e 's#/usr/local/RepeatScout-1.0.6#/software/RepeatScout-1.0.6#g' RepModelConfig.pm && \
 	sed -i -e 's#/usr/local/rmblast/bin#/software/rmblast-2.9.0-p2/bin#g' RepModelConfig.pm && \
 	sed -i -e 's#/usr/local/RepeatMasker#/software/RepeatMasker#g' RepModelConfig.pm && \
 	sed -i -e 's#/usr/local/bin/trf409.linux64#/usr/local/bin/trf#g' RepModelConfig.pm
 	#cp RepModelConfig.pm.tmpl RepModelConfig.pm

#Repeatmasker
# fixing the perl paths is from: https://github.com/chrishah/maker-docker/blob/master/repeatmasker/Dockerfile
RUN wget http://www.repeatmasker.org/RepeatMasker/RepeatMasker-4.1.3-p1.tar.gz && \
	tar xvfz RepeatMasker-4.1.3-p1.tar.gz 
RUN cd RepeatMasker && ./configure -trf_prgm $(which trf) -default_search_engine ncbi -rmblast_dir /software/rmblast-2.9.0-p2/bin

#	&& \
#	sed -i -e 's#TRF_PRGM = ""#TRF_PRGM = \"/usr/local/bin/trf\"#g' RepeatMasker/RepeatMaskerConfig.tmpl && \
#	sed -i -e 's#DEFAULT_SEARCH_ENGINE = \"crossmatch\"#DEFAULT_SEARCH_ENGINE = \"ncbi\"#g' RepeatMasker/RepeatMaskerConfig.tmpl && \
#	sed -i -e 's#RMBLAST_DIR   = \"/usr/local/rmblast\"#RMBLAST_DIR   = \"/software/rmblast-2.9.0-p2/bin\"#g' RepeatMasker/RepeatMaskerConfig.tmpl && \
#	sed -i -e 's#REPEATMASKER_DIR          = \"\$FindBin::RealBin\"#REPEATMASKER_DIR          = \"/software/RepeatMasker\"#g' RepeatMasker/RepeatMaskerConfig.tmpl && \
#	cp RepeatMasker/RepeatMaskerConfig.tmpl RepeatMasker/RepeatMaskerConfig.pm && \
#	perl -i -0pe 's/^#\!.*perl.*/#\!\/usr\/bin\/env perl/g' \
#	RepeatMasker/RepeatMasker \
#	RepeatMasker/DateRepeats \
#	RepeatMasker/ProcessRepeats \
#	RepeatMasker/RepeatProteinMask \
#	RepeatMasker/DupMasker \
#	RepeatMasker/util/queryRepeatDatabase.pl \
#	RepeatMasker/util/queryTaxonomyDatabase.pl \
#	RepeatMasker/util/rmOutToGFF3.pl \
#	RepeatMasker/util/rmToUCSCTables.pl
#

RUN apt-get install -y cpanminus
RUN cpanm 

ENV RMBLAST_DIR=/software/rmblast-2.9.0-p2/bin
ENV RECON_DIR=/software/RECON-1.08/bin
ENV PATH="/software/RepeatModeler-2.0.3:$PATH"

WORKDIR /data

CMD ["/bin/bash"]
