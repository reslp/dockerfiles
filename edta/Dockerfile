FROM reslp/mamba:0.7.0
MAINTAINER <philipp.resl@uni-graz.at>

RUN mamba install -y -c conda-forge -c bioconda edta=1.9.6 python=3.6 tensorflow=1.14 repeatmasker=4.0.9 'h5py<3'
RUN mamba install -y rmblast=2.9.0
ENV TERM=xterm
RUN cd /opt/conda/share/RepeatMasker/ && perl ./configure -trfbin=`which trf` -hmmerbin=`which nhmmscan` -rmblastbin=`which rmblastn`
#RUN cd /opt/conda/share/RepeatMasker/ && perl ./configure -trf_prgm=`which trf` -hmmer_dir=`which nhmmscan` -rmblast_dir=`which rmblastn`
#RUN wget https://raw.githubusercontent.com/chrishah/maker-docker/master/repeatmasker/to_include/rebuild && cp rebuild /opt/conda/share/RepeatMasker/rebuild
#ENV REPEATMASKER_DIR=/opt/conda/share/RepeatMasker
#RUN perl /opt/conda/share/RepeatMasker/rebuild

#change shebang line so that this works in singularity (-login shells lack PATH)
RUN sed -ie 's$#!\/bin\/bash -login$#!\/bin\/bash$' /opt/conda/share/EDTA/util/run_helitron_scanner.sh

#lower maxmemory for helitronscanner to 10GB
RUN sed -ie 's$MEMGB=300$MEMGB=50$' /opt/conda/share/EDTA/util/run_helitron_scanner.sh

WORKDIR /data

CMD ["EDTA.pl"]
