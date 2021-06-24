FROM reslp/mamba:0.7.0
MAINTAINER <philipp.resl@uni-graz.at>

#the following to lines are used to install the legacy package from bioconda
RUN mamba install -y -c conda-forge -c bioconda edta=1.9.6 python=3.6 tensorflow=1.14 repeatmasker=4.0.9 'h5py<3'
RUN mamba install -y rmblast=2.9.0

# the next three lines install EDTA from Github
WORKDIR /software
RUN git clone https://github.com/oushujun/EDTA.git
#RUN cd EDTA && mamba install -y $(tail -n +7 EDTA.yml | sed 's/ - //' | tr "\n" " ")
RUN cp -rf /software/EDTA /opt/conda/share/

ENV TERM=xterm
#RUN cd /opt/conda/share/RepeatMasker/ && perl ./configure -trf_prgm=`which trf` -hmmer_dir=`dirname $(which nhmmscan)` -rmblast_dir=`dirname $(which rmblastn)`
RUN cd /opt/conda/share/RepeatMasker/ && perl ./configure -trf_prgm=`which trf` -hmmerbin=`which nhmmscan` -rmblastbin=`which rmblastn`

#change shebang line so that this works in singularity (-login shells lack PATH)
RUN sed -ie 's$#!\/bin\/bash -login$#!\/bin\/bash$' /opt/conda/share/EDTA/util/run_helitron_scanner.sh

#lower maxmemory for helitronscanner to 10GB
RUN sed -ie 's$MEMGB=300$MEMGB=50$' /opt/conda/share/EDTA/util/run_helitron_scanner.sh

WORKDIR /data

CMD ["EDTA.pl"]
