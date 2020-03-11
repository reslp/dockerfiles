FROM ubuntu:18.04
MAINTAINER <philipp.resl@uni-graz.at>

RUN conda install -c bioconda maxbin2=2.2.6 && \
	conda clean -a -y

CMD ["run_MaxBin.pl"]