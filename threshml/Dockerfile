FROM rocker/r-base:3.6.3

RUN apt-get update && apt-get install -y wget

RUN install2.r ape

RUN R -e "install.packages('http://www.phytools.org/Rphylip/Rphylip_0.1-27.tar.gz',type='source',repos=NULL)"

RUN wget http://evolution.gs.washington.edu/phylip/download/threshml/threshml.zip && \
	unzip threshml.zip && \
	mv threshml/threshml.linuxi64 /usr/bin/threshml

CMD ["R"]
