FROM continuumio/miniconda3:4.8.2
MAINTAINER philipp.resl@uni-graz.at>

RUN conda config --append channels bioconda && \
	conda config --append channels conda-forge && \
	conda config --append channels anaconda && \
	conda config --set notify_outdated_conda false && \
	conda install -c anaconda requests beautifulsoup4 && \
	conda clean -a -y
RUN apt-get install -y wget && \
	wget https://raw.githubusercontent.com/reslp/genomics/master/scrape_cazy.py && \
	mv scrape_cazy.py /usr/local/bin && \
	chmod +x /usr/local/bin/scrape_cazy.py
WORKDIR /data
CMD ["scrape_cazy.py","-h"]