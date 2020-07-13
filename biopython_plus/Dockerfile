FROM continuumio/miniconda3:4.7.12
MAINTAINER <philipp.resl@uni-graz.at>

RUN conda install -c anaconda biopython=1.77 && \
	conda install -c anaconda scipy

CMD ["python", "-v"]