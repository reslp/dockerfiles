FROM continuumio/miniconda3:4.7.12
MAINTAINER <philipp.resl@uni-graz.at>

RUN conda install -c anaconda biopython=1.77 && \
	conda install -c anaconda scipy && \
	conda install -c anaconda matplotlib=3.3.1 && \
	conda install -c anaconda pandas=1.1.1 && \
	conda install -c anaconda seaborn=0.10.1

CMD ["python", "-v"]
