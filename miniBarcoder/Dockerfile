FROM continuumio/miniconda:4.7.10

WORKDIR /software

RUN conda install mamba -c conda-forge

RUN conda config --add channels bioconda
RUN conda install python=2.7 mafft racon graphmap=0.5.2 blast seqtk git fasta3
RUN conda install -c anaconda biopython 

#ENV PATH /opt/conda/envs/mbarcode/bin:$PATH
#ENV CONDA_DEFAULT_ENV mbarcode
#RUN /bin/bash -c "source activate mbarcode"

RUN git clone https://github.com/asrivathsan/miniBarcoder.git
RUN cd miniBarcoder && \
	python setup.py install

WORKDIR /data

CMD ["/bin/bash"]