FROM continuumio/miniconda
SHELL ["/bin/bash", "-c"]

RUN . /opt/conda/etc/profile.d/conda.sh
RUN conda config --set notify_outdated_conda false
RUN conda install -y -c bioconda trimmomatic=0.38

CMD ["trimmomatic"]

