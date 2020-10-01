FROM continuumio/miniconda
SHELL ["/bin/bash", "-c"]

RUN . /opt/conda/etc/profile.d/conda.sh
#RUN source $HOME/.bashrc
RUN conda config --set notify_outdated_conda false
RUN conda install -y -c bioconda spades=3.13.0

WORKDIR /data

CMD ["spades.py", "-h"]

