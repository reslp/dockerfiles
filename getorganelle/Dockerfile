FROM continuumio/miniconda3:4.8.2

RUN conda install -c bioconda getorganelle=1.7.1
RUN get_organelle_config.py -a animal_mt,embplant_mt,embplant_nr,embplant_pt,fungus_mt,other_pt

CMD [get_organelle_from_reads.py]

