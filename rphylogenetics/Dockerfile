FROM reslp/rmarkdown:4.0.3
MAINTAINER <philipp.resl@uni-graz.at>

#install some packages for phylogenetics, statistics and plotting
RUN install2.r -e ape phytools wesanderson RColorBrewer
#phytools geiger brms wesanderson RColorBrewer 
RUN install2.r -e BiocManager
RUN R -e 'BiocManager::install("ggtree")'
RUN install2.r -e cowplot
WORKDIR /data

CMD ["R"]
