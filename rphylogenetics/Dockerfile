FROM r-base:4.0.3
MAINTAINER <philipp.resl@uni-graz.at>

RUN apt-get update && apt-get install -y libxml2-dev libcurl4-openssl-dev libssl-dev
#install some packages for phylogenetics, statistics and plotting
RUN install2.r -e tidyverse dplyr readr ggplot2 viridis
RUN install2.r -e plotly plyr reshape2
RUN install2.r -e ggpubr cowplot patchwork
RUN install2.r -e ape phytools wesanderson RColorBrewer
#phytools geiger brms wesanderson RColorBrewer 
RUN install2.r -e BiocManager
RUN R -e 'BiocManager::install("ggtree")'
WORKDIR /data

CMD ["R"]
