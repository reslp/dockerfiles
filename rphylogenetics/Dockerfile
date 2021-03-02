FROM reslp/rmarkdown:4.0.3
MAINTAINER <philipp.resl@uni-graz.at>

#install some packages for phylogenetics, statistics and plotting
RUN install2.r -e ape phytools wesanderson RColorBrewer
#phytools geiger brms wesanderson RColorBrewer 

WORKDIR /data

CMD ["R"]
