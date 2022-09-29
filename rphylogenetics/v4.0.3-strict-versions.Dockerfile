FROM r-base:4.0.3
MAINTAINER <philipp.resl@uni-graz.at>

RUN apt-get update && apt-get install -y libxml2-dev libcurl4-openssl-dev libssl-dev build-essential
#install some packages for phylogenetics, statistics and plotting
RUN apt install -y libtiff-dev libfontconfig1-dev libharfbuzz-dev libfribidi-dev
RUN apt clean
RUN install2.r -e devtools
#RUN install2.r -e tidyverse dplyr readr ggplot2 viridis
RUN R -e 'require(devtools); install_version("tidyverse", version = "1.3.1", repos = "http://cran.us.r-project.org")'
RUN R -e 'require(devtools); install_version("viridis", version = "0.6.1", repos = "http://cran.us.r-project.org")'
#RUN install2.r -e plotly plyr reshape2
RUN R -e 'require(devtools); install_version("plotly", version = "4.9.3", repos = "http://cran.us.r-project.org")'
RUN R -e 'require(devtools); install_version("plyr", version = "1.8.6", repos = "http://cran.us.r-project.org")'
RUN R -e 'require(devtools); install_version("reshape2", version = "1.4.4", repos = "http://cran.us.r-project.org")'
#RUN install2.r -e ggpubr cowplot patchwork
RUN R -e 'require(devtools); install_version("ggpubr", version = "0.4.0", repos = "http://cran.us.r-project.org")'
RUN R -e 'require(devtools); install_version("cowplot", version = "1.1.1", repos = "http://cran.us.r-project.org")'
RUN R -e 'require(devtools); install_version("patchwork", version = "1.1.1", repos = "http://cran.us.r-project.org")'
#RUN install2.r -e ape phytools wesanderson RColorBrewer
#phytools geiger brms wesanderson RColorBrewer 
RUN R -e 'require(devtools); install_version("ape", version = "5.5", repos = "http://cran.us.r-project.org")'
RUN R -e 'require(devtools); install_version("wesanderson", version = "0.3.6", repos = "http://cran.us.r-project.org")'
RUN R -e 'require(devtools); install_version("RColorBrewer", version = "1.1-2", repos = "http://cran.us.r-project.org")'
#RUN install2.r -e ggrepel
RUN R -e 'require(devtools); install_version("ggrepel", version = "0.9.1", repos = "http://cran.us.r-project.org")'
#RUN install2.r -e BiocManager
#RUN R -e 'require(devtools); install_version("phytools", version = "0.7-80", repos = "http://cran.us.r-project.org")'
RUN R -e 'require(devtools); install_version("phangorn", version = "2.7.0", repos = "http://cran.us.r-project.org")'
RUN R -e 'require(devtools); install_version("phytools", version = "0.7-70", repos = "http://cran.us.r-project.org")'
RUN R -e 'require(devtools); install_version("pheatmap", version = "1.0.12", repos = "http://cran.us.r-project.org")'

RUN install2.r -e ggfun BiocManager
Run R -e 'devtools::install_version("rvcheck", version = "0.1.8", repos = "http://cran.us.r-project.org")'
RUN R -e 'BiocManager::install("ggtree")'

RUN R -e 'require(devtools); install_version("ggnewscale", version = "0.4.7", repos = "http://cran.us.r-project.org")'
WORKDIR /data

CMD ["R"]
