FROM r-base:4.0.3


#RUN useradd docker \
#	&& mkdir /home/docker \
#	&& chown docker:docker /home/docker \
#	&& addgroup docker staff

#RUN echo "deb http://http.debian.net/debian sid main" > /etc/apt/sources.list.d/debian-unstable.list \
#        && echo 'APT::Default-Release "testing";' > /etc/apt/apt.conf.d/default

RUN apt-get update && apt-get install -y pandoc littler r-cran-littler
#USER docker
#RUN ln -s /usr/lib/R/site-library/litter/examples/install.r /usr/local/bin/install.r \
#	&& ln -s /usr/lib/R/site-library/litter/examples/install2.r /usr/local/bin/install2.r \
#	&& ln -s /usr/lib/R/site-library/litter/examples/installGithub.r /usr/local/bin/installGithub.r \
#	&& ln -s /usr/lib/R/site-library/litter/examples/testInstalled.r /usr/local/bin/testInstalled.r

#install some packages relevant for rmarkdown
RUN apt-get install -y libxml2-dev libcurl4-openssl-dev libssl-dev
RUN install2.r -e tinytex tidyverse kableExtra dplyr tidyverse readr ggplot2 viridis
RUN R -e "tinytex::install_tinytex()"	
RUN install2.r -e plotly plyr reshape2
RUN install2.r -e ggpubr cowplot patchwork
RUN install2.r -e flexdashboard
RUN install2.r -e formattable
#RUN apt-get install -y texlive-latex-extra
#RUN tlmgr init-usertree
#RUN tlmgr update --self && tlmgr install framed
#RUN apt-get install -y texlive-xetex
#RUN R -e "tinytex::install_tinytex()"
#RUN chmod a+w /usr/local/texlive
#RUN chmod a+w /usr/local/texlive/tlpkg

WORKDIR /data

CMD ["R"]
