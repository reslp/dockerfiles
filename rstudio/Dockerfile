FROM rocker/rstudio:4.0.0-ubuntu18.04

RUN apt-get install -y libxml2-dev
#RUN R -e "install.packages(c('knitr', 'kableExtra', 'tinytex'), dependencies=TRUE, repos='http://cran.rstudio.com/')"
RUN install2.r knitr \
	kableExtra \
	tinytex

CMD ["/init"]
