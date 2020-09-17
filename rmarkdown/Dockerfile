FROM r-base:3.6.3

RUN apt-get update && apt-get install -y pandoc texlive-latex-recommended
RUN R -e "install.packages(c('knitr', 'rmarkdown', method='wget'))"
RUN R -e "install.packages('tinytex', method='wget')"
RUN R -e "install.packages('ggplot2', method='wget')"
#RUN apt-get install -y texlive-latex-extra
RUN tlmgr init-usertree
RUN tlmgr update --self && tlmgr install framed

WORKDIR /data

CMD ["R"]