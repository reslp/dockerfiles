FROM julia:1.6.3-buster
MAINTAINER <philipp.resl@uni-graz.at>

# first add some BioJulia packages:
RUN julia -e 'using Pkg; Pkg.add(["BioSequences", "FASTX", "BioCore", "BioAlignments", "BioSymbols", "BioServices", "GFF3", "BED", "BioFetch", "GenomicAnnotations"])'
# second install tools for plotting:
RUN julia -e 'using Pkg; Pkg.add(["Plots", "PyPlot", "Makie", "Gadfly"])'
# install some reporting tools:
RUN julia -e 'using Pkg; Pkg.add(["Weave"])'
# install some data tools:
RUN julia -e 'using Pkg; Pkg.add(["DataFrames"])'

WORKDIR /data

RUN ["julia"]
