FROM reslp/mamba:0.7.0
MAINTAINER <philipp.resl@uni-graz.at>

RUN mamba install psutil scipy matplotlib

WORKDIR /software

#not sure if this command works always:
RUN wget --header 'Host: files.dtu.dk' --user-agent 'Mozilla/5.0 (X11; Ubuntu; Linux x86_64; rv:84.0) Gecko/20100101 Firefox/84.0' --header 'Accept: text/html,application/xhtml+xml,application/xml;q=0.9,image/webp,*/*;q=0.8' --header 'Accept-Language: en-US,en;q=0.5' --referer 'https://files.dtu.dk/userportal/' --header 'Cookie: JSESSIONID=B020F421C63E86507B8B3EB807E65819; BIGipServer~BIT-DSG~pool-files.dtu.dk=549850890.20480.0000' --header 'Upgrade-Insecure-Requests: 1' 'https://files.dtu.dk/fss/public/link/public/stream/read/CUPP_v1.0.14.tar.gz?linkToken=hLin6ni4p-SWuKfp&itemName=CUPP_program' --output-document 'CUPP_v1.0.14.tar.gz'

RUN tar xvfz CUPP_v1.0.14.tar.gz
RUN chmod +x CUPP_v1.0.14/CUPPclustering.py && chmod +x CUPP_v1.0.14/CUPPprediction.py
RUN sed -i 's$#!/usr/bin/python$#!/usr/bin/env python$' CUPP_v1.0.14/CUPPclustering.py && sed -i 's$#!/usr/bin/python$#!/usr/bin/env python$' CUPP_v1.0.14/CUPPprediction.py
ENV PATH="$PATH:/software/CUPP_v1.0.14"
