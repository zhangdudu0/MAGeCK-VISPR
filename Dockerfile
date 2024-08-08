FROM centos:centos7

# install tools
RUN yum -y install git 
RUN yum -y install wget
RUN yum -y groupinstall 'Development Tools'

# install conda
RUN wget \
    https://repo.anaconda.com/miniconda/Miniconda3-py37_4.8.3-Linux-x86_64.sh \
    && mkdir /root/.conda \
    && bash Miniconda3-py37_4.8.3-Linux-x86_64.sh -b \
    && rm -f Miniconda3-py37_4.8.3-Linux-x86_64.sh 

# set up conda
ENV PATH="/root/miniconda3/bin:${PATH}"
ARG PATH="/root/miniconda3/bin:${PATH}"
RUN conda config --add channels bioconda
RUN conda config --add channels conda-forge

# install snakemake and mageck-vispr
RUN conda install -c bioconda -c conda-forge -y mamba
#RUN mamba install -c bioconda -c conda-forge  -y snakemake
#RUN mamba install -c bioconda -c conda-forge -y mageck-vispr
RUN mamba install -c bioconda -c conda-forge  -y fastqc snakemake cutadapt jinja2 numpy scipy


# install mageck and vispr

# mageck
RUN mkdir /mageck
WORKDIR /mageck
RUN git clone https://davidliwei@bitbucket.org/liulab/mageck.git
WORKDIR /mageck/mageck
RUN python setup.py install

# vispr
RUN mkdir /vispr
WORKDIR /vispr
RUN git clone https://davidliwei@bitbucket.org/liulab/vispr.git  
WORKDIR /vispr/vispr
RUN python setup.py install

# install mageck-vispr from source
COPY . /src
WORKDIR /src
# install mageck-vispr
RUN python setup.py install

WORKDIR /
