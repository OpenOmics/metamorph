# Base image for metaphlan-utils,
# noble base image for Ubuntu 24.04
FROM ubuntu:24.04

# Depedencies of metaphlan-utils:
# • https://github.com/biobakery/MetaPhlAn/blob/master/bioconda_recipe/meta.yaml
# @TLDR:
#    • python3                 apt-get, installs: 3.12.3
#    • python3-pip             apt-get, installs: 24.0
#    • python3-dendropy        apt-get, installs: 4.6.1
#    • python3-numpy           apt-get, installs: 1.26.4
#    • python3-biom-format     apt-get, installs: 2.1.15
#    • python3-matplotlib      apt-get, installs: 3.6.3
#    • python3-biopython       apt-get, installs: 1.83
#    • python3-pandas          apt-get, installs: 2.1.4
#    • python3-scipy           apt-get, installs: 1.11.4
#    • python3-requests        apt-get, installs: 2.31.0
#    • python3-pysam           apt-get, installs: 0.22.0
#    • r-base                  apt-get, installs: 4.3.3
#    • r-cran-optparse         apt-get, installs: 1.7.4
#    • r-cran-ape              apt-get, installs: 5.7
#    • r-cran-biocmanager      apt-get, installs: 1.30.22
#    • r-bioc-microbiome       bioconductor, installs: 1.30.0
#    • r-cran-rbiom            cran, installs: 2.2.0
#    • r-cran-compositions     cran, installs: 2.0-8

LABEL maintainer=kuhnsa@nih.gov

############### INIT ################
# Create Container filesystem specific 
# working directory and opt directories
# to avoid collisions with the host's
# filesystem, i.e. /opt and /data
RUN mkdir -p /opt2 && mkdir -p /data2
WORKDIR /opt2 

# Set time zone to US east coast 
ENV TZ=America/New_York
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime \
        && echo $TZ > /etc/timezone

############### SETUP ################
# This section installs system packages 
# required for your project. If you need 
# extra system packages add them here.
RUN apt-get update \
    && apt-get -y upgrade \
    && DEBIAN_FRONTEND=noninteractive apt-get install -y \
        build-essential \
        bzip2 \
        cmake \
        curl \
        gawk \
        git \
        gzip \
        locales \
        r-base \
        python3 \
        python3-pip \
        python3-dendropy \
        python3-numpy \
        python3-biom-format \
        python3-matplotlib \
        python3-biopython \
        python3-pandas \
        python3-scipy \
        python3-requests \
        python3-pysam \
        r-base \
        r-cran-optparse \
        r-cran-ape \
        r-cran-biocmanager \
        samtools \
        unzip \
        wget \
        zlib1g-dev \
    && apt-get clean && apt-get purge \
    && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

# Set the locale
RUN localedef -i en_US -f UTF-8 en_US.UTF-8

# Make python3 the default interpreter
# and install Python Packages 
RUN ln -sf /usr/bin/python3 /usr/bin/python

############### INSTALL ################
# Install GNU which, increases compatiability
# with which distros included with Debian/Rocky
# linux installations, this version of which
# includes extra options that the ubuntu version 
# does not have. The gnu verison is backwards 
# compatiable with the ubunutu version.
RUN mkdir -p /opt2/ \
        && wget --no-check-certificate -O /opt2/which.tar.gz https://ftp.gnu.org/gnu/which/which-2.21.tar.gz \
        && tar -xvf /opt2/which.tar.gz -C /opt2/ \
        && rm /opt2/which.tar.gz \
        && cd /opt2/which-2.21/ \
        && ./configure --prefix=/opt2/which-2.21 \
        && make \
        && make install

# Install remaining R packages which
# are not available via apt-get
RUN cd /opt2 \
    && Rscript -e 'install.packages(c("rbiom", "compositions"), Ncpus = 4, repos="http://cran.r-project.org")' \
    && Rscript -e 'BiocManager::install(c("microbiome"), Ncpus = 4, update = FALSE, ask = FALSE)'

# Install calculate diversity script
RUN cd /opt2 \
    && mkdir -p /opt2/metaphlan_utils \
    && wget --no-check-certificate -O /opt2/metaphlan_utils/calculate_diversity.R https://raw.githubusercontent.com/biobakery/MetaPhlAn/refs/heads/master/metaphlan/utils/calculate_diversity.R \
    && chmod -R a+rx /opt2/metaphlan_utils

################ POST #################
# Add Dockerfile and export important 
# environment variables
ADD Dockerfile /opt2/metaphlan-utils_v0-1-0.dockerfile
RUN chmod -R a+rX /opt2
ENV PATH="/opt2/which-2.21/bin:/opt2/metaphlan_utils:${PATH}"
WORKDIR /data2