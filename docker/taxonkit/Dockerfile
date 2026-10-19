# Base image for taxonkit,
# noble base image for Ubuntu 24.04
FROM ubuntu:24.04

# TaxonKit: https://github.com/shenwei356/taxonkit
# A practical and efficient NCBI taxonomy toolkit, used here to expand
# a taxID into its full lineage (domain/phylum/class/order/family/genus/
# species). It ships as a single static Go binary with no runtime deps.
# @TLDR:
#    • taxonkit    github release, installs: v0.20.0
ENV TAXONKIT_VERSION=0.20.0

LABEL maintainer=rroutsong@nih.gov

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
        ca-certificates \
        curl \
        gzip \
        locales \
        tar \
    && apt-get clean && apt-get purge \
    && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

# Set the locale
RUN localedef -i en_US -f UTF-8 en_US.UTF-8

############### INSTALL ################
# Install taxonkit from its upstream
# GitHub release (statically-linked
# Go binary, no build step needed)
RUN mkdir -p /opt2/taxonkit \
    && curl -fsSL -o /opt2/taxonkit.tar.gz \
        "https://github.com/shenwei356/taxonkit/releases/download/v${TAXONKIT_VERSION}/taxonkit_linux_amd64.tar.gz" \
    && tar -zxvf /opt2/taxonkit.tar.gz -C /opt2/taxonkit \
    && rm /opt2/taxonkit.tar.gz \
    && chmod a+rx /opt2/taxonkit/taxonkit

# NOTE: taxonkit needs an NCBI-taxonomy-style taxdump (names.dmp,
# nodes.dmp, delnodes.dmp, merged.dmp) to do anything useful. That
# database is NOT baked into this image -- like this pipeline's other
# reference databases (see config/resources.json), it's bind-mounted
# in at runtime and pointed to via TAXONKIT_DB rather than the
# $HOME/.taxonkit default, so the image stays small and the database
# stays independently updatable.
ENV TAXONKIT_DB=/data2/NCBI_TAX_DB

################ POST #################
# Add Dockerfile and export important
# environment variables
ADD Dockerfile /opt2/taxonkit_v${TAXONKIT_VERSION}.dockerfile
RUN chmod -R a+rX /opt2
ENV PATH="/opt2/taxonkit:${PATH}"
WORKDIR /data2
