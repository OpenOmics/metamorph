FROM mambaorg/micromamba 
# $MAMBA_USER from inherited image

# install sudo, update apt pkgs
USER root
RUN apt-get update; apt-get install -y sudo vim
RUN echo '%sudo ALL=(ALL) NOPASSWD:ALL' >> /etc/sudoers
RUN usermod -a -G sudo $MAMBA_USER

# back to mamba user
USER $MAMBA_USER

# update chanels
RUN micromamba update -y -n base -c conda-forge; micromamba update -y -n base -c bioconda; micromamba update -y -n base -c defaults

# setup up metagenome environment
COPY env.yaml /env.yaml
RUN micromamba create -f /env.yaml
ENV ENV_NAME=metagenome