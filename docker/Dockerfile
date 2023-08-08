FROM mambaorg/micromamba 
# $MAMBA_USER from inherited image

# install sudo, update apt pkgs
USER root
RUN apt-get update; apt-get install -y sudo vim
RUN echo '%sudo ALL=(ALL) NOPASSWD:ALL' >> /etc/sudoers
RUN usermod -a -G sudo $MAMBA_USER

# back to mamba user
USER $MAMBA_USER

# setup up metagenome environment
COPY assets /assets
ENV ENV_NAME=metagenome
#RUN micromamba create -f /assets/base.yaml
RUN micromamba create -f /assets/metawrap.yaml
#RUN micromamba run -n metawrap_env -c ursky -c bioconda metawrap-mg=1.2.1
