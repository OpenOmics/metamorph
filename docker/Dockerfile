FROM mambaorg/micromamba

# install sudo, update system pkgs
USER root
RUN apt-get update; apt-get upgrade; apt-get install -y sudo vim debianutils bash
RUN echo '%sudo ALL=(ALL) NOPASSWD:ALL' >> /etc/sudoers
RUN echo 'routsongrm ALL=(ALL) NOPASSWD:ALL' >> /etc/sudoers

# $MAMBA_USER from inherited image
RUN usermod -a -G sudo $MAMBA_USER
RUN passwd -d $MAMBA_USER

# back to mamba user
USER $MAMBA_USER

# setup up metagenome environment
COPY assets /assets
