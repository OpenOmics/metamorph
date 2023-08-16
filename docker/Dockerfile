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
RUN micromamba create -f /assets/base.yaml

# setup metawrap environment
RUN micromamba create -n metawrap_env
RUN micromamba install -y --only-deps -n metawrap_env -c defaults -c conda-forge -c bioconda -c ursky metawrap-mg
RUN micromamba install -y -n metawrap_env -c ursky metawrap-mg=1.3.2
RUN micromamba install -y -n metawrap_env -c defaults -c conda-forge -c bioconda -c ursky biopython blas=2.5 blast=2.6.0 bmtagger bowtie2 bwa checkm-genome fastqc kraken=1.1 kraken2 krona=2.7 matplotlib maxbin2 megahit metabat2 pandas prokka quast r-ggplot2 r-recommended salmon samtools=1.9 seaborn spades trim-galore
RUN micromamba install -y -n metawrap_env -c ursky concoct=1.0.0 pplacer
