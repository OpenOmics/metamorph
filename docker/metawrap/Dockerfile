FROM condaforge/miniforge3:latest
ENV TZ="America/Costa_Rica"
RUN DEBIAN_FRONTEND=noninteractive
RUN ln -sf /usr/share/zoneinfo/America/Costa_Rica /etc/localtime
RUN apt-get update; apt-get install -y -qq curl build-essential vim dos2unix bash python2.7 cmake zlib1g-dev
RUN mamba create -y -n metawrap-env
RUN conda config --add channels defaults; conda config --add channels conda-forge; \
    conda config --add channels bioconda; conda config --add channels ursky; \
    conda config --add channels biobakery
RUN mamba install -y --only-deps -c ursky -n metawrap-env metawrap-mg==1.3.2
RUN mamba install -y -n metawrap-env -c bioconda spades=3.15.5
RUN cd /home; git clone https://github.com/rroutsong/metaWRAP.git; chmod -R 777 metaWRAP
RUN mamba create -n checkm python=3.8
RUN mamba install -y -n checkm -c bioconda numpy matplotlib pysam hmmer prodigal pplacer fastani samtools bowtie2 trnascan-se gunc prokka gtdbtk=2.3.2
RUN mamba run -n checkm pip3 install checkm-genome drep
RUN mamba create -n bb3 python=3.7 
RUN mamba install -y -n bb3 -c biobakery humann
RUN mkdir /install; cd /install; wget https://carlowood.github.io/which/which-2.21.tar.gz; tar xvf which-2.21.tar.gz
RUN cd /install/which-2.21; ./configure; make && make install
RUN rm /usr/bin/which; ln -s /usr/local/bin/which /usr/bin/which
RUN cd /home; git clone https://github.com/dutilh/CAT.git; chmod +x CAT_pack/CAT_pack; ln -s /home/CAT/CAT_pack/CAT_pack /usr/bin/CAT; chmod +x /usr/bin/CAT 
RUN cd /install; wget https://github.com/hyattpd/Prodigal/archive/refs/tags/v2.6.3.tar.gz; tar xvf v2.6.3.tar.gz
RUN cd /install/Prodigal-2.6.3; make install INSTALLDIR=/usr/bin
RUN cd /install; wget https://github.com/bbuchfink/diamond/archive/refs/tags/v2.1.9.tar.gz; tar xvf v2.1.9.tar.gz
RUN cd /install/diamond-2.1.9; mkdir bin; cd bin; cmake ..; make -j4; ln -s /install/diamond-2.1.9/bin/diamond /usr/bin/diamond
RUN cd /home; wget https://sourceforge.net/projects/bbmap/files/BBMap_39.06.tar.gz; tar xvf BBMap_39.06.tar.gz 
ENV PATH="/home/metaWRAP/bin:/opt/conda/envs/metawrap-env/bin:/home/bbmap:$PATH"
RUN cd /install; mkdir mash; cd mash; wget https://github.com/marbl/Mash/releases/download/v2.3/mash-Linux64-v2.3.tar; tar xvf mash-Linux64-v2.3.tar
RUN ln -s /install/mash/mash-Linux64-v2.3/mash /usr/bin/mash
RUN cd /install; mkdir mummer; cd mummer; wget https://github.com/mummer4/mummer/releases/download/v4.0.0rc1/mummer-4.0.0rc1.tar.gz; tar xvf mummer-4.0.0rc1.tar.gz
RUN cd /install/mummer/mummer-4.0.0rc1; ./configure --prefix=/usr; make && make install
COPY docker/metawrap/config-metawrap /home/metaWRAP/bin/config-metawrap
COPY docker/metawrap/Dockerfile /Dockerfile
RUN md5sum Dockerfile > /Dockerfile.md5
ADD docker/metawrap/mw /home/metaWRAP/bin/mw
RUN dos2unix /home/metaWRAP/bin/config-metawrap
RUN dos2unix /home/metaWRAP/bin/mw
RUN chmod a+rx /home/metaWRAP/bin/config-metawrap /home/metaWRAP/bin/mw
ENV BASH_ENV="/etc/bash.bashrc"
RUN echo ". /opt/conda/etc/profile.d/conda.sh && conda activate metawrap-env" >> /etc/bash.bashrc
RUN echo 'export PATH="/home/metaWRAP/bin:/opt/conda/envs/metawrap-env/bin:$PATH"' >> /etc/bash.bashrc
ENTRYPOINT ["/bin/bash"]
