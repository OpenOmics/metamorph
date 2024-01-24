FROM condaforge/miniforge3:latest
RUN apt-get update; apt-get install -y -qq curl vim dos2unix bash python2.7 gnu-which
RUN mamba create -y -n metawrap-env
RUN conda config --add channels defaults; conda config --add channels conda-forge; \
    conda config --add channels bioconda;conda config --add channels ursky
RUN mamba install -y --only-deps -c ursky -n metawrap-env metawrap-mg==1.3.2
RUN cd /home; git clone https://github.com/bxlab/metaWRAP.git; chmod -R 777 metaWRAP
RUN mamba run -n metawrap-env pip3 install drep
ENV PATH="/home/metaWRAP/bin:$PATH"
COPY docker/metawrap/config-metawrap /home/metaWRAP/bin/config-metawrap
COPY docker/metawrap/Dockerfile /Dockerfile
RUN md5sum Dockerfile > /Dockerfile.md5
ADD docker/metawrap/mw /home/metaWRAP/bin/mw
RUN dos2unix /home/metaWRAP/bin/config-metawrap
RUN dos2unix /home/metaWRAP/bin/mw
RUN chmod a+rx /home/metaWRAP/bin/config-metawrap /home/metaWRAP/bin/mw
RUN echo ". /opt/conda/etc/profile.d/conda.sh && conda activate metawrap-env" >> /etc/bash.bashrc
RUN echo 'export PATH="/home/metaWRAP/bin:$PATH"' >> /etc/bash.bashrc
ENTRYPOINT ["/bin/bash"]