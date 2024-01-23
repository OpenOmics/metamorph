FROM condaforge/miniforge3:latest
RUN apt-get update; apt-get install -y -qq curl vim dos2unix bash
RUN mamba create -y -n metawrap-env
RUN conda config --add channels defaults; conda config --add channels conda-forge; \
    conda config --add channels bioconda;conda config --add channels ursky
RUN mamba install -y --only-deps -c ursky -n metawrap-env metawrap-mg==1.3.2
RUN cd /home; git clone https://github.com/bxlab/metaWRAP.git; chmod -R 777 metaWRAP
ENV PATH="/home/metaWRAP/bin:$PATH"
COPY docker/metawrap/config-metawrap /home/metaWRAP/bin/config-metawrap
ADD docker/metawrap/mw /home/metaWRAP/bin/mw
RUN dos2unix /home/metaWRAP/bin/config-metawrap
RUN chmod a+rx /home/metaWRAP/bin/config-metawrap /home/metaWRAP/bin/mw
RUN mamba run -n metawrap-env pip3 install drep
ENV BASH_ENV=/etc/skel/.bashrc
ENTRYPOINT ["/bin/bash"]
