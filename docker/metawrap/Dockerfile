FROM rroutsong/metamorph_base:latest
# setup metawrap environment
RUN micromamba create -n metawrap_env
ENV ENV_NAME=metawrap_env
RUN micromamba install -y --only-deps -n metawrap_env -c defaults -c conda-forge -c bioconda -c ursky metawrap-mg
RUN micromamba install -y -n metawrap_env -c ursky metawrap-mg=1.3.2
RUN micromamba install -y -n metawrap_env -c defaults -c conda-forge -c bioconda -c ursky biopython blas=2.5 blast=2.6.0 bmtagger bowtie2 multiqc bwa checkm-genome fastqc kraken=1.1 kraken2 krona=2.7 matplotlib maxbin2 megahit metabat2 pandas prokka quast r-ggplot2 r-recommended salmon samtools=1.9 seaborn spades trim-galore
RUN micromamba install -y -n metawrap_env -c ursky concoct=1.0.0 pplacer