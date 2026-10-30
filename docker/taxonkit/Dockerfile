FROM mambaorg/micromamba:2.10.0

# Bioconda supplies TaxonKit and phyloseq; conda-forge supplies their runtimes.
RUN micromamba install --yes --name base --strict-channel-priority \
      --channel conda-forge --channel bioconda \
      python=3.12 numpy r-base r-data.table bioconductor-phyloseq 'taxonkit>=0.20' \
    && micromamba clean --all --yes

WORKDIR /opt/taxonomy-tools
COPY --chown=$MAMBA_USER:$MAMBA_USER \
    expand_taxonomy.py aggregate_phyloseq.R test_expand_taxonomy.py ./

# Activate the environment during build and fail if dependencies/tests fail.
ARG MAMBA_DOCKERFILE_ACTIVATE=1
RUN Rscript aggregate_phyloseq.R --check \
    && taxonkit version \
    && python -m unittest -v test_expand_taxonomy.py

ENV PYTHONUNBUFFERED=1
WORKDIR /data

# Preserve micromamba's activation wrapper when setting our own entrypoint.
ENTRYPOINT ["/usr/local/bin/_entrypoint.sh", "python", "/opt/taxonomy-tools/expand_taxonomy.py"]
CMD ["--help"]
