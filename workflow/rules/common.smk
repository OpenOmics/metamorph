# ~~~~~~~~~~
# Common utility snakemake utility
# ~~~~~~~~~~
rule concatenate_reads:
    input:
        to_cat=expand(join(workpath, "{name}.R2_readqc.fastq.gz"), name=samples),
    output:
        cohort_reads=join(workpath, "cohort.reads.fastq.gz"),
    shell:
        """
            cat {input.to_cat} > {output.cohort_reads}
        """
