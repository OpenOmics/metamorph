# ~~~~~~~~~~
# Metawrap metagenome assembly and analysis rules
# ~~~~~~~~~~
from os.path import join
from itertools import chain

# ~~~~~~~~~~
# Constants and paths
# ~~~~~~~~~~
workpath                            = config["project"]["workpath"]
datapath                            = config["project"]["datapath"]
rna_coasm                           = config["options"]["rnacoa"]
rna_sample_stems                    = config["rna"]


rule decompress_rna_reads:
    input:
        this_compressed_read        = start_dc
    output:
        this_uncompressed_read      = join(workpath, "{stem}.fastq"),
    params:
        rname                       = "decompress_rna_reads",
    shell:
        """
        """


rule concat_rna_reads:
    input:
    output:
    params:
        rname                       = "concat_rna_reads",
    shell:
        """
        """


rule rna_read_qc:
    input:
    output:
    params:
        rname                       = "rna_read_qc",
    shell:
        """
        mw read_qc \
        -1 $(RAWDATA.RNA)/$${sample}_1.fastq \
        -2 $(RAWDATA.RNA)/$${sample}_2.fastq \
        --skip-trimming --skip-pre-qc-report --skip-post-qc-report \
        -x hg38 \
        -t {threads} \
        -o $(ANALYSIS)/READ_QC_RNA/$${sample}
        fastqc -o $(ANALYSIS)/FASTQC_RNA/PRE -t {threads} -f fastq $(RAWDATA.RNA)/$${sample}_1.fastq $(RAWDATA.RNA)/$${sample}_2.fastq
		fastqc -o $(ANALYSIS)/FASTQC_RNA/POST -t {threads} -f fastq $(ANALYSIS)/READ_QC_RNA/$${sample}_1.fastq $(ANALYSIS)/READ_QC_RNA/$${sample}_2.fastq
        """

rule map_rna_to_metagenome:
    input:
    output:
    params:
        rname                       = "map_rna_to_metagenome",
    shell:
        """
        humann --threads 16 --input $(ANALYSIS)/READ_QC_RNA/$${sample}_concat.fastq --remove-temp-output --input-format fastq --output-basename $${sample} --output ./
        """
    