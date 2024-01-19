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

rule concat

rule map_rna_to_metagenome:
    input:
        concat_rna_read             = 
    output:
    params:
    shell:
        """
        humann --threads 16 --input $(ANALYSIS)/READ_QC_RNA/$${sample}_concat.fastq --remove-temp-output --input-format fastq --output-basename $${sample} --output ./
        """
    