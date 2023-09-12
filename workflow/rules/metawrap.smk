# ~~~~~~~~~~
# Metawrap metagenome assembly and analysis rules
# ~~~~~~~~~~
from os.path import join
from itertools import chain


# ~~~~~~~~~~
# Constants and paths
# ~~~~~~~~~~ 
default_threads            = cluster["__default__"]['threads']
default_memory             = cluster["__default__"]['mem']
top_readqc_dir             = join(workpath, "metawrap_readqc")
top_assembly_dir           = join(workpath, "metawrap_assembly")
top_tax_dir                = join(workpath, "metawrap_taxonomy")
top_binning_dir            = join(workpath, "metawrap_binning")
metawrap_container         = config["images"]["metawrap"]


rule metawrap_read_qc:
    """
        Metawrap read quality control rule for producing high quality reads for assembly.

        Quality-control step accomplishes three tasks:
            - Trims adapter sequences from input read pairs (.fastq.gz)
            - Quality filters from input read pairs
            - Removes human contamination from input read pairs
        @Container requirements
            - [MetaWrap](https://github.com/bxlab/metaWRAP) orchestrates execution of 
              these tools in a psuedo-pipeline:
                - [FastQC](https://github.com/s-andrews/FastQC)
                - [bmtagger](ftp://ftp.ncbi.nlm.nih.gov/pub/agarwala/bmtagger/)
                - [Trim-Galore](https://github.com/FelixKrueger/TrimGalore)
                - [Cutadapt](https://github.com/marcelm/cutadapt)
        
        @Environment specifications
            - minimum 16gb memory
                - will load bmtagger index in memory (8gb)
        @Input:
            - Raw fastq.gz reads (R1 & R2 per-sample) from the instrument
            
        @Outputs:
            - Trimmed, quality filtered reads (R1 & R2)
            - FastQC html report and zip file on trimmed data
    """
    input:
        R1                  = join(workpath, "{name}.R1.fastq.gz"),
        R2                  = join(workpath, "{name}.R2.fastq.gz"),
    output:
        R1_bmtagger_report  = join(top_readqc_dir, "{name}", "{name}.R1_bmtagger_report.html"),
        R2_bmtagger_report  = join(top_readqc_dir, "{name}", "{name}.R2_bmtagger_report.html"),
        R1_fastqc_report    = join(top_readqc_dir, "{name}", "{name}.R1_fastqc_report.html"),
        R2_fastqc_report    = join(top_readqc_dir, "{name}", "{name}.R2_fastqc_report.html"),
        R1_qc_reads         = join(workpath, "{name}.R1_readqc.fastq"),
        R2_qc_reads         = join(workpath, "{name}.R2_readqc.fastq"),
    params:
        readqc_dir          = join(top_readqc_dir, "{name}"),
        R1_mw_named         = join(top_readqc_dir, "{name}", "{name}_1.fastq"),
        R2_mw_named         = join(top_readqc_dir, "{name}", "{name}_2.fastq"),
        
    singularity: metawrap_container,
    threads: int(cluster["metawrap_genome_assembly"].get('threads', default_threads)),
    shell: 
        """
            # Setting up inputs to metawrap,
            # needs uncompressed files as input,
            # input files also need to end with
            # the following ext: _[1,2].fastq
            gunzip -c {input.R1} > {params.R1_mw_named}
            gunzip -c {input.R2} > {params.R2_mw_named}
            
            # Running metawraps read_qc module
            metawrap read_qc -1 {params.R1_mw_named} -2 {params.R2_mw_named} -t {threads} -o {params.readqc_dir}
            
            # Rename R1 output files from metawrap
            mv {params.readqc_dir}/final_pure_reads_1.fastq {params.R1_qc_reads}
            rm -f {params.R1_mw_named}
            cp {params.readqc_dir}/post-QC_report/final_pure_reads_1_fastqc.html {output.R1_bmtagger_report}
            cp {params.readqc_dir}/pre-QC_report/{params.R1_mw_named}_fastqc.html {output.R1_fastqc_report}

            # Rename R2 output files from metawrap
            mv {params.readqc_dir}/final_pure_reads_2.fastq {params.R2_qc_reads}
            rm -f {params.R2_mw_named}
            cp {params.readqc_dir}/post-QC_report/final_pure_reads_2_fastqc.html {output.R2_bmtagger_report}
            cp {params.readqc_dir}/pre-QC_report/{params.R2_mw_named}_fastqc.html {output.R2_fastqc_report}
        """


rule metawrap_genome_assembly:
    """
        TODO
        @Input:
            Clean trimmed fastq.gz reads (R1 & R2 per sample)
        @Output:
            Megahit assembled contigs and reports
            Metaspades assembled contigs and reports
            Ensemble assembled contigs and reports
    """
    input:
        R1                          = join(workpath, "{name}.R1_readqc.fastq"),
        R2                          = join(workpath, "{name}.R2_readqc.fastq"),
    output:
        # megahit outputs
        megahit_assembly            = join(top_assembly_dir, "{name}", "megahit", "final.contigs.fa"),
        megahit_longcontigs         = join(top_assembly_dir, "{name}", "megahit", "long.contigs.fa"),
        megahit_log                 = join(top_assembly_dir, "{name}", "megahit", "log"),
        # metaspades outsputs
        metaspades_assembly         = join(top_assembly_dir, "{name}", "metaspades", "contigs.fasta"),
        metaspades_graph            = join(top_assembly_dir, "{name}", "metaspades", "assembly_graph.fastg"),
        metaspades_longscaffolds    = join(top_assembly_dir, "{name}", "metaspades", "long_scaffolds.fasta"),
        metaspades_scaffolds        = join(top_assembly_dir, "{name}", "metaspades", "scaffolds.fasta"),
        metaspades_cor_readsr1      = join(top_assembly_dir, "{name}", "metaspades", "{name}_1.fastq.00.0_0.cor.fastq.gz"),
        metaspades_cor_readsr2      = join(top_assembly_dir, "{name}", "metaspades", "corrected", "{name}_2.fastq.00.0_0.cor.fastq.gz"),
        # ensemble outputs
        final_assembly              = join(top_assembly_dir, "{name}", "final_assembly.fasta"),
        final_assembly_report       = join(top_assembly_dir, "{name}", "assembly_report.html"),
    singularity: metawrap_container,
    params:
        assembly_dir                = join(top_assembly_dir, "{name}"),
        memlimit                    = cluster["metawrap_genome_assembly"].get('mem', default_memory),
        contig_min_len              = "1000",
        assembly_R1                 = join(workpath, "{name}_1.fastq"),
        assembly_R2                 = join(workpath, "{name}_2.fastq"),
    threads: int(cluster["metawrap_genome_assembly"].get('threads', default_threads)),
    shell:
        """
            # link to the file names metawrap expects
            ln -s {input.R1} {params.assembly_R1}
            ln -s {input.R2} {params.assembly_R2}
            # run genome assembler
            metawrap assembly --megahit --metaspades -m {params.memlimit} -t {threads} -l {params.contig_min_len} -1 {params.assembly_R1} -2 {params.assembly_R2} -o {params.assembly_dir}
        """


rule metawrap_tax_classification:
    """
        TODO: docstring
    """
    input:
        final_assembly              = join(top_assembly_dir, "{name}", "final_assembly.fasta"),
        R1                          = expand(join(workpath, "{name}.R1_readqc.fastq.gz")),
        R2                          = expand(join(workpath, "{name}.R2_readqc.fastq.gz")),
        reads                       = list(chain(*zip(R1, R2)))
    params:
        tax_subsample               = str(int(1e6)),
        tax_dir                     = join(top_tax_dir, "{name}"),
    singularity: metawrap_container
    threads: cluster["metawrap_tax_classification"].get("threads", default_threads),
    shell:
        """
            metawrap kraken2 -t {threads} -s {params.tax_subsample} -o {params.tax_dir} {input.final_assembly} {input.reads}
        """

"""
rule metawrap_tax_viz:
    Intended for krona viz, but not sure if necessary
"""

rule metawrap_binning_rename:
    """
        TODO: docstring
    """
    input:
        R1_from_qc                          = join(workpath, "{name}.R{pair}_readqc.fastq"),
        R2_from_qc                          = join(workpath, "{name}.R{pair}_readqc.fastq"),
    output:
        R1_bin_name                         = join(workpath, "{name}_R{pair}.fastq")
        R2_bin_name                         = join(workpath, "{name}_R{pair}.fastq")
    shell:
        """
            ln -s {input.R1_from_qc} {input.R1_bin_name}
            ln -s {input.R2_from_qc} {input.R2_bin_name}
        """


rule metawrap_assembly_binning:
    """
        TODO: docstring
    """
    input:
        assembly                    = join(top_assembly_dir, "{name}", "final_assembly.fasta"),
        
    params:
        bin_mem                     = int(cluster.get("metawrap_assembly_binning", default_memory)),
        bin_dir                     = join(workpath, "metawrap_binning", "{name}"),
    singularity: metawrap_container
    threads: cluster["metawrap_tax_classification"].get("threads", default_threads),
    shell:
        """
            metawrap binning --metabat2 --maxbin2 --concoct --interleaved -m {params.bin_mem} -t {threads} -a {input.assembly} -o $(ANALYSIS)/BINNING/$${sample}_bins $(ANALYSIS)/CLEAN_READS/$${sample}_*.fastq
        """
