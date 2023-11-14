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
top_refine_dir             = join(workpath, "metawrap_refine_bins")
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
        R1                  = expand(join(workpath, "{name}.R1.fastq.gz"), name=samples),
        R2                  = expand(join(workpath, "{name}.R2.fastq.gz"), name=samples),
    output:
        R1_bmtagger_report  = expand(join(top_readqc_dir, "{name}", "{name}.R1_bmtagger_report.html"), name=samples),
        R2_bmtagger_report  = expand(join(top_readqc_dir, "{name}", "{name}.R2_bmtagger_report.html"), name=samples),
        R1_fastqc_report    = expand(join(top_readqc_dir, "{name}", "{name}.R1_fastqc_report.html"), name=samples),
        R2_fastqc_report    = expand(join(top_readqc_dir, "{name}", "{name}.R2_fastqc_report.html"), name=samples),
        R1_qc_reads         = expand(join(workpath, "{name}", "{name}.R1_readqc.fastq"), name=samples),
        R2_qc_reads         = expand(join(workpath, "{name}", "{name}.R2_readqc.fastq"), name=samples),
        readqc_dir          = expand(join(top_readqc_dir, "{name}"), name=samples),
    params:
        R1_mw_named         = expand(join(top_readqc_dir, "{name}", "{name}_1.fastq"), name=samples),
        R2_mw_named         = expand(join(top_readqc_dir, "{name}", "{name}_2.fastq"), name=samples),
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
            metawrap read_qc -1 {params.R1_mw_named} -2 {params.R2_mw_named} -t {threads} -o {output.readqc_dir}
            
            # Rename R1 output files from metawrap
            mv {output.readqc_dir}/final_pure_reads_1.fastq {output.R1_qc_reads}
            rm -f {params.R1_mw_named}
            cp {output.readqc_dir}/post-QC_report/final_pure_reads_1_fastqc.html {output.R1_bmtagger_report}
            cp {output.readqc_dir}/pre-QC_report/{params.R1_mw_named}_fastqc.html {output.R1_fastqc_report}

            # Rename R2 output files from metawrap
            mv {output.readqc_dir}/final_pure_reads_2.fastq {output.R2_qc_reads}
            rm -f {params.R2_mw_named}
            cp {output.readqc_dir}/post-QC_report/final_pure_reads_2_fastqc.html {output.R2_bmtagger_report}
            cp {output.readqc_dir}/pre-QC_report/{params.R2_mw_named}_fastqc.html {output.R2_fastqc_report}
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
        R1                          = expand(join(workpath, "{name}", "{name}.R1_readqc.fastq"), name=samples),
        R2                          = expand(join(workpath, "{name}", "{name}.R2_readqc.fastq"), name=samples),
    output:
        # megahit outputs
        megahit_assembly            = expand(join(top_assembly_dir, "{name}", "megahit", "final.contigs.fa"), name=samples),
        megahit_longcontigs         = expand(join(top_assembly_dir, "{name}", "megahit", "long.contigs.fa"), name=samples),
        megahit_log                 = expand(join(top_assembly_dir, "{name}", "megahit", "log"), name=samples),
        # metaspades outsputs
        metaspades_assembly         = expand(join(top_assembly_dir, "{name}", "metaspades", "contigs.fasta"), name=samples),
        metaspades_graph            = expand(join(top_assembly_dir, "{name}", "metaspades", "assembly_graph.fastg"), name=samples),
        metaspades_longscaffolds    = expand(join(top_assembly_dir, "{name}", "metaspades", "long_scaffolds.fasta"), name=samples),
        metaspades_scaffolds        = expand(join(top_assembly_dir, "{name}", "metaspades", "scaffolds.fasta"), name=samples),
        metaspades_cor_readsr1      = expand(join(top_assembly_dir, "{name}", "metaspades", "{name}_1.fastq.00.0_0.cor.fastq.gz"), name=samples),
        metaspades_cor_readsr2      = expand(join(top_assembly_dir, "{name}", "metaspades", "corrected", "{name}_2.fastq.00.0_0.cor.fastq.gz"), name=samples),
        # ensemble outputs
        final_assembly              = expand(join(top_assembly_dir, "{name}", "final_assembly.fasta"), name=samples),
        final_assembly_report       = expand(join(top_assembly_dir, "{name}", "assembly_report.html"), name=samples),
        assembly_R1                 = expand(join(top_assembly_dir, "{name}", "{name}_1.fastq"), name=samples),
        assembly_R2                 = expand(join(top_assembly_dir, "{name}", "{name}_2.fastq"), name=samples),
        assembly_dir                = expand(join(top_assembly_dir, "{name}"), name=samples),
    singularity: metawrap_container,
    params:
        memlimit                    = cluster["metawrap_genome_assembly"].get('mem', default_memory),
        contig_min_len              = "1000",
    threads: int(cluster["metawrap_genome_assembly"].get('threads', default_threads)),
    shell:
        """
            # link to the file names metawrap expects
            ln -s {input.R1} {output.assembly_R1}
            ln -s {input.R2} {output.assembly_R2}
            # run genome assembler
            metawrap assembly --megahit --metaspades -m {params.memlimit} -t {threads} -l {params.contig_min_len} -1 {output.assembly_R1} -2 {output.assembly_R2} -o {output.assembly_dir}
        """


rule metawrap_tax_classification:
    """
        TODO: docstring
    """
    input:
        reads                       = expand(join(workpath, "{name}", "{name}.R{pair}_readqc.fastq"), name=samples, pair=['1', '2']),
        final_assembly              = expand(join(top_assembly_dir, "{name}", "final_assembly.fasta"), name=samples),
    output:
        krak2_asm                   = expand(join(top_tax_dir, "{name}", "final_assembly.krak2"), name=samples),
        kraken2_asm                 = expand(join(top_tax_dir, "{name}", "final_assembly.kraken2"), name=samples),
        krona_asm                   = expand(join(top_tax_dir, "{name}", "final_assembly.krona"), name=samples),
        kronagram                   = expand(join(top_tax_dir, "{name}", "kronagram.html"), name=samples),
        tax_dir                     = expand(join(top_tax_dir, "{name}"), name=samples),
    params:
        tax_subsample               = str(int(1e6)),
    singularity: metawrap_container,
    threads: cluster["metawrap_tax_classification"].get("threads", default_threads),
    shell:
        """
            metawrap kraken2 -t {threads} -s {params.tax_subsample} -o {output.tax_dir} {input.final_assembly} {input.reads}
        """


rule metawrap_setup_binning:
    input:
        R1_from_qc                  = expand(join(workpath, "{name}", "{name}.R1_readqc.fastq"), name=samples),
        R2_from_qc                  = expand(join(workpath, "{name}", "{name}.R2_readqc.fastq"), name=samples),
    output:
        R1_bin_name                 = expand(join(workpath, "{name}", "{name}_{pair}.fastq"), name=samples, pair=['1']),
        R2_bin_name                 = expand(join(workpath, "{name}", "{name}_{pair}.fastq"), name=samples, pair=['2']),
    shell:
        """
            ln -s {input.R1_from_qc} {output.R1_bin_name}
            ln -s {input.R2_from_qc} {output.R2_bin_name}
        """
        

rule metawrap_binning:
    """
        TODO: docstring
    """
    input:
        paired_reads                = expand(join(workpath, "{name}", "{name}_{pair}.fastq"), name=samples, pair=['1', '2']),
        assembly                    = expand(join(top_assembly_dir, "{name}", "final_assembly.fasta"), name=samples),
    output:
        bin_dir                     = expand(join(top_binning_dir, "{name}"), name=samples),
        refine_dir                  = expand(join(top_refine_dir, "{name}"),  name=samples),
    params:
        bin_mem                     = cluster.get("metawrap_assembly_binning", default_memory),
        # minimum percentage completions of bins
        min_perc_complete           = "50",
        # maximum percentage of contamination
        max_perc_contam             = "5",
    singularity: metawrap_container,
    threads: cluster["metawrap_tax_classification"].get("threads", default_threads),
    shell:
        """
            metawrap binning --metabat2 --maxbin2 --concoct --interleaved -m {params.bin_mem} -t {threads} -a {input.assembly} -o {output.bin_dir} {input.paired_reads}
            metawrap bin_refinement -o {output.refine_dir} -t {threads} -A {output.bin_dir}/metabat2_bins -B {output.bin_dir}/maxbin2_bins -c {params.min_perc_complete} -x {params.max_perc_contam}
        """
