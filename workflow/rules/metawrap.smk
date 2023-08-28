# ~~~~~~~~~~
# Metawrap metagenome assembly and analysis rules
# ~~~~~~~~~~
default_threads            = cluster["__default__"]['threads']
default_memory             = cluster["__defualt__"]['mem']
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
              these tools in psuedo-pipeline:
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
        R1_bmtagger_report  = join(workpath, "metawrap_readqc", "{name}", "{name}.R1_bmtagger_report.html"),
        R2_bmtagger_report  = join(workpath, "metawrap_readqc", "{name}", "{name}.R2_bmtagger_report.html"),
        R1_fastqc_report    = join(workpath, "metawrap_readqc", "{name}", "{name}.R1_fastqc_report.html"),
        R2_fastqc_report    = join(workpath, "metawrap_readqc", "{name}", "{name}.R2_fastqc_report.html"),
        R1_qc_reads         = join(workpath, "{name}.R1_readqc.fastq.gz"),
        R2_qc_reads         = join(workpath, "{name}.R2_readqc.fastq.gz"),
    params:
        readqc_dir          = join(workpath, "metawrap_readqc", "{name}"),
        R1_mw_named         = join(workpath, "{name}_1.fastq"),
        R2_mw_named         = join(workpath, "{name}_2.fastq"),
        intermediate_qc_R1  = join(workpath, "{name}.R1_readqc.fastq"),
        intermediate_qc_R2  = join(workpath, "{name}.R2_readqc.fastq"),
    container: metawrap_container,
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
            mv {params.readqc_dir}/final_pure_reads_1.fastq {params.intermediate_qc_R1}
            gzip {params.intermediate_qc_R1}
            rm -f {params.R1_mw_named}
            cp {params.readqc_dir}/post-QC_report/final_pure_reads_1_fastqc.html {output.R1_bmtagger_report}
            cp {params.readqc_dir}/pre-QC_report/{params.R1_mw_named}_fastqc.html {output.R1_fastqc_report}

            # Rename R2 output files from metawrap
            mv {params.readqc_dir}/final_pure_reads_2.fastq {params.intermediate_qc_R2}
            gzip {params.intermediate_qc_R2}
            rm -f {params.R2_mw_named}
            cp {params.readqc_dir}/post-QC_report/final_pure_reads_2_fastqc.html {output.R2_bmtagger_report}
            cp {params.readqc_dir}/pre-QC_report/{params.R2_mw_named}_fastqc.html {output.R2_fastqc_report}
        """


rule metawrap_genome_assembly:
    """
        @Input:
            Clean trimmed fastq.gz reads (R1 & R2 per sample)
        @Output:
            Megahit assembled contigs and reports
            Metaspades assembled contigs and reports
            Ensemble assembled contigs and reports
    """
    input:
        R1                          = join(workpath, "{name}.R1_readqc.fastq.gz"),
        R2                          = join(workpath, "{name}.R2_readqc.fastq.gz"),
    output:
        assembly_dir                = join(workpath, "metawrap_assembly", "{name}"),
        # megahit outputs
        megahit_dir                 = join(assembly_dir, "{name}", "megahit"),
        megahit_assembly            = join(megahit_dir, "final.contigs.fa"),
        megahit_longcontigs         = join(megahit_dir, "long.contigs.fa"),
        megahit_log                 = join(megahit_dir, "log"),
        # metaspades outsputs
        metaspades_dir              = join(workpath, "metawrap_assembly", "{name}", "metaspades"),
        metaspades_assembly         = join(metaspades_dir, "contigs.fasta"),
        metaspades_graph            = join(metaspades_dir, "assembly_graph.fastg"),
        metaspades_longscaffolds    = join(metaspades_dir, "long_scaffolds.fasta"),
        metaspades_scaffolds        = join(metaspades_dir, "scaffolds.fasta"),
        metaspades_cor_readsr1      = join(metaspades_dir, "{name}_1.fastq.00.0_0.cor.fastq.gz"),
        metaspades_cor_readsr2      = join(metaspades_dir, "corrected", "{name}_2.fastq.00.0_0.cor.fastq.gz"),
        # ensemble outputs
        final_assembly              = join(assembly_dir, "{name}", "final_assembly.fasta"),
        final_assembly_report       = join(assembly_dir, "{name}", "assembly_report.html"),
    container: metawrap_container,
    params:
        memlimit                = cluster["metawrap_genome_assembly"].get('memory', default_memory),
        contig_min_len          = "1000",
        assembly_R1             = join(workpath, "{name}_1.fastq.gz"),
        assembly_R2             = join(workpath, "{name}_2.fastq.gz"),
    threads: int(cluster["metawrap_genome_assembly"].get('threads', default_threads)),
    shell:
        """
            ln {input.R1} {params.assembly_R1}
            ln {input.R2} {params.assembly_R2}
            metawrap assembly --megahit --metaspades -m {params.memlimit} -t {threads} -l {params.contig_min_len} -1 {params.assembly_R1} -2 {params.assembly_R2} -o {output.assembly_dir}
        """
