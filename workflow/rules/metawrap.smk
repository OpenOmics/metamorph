rule metawrap_read_qc:
    """
        Quality-control step accomplishes three tasks:
            - Trims adapter sequences from input read pairs (.fastq.gz)
            - Quality filters from input read pairs
            - Removes human contamination from input read pairs
        Tools used:
            - [MetaWrap](https://github.com/bxlab/metaWRAP) orchestrates execution of 
              these tools in psuedo-pipeline:
                - [FastQC](https://github.com/s-andrews/FastQC)
                - [bmtagger](ftp://ftp.ncbi.nlm.nih.gov/pub/agarwala/bmtagger/)
                - [Trim-Galore](https://github.com/FelixKrueger/TrimGalore)
                - [Cutadapt](https://github.com/marcelm/cutadapt)

        @Input:
            Trimmed FastQ file (scatter)
        @Output:
            FastQC html report and zip file on trimmed data
        
        requires minimum of ~16gb memory for loading bmtagger index in memory (8gb)
    """
    input:
        R1 = join(workpath, '{name}.R1.fastq.gz'),
        R2 = join(workpath, '{name}.R2.fastq.gz'),
    output:
        R1_bmtagger_report  = join(workpath, 'metawrap_readqc', '{name}', '{name}.R1_bmtagger_report.html'),
        R2_bmtagger_report  = join(workpath, 'metawrap_readqc', '{name}', '{name}.R2_bmtagger_report.html'),
        R1_fastqc_report    = join(workpath, 'metawrap_readqc', '{name}', '{name}.R1_fastqc_report.html'),
        R2_fastqc_report    = join(workpath, 'metawrap_readqc', '{name}', '{name}.R2_fastqc_report.html'),
        R1_qc_reads         = join(workpath, '{name}.R1_readqc.fastq.gz'),
        R2_qc_reads         = join(workpath, '{name}.R2_readqc.fastq.gz'),
    params:
        readqc_dir          = join(workpath, 'metawrap_readqc', "{name}"),
        R1_mw_named         = join(workpath, '{name}_1.fastq'),
        R2_mw_named         = join(workpath, '{name}_2.fastq'),
        intermediate_qc_R1  = join(workpath, '{name}.R1_readqc.fastq'),
        intermediate_qc_R2  = join(workpath, '{name}.R2_readqc.fastq'),
    container: config['images']['metawrap'],
    threads: int(allocated("threads", "metawrap_read_qc", cluster)),
    shell: 
        """
        # Setting up inputs to metawrap,
        # needs uncompressed files as input,
        # input files also need to end with
        # the following ext: _[1,2].fastq
        gunzip -c {input.R1} > {params.R1_mw_named}
        gunzip -c {input.R2} > {params.R2_mw_named}
        
        # Running metawraps read_qc module
        metawrap read_qc \\
            -1 {params.R1_mw_named} \\
            -2 {params.R2_mw_named} \\
            -t {threads} \\
            -o {params.readqc_dir}
        
        # Rename R1 output files from metawrap
        mv {params.readqc_dir}/final_pure_reads_1.fastq {params.intermediate_qc_R1}
        gzip {params.intermediate_qc_R1}
        rm -f {params.R1_mw_named}
        cp {params.readqc_dir}/post-QC_report/final_pure_reads_1_fastqc.html {output.R1_bmtagger_report}
        cp {params.readqc_dir}/pre-QC_report/{params.R1_mw_named}_fastqc.html {output.R1_fastqc_report}

        # Rename R2 output files from metawrap
        mv {params.readqc_dir}/final_pure_reads_2.fastq  {params.intermediate_qc_R2}
        gzip {params.intermediate_qc_R2}
        rm -f {params.R2_mw_named}
        cp {params.readqc_dir}/post-QC_report/final_pure_reads_2_fastqc.html {output.R2_bmtagger_report}
        cp {params.readqc_dir}/pre-QC_report/{params.R2_mw_named}_fastqc.html {output.R2_fastqc_report}
        """