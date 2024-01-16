# ~~~~~~~~~~
# Metawrap metagenome assembly and analysis rules
# ~~~~~~~~~~
from os.path import join
from itertools import chain


# ~~~~~~~~~~
# Constants and paths
# ~~~~~~~~~~
workpath                   = config["project"]["workpath"]
datapath                   = config["project"]["datapath"]
default_threads            = cluster["__default__"]['threads']
default_memory             = cluster["__default__"]['mem']
top_log_dir                = join(workpath, config['project']['id'], "logs")
top_readqc_dir             = join(workpath, config['project']['id'], "metawrap_qc")
top_assembly_dir           = join(workpath, config['project']['id'], "metawrap_assembly")
top_tax_dir                = join(workpath, config['project']['id'], "metawrap_kmer")
top_binning_dir            = join(workpath, config['project']['id'], "metawrap_binning")
top_refine_dir             = join(workpath, config['project']['id'], "metawrap_bin_refine")
metawrap_container         = config["containers"]["metawrap"]


"""
    1. read qc
    2. assembly
    3. binning
    4. bin refinement
    5. depreplicate bins
    6. annotate bins
    7. index depreplicated genomes
    8. align DNA to assembly
""" 


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
        R1                  = join(datapath, "{name}_R1.fastq.gz"),
        R2                  = join(datapath, "{name}_R2.fastq.gz"),
    output:
        R1_bmtagger_report  = join(top_readqc_dir, "{name}", "{name}.R1_bmtagger_report.html"),
        R2_bmtagger_report  = join(top_readqc_dir, "{name}", "{name}.R2_bmtagger_report.html"),
        R1_fastqc_report    = join(top_readqc_dir, "{name}", "{name}.R1_fastqc_report.html"),
        R2_fastqc_report    = join(top_readqc_dir, "{name}", "{name}.R2_fastqc_report.html"),
        R1_qc_reads         = join(workpath, "{name}", "{name}.R1_readqc.fastq"),
        R2_qc_reads         = join(workpath, "{name}", "{name}.R2_readqc.fastq"),
    params:
        rname               = "metawrap_read_qc",
        this_qc_dir         = join(top_readqc_dir, "{name}"),
        R1_mw_named         = join(top_readqc_dir, "{name}", "{name}_1.fastq"),
        R2_mw_named         = join(top_readqc_dir, "{name}", "{name}_2.fastq"),
    containerized: metawrap_container,
    log: join(top_log_dir, "read_qc", "{name}")
    threads: int(cluster["metawrap_genome_assembly"].get('threads', default_threads)),
    shell: 
        """
        gunzip -c {input.R1} > {params.R1_mw_named}
        gunzip -c {input.R2} > {params.R2_mw_named}
        metawrap read_qc -1 {params.R1_mw_named} -2 {params.R2_mw_named} -t {threads} -o {params.this_qc_dir}
        mv {params.this_qc_dir}/final_pure_reads_1.fastq {output.R1_qc_reads}
        rm -f {params.R1_mw_named}
        cp {params.this_qc_dir}/post-QC_report/final_pure_reads_1_fastqc.html {output.R1_bmtagger_report}
        cp {params.this_qc_dir}/pre-QC_report/{params.R1_mw_named}_fastqc.html {output.R1_fastqc_report}
        mv {params.this_qc_dir}/final_pure_reads_2.fastq {output.R2_qc_reads}
        rm -f {params.R2_mw_named}
        cp {params.this_qc_dir}/post-QC_report/final_pure_reads_2_fastqc.html {output.R2_bmtagger_report}
        cp {params.this_qc_dir}/pre-QC_report/{params.R2_mw_named}_fastqc.html {output.R2_fastqc_report}
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
        rname                       = "metawrap_genome_assembly",
        memlimit                    = cluster["metawrap_genome_assembly"].get('mem', default_memory),
        contig_min_len              = "1000",
    threads: int(cluster["metawrap_genome_assembly"].get('threads', default_threads)),
    shell:
        """
            # link to the file names metawrap expects
            ln -s {input.R1} {output.assembly_R1}
            ln -s {input.R2} {output.assembly_R2}
            # run genome assembler
            metawrap assembly \
            --megahit \
            --metaspades \
            -m {params.memlimit} \
            -t {threads} \
            -l {params.contig_min_len} \
            -1 {output.assembly_R1} \
            -2 {output.assembly_R2} \
            -o {output.assembly_dir}
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
        rname                       = "metawrap_tax_classification",
        tax_subsample               = str(int(1e6)),
    singularity: metawrap_container,
    threads: cluster["metawrap_tax_classification"].get("threads", default_threads),
    shell:
        """
            metawrap kraken2 \
            -t {threads} \
            -s {params.tax_subsample} \
            -o {output.tax_dir} \
            {input.final_assembly} \
            {input.reads}
        """


rule metawrap_setup_binning:
    input:
        R1_from_qc                  = join(workpath, "{name}", "{name}.R1_readqc.fastq"),
        R2_from_qc                  = join(workpath, "{name}", "{name}.R2_readqc.fastq"),
    output:
        R1_bin_name                 = join(workpath, "{name}", "{name}_1.fastq"),
        R2_bin_name                 = join(workpath, "{name}", "{name}_2.fastq"),
    params:
        rname                       = "metawrap_setup_binning",
    shell:
        """
        ln -s {input.R1_from_qc} {output.R1_bin_name}
        ln -s {input.R2_from_qc} {output.R2_bin_name}
        """


rule metawrap_binning:
    input:
        r1_read                     = join(workpath, "{name}", "{name}_1.fastq"),
        r2_read                     = join(workpath, "{name}", "{name}_2.fastq"),
        assembly                    = join(top_assembly_dir, "{name}", "final_assembly.fasta"),
    output:
        maxbin_bins                 = directory(join(top_binning_dir, "{name}", "maxbin2_bins")),
        metabat2_bins               = directory(join(top_binning_dir, "{name}", "metabat2_bins")),
        metawrap_bins               = directory(join(top_binning_dir, "{name}", "metawrap_50_5_bins")),
        maxbin_contigs              = join(top_binning_dir, "{name}", "maxbin2_bins.contigs"),
        maxbin_stats                = join(top_binning_dir, "{name}", "maxbin2_bins.stats"),
        metabat2_contigs            = join(top_binning_dir, "{name}", "metabat2_bins.contigs"),
        metabat2_stats              = join(top_binning_dir, "{name}", "metabat2_bins.stats"),
        metawrap_contigs            = join(top_binning_dir, "{name}", "metawrap_50_5_bins.contigs"),
        metawrap_stats              = join(top_binning_dir, "{name}", "metawrap_50_5_bins.stats"),
        bin_figure                  = join(top_binning_dir, "{name}", "figures", "binning_results.png"),
    params:
        rname                       = "metawrap_binning",
        bin_dir                     = join(top_binning_dir, "{name}"),
        refine_dir                  = join(top_refine_dir, "{name}"),
        bin_mem                     = cluster.get("metawrap_assembly_binning", default_memory),
        min_perc_complete           = "50",
        max_perc_contam             = "5",
    singularity: metawrap_container,
    threads: cluster["metawrap_tax_classification"].get("threads", default_threads),
    shell:
        """
        # metawrap binning
        metawrap binning \
        --metabat2 --maxbin2 --concoct \
        -m {params.bin_mem} \
        -t {threads} \
        -a {input.assembly} \
        -o {params.bin_dir} \
        {input.r1_read} {input.r2_read}
        # metawrap bin refinement
        metawrap bin_refinement \
        -o {params.refine_dir} \
        -t {threads} \
        -A {params.bin_dir}/metabat2_bins \
        -B {params.bin_dir}/maxbin2_bins \
        -c {params.min_perc_complete} \
        -x {params.max_perc_contam}
        """


rule derep_bins:
    input:
        maxbin_contigs              = join(top_binning_dir, "{name}", "maxbin2_bins.contigs"),
        maxbin_stats                = join(top_binning_dir, "{name}", "maxbin2_bins.stats"),
        metabat2_contigs            = join(top_binning_dir, "{name}", "metabat2_bins.contigs"),
        metabat2_stats              = join(top_binning_dir, "{name}", "metabat2_bins.stats"),
        metawrap_contigs            = join(top_binning_dir, "{name}", "metawrap_50_5_bins.contigs"),
        metawrap_stats              = join(top_binning_dir, "{name}", "metawrap_50_5_bins.stats"),
    output:
        dereplicated_bins           = directory(join(top_refine_dir, "{name}", "dereplicated_bins")),
        search_rep_bc               = join(top_binning_dir, "{name}", "maxbin2_bins"),
    singularity: metawrap_container,
    threads: 32
    params:
        rname                       = "derep_bins",
        sid                         = "{name}",
        maxbin_bins                 = join(top_binning_dir, "{name}", "maxbin2_bins"),
        metabat2_bins               = join(top_binning_dir, "{name}", "metabat2_bins"),
        metawrap_bins               = join(top_binning_dir, "{name}", "metawrap_50_5_bins"),
        # -l LENGTH: Minimum genome length (default: 50000)
        minimum_genome_length = "1000",
        # -pa[--P_ani] P_ANI: ANI threshold to form primary (MASH) clusters (default: 0.9)
        ani_primary_threshold = "0.9",
        # -sa[--S_ani] S_ANI: ANI threshold to form secondary clusters (default: 0.95)
        ani_secondary_threshold = "0.95",
        # -nc[--cov_thresh] COV_THRESH: Minmum level of overlap between genomes when doing secondary comparisons (default: 0.1)
        min_overlap = "0.1",
        # -cm[--coverage_method] {total,larger}  {total,larger}: Method to calculate coverage of an alignment
        coverage_method = 'larger',
    shell:
        """
		sed -i 's/^bin./{name}_bin./g' {input.maxbin_stats} && \
		sed -i 's/^bin./{name}_bin./g' {input.metabat2_stats} && \
		sed -i 's/^bin./{name}_bin./g' {input.metawrap_stats} && \
        touch {output.search_rep_bc}
        dRep dereplicate \
        -p {threads} \
        -l {params.minimum_genome_length} \
        -pa {params.ani_primary_threshold} \
        -sa {params.ani_secondary_threshold} \
        -nc {params.min_overlap} \
        -cm {params.coverage_method} \
        -g {params.metawrap_bins}/* \
        {output.dereplicated_bins}
        """
