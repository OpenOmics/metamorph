# ~~~~~~~~~~
# Metawrap metagenome assembly and analysis rules
# ~~~~~~~~~~
from os import listdir
from os.path import join, basename
from itertools import chain
from uuid import uuid4


# ~~~~~~~~~~
# Constants and paths
# ~~~~~~~~~~
# resource defaults
default_threads            = cluster["__default__"]['threads']
default_memory             = cluster["__default__"]['mem']

# directories
workpath                   = config["project"]["workpath"]
datapath                   = config["project"]["datapath"]
samples                    = config["samples"] if not coassemble else ['concatenated']
top_log_dir                = join(workpath, "logfiles")
top_readqc_dir             = join(workpath, config['project']['id'], "metawrap_read_qc")
top_trim_dir               = join(workpath, config['project']['id'], "trimmed_reads")
top_assembly_dir           = join(workpath, config['project']['id'], "metawrap_assembly")
top_tax_dir                = join(workpath, config['project']['id'], "metawrap_kmer")
top_binning_dir            = join(workpath, config['project']['id'], "metawrap_binning")
top_refine_dir             = join(workpath, config['project']['id'], "metawrap_bin_refinement")

# workflow flags
metawrap_container         = config["containers"]["metawrap"]
pairedness                 = list(range(1, config['project']['nends']+1))
mem2int                    = lambda x: int(str(x).lower().replace('gb', '').replace('g', ''))

"""
    Step-wise pipeline outline:
        0. concat or no concat
        1. read qc
        2. assembly
        3. binning
        4. bin refinement
        5. depreplicate bins
        6. annotate bins
        7. index depreplicated genomes
        8. align DNA to assembly
"""

rule concat_reads:
    input:
        all_r1_reads                = expand(join(workpath, "dna", "{sid}_R1.fastq.gz"), sid=config['samples'] if config['coassembly'] else []),
        all_r2_reads                = expand(join(workpath, "dna", "{sid}_R2.fastq.gz"), sid=config['samples'] if config['coassembly'] else []),
    output:
        big_compressed_read_r1      = join(workpath, "dna", "concatenated_R1.fastq.gz"),
        big_compressed_read_r2      = join(workpath, "dna", "concatenated_R2.fastq.gz"),
        big_read1_hash              = join(workpath, "dna", "concatenated_R1.md5"),
        big_read2_hash              = join(workpath, "dna", "concatenated_R2.md5"),
    params:
        big_read_r1                 = join(workpath, "dna", "concatenated_R1.fastq"),
        big_read_r2                 = join(workpath, "dna", "concatenated_R2.fastq"),
        rname                       = "concat_reads",
        input_dir                   = join(workpath, "dna"),
    threads: int(cluster["concat_reads"].get('threads', default_threads)),
    shell: 
        """
        shopt -s extglob
        # concat r1
        for fastq in {params.input_dir}/*R1*f?(ast)q*; do
            ext=$(echo "${{fastq: -2}}" | tr '[:upper:]' '[:lower:]')
            if [[ "$ext" == "gz" ]]; then
                zcat $fastq >> {params.big_read_r1}
            else
                cat $fastq >> {params.big_read_r1}
            fi;
        done

        # concat r2
        for fastq in {params.input_dir}/*R2*f?(ast)q*; do 
            ext=$(echo "${{fastq: -2}}" | tr '[:upper:]' '[:lower:]')
            if [[ "$ext" == "gz" ]]; then
                zcat $fastq > {params.big_read_r2}
            else
                cat $fastq > {params.big_read_r2}
            fi;
        done
        shopt -u extglob

        pigz -9 -p 28 -c {params.big_read_r1} > {output.big_compressed_read_r1}
        pigz -9 -p 28 -c {params.big_read_r2} > {output.big_compressed_read_r2}
        md5sum {output.big_compressed_read_r1} > {output.big_read1_hash}
        md5sum {output.big_compressed_read_r2} > {output.big_read2_hash}
        rm {params.big_read_r1} {params.big_read_r2}
        """


rule metawrap_read_qc:
    """
    The metaWRAP::Read_qc module is meant to pre-process raw Illumina sequencing reads in preparation for assembly and alignment. 
    The raw reads are trimmed based on adapted content and PHRED scored with the default setting of Trim-galore, ensuring that only high-quality 
    sequences are left. Then reads are then aligned to the host genome (e.g. human) with bmtagger, and any host reads are removed from the 
    metagenomic data to remove host contamination. Read pairs where only one read was aligned to the host genome are also removed. 
    Finally, FASTQC is used to generate quality reports of the raw and final read sets in order to assess read quality improvement. 
    The users have control over which of the above features they wish to use.

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
        - Raw fastq reads (R1 & R2 per-sample) from the instrument
        
    @Outputs:
        - Trimmed, quality filtered reads (R1 & R2)
        - FastQC html report and zip file on trimmed data
    """
    input:
        R1                  = join(workpath, "dna", "{name}_R1.fastq.gz"),
        R2                  = join(workpath, "dna", "{name}_R2.fastq.gz"),
    output:
        R1_pretrim_report   = join(top_readqc_dir, "{name}", "{name}_R1_pretrim_report.html"),
        R2_pretrim_report   = join(top_readqc_dir, "{name}", "{name}_R2_pretrim_report.html"),
        R1_postrim_report   = join(top_readqc_dir, "{name}", "{name}_R1_postrim_report.html"),
        R2_postrim_report   = join(top_readqc_dir, "{name}", "{name}_R2_postrim_report.html"),
        R1_trimmed          = join(top_trim_dir, "{name}", "{name}_R1_trimmed.fastq"),
        R2_trimmed          = join(top_trim_dir, "{name}", "{name}_R2_trimmed.fastq"),
        R1_trimmed_gz       = join(top_trim_dir, "{name}", "{name}_R1_trimmed.fastq.gz"),
        R2_trimmed_gz       = join(top_trim_dir, "{name}", "{name}_R2_trimmed.fastq.gz"),
    params:
        rname               = "metawrap_read_qc",
        sid                 = "{name}",
        this_qc_dir         = join(top_readqc_dir, "{name}"),
        trim_out            = join(top_trim_dir, "{name}"),
        tmp_safe_dir        = join(config['options']['tmp_dir'], 'read_qc'),
        tmpr1               = lambda _, output, input: join(config['options']['tmp_dir'], 'read_qc', str(basename(str(input.R1))).replace('_R1.', '_1.').replace('.gz', '')),
        tmpr2               = lambda _, output, input: join(config['options']['tmp_dir'], 'read_qc', str(basename(str(input.R2))).replace('_R2.', '_2.').replace('.gz', '')),
    containerized: metawrap_container,
    threads: int(cluster["metawrap_genome_assembly"].get('threads', default_threads)),
    shell: 
        """
            # safe temp directory
            if [ ! -d "{params.tmp_safe_dir}" ]; then mkdir -p "{params.tmp_safe_dir}"; fi
            tmp=$(mktemp -d -p "{params.tmp_safe_dir}")
            trap 'rm -rf "{params.tmp_safe_dir}"' EXIT

            # uncompress to lscratch
            rone="{input.R1}"
            ext=$(echo "${{rone: -2}}" | tr '[:upper:]' '[:lower:]')
            if [[ "$ext" == "gz" ]]; then
                zcat {input.R1} > {params.tmpr1}
            else
                ln -s {input.R1} {params.tmpr1}
            fi;
            rtwo="{input.R2}"
            ext=$(echo "${{rtwo: -2}}" | tr '[:upper:]' '[:lower:]')
            if [[ "$ext" == "gz" ]]; then
                zcat {input.R2} > {params.tmpr2}
            else
                ln -s {input.R2} {params.tmpr2}
            fi;

            # read quality control, host removal
            # TODO: add support for mouse reads (mm10 genome prefix, "-x mm10")
            mw read_qc -1 {params.tmpr1} -2 {params.tmpr2} -t {threads} -o {params.this_qc_dir}

            # collate fastq outputs to facilitate workflow, compress
            ln -s {params.this_qc_dir}/final_pure_reads_1.fastq {params.trim_out}/{params.sid}_R1_trimmed.fastq
            ln -s {params.this_qc_dir}/final_pure_reads_2.fastq {params.trim_out}/{params.sid}_R2_trimmed.fastq
            pigz -9 -p {threads} -c {params.this_qc_dir}/final_pure_reads_1.fastq  > {params.trim_out}/{params.sid}_R1_trimmed.fastq.gz
            pigz -9 -p {threads} -c {params.this_qc_dir}/final_pure_reads_2.fastq > {params.trim_out}/{params.sid}_R2_trimmed.fastq.gz
            
            # collate outputs to facilitate
            mkdir -p {params.this_qc_dir}/dna
            ln -s {params.this_qc_dir}/post-QC_report/final_pure_reads_1_fastqc.html {params.this_qc_dir}/{params.sid}_R1_postrim_report.html
            ln -s {params.this_qc_dir}/post-QC_report/final_pure_reads_2_fastqc.html {params.this_qc_dir}/{params.sid}_R2_postrim_report.html
            ln -s {params.this_qc_dir}/pre-QC_report/{params.sid}_1_fastqc.html {params.this_qc_dir}/{params.sid}_R1_pretrim_report.html
            ln -s {params.this_qc_dir}/pre-QC_report/{params.sid}_2_fastqc.html {params.this_qc_dir}/{params.sid}_R2_pretrim_report.html
        """


rule metawrap_genome_assembly:
    """
        The metaWRAP::Assembly module allows the user to assemble a set of metagenomic reads with either 
        metaSPAdes or MegaHit (default). While metaSPAdes results in a superior assembly in most samples, 
        MegaHit is much more memory efficient, faster, and scales well with large datasets. 
        In addition to simplifying parameter selection for the user, this module also sorts and formats 
        the MegaHit assembly in a way that makes it easier to inspect. The contigs are sorted by length 
        and their naming is changed to resemble that of SPAdes, including the contig ID, length, and coverage. 
        Finally, short scaffolds are discarded (<1000bp), and an assembly report is generated with QUAST.

        @Input:
            Clean trimmed fastq reads (R1 & R2 per sample)

        @Output:
            Megahit assembled contigs and reports
            Metaspades assembled contigs and reports
            Ensemble assembled contigs and reports
    """
    input:
        R1                          = join(top_trim_dir, "{name}", "{name}_R1_trimmed.fastq"),
        R2                          = join(top_trim_dir, "{name}", "{name}_R2_trimmed.fastq"),
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
        metaspades_cor_readsr1      = expand(join(top_assembly_dir, "{name}", "metaspades", "corrected", "{name}_1.00.0_0.cor.fastq.gz"), name=samples),
        metaspades_cor_readsr2      = expand(join(top_assembly_dir, "{name}", "metaspades", "corrected", "{name}_2.00.0_0.cor.fastq.gz"), name=samples),
        # ensemble outputs
        final_assembly              = expand(join(top_assembly_dir, "{name}", "final_assembly.fasta"), name=samples),
        final_assembly_report       = expand(join(top_assembly_dir, "{name}", "assembly_report.html"), name=samples),
        assembly_R1                 = expand(join(top_assembly_dir, "{name}", "{name}_1.fastq"), name=samples),
        assembly_R2                 = expand(join(top_assembly_dir, "{name}", "{name}_2.fastq"), name=samples),
    singularity: metawrap_container,
    params:
        rname                       = "metawrap_genome_assembly",
        memlimit                    = mem2int(cluster["metawrap_genome_assembly"].get('mem', default_memory)),
        mh_dir                      = expand(join(top_assembly_dir, "{name}", "megahit"), name=samples),
        assembly_dir                = expand(join(top_assembly_dir, "{name}"), name=samples),
        contig_min_len              = "1000",
    threads: int(cluster["metawrap_genome_assembly"].get('threads', default_threads)),
    shell:
        """
            # remove empty directories by snakemake, to prevent metawrap error
            rm -rf {params.mh_dir:q}
            # link to the file names metawrap expects
            ln -s {input.R1} {output.assembly_R1}
            ln -s {input.R2} {output.assembly_R2}
            # run genome assembler
            mw assembly \
            --megahit \
            --metaspades \
            -m {params.memlimit} \
            -t {threads} \
            -l {params.contig_min_len} \
            -1 {output.assembly_R1} \
            -2 {output.assembly_R2} \
            -o {params.assembly_dir}
        """


rule metawrap_tax_classification:
    """
        The metaWRAP::Taxonomic Classification module takes in any number of fastq or fasta files, classifies the contained sequences 
        with KRAKEN, and reports the taxonomy distribution in a kronagram using KronaTools. If the sequences passed to the module 
        belong to an assembly and follow the contig naming convention of the Assembly module, the taxonomy of each contig is weighted based on 
        its length and coverage [weight=coverage*length]. The classifications of the sequences are then summarized in a format that 
        KronaTools' ktImportText function recognizes, and a final kronagram in html format is made with all the samples.

        @Input:
            - clean & trimmed reads (R1 & R2)
            - ensemble genome assembly

        @Output:
            - kraken2 kmer classification reports and tabular data outputs
            - krona tabular outputs
            - krona plot (interactive circular pie charts) of classified taxonomies 

    """
    input:
        r1                          = join(top_trim_dir, "{name}", "{name}_R1_trimmed.fastq"), 
        r2                          = join(top_trim_dir, "{name}", "{name}_R2_trimmed.fastq"),
        final_assembly              = join(top_assembly_dir, "{name}", "final_assembly.fasta"),
    output:
        krak2_asm                   = expand(join(top_tax_dir, "{name}", "final_assembly.krak2"), name=samples),
        kraken2_asm                 = expand(join(top_tax_dir, "{name}", "final_assembly.kraken2"), name=samples),
        krona_asm                   = expand(join(top_tax_dir, "{name}", "final_assembly.krona"), name=samples),
        kronagram                   = expand(join(top_tax_dir, "{name}", "kronagram.html"), name=samples),
    params:
        tax_dir                     = expand(join(top_tax_dir, "{name}"), name=samples),
        rname                       = "metawrap_tax_classification",
        tax_subsample               = str(int(1e6)),
        reads                       = lambda _, output, input: ' '.join([input.R1, input.R2]),
    singularity: metawrap_container,
    threads: int(cluster["metawrap_tax_classification"].get("threads", default_threads)),
    shell:
        """
            mkdir -p """+top_tax_dir+"""
            mw kraken2 \
            -t {threads} \
            -s {params.tax_subsample} \
            -o {params.tax_dir} \
            {input.final_assembly} \
            {params.reads}
        """


rule metawrap_binning:
    """
        The metaWRAP::Binning module is meant to be a convenient wrapper around three metagenomic binning software: MaxBin2, metaBAT2, and CONCOCT.

        First the metagenomic assembly is indexed with bwa-index, and then paired end reads from any number of samples are aligned to it. 
        The alignments are sorted and compressed with samtools, and library insert size statistics are also gathered at the same time 
        (insert size average and standard deviation). 

        metaBAT2's `jgi_summarize_bam_contig_depths` function is used to generate contig adundance table, and it is then converted 
        into the correct format for each of the three binners to take as input. After MaxBin2, metaBAT2, and CONCOCT finish binning 
        the contigs with default settings (the user can specify which software he wants to bin with), the final bins folders are 
        created with formatted bin fasta files for easy inspection. 
        
        Optionally, the user can chose to immediately run CheckM on the bins to determine the success of the binning. 
        
        CheckM's `lineage_wf` function is used to predict essential genes and estimate the completion and 
        contamination of each bin, and a custom script is used to generate an easy to view report on each bin set.

        Orchestrates execution of this ensemble of metagenomic binning software:
            - MaxBin2
            - metaBAT2
            - CONCOCT

        Stepwise algorithm for metawrap genome binning:
            1. Alignment
                a. metagenomic assembly is indexed with bwa-index
                b. paired end reads from any number of samples are aligned
            2. Collate, sort, stage alignment outputs 
                a. Alignments are sorted and compressed with samtools, and library insert size statistics are also 
                   gathered at the same time (insert size average and standard deviation). 
                b. metaBAT2's `jgi_summarize_bam_contig_depths` function is used to generate contig adundance table
                    i. It is converted into the correct format for each of the three binners to take as input. 
            3. Ensemble binning 
                a. MaxBin2, metaBAT2, CONCOCT are run with default settings
                b. Final bins directories are created with formatted bin fasta files for easy inspection. 
            4. CheckM is run on the bin output from step 3 to determine the success of the binning. 
                a. CheckM's `lineage_wf` function is used to predict essential genes and estimate the completion and contamination of each bin.
                b. Outputs are formatted and collected for better viewing.

        @Input:
            Clean trimmed fastq.gz reads (R1 & R2 per sample)

        @Output:
            Megahit assembled draft-genomes (bins) and reports
            Metaspades assembled draft-genomes (bins) and reports
            Ensemble assembled draft-genomes (bins) and reports

    """
    input:
        R1                          = join(top_trim_dir, "{name}", "{name}_R1_trimmed.fastq"),
        R2                          = join(top_trim_dir, "{name}", "{name}_R2_trimmed.fastq"),
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
        bin_parent_dir              = top_binning_dir,
        bin_dir                     = join(top_binning_dir, "{name}"),
        bin_refine_dir              = join(top_binning_dir, "{name}"),
        bin_mem                     = mem2int(cluster['metawrap_binning'].get("mem", default_memory)),
        mw_trim_linker_R1           = join(top_trim_dir, "{name}", "{name}_1.fastq"),
        mw_trim_linker_R2           = join(top_trim_dir, "{name}", "{name}_2.fastq"),
        min_perc_complete           = "50",
        max_perc_contam             = "5",
    singularity: metawrap_container,
    threads: int(cluster["metawrap_binning"].get("threads", default_threads)),
    shell:
        """
        # set checkm data
        export CHECKM_DATA_PATH="/data2/CHECKM_DB"

        # make base dir if not exists
        mkdir -p {params.bin_parent_dir}
        if [ -d "{params.bin_dir}" ]; then rm -rf {params.bin_dir}; fi

        # setup links for metawrap input
        [[ -f "{params.mw_trim_linker_R1}" ]] || ln -s {input.R1} {params.mw_trim_linker_R1}
        [[ -f "{params.mw_trim_linker_R2}" ]] || ln -s {input.R2} {params.mw_trim_linker_R2}

        # metawrap binning
        mw binning \
        --metabat2 --maxbin2 --concoct \
        -m {params.bin_mem} \
        -t {threads} \
        -a {input.assembly} \
        -o {params.bin_dir} \
        {params.mw_trim_linker_R1} {params.mw_trim_linker_R2}

        # metawrap bin refinement
        mw bin_refinement \
        -o {params.bin_dir} \
        -t {threads} \
        -A {params.bin_dir}/metabat2_bins \
        -B {params.bin_dir}/maxbin2_bins \
        -c {params.min_perc_complete} \
        -x {params.max_perc_contam}
        """


rule bin_stats:
    input:
        maxbin_bins                 = directory(join(top_binning_dir, "{name}", "maxbin2_bins")),
        metabat2_bins               = directory(join(top_binning_dir, "{name}", "metabat2_bins")),
        metawrap_bins               = directory(join(top_binning_dir, "{name}", "metawrap_50_5_bins")),
        maxbin_stats                = join(top_binning_dir, "{name}", "maxbin2_bins.stats"),
        metabat2_stats              = join(top_binning_dir, "{name}", "metabat2_bins.stats"),
        metawrap_stats              = join(top_binning_dir, "{name}", "metawrap_50_5_bins.stats"),
    output:
        this_refine_summary         = join(top_refine_dir, "{name}", "RefinedBins_summmary.txt"),
        named_stats_maxbin2         = join(top_refine_dir, "{name}", "named_maxbin2_bins.stats"),
        named_stats_metabat2        = join(top_refine_dir, "{name}", "named_metabat2_bins.stats"),
        named_stats_metawrap        = join(top_refine_dir, "{name}", "named_metawrap_bins.stats"),
    params:
        sid                         = "{name}",
        this_bin_dir                = join(top_refine_dir, "{name}"),
        # count number of fasta files in the bin folders to get count of bins
        metabat2_num_bins           = lambda wc, _out, _in: str(len([fn for fn in os.listdir(_in.metabat2_bins) if "unbinned" not in fn.lower()])),
        maxbin_num_bins             = lambda wc, _out, _in: str(len([fn for fn in os.listdir(_in.maxbin_bins) if "unbinned" not in fn.lower()])),
        metawrap_num_bins           = lambda wc, _out, _in: str(len([fn for fn in os.listdir(_in.metawrap_bins) if "unbinned" not in fn.lower()])),
    shell:
        """
        # count cumulative lines starting with `>`
        metabat2_contigs=$(cat {input.metabat2_bins}/*fa | grep -c "^>")
        maxbin_configs=$(cat {input.maxbin_bins}/*fa | grep -c "^>")
        metawrap_contigs=$(cat {input.metawrap_bins}/*fa | grep -c "^>")
        echo "SampleID\tmetabat2_bins\tmaxbin2_bins\tmetaWRAP_50_5_bins\tmetabat2_contigs\tmaxbin2_contigs\tmetaWRAP_50_5_contigs" > {output.this_refine_summary}
        echo "{params.sid}\t{params.metabat2_num_bins}\t{params.maxbin_num_bins}\t{params.metawrap_num_bins}\t$metabat2_contigs\t$maxbin_configs\t$metawrap_contigs"

        # name contigs with SID
        cat {input.maxbin_stats} | sed 's/^bin./{name}_bin./g' > {output.named_stats_maxbin2}
		cat {input.metabat2_stats} | sed 's/^bin./{name}_bin./g' > {output.named_stats_metabat2}
		cat {input.metawrap_stats} | sed 's/^bin./{name}_bin./g' > {output.named_stats_metawrap}
        """


rule cumulative_bin_stats:
    input:
        this_refine_summary         = expand(join(top_binning_dir, "{name}", "RefinedBins_summmary.txt"), name=samples),
    output:
        cumulative_bin_summary      = join(top_refine_dir, "{name}", "RefinedBins_summmary.txt"),
        cumulative_maxbin_stats     = join(top_refine_dir, "{name}", "cumulative_stats_maxbin.txt"),
        cumulative_metabat2_stats   = join(top_refine_dir, "{name}", "cumulative_stats_metabat2.txt"),
        cumulative_metawrap_stats   = join(top_refine_dir, "{name}", "cumulative_stats_metawrap.txt"),
    params:
        bin_dir                     = top_binning_dir
    shell:
        """
        # generate cumulative binning summary
        echo "SampleID\tmetabat2_bins\tmaxbin2_bins\tmetaWRAP_50_5_bins\tmetabat2_contigs\tmaxbin2_contigs\tmetaWRAP_50_5_contigs" > {output.cumulative_bin_summary}
        for report in {params.bin_dir}/*/RefinedBins_summmary.txt; do 
            tail -n+2 $report >> {output.cumulative_bin_summary}
            echo "tail -n+2 $report >> {output.cumulative_bin_summary}"
        done

        # generate cumulative statistic report for binning
        for report in {params.bin_dir}/*/named_maxbin2_bins.stats; do
            cat $report >> {output.cumulative_maxbin_stats}
        done
        for report in {params.bin_dir}/*/named_maxbin2_bins.stats; do
            cat $report >> {output.cumulative_metabat2_stats}
        done
        for report in {params.bin_dir}/*/named_maxbin2_bins.stats; do
            cat $report >> {output.cumulative_metawrap_stats}
        done

        if $(wc -l < {output.cumulative_bin_summary}) -le 1; then
            echo "No information in bin refinement summary!" 1>&2
            exit 1
        fi
        """


rule derep_bins:
    """
        dRep is a step which further refines draft-quality genomes (bins) by using a 
        fast, inaccurate estimation of genome distance, and a slow, accurate measure 
        of average nucleotide identity.

        @Input:
            maxbin2 ssembly bins, contigs, stat summaries
            metabat2 ssembly bins, contigs, stat summaries
            metawrap ssembly bins, contigs, stat summaries

        @Output:
            directory of consensus ensemble bins (deterministic output)

    """
    input:
        maxbin_contigs              = join(top_binning_dir, "{name}", "maxbin2_bins.contigs"),
        maxbin_stats                = join(top_binning_dir, "{name}", "maxbin2_bins.stats"),
        metabat2_contigs            = join(top_binning_dir, "{name}", "metabat2_bins.contigs"),
        metabat2_stats              = join(top_binning_dir, "{name}", "metabat2_bins.stats"),
        metawrap_contigs            = join(top_binning_dir, "{name}", "metawrap_50_5_bins.contigs"),
        metawrap_stats              = join(top_binning_dir, "{name}", "metawrap_50_5_bins.stats"),
    output:
        dereplicated_bins           = directory(join(top_refine_dir, "{name}", "dereplicated_bins")),
    singularity: metawrap_container,
    threads: int(cluster["derep_bins"].get("threads", default_threads)),
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
        mkdir -p {output.dereplicated_bins}
        dRep dereplicate \
        -g $(ls {params.metawrap_bins}/* | tr '\\n' ' ') \
        -p {threads} \
        -l {params.minimum_genome_length} \
        -pa {params.ani_primary_threshold} \
        -sa {params.ani_secondary_threshold} \
        -nc {params.min_overlap} \
        -cm {params.coverage_method} \
        {output.dereplicated_bins}
        """
