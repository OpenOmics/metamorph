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
samples                    = config["samples"]
top_log_dir                = join(workpath, "logfiles")
top_readqc_dir             = join(workpath, config['project']['id'], "metawrap_read_qc_skipBmtagger")
top_trim_dir               = join(workpath, config['project']['id'], "trimmed_reads")
top_assembly_dir           = join(workpath, config['project']['id'], "metawrap_assembly")
top_tax_dir                = join(workpath, config['project']['id'], "metawrap_kmer")
top_binning_dir            = join(workpath, config['project']['id'], "metawrap_binning")
top_refine_dir             = join(workpath, config['project']['id'], "metawrap_bin_refine")
top_mags_dir               = join(workpath, config['project']['id'], "mags")
top_mapping_dir            = join(workpath, config['project']['id'], "humann3_dna")
top_centrifuger_dir        = join(workpath, config['project']['id'], "centrifuger_dna")
top_resources_dir          = join(workpath, config['project']['id'], "resources")
top_qc_summary_dir         = join(workpath, config['project']['id'], "qc_summary")

# Generic, sample-agnostic reconciliation thresholds for dna_sample_qc_summary
# (Section below) -- not tuned to any particular organism or dataset. A
# sample below this fraction of its raw read pairs still present after
# dehosting+trimming is flagged regardless of why; a sample above this host-
# aligned fraction is flagged regardless of what the host reference was.
qc_min_overall_retained_pct = 50.0
qc_max_host_pct              = 20.0

# user-injected additional host/contaminant genome(s) for dehosting (--host-genome);
# empty list (the default) means dehosting behaves exactly as it always has
host_genome_paths          = config['options'].get('host_genome', [])
injected_host_idx_dir      = join(top_resources_dir, "injected_host_genome")
injected_host_idx_prefix   = join(injected_host_idx_dir, "injected_host")
# bowtie2_dehost's own output is the pipeline's final dehosted fastq only when
# there's no second screening pass to chain after it; otherwise it's an
# intermediate that screen_injected_host_dna consumes and finalizes
dna_dehost_tag             = "_dehost" if not host_genome_paths else "_dehost_prehostscreen"

# workflow flags
metawrap_container          = config["containers"]["metawrap"]
bowtie2_samtools_container  = config["images"]["bowtie2_samtools"]
metaphlan_utils_container   = config["containers"]["metaphlan_utils"]
centrifuger_sylph_container = config["containers"]["centrifuger_sylph"]
pairedness                  = list(range(1, config['project']['nends']+1))
mem2int                     = lambda x: int(str(x).lower().replace('gb', '').replace('g', ''))
assemblers                  = config["options"]["assemblers"]
use_megahit                 = 'megahit' in assemblers
use_metaspades              = 'metaspades' in assemblers
assembler_flags             = {
                                ('metaspades',): "--metaspades ",
                                ('megahit', 'metaspades'):  "--megahit --metaspades ",
                                ('megahit',): "--megahit "
                              }
"""
    Step-wise pipeline outline:
        1. read qc
        2. assembly
        3. binning
        4. bin refinement
        5. depreplicate bins
        6. annotate bins
        7. index depreplicated genomes
        8. align DNA to assembly
"""

rule metawrap_read_qc_skipBmtagger:
    """
    The metaWRAP::Read_qc module is meant to pre-process raw Illumina sequencing reads in preparation for assembly and alignment. 
    The raw reads are trimmed based on adapted content and PHRED scored with the default setting of Trim-galore, ensuring that only high-quality 
    sequences are left. Finally, FASTQC is used to generate quality reports of the raw and final read sets in order to assess read quality improvement. 
    The users have control over which of the above features they wish to use.

    Quality-control step accomplishes two tasks:
        - Trims adapter sequences from input read pairs (.fastq.gz)
        - Quality filters from input read pairs


    @Container requirements
        - [MetaWrap](https://github.com/bxlab/metaWRAP) orchestrates execution of 
            these tools in a psuedo-pipeline:
            - [FastQC](https://github.com/s-andrews/FastQC)
            - [Trim-Galore](https://github.com/FelixKrueger/TrimGalore)
            - [Cutadapt](https://github.com/marcelm/cutadapt)
    
    @Environment specifications
        - minimum 16gb memory


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
        R1_trimmed_gz       = join(top_trim_dir, "{name}", "{name}_R1_trimmed.fastq.gz"),
        R2_trimmed_gz       = join(top_trim_dir, "{name}", "{name}_R2_trimmed.fastq.gz"),
    params:
        rname               = "metawrap_read_qc_skipBmtagger",
        sid                 = "{name}",
        this_qc_dir         = join(top_readqc_dir, "{name}"),
        trim_out            = join(top_trim_dir, "{name}"),
        tmp_safe_dir        = join(config['options']['tmp_dir'], 'read_qc', "{name}"),
        tmpr1               = lambda _, output, input: join(config['options']['tmp_dir'], 'read_qc', _.name, str(basename(str(input.R1))).replace('_R1.', '_1.').replace('.gz', '')),
        tmpr2               = lambda _, output, input: join(config['options']['tmp_dir'], 'read_qc', _.name, str(basename(str(input.R2))).replace('_R2.', '_2.').replace('.gz', '')),
    containerized: metawrap_container,
    threads: int(cluster["metawrap_read_qc_skipBmtagger"].get('threads', default_threads)),
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

            # read quality control
            mw read_qc -1 {params.tmpr1} -2 {params.tmpr2} -t {threads} -o {params.this_qc_dir} --skip-bmtagger

            # collate fastq outputs to facilitate workflow, compress
            pigz -9 -p {threads} -c {params.this_qc_dir}/final_pure_reads_1.fastq  > {params.trim_out}/{params.sid}_R1_trimmed.fastq.gz
            pigz -9 -p {threads} -c {params.this_qc_dir}/final_pure_reads_2.fastq > {params.trim_out}/{params.sid}_R2_trimmed.fastq.gz

            # remove uncompressed trimmed fastqs now that gzipped copies exist
            rm -f {params.this_qc_dir}/final_pure_reads_1.fastq {params.this_qc_dir}/final_pure_reads_2.fastq

            # collate outputs to facilitate
            mkdir -p {params.this_qc_dir}/dna
            ln -s {params.this_qc_dir}/post-QC_report/final_pure_reads_1_fastqc.html {params.this_qc_dir}/{params.sid}_R1_postrim_report.html
            ln -s {params.this_qc_dir}/post-QC_report/final_pure_reads_2_fastqc.html {params.this_qc_dir}/{params.sid}_R2_postrim_report.html
            ln -s {params.this_qc_dir}/pre-QC_report/{params.sid}_1_fastqc.html {params.this_qc_dir}/{params.sid}_R1_pretrim_report.html
            ln -s {params.this_qc_dir}/pre-QC_report/{params.sid}_2_fastqc.html {params.this_qc_dir}/{params.sid}_R2_pretrim_report.html
        """

rule bowtie2_dehost:
    """
        This module utilizes bowtie2 to map trimmed reads to hg38, and select the untrimmed reads from the bam file as clean read.
    """
    input:
        R1                          = join(top_trim_dir, "{name}", "{name}_R1_trimmed.fastq.gz"),
        R2                          = join(top_trim_dir, "{name}", "{name}_R2_trimmed.fastq.gz"),
    output:
        R1_dehost                   = join(top_trim_dir, "{name}", "{name}_R1" + dna_dehost_tag + ".fastq.gz"),
        R2_dehost                   = join(top_trim_dir, "{name}", "{name}_R2" + dna_dehost_tag + ".fastq.gz"),
        dehost_summary              = join(top_trim_dir, "{name}", "{name}_hg38_dehost_summary.txt"),
    params:
        rname                       = "bowtie2_dehost",
        sid                         = "{name}",
        tmp_safe_dir                = join(config['options']['tmp_dir'], 'read_qc', "{name}"),
        bowtie2sam                  = join(config['options']['tmp_dir'], 'read_qc', "{name}", "{name}.mapped_and_unmapped.sam"),
        bowtie2bam                  = join(config['options']['tmp_dir'], 'read_qc', "{name}", "{name}.mapped_and_unmapped.bam"),
        bowtie2bamUnmapped          = join(config['options']['tmp_dir'], 'read_qc', "{name}", "{name}.unmapped.bam"),
        bowtie2bamUnmappedSorted    = join(config['options']['tmp_dir'], 'read_qc', "{name}", "{name}.unmappedSorted.bam"),
        hg38_bowtie2_idx_prefix     = "/data2/bowtie2/hg38",
    containerized: bowtie2_samtools_container,
    threads: int(cluster["bowtie2_dehost"].get('threads', default_threads)),
    shell:
        """
            # safe temp directory
            if [ ! -d "{params.tmp_safe_dir}" ]; then mkdir -p "{params.tmp_safe_dir}"; fi
            tmp=$(mktemp -d -p "{params.tmp_safe_dir}")
            trap 'rm -rf "{params.tmp_safe_dir}"' EXIT

            #run bowtie to map reads to hg38
            bowtie2 -p {threads} -x {params.hg38_bowtie2_idx_prefix} \
            -1 {input.R1} \
            -2 {input.R2} \
            -S {params.bowtie2sam}

            # run samtools to convert sam to bam
            samtools view -@ {threads} -bS {params.bowtie2sam} > {params.bowtie2bam}

            # run samtools to filter out the reads that none of the pair mapped to hg38
            samtools view -@ {threads} -b -f 12 -F 256 {params.bowtie2bam} > {params.bowtie2bamUnmapped}

            # split paired-end reads into separated fastq files .._R1 .._R2
            ## sort bam file by read name ( -n ) to have paired reads next to each other
            samtools sort -@ {threads} -n {params.bowtie2bamUnmapped} -o {params.bowtie2bamUnmappedSorted}
            ## Convert to fastq
            samtools fastq -@ {threads} {params.bowtie2bamUnmappedSorted} \
            -1 {output.R1_dehost} \
            -2 {output.R2_dehost} \
            -0 /dev/null -s /dev/null -n

            # summarize the filter: read-pair counts straight from the BAMs
            # already produced above, rather than re-decompressing fastqs
            total_pairs=$(( $(samtools view -c -@ {threads} {params.bowtie2bam}) / 2 ))
            retained_pairs=$(( $(samtools view -c -@ {threads} {params.bowtie2bamUnmappedSorted}) / 2 ))
            removed_pairs=$(( total_pairs - retained_pairs ))
            pct_removed=$(awk -v r="$removed_pairs" -v t="$total_pairs" 'BEGIN{{ printf "%.4f", (t>0)? r/t*100 : 0 }}')
            pct_retained=$(awk -v r="$retained_pairs" -v t="$total_pairs" 'BEGIN{{ printf "%.4f", (t>0)? r/t*100 : 0 }}')
            {{
                echo "metamorph dehosting summary"
                echo "============================"
                echo "Sample:                  {params.sid}"
                echo "Filter stage:            default host reference (hg38)"
                echo "Reference index:         {params.hg38_bowtie2_idx_prefix}"
                echo "Program:                 bowtie2 + samtools"
                echo "Retention criterion:     read pair kept only if BOTH mates are unmapped to this reference (samtools view -f 12 -F 256)"
                echo ""
                echo "Input read pairs:        ${{total_pairs}}"
                echo "Removed (mapped):        ${{removed_pairs}} (${{pct_removed}}%)"
                echo "Retained (unmapped):     ${{retained_pairs}} (${{pct_retained}}%)"
            }} > {output.dehost_summary}
        """


if host_genome_paths:

    rule build_injected_host_index:
        """
            Builds one combined bowtie2 index from every FASTA passed via
            --host-genome, so the dehosting steps can screen reads against it
            in addition to the default host reference. Only exists in the DAG
            when --host-genome was actually given.
        """
        input:
            genomes                 = host_genome_paths,
        output:
            index_complete           = touch(join(injected_host_idx_dir, "INJECTED_HOST_INDEX_COMPLETE")),
        params:
            rname                    = "build_injected_host_index",
            idx_dir                  = injected_host_idx_dir,
            idx_prefix               = injected_host_idx_prefix,
            combined_fa              = join(injected_host_idx_dir, "injected_host_combined.fa"),
            genome_files              = host_genome_paths,
        containerized: bowtie2_samtools_container,
        threads: int(cluster["build_injected_host_index"].get('threads', default_threads)),
        shell:
            """
                mkdir -p {params.idx_dir}

                # Any one (or combination) of the injected FASTA(s) may reuse
                # generic, assembler-default contig names (e.g. "ptg000001l")
                # that collide with each other -- seen directly with a
                # multi-genome FASTA built from two separate assemblies that
                # both used the same naming convention. bowtie2-build itself
                # doesn't care, but a SAM/BAM header requires globally unique
                # reference names, so samtools fails downstream otherwise.
                # Rewrite every sequence to a guaranteed-unique ID up front,
                # keeping the original name as a suffix for traceability.
                for g in {params.genome_files}; do
                    case "$g" in
                        *.gz) zcat "$g" ;;
                        *)    cat  "$g" ;;
                    esac
                done | awk '/^>/ {{n++; print ">injected_seq" n "_" substr($0,2); next}} {{print}}' > {params.combined_fa}

                bowtie2-build --threads {threads} {params.combined_fa} {params.idx_prefix}
            """


    rule screen_injected_host_dna:
        """
            Second dehosting pass: screens reads already cleared of the
            default host reference (bowtie2_dehost) against every genome
            supplied via --host-genome, using the same both-mates-unmapped
            criterion. Only exists in the DAG when --host-genome was given;
            this rule's output is the pipeline's real, final dehosted fastq
            in that case.
        """
        input:
            R1                      = join(top_trim_dir, "{name}", "{name}_R1_dehost_prehostscreen.fastq.gz"),
            R2                      = join(top_trim_dir, "{name}", "{name}_R2_dehost_prehostscreen.fastq.gz"),
            index_complete           = join(injected_host_idx_dir, "INJECTED_HOST_INDEX_COMPLETE"),
        output:
            R1_dehost                = join(top_trim_dir, "{name}", "{name}_R1_dehost.fastq.gz"),
            R2_dehost                = join(top_trim_dir, "{name}", "{name}_R2_dehost.fastq.gz"),
            dehost_summary           = join(top_trim_dir, "{name}", "{name}_injected_host_dehost_summary.txt"),
        params:
            rname                    = "screen_injected_host_dna",
            sid                      = "{name}",
            tmp_safe_dir             = join(config['options']['tmp_dir'], 'host_screen', "{name}"),
            bowtie2sam               = join(config['options']['tmp_dir'], 'host_screen', "{name}", "{name}.injected_host.sam"),
            bowtie2bam               = join(config['options']['tmp_dir'], 'host_screen', "{name}", "{name}.injected_host.bam"),
            bowtie2bamUnmapped       = join(config['options']['tmp_dir'], 'host_screen', "{name}", "{name}.injected_host.unmapped.bam"),
            bowtie2bamUnmappedSorted = join(config['options']['tmp_dir'], 'host_screen', "{name}", "{name}.injected_host.unmappedSorted.bam"),
            injected_host_idx_prefix = injected_host_idx_prefix,
            injected_genomes         = " ".join(host_genome_paths),
        containerized: bowtie2_samtools_container,
        threads: int(cluster["screen_injected_host_dna"].get('threads', default_threads)),
        shell:
            """
                if [ ! -d "{params.tmp_safe_dir}" ]; then mkdir -p "{params.tmp_safe_dir}"; fi
                tmp=$(mktemp -d -p "{params.tmp_safe_dir}")
                trap 'rm -rf "{params.tmp_safe_dir}"' EXIT

                # run bowtie against every user-injected host/contaminant genome
                bowtie2 -p {threads} -x {params.injected_host_idx_prefix} \
                -1 {input.R1} \
                -2 {input.R2} \
                -S {params.bowtie2sam}

                samtools view -@ {threads} -bS {params.bowtie2sam} > {params.bowtie2bam}

                # keep only pairs where neither mate mapped to any injected genome
                samtools view -@ {threads} -b -f 12 -F 256 {params.bowtie2bam} > {params.bowtie2bamUnmapped}

                samtools sort -@ {threads} -n {params.bowtie2bamUnmapped} -o {params.bowtie2bamUnmappedSorted}
                samtools fastq -@ {threads} {params.bowtie2bamUnmappedSorted} \
                -1 {output.R1_dehost} \
                -2 {output.R2_dehost} \
                -0 /dev/null -s /dev/null -n

                # summarize the filter: read-pair counts straight from the BAMs
                # already produced above, rather than re-decompressing fastqs
                total_pairs=$(( $(samtools view -c -@ {threads} {params.bowtie2bam}) / 2 ))
                retained_pairs=$(( $(samtools view -c -@ {threads} {params.bowtie2bamUnmappedSorted}) / 2 ))
                removed_pairs=$(( total_pairs - retained_pairs ))
                pct_removed=$(awk -v r="$removed_pairs" -v t="$total_pairs" 'BEGIN{{ printf "%.4f", (t>0)? r/t*100 : 0 }}')
                pct_retained=$(awk -v r="$retained_pairs" -v t="$total_pairs" 'BEGIN{{ printf "%.4f", (t>0)? r/t*100 : 0 }}')
                {{
                    echo "metamorph dehosting summary"
                    echo "============================"
                    echo "Sample:                  {params.sid}"
                    echo "Filter stage:            user-injected host/contaminant genome(s) (--host-genome)"
                    echo "Injected genome FASTA(s): {params.injected_genomes}"
                    echo "Combined index:          {params.injected_host_idx_prefix}"
                    echo "Program:                 bowtie2 + samtools"
                    echo "Retention criterion:     read pair kept only if BOTH mates are unmapped to this combined reference (samtools view -f 12 -F 256)"
                    echo "Input to this stage:     output of the default hg38 screen (bowtie2_dehost); see {params.sid}_hg38_dehost_summary.txt"
                    echo ""
                    echo "Input read pairs:        ${{total_pairs}}"
                    echo "Removed (mapped):        ${{removed_pairs}} (${{pct_removed}}%)"
                    echo "Retained (unmapped):     ${{retained_pairs}} (${{pct_retained}}%)"
                }} > {output.dehost_summary}
            """


rule dna_centrifuger:
    """
    Data processing step to classify taxonomic composition of the host removed reads.
    For more information about centrifuger, please visit: 
       - https://github.com/mourisl/centrifuger
    @Environment specifications
        - Minimum 200 GB memory, centrifuger index is around 167 GB
    @Input:
        - Dehost fastq files (scatter-per-sample)
    @Outputs:
        - Centrifuge classification file
    """
    input:
        R1_dehost                   = join(top_trim_dir, "{name}", "{name}_R1_dehost.fastq.gz"),
        R2_dehost                   = join(top_trim_dir, "{name}", "{name}_R2_dehost.fastq.gz"),
    output:
        classification              = join(top_centrifuger_dir, "{name}_centrifuger_classification.tsv"),
        centrifuger_quant           = join(top_centrifuger_dir, "{name}_centrifuger_quantification_report.tsv"),
        metaphlan_quant             = join(top_centrifuger_dir, "{name}_centrifuger_quantification_report_mpl.tsv"),
    params:
        rname                       = "dna_centrifuger",
        sid                         = "{name}",
        centrifuger_idx_prefix      = "/data2/centrifuger/gtdb_r226+refseq_hvfc/cfr_gtdb_r226+refseq_hvfc",
    container: centrifuger_sylph_container,
    threads: int(cluster["dna_centrifuger"].get('threads', default_threads)),
    shell: 
        """
        # Runs centrifuger on gDNA to classify taxonomic
        # composition of host removed reads
        centrifuger \\
            -t {threads} \\
            -x {params.centrifuger_idx_prefix} \\
            -1 {input.R1_dehost} \\
            -2 {input.R2_dehost}  \\
        > {output.classification}
        # Runs centrifuger quant on gDNA to quantify
        # taxonomic composition of host removed reads
        # Centrifuger output file format
        centrifuger-quant \\
            -x {params.centrifuger_idx_prefix} \\
            -c {output.classification} \\
            --output-format 0 \\
        > {output.centrifuger_quant}
        # Metaphlan output file format
        centrifuger-quant \\
            -x {params.centrifuger_idx_prefix} \\
            -c {output.classification} \\
            --output-format 1 \\
        > {output.metaphlan_quant}
        """

rule dna_sample_qc_summary:
    """
    Joins measurements already produced by the read-qc, dehosting, and
    classification rules above into one QC/provenance row per sample:
    raw/trimmed read counts and FastQC module statuses, dehosting pair
    counts and host-alignment rate, and Centrifuger's classified fraction
    and domain-level breakdown -- plus reconciliation flags for the ways
    this join could silently be wrong (R1/R2 desync, a count that rose
    instead of fell between stages, the classifier's actual input not
    matching what dehosting reported producing, a missing/empty file,
    unusually high host content, or unusually low overall retention).
    See workflow/scripts/sample_qc_summary.py for the actual logic; this
    rule only locates the right files and counts the two numbers that
    script deliberately doesn't trust from any single upstream file
    (the dehosted FASTQs' own R1/R2 pair counts).
    """
    input:
        fastqc_pretrim_r1              = join(top_readqc_dir, "{name}", "{name}_R1_pretrim_report.html"),
        fastqc_pretrim_r2               = join(top_readqc_dir, "{name}", "{name}_R2_pretrim_report.html"),
        fastqc_posttrim_r1              = join(top_readqc_dir, "{name}", "{name}_R1_postrim_report.html"),
        fastqc_posttrim_r2               = join(top_readqc_dir, "{name}", "{name}_R2_postrim_report.html"),
        hg38_dehost_summary             = join(top_trim_dir, "{name}", "{name}_hg38_dehost_summary.txt"),
        injected_host_dehost_summary    = (join(top_trim_dir, "{name}", "{name}_injected_host_dehost_summary.txt") if host_genome_paths else []),
        R1_dehost                       = join(top_trim_dir, "{name}", "{name}_R1_dehost.fastq.gz"),
        R2_dehost                       = join(top_trim_dir, "{name}", "{name}_R2_dehost.fastq.gz"),
        classification                  = join(top_centrifuger_dir, "{name}_centrifuger_classification.tsv"),
        centrifuger_quant               = join(top_centrifuger_dir, "{name}_centrifuger_quantification_report.tsv"),
    output:
        qc_summary                      = join(top_qc_summary_dir, "{name}_qc_summary.tsv"),
    params:
        rname                           = "dna_sample_qc_summary",
        sid                             = "{name}",
        workflow_mode                   = config["options"]["workflow"],
        pipeline_version                = config["project"]["version"],
        pipeline_git_commit             = config["project"]["git_commit_hash"],
        max_host_pct                    = qc_max_host_pct,
        min_overall_retained_pct        = qc_min_overall_retained_pct,
        injected_summary_arg            = lambda wc: (
            "--injected-host-dehost-summary " + join(top_trim_dir, wc.name, wc.name + "_injected_host_dehost_summary.txt")
            if host_genome_paths else ""
        ),
    containerized: metawrap_container,
    threads: int(cluster["dna_sample_qc_summary"].get("threads", default_threads)),
    shell:
        """
        r1_dehost_count=$(( $(zcat {input.R1_dehost} | wc -l) / 4 ))
        r2_dehost_count=$(( $(zcat {input.R2_dehost} | wc -l) / 4 ))

        python3 workflow/scripts/sample_qc_summary.py \\
            --sample-id {params.sid} \\
            --library-type DNA \\
            --workflow {params.workflow_mode} \\
            --pipeline-version {params.pipeline_version} \\
            --pipeline-git-commit {params.pipeline_git_commit} \\
            --fastqc-pretrim-r1 {input.fastqc_pretrim_r1} \\
            --fastqc-pretrim-r2 {input.fastqc_pretrim_r2} \\
            --fastqc-posttrim-r1 {input.fastqc_posttrim_r1} \\
            --fastqc-posttrim-r2 {input.fastqc_posttrim_r2} \\
            --hg38-dehost-summary {input.hg38_dehost_summary} \\
            {params.injected_summary_arg} \\
            --classification-tsv {input.classification} \\
            --quant-report {input.centrifuger_quant} \\
            --dehosted-r1-count ${{r1_dehost_count}} \\
            --dehosted-r2-count ${{r2_dehost_count}} \\
            --max-host-pct {params.max_host_pct} \\
            --min-overall-retained-pct {params.min_overall_retained_pct} \\
            --output {output.qc_summary}
        """

rule dna_qc_summary_merge:
    """
    Concatenates every sample's dna_sample_qc_summary row into one
    cohort-level TSV with a single header -- the per-sample files stay
    around too (same pattern as the per-sample vs. merged bugs-list
    tables elsewhere in this workflow), this just adds the one-file view.
    """
    input:
        expand(join(top_qc_summary_dir, "{name}_qc_summary.tsv"), name=samples),
    output:
        qc_summary_merged                = join(workpath, config['project']['id'], "qc_summary.tsv"),
    params:
        rname                            = "dna_qc_summary_merge",
    containerized: metawrap_container,
    threads: int(cluster["dna_qc_summary_merge"].get("threads", default_threads)),
    shell:
        """
        head -n 1 {input[0]} > {output.qc_summary_merged}
        for f in {input}; do
            tail -n +2 "$f" >> {output.qc_summary_merged}
        done
        """

rule dna_humann_classify:
    input:
        R1                  = join(top_trim_dir, "{name}", "{name}_R1_dehost.fastq.gz"),
        R2                  = join(top_trim_dir, "{name}", "{name}_R2_dehost.fastq.gz"),
    output:
        hm3_gene_fam        = join(top_mapping_dir, '{name}_genefamilies.tsv'),
        hm3_path_abd        = join(top_mapping_dir, '{name}_pathabundance.tsv'),
        hm3_path_cov        = join(top_mapping_dir, '{name}_pathcoverage.tsv'),
        mpl4_bugs_list      = join(top_mapping_dir, '{name}_bugs_list.tsv'),
        humann_log          = join(top_mapping_dir, '{name}_humann3.log'),
        humann_config       = join(top_mapping_dir, '{name}_humann3.conf'),
    params:
        rname               = "dna_humann_classify",
        sid                 = "{name}",
        tmp_safe_dir        = join(config['options']['tmp_dir'], 'dna_map'),
        tmpread             = join(config['options']['tmp_dir'], 'dna_map', "{name}_concat.fastq.gz"),
        hm3_map_dir         = top_mapping_dir,
        uniref_db           = "/data2/uniref",      # from <root>/config/resources.json
        chocophlan_db       = "/data2/chocophlan",  # from <root>/config/resources.json
        util_map_db         = "/data2/um",          # from <root>/config/resources.json
        metaphlan_db        = "/data2/metaphlan",   # from <root>/config/resources.json
        deep_mode           = "\\\n        --bypass-translated-search" if not config["options"]["shallow_profile"] else "",
    threads: int(cluster["dna_humann_classify"].get('threads', default_threads)),
    containerized: metawrap_container,
    shell:
        """
        . /opt/conda/etc/profile.d/conda.sh && conda activate bb3
        # safe temp directory
        if [ ! -d "{params.tmp_safe_dir}" ]; then mkdir -p "{params.tmp_safe_dir}"; fi
        tmp=$(mktemp -d -p "{params.tmp_safe_dir}")
        trap 'rm -rf "{params.tmp_safe_dir}"' EXIT

        # human/metaphlan configuration
        humann_config --print > {output.humann_config}
        export DEFAULT_DB_FOLDER={params.metaphlan_db}

        cat {input.R1} {input.R2} > {params.tmpread}
        humann \
        --threads {threads} \
        --input {params.tmpread} \
        --input-format fastq.gz {params.deep_mode} \
        --metaphlan-options "--bowtie2db {params.metaphlan_db} --nproc {threads} --offline" \
        --output-basename {params.sid} \
        --log-level DEBUG \
        --o-log {output.humann_log} \
        --output {params.hm3_map_dir}

        mv {params.hm3_map_dir}/{params.sid}_humann_temp/{params.sid}_metaphlan_bugs_list.tsv {output.mpl4_bugs_list}
        rm -r {params.hm3_map_dir}/{params.sid}_humann_temp
        """

rule dna_humann_diversity_calculation:
    """
    This step calculates the following diversity metrics from the merged buglist:
            1. observed species
            2. Shannon index
            3. Jaccard distance
            4. Bray-curtis distance
            5. Weighted and unweighted unifrac distance
    """
    input:
        bugs_list           = join(top_mapping_dir, 'merged_bugs_list.tsv'),
    output:
        shannon             = join(top_mapping_dir, 'diversity_analysis', 'merged_bugs_list_shannon.tsv'),
        richness            = join(top_mapping_dir, 'diversity_analysis', 'merged_bugs_list_richness.tsv'),
        jaccard             = join(top_mapping_dir, 'diversity_analysis', 'merged_bugs_list_jaccard.tsv'),
        wunif               = join(top_mapping_dir, 'diversity_analysis', 'merged_bugs_list_weighted-unifrac.tsv'),
        unwunif             = join(top_mapping_dir, 'diversity_analysis', 'merged_bugs_list_unweighted-unifrac.tsv'),
    params:
        rname               = "dna_humann_diversity",
        div_dir             = join("metagenome_results/humann3_dna", "diversity_analysis"),
        wunif_log           = join("metagenome_results/humann3_dna", "diversity_analysis", "merged_bugs_list_weighted-unifrac.log"),
        unwunif_log         = join("metagenome_results/humann3_dna", "diversity_analysis", "merged_bugs_list_unweighted-unifrac.log"),
        tree                = "/data2/mpa_nwk_tree/mpa_vJan21_CHOCOPhlAnSGB_202103.nwk",
    containerized: metaphlan_utils_container,
    threads: int(cluster["dna_humann_diversity_calculation"].get('threads', default_threads)),
    shell:
        """
        calculate_diversity.R -f {input.bugs_list} -o {params.div_dir} -d beta -m jaccard -s t__
     	calculate_diversity.R -f {input.bugs_list} -o {params.div_dir} -d alpha -m richness -s t__
    	calculate_diversity.R -f {input.bugs_list} -o {params.div_dir} -d alpha -m shannon -s t__
    	calculate_diversity.R -f {input.bugs_list} -o {params.div_dir} -t {params.tree} -m weighted-unifrac -s t__ > {params.wunif_log}
    	calculate_diversity.R -f {input.bugs_list} -o {params.div_dir} -t {params.tree} -m unweighted-unifrac -s t__ > {params.unwunif_log}
        """

rule dna_humann_summarize:
    """
        This step merges all per-sample tables as a single table for each type of input:
     	1. Merges the buglist files generated from metaphlan4 to a single table with merge_metaphlan_tables.py
        2.  
            a) calculate CPM for path abundance.tsv (with humann_renorm_table)
            b) separate the CPM table to one stratified and one unstratified table (with humann_split_stratified_table)
            c) merge per-sample stratified CPM tables as a single table with (humann_join_tables)
            d) do the same thing for the unstratified per-sample cpm tables
    """
    input:
        bugs_list           = expand(join(top_mapping_dir, '{name}_bugs_list.tsv'), name=samples),
        path_abd            = expand(join(top_mapping_dir, '{name}_pathabundance.tsv'), name=samples),
  
    output:
        mer_bugs_list       = join(top_mapping_dir, 'merged_bugs_list.tsv'),
        mer_path_abd_str    = join(top_mapping_dir, 'merged_pathabundance.cpm.stratified.tsv'),
        mer_path_abd_unstr  = join(top_mapping_dir, 'merged_pathabundance.cpm.unstratified.tsv'),
    params:
        rname               = "dna_humann_summarize",
        hm3_map_dir         = top_mapping_dir,
        stratified_dir      = join(top_mapping_dir,'stratified'),
        unstratified_dir    = join(top_mapping_dir,'unstratified'),
    containerized: metawrap_container,
    threads: int(cluster["dna_humann_summarize"].get('threads', default_threads)),
    shell:
        """
        . /opt/conda/etc/profile.d/conda.sh && conda activate bb3

        # merge bugs list
        merge_metaphlan_tables.py {input.bugs_list} > {output.mer_bugs_list}

        # STEP 1: convert RPK to CPM
        for file in {input.path_abd}; do 
            humann_renorm_table \
            --input ${{file}} \
            --output ${{file%.*}}-cpm.tsv \
            --units cpm --update-snames; 
        done
        # STEP 2
        # split stratified table
        for file in {params.hm3_map_dir}/*pathabundance-cpm.tsv; do 
            humann_split_stratified_table \
            --input ${{file}} \
            --output {params.hm3_map_dir}; 
        done


        # STEP 3: place the cpm tables into the dirs they belong to
        mkdir -p {params.stratified_dir}
        mkdir -p {params.unstratified_dir}
        mv {params.hm3_map_dir}/*_pathabundance-cpm_stratified.tsv {params.stratified_dir}
        mv {params.hm3_map_dir}/*_pathabundance-cpm_unstratified.tsv {params.unstratified_dir}

        # STEP 4:join all the tables of the same type
        humann_join_tables \
        --input {params.stratified_dir} \
        --output {output.mer_path_abd_str}
        humann_join_tables \
        --input {params.unstratified_dir} \
        --output {output.mer_path_abd_unstr}

        """

rule dna_decompress_dehost_reads:
    """ 
       All steps in the assembly option involves using the decompressed dehost reads, so this step decompressed teh dehost reads to make them ready to use.
    """
    input:
        R1                  = join(top_trim_dir, "{name}", "{name}_R1_dehost.fastq.gz"),
        R2                  = join(top_trim_dir, "{name}", "{name}_R2_dehost.fastq.gz"),    
    output:
        R1                  = join(top_trim_dir, "{name}", "{name}_R1_dehost.fastq"),
        R2                  = join(top_trim_dir, "{name}", "{name}_R2_dehost.fastq"), 
    params:
        rname               = "dna_decompress_dehost_reads",
        sid                 = "{name}",
    threads: int(cluster["dna_decompress_dehost_reads"].get('threads', default_threads)),
    shell:
        """
        zcat {input.R1} > {output.R1}
        zcat {input.R2} > {output.R2}
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
            Dehosted reads generated from bowtie2_dehost (R1 & R2 per sample)

        @Output:
            Megahit assembled contigs and reports
            Metaspades assembled contigs and reports
            Ensemble assembled contigs and reports
    """
    input:
        R1                  = join(top_trim_dir, "{name}", "{name}_R1_dehost.fastq"),
        R2                  = join(top_trim_dir, "{name}", "{name}_R2_dehost.fastq"),    
    output:
        # megahit outputs
        megahit_assembly            = join(top_assembly_dir, "{name}", "megahit", "final.contigs.fa") if use_megahit else [],
        # megahit_longcontigs         = join(top_assembly_dir, "{name}", "megahit", "long.contigs.fa"),
        megahit_log                 = join(top_assembly_dir, "{name}", "megahit", "log") if use_megahit else [],
        # metaspades outputs
        metaspades_assembly         = join(top_assembly_dir, "{name}", "metaspades", "contigs.fasta") if use_metaspades else [],
        metaspades_graph            = join(top_assembly_dir, "{name}", "metaspades", "assembly_graph.fastg") if use_metaspades else [],
        metaspades_longscaffolds    = join(top_assembly_dir, "{name}", "metaspades", "long_scaffolds.fasta") if use_metaspades else [],
        metaspades_scaffolds        = join(top_assembly_dir, "{name}", "metaspades", "scaffolds.fasta") if use_metaspades else [],
        metaspades_cor_reads        = directory(join(top_assembly_dir, "{name}", "metaspades", "corrected")) if use_metaspades else [],
        # ensemble outputs
        final_assembly              = join(top_assembly_dir, "{name}", "final_assembly.fasta"),
        final_assembly_report       = join(top_assembly_dir, "{name}", "assembly_report.html"),
        assembly_R1                 = join(top_assembly_dir, "{name}", "{name}_1.fastq"),
        assembly_R2                 = join(top_assembly_dir, "{name}", "{name}_2.fastq"),
    singularity: metawrap_container,
    params:
        rname                       = "metawrap_genome_assembly",
        memlimit                    = mem2int(cluster["metawrap_genome_assembly"].get('mem', default_memory)),
        mh_dir                      = join(top_assembly_dir, "{name}", "megahit"),
        assembly_dir                = join(top_assembly_dir, "{name}"),
        contig_min_len              = "1000",
        assembler_flag              = assembler_flags[tuple(assemblers)]
    threads: int(cluster["metawrap_genome_assembly"].get('threads', default_threads)),
    shell:
        """
            if [ -d "{params.assembly_dir}" ]; then rm -rf "{params.assembly_dir}"; fi      
            mkdir -p "{params.assembly_dir}"
            if [ -d {params.mh_dir:q} ]; then rm -rf {params.mh_dir:q}; fi
            # link to the file names metawrap expects
            ln -s {input.R1} {output.assembly_R1}
            ln -s {input.R2} {output.assembly_R2}
            # run genome assembler
            mw assembly {params.assembler_flag}\\
                -m {params.memlimit} \\
                -t {threads} \\
                -l {params.contig_min_len} \\
                -1 {output.assembly_R1} \\
                -2 {output.assembly_R2} \\
                -o {params.assembly_dir}
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
            - metaBAT2cd
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
            Dehosted reads generated from bowtie2_dehost (R1 & R2 per sample)

        @Output:
            Megahit assembled draft-genomes (bins) and reports
            Metaspades assembled draft-genomes (bins) and reports
            Ensemble assembled draft-genomes (bins) and reports

    """
    input:
        R1                          = join(top_trim_dir, "{name}", "{name}_R1_dehost.fastq"),
        R2                          = join(top_trim_dir, "{name}", "{name}_R2_dehost.fastq"),   
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
        bin_breadcrumb              = join(top_binning_dir, "{name}", "{name}_BINNING_COMPLETE"),
        bin_figure                  = join(top_binning_dir, "{name}", "figures", "binning_results.png"),
    params:
        rname                       = "metawrap_binning",
        bin_parent_dir              = top_binning_dir,
        bin_dir                     = join(top_binning_dir, "{name}"),
        bin_refine_dir              = join(top_binning_dir, "{name}"),
        bin_mem                     = mem2int(cluster['metawrap_binning'].get("mem", default_memory)),
        mw_trim_linker_R1           = join(top_trim_dir, "{name}", "{name}_1.fastq"),
        mw_trim_linker_R2           = join(top_trim_dir, "{name}", "{name}_2.fastq"),
        tmp_bin_dir                 = join(config['options']['tmp_dir'], 'bin', "{name}"),
        tmp_binref_dir              = join(config['options']['tmp_dir'], 'bin_rf', "{name}"),
        min_perc_complete           = "50",
        max_perc_contam             = "5",
    singularity: metawrap_container,
    threads: int(cluster["metawrap_binning"].get("threads", default_threads)),
    shell:
        """
        # safe temp directory
        if [ ! -d "{params.tmp_bin_dir}" ]; then mkdir -p "{params.tmp_bin_dir}"; fi
        if [ ! -d "{params.tmp_binref_dir}" ]; then mkdir -p "{params.tmp_binref_dir}"; fi
        tmp_bin=$(mktemp -d -p "{params.tmp_bin_dir}")
        tmp_bin_ref=$(mktemp -d -p "{params.tmp_binref_dir}")
        trap 'rm -rf "{params.tmp_bin_dir}" "{params.tmp_binref_dir}"' EXIT

        # set checkm data
        checkm data setRoot /data2/CHECKM_DB
        export CHECKM_DATA_PATH="/data2/CHECKM_DB"

        # make base dir if not exists
        mkdir -p {params.bin_parent_dir}
        if [ -d "{params.bin_dir}" ]; then rm -rf {params.bin_dir}; fi

        # setup links for metawrap input
        rm -f "{params.mw_trim_linker_R1}"
        ln -s {input.R1} {params.mw_trim_linker_R1}
        rm -f "{params.mw_trim_linker_R2}"
        ln -s {input.R2} {params.mw_trim_linker_R2}

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
            -o {params.tmp_binref_dir} \
            -m {params.bin_mem} \
            -t {threads} \
            -A {params.bin_dir}/metabat2_bins \
            -B {params.bin_dir}/maxbin2_bins \
            -c {params.min_perc_complete} \
            -x {params.max_perc_contam}
        
        cp -r {params.tmp_binref_dir}/* {params.bin_dir}
        bold=$(tput bold)
        normal=$(tput sgr0)
        for this_bin in `find {params.bin_dir} -type f -regex ".*/bin.[0-9][0-9]?[0-9]?.fa"`; do
            fn=$(basename $this_bin)
            to=$(dirname $this_bin)/{wildcards.name}_$fn
            echo "relabeling ${{bold}}$fn${{normal}} to ${{bold}}{wildcards.name}_$fn${{normal}}"
            mv $this_bin $to
        done
        
        touch {output.bin_breadcrumb}
        chmod -R 775 {params.bin_dir}
        """


rule bin_stats:
    input:
        maxbin_bins                 = join(top_binning_dir, "{name}", "maxbin2_bins"),
        metabat2_bins               = join(top_binning_dir, "{name}", "metabat2_bins"),
        metawrap_bins               = join(top_binning_dir, "{name}", "metawrap_50_5_bins"),
        maxbin_stats                = join(top_binning_dir, "{name}", "maxbin2_bins.stats"),
        metabat2_stats              = join(top_binning_dir, "{name}", "metabat2_bins.stats"),
        metawrap_stats              = join(top_binning_dir, "{name}", "metawrap_50_5_bins.stats"),
    output:
        this_refine_summary         = join(top_refine_dir, "{name}", "RefinedBins_summmary.txt"),
        named_stats_maxbin2         = join(top_refine_dir, "{name}", "named_maxbin2_bins.stats"),
        named_stats_metabat2        = join(top_refine_dir, "{name}", "named_metabat2_bins.stats"),
        named_stats_metawrap        = join(top_refine_dir, "{name}", "named_metawrap_bins.stats"),
    params:
        rname                       = "bin_stats",
        sid                         = "{name}",
        this_bin_dir                = join(top_refine_dir, "{name}"),
        # count number of fasta files in the bin folders to get count of bins
        metabat2_num_bins           = lambda wc, input: str(len([fn for fn in os.listdir(input.metabat2_bins) if "unbinned" not in fn.lower()])),
        maxbin_num_bins             = lambda wc, input: str(len([fn for fn in os.listdir(input.maxbin_bins) if "unbinned" not in fn.lower()])),
        metawrap_num_bins           = lambda wc, input: str(len([fn for fn in os.listdir(input.metawrap_bins) if "unbinned" not in fn.lower()])),
    shell:
        """
        metabat2_contigs=$(cat {input.metabat2_bins}/*fa | grep -c "^>")
        maxbin_configs=$(cat {input.maxbin_bins}/*fa | grep -c "^>")
        metawrap_contigs=$(cat {input.metawrap_bins}/*fa | grep -c "^>")
        echo "SampleID\tmetabat2_bins\tmaxbin2_bins\tmetaWRAP_50_5_bins\tmetabat2_contigs\tmaxbin2_contigs\tmetaWRAP_50_5_contigs" > {output.this_refine_summary}
        echo "{params.sid}\t{params.metabat2_num_bins}\t{params.maxbin_num_bins}\t{params.metawrap_num_bins}\t$metabat2_contigs\t$maxbin_configs\t$metawrap_contigs" >> {output.this_refine_summary}
        
        # name contigs with SID
        cat {input.maxbin_stats} | sed 's/^bin./{params.sid}_bin./g' > {output.named_stats_maxbin2}
        sed -i '1 s/^.*$/genome\tcompleteness\tcontamination\tGC\tlineage\tN50\tsize\tbinner/' {output.named_stats_maxbin2}
        cat {input.metabat2_stats} | sed 's/^bin./{params.sid}_bin./g' > {output.named_stats_metabat2}
        sed -i '1 s/^.*$/genome\tcompleteness\tcontamination\tGC\tlineage\tN50\tsize\tbinner/' {output.named_stats_metabat2}
        cat {input.metawrap_stats} | sed 's/^bin./{params.sid}_bin./g' > {output.named_stats_metawrap}
        sed -i '1 s/^.*$/genome\tcompleteness\tcontamination\tGC\tlineage\tN50\tsize\tbinner/' {output.named_stats_metawrap}
        """


rule cumulative_bin_stats:
    input:
        this_refine_summary         = expand(join(top_refine_dir, "{name}", "RefinedBins_summmary.txt"), name=samples),
        named_stats_maxbin2         = expand(join(top_refine_dir, "{name}", "named_maxbin2_bins.stats"), name=samples),
        named_stats_metabat2        = expand(join(top_refine_dir, "{name}", "named_metabat2_bins.stats"), name=samples),
        named_stats_metawrap        = expand(join(top_refine_dir, "{name}", "named_metawrap_bins.stats"), name=samples),
    output:
        cumulative_bin_summary      = join(top_refine_dir, "RefinedBins_summmary.txt"),
        cumulative_maxbin_stats     = join(top_refine_dir, "cumulative_stats_maxbin.txt"),
        cumulative_metabat2_stats   = join(top_refine_dir, "cumulative_stats_metabat2.txt"),
        cumulative_metawrap_stats   = join(top_refine_dir, "cumulative_stats_metawrap.txt"),
    params:
        rname                       = "cumulative_bin_stats",
        refine_dir                  = top_refine_dir
    shell:
        """
        # generate cumulative binning summary
        echo "SampleID\tmetabat2_bins\tmaxbin2_bins\tmetaWRAP_50_5_bins\tmetabat2_contigs\tmaxbin2_contigs\tmetaWRAP_50_5_contigs" > {output.cumulative_bin_summary}
        for report in `ls {params.refine_dir}/*/RefinedBins_summmary.txt`; do 
            tail -n+2 $report >> {output.cumulative_bin_summary}
        done

        # generate cumulative statistic report for binning
        FNs="genome\tcompleteness\tcontamination\tGC\tlineage\tN50\tsize\tbinner"
        echo $FNs > {output.cumulative_metabat2_stats}
        for report in `ls {params.refine_dir}/*/named_maxbin2_bins.stats`; do
            tail -n+2 $report  >> {output.cumulative_maxbin_stats}
        done
        echo $FNs > {output.cumulative_metabat2_stats}
        for report in `ls {params.refine_dir}/*/named_maxbin2_bins.stats`; do
            tail -n+2 $report  >> {output.cumulative_metabat2_stats}
        done
        echo $FNs > {output.cumulative_metawrap_stats}
        for report in `ls {params.refine_dir}/*/named_maxbin2_bins.stats`; do
            tail -n+2 $report >> {output.cumulative_metawrap_stats}
        done

        if $(wc -l < {output.cumulative_bin_summary}) -le 1; then
            echo "No information in bin refinement summary!" 1>&2
            exit 1
        fi
        """


rule prep_genome_info:
    input:
        expand(join(top_refine_dir, "{name}", "named_metawrap_bins.stats"), name=samples),
    output:
        join(top_refine_dir, "genomeInfo.csv")
    localrule: True
    run:
        from csv import DictReader, DictWriter
        combined_info = []
        for sample_info in input:
            with open(sample_info, 'r') as this_csv:
                reader = DictReader(this_csv, delimiter="\t")
                for row in reader:
                    row['genome'] = row['genome'] + '.fa'
                    combined_info.append(dict(
                        genome=row['genome'], 
                        completeness=row['completeness'], 
                        contamination=row['contamination']
                    ))
        with open(output[0], 'w') as genome_info_out:
            wrt = DictWriter(genome_info_out, ['genome', 'completeness', 'contamination'], delimiter=",")
            wrt.writeheader()
            wrt.writerows(combined_info)


rule derep_bins:
    """
        dRep is a step which further refines draft-quality genomes (bins) by using a 
        fast, inaccurate estimation of genome distance, and a slow, accurate measure 
        of average nucleotide identity.

        @Input:
            maxbin2 assembly bins, contigs, stat summaries
            metabat2 assembly bins, contigs, stat summaries
            metawrap assembly bins, contigs, stat summaries

        @Output:
            directory of consensus ensemble bins (deterministic output)
    """
    input:
        bins                        = expand(join(top_binning_dir, "{name}", "{name}_BINNING_COMPLETE"), name=samples),
        ginfo                       = join(top_refine_dir, "genomeInfo.csv"),
    output:
        derep_genome_info           = join(top_refine_dir, "dRep", "data_tables", "Widb.csv"),
        derep_winning_figure        = join(top_refine_dir, "dRep", "figures", "Winning_genomes.pdf"),
        derep_args                  = join(top_refine_dir, "dRep", "log", "cluster_arguments.json"),
        derep_genome_dir            = directory(join(top_refine_dir, "dRep", "dereplicated_genomes")),
        derep_dir                   = directory(join(top_refine_dir, "dRep")),
    singularity: metawrap_container,
    threads: int(cluster["derep_bins"].get("threads", default_threads)),
    params:
        rname                       = "derep_bins",
        tmpdir                      = config['options']['tmp_dir'],
        bindir                      = join(top_binning_dir),
        outto                       = join(top_refine_dir, "dRep"),
        metawrap_dir_name           = "metawrap_50_5_bins",
        # -l LENGTH: Minimum genome length (default: 50000)
        minimum_genome_length       = "10000",
        # -pa[--P_ani] P_ANI: ANI threshold to form primary (MASH) clusters (default: 0.9)
        ani_primary_threshold       = "0.9",
        # -sa[--S_ani] S_ANI: ANI threshold to form secondary clusters (default: 0.95)
        ani_secondary_threshold     = "0.95",
        # -nc[--cov_thresh] COV_THRESH: Minmum level of overlap between genomes when doing secondary comparisons (default: 0.1)
        min_overlap                 = "0.5",
        # -cm[--coverage_method] {total,larger}  {total,larger}: Method to calculate coverage of an alignment
        coverage_method             = 'larger',
        # -comp COMPLETENESS, --completeness COMPLETENESS: Minimum genome completeness (default: 75)
        completeness                = "50",
        # --S_algorithm {fastANI,ANImf,ANIn,goANI,gANI}: Algorithm for secondary clustering comaprisons:
        second_cluster_algo         = "fastANI"
    shell:
        """
        # dRep's CheckM calls (prodigal + hmmfetch marker-gene extraction, run with
        # {threads} parallel workers across every bin) do a lot of small-file I/O.
        # Point dRep's actual working/output directory at node-local lscratch (this
        # rule reserves it via cluster.json's gres) instead of the shared network
        # --tmp-dir: it's faster, and it stops that many-small-files churn from
        # contending with every other job on the shared filesystem. Falls back to
        # the configured --tmp-dir when not running under Slurm (e.g. local testing).
        if [ -n "${{SLURM_JOB_ID:-}}" ] && [ -d "/lscratch/${{SLURM_JOB_ID}}" ]; then
            DREP_WORKROOT="/lscratch/${{SLURM_JOB_ID}}"
        else
            DREP_WORKROOT="{params.tmpdir}"
        fi
        if [ ! -d "$DREP_WORKROOT" ]; then mkdir -p "$DREP_WORKROOT"; fi
        tmp=$(mktemp -d -p "$DREP_WORKROOT")
        drep_work=$(mktemp -d -p "$DREP_WORKROOT")
        # dRep's CheckM calls use Python's multiprocessing.Manager(), which binds a
        # Unix domain socket under $TMPDIR; AF_UNIX socket paths have a hard ~108-byte
        # kernel limit. $DREP_WORKROOT can be arbitrarily long, so point TMPDIR at a
        # short-named symlink (under the system tmp directory, not $DREP_WORKROOT)
        # that resolves to the real tmp dir instead of using the long path directly.
        short_tmp=$(mktemp -u -t drep_ckm.XXXXXX)
        ln -s "${{tmp}}" "${{short_tmp}}"
        export TMPDIR=${{short_tmp}}
        trap 'rm -rf "${{tmp}}" "${{drep_work}}"; rm -f "${{short_tmp}}"' EXIT

        # activate conda environment, initialize checkm, check deps
        . /opt/conda/etc/profile.d/conda.sh && conda activate checkm
        checkm data setRoot /data2/CHECKM_DB
        export CHECKM_DATA_PATH="/data2/CHECKM_DB"
        dRep check_dependencies

        # run drep
        DREP_BINS=$(ls {params.bindir}/*/{params.metawrap_dir_name}/*.fa | tr '\\n' ' ')
        NUM_BINS=$(ls 2>/dev/null -Ubad1 -- {params.bindir}/*/{params.metawrap_dir_name}/*.fa | wc -l)
        # ceiling division so --primary_chunksize is never 0 (dRep raises
        # "ValueError: range() arg 3 must not be zero" on a 0 chunksize,
        # which floor division produces whenever there are fewer than 3 bins)
        NUM_CHUNK=$(( (NUM_BINS + 2) / 3 ))
        NUM_CHUNK=$(echo $NUM_CHUNK | awk '{{print int($1+0.5)}}')
        mkdir -p {params.outto}

        echo "drep bins: ${{DREP_BINS}}"
        echo "num bins: ${{NUM_BINS}}"
        echo "dRep scratch work dir: ${{drep_work}}"

        # threshold informed by NIDDK-5: /data/NIDDK_IDSS/projects/NIDDK-5_metamorph/genomeInfo.csv
        if [[ ${{NUM_BINS}} -ge 2000 ]]; then
            echo "\033[0;36mLARGE bin set detected, using flags '--multiround_primary_clustering' and '--genomeInfo'\033[0m"
            dRep dereplicate -d \\
                -g ${{DREP_BINS}} \\
                -p {threads} \\
                -l {params.minimum_genome_length} \\
                -pa {params.ani_primary_threshold} \\
                -sa {params.ani_secondary_threshold} \\
                -nc {params.min_overlap} \\
                -cm {params.coverage_method} \\
                --genomeInfo {input.ginfo} \\
                --S_algorithm {params.second_cluster_algo} \\
                -comp {params.completeness} \\
                --multiround_primary_clustering \\
                --primary_chunksize $NUM_CHUNK \\
                --run_tertiary_clustering \\
                ${{drep_work}}
        else
            echo "\033[0;36mSMALL bin set detected, **NOT** using flags '--multiround_primary_clustering', '--genomeInfo', or '--run_tertiary_clustering'\033[0m"
            # --run_tertiary_clustering re-clusters dRep's winning/representative
            # genomes as an extra consistency check. With a small, closely-related
            # bin set it's common for every genome to collapse into a single
            # winning representative, which leaves tertiary clustering with an
            # empty distance matrix and crashes dRep ("ValueError: The number of
            # observations cannot be determined on an empty distance matrix").
            # Skip it here; it's only meaningful with enough winners to compare.
            dRep dereplicate -d \\
                -g ${{DREP_BINS}} \\
                -p {threads} \\
                -l {params.minimum_genome_length} \\
                -pa {params.ani_primary_threshold} \\
                -sa {params.ani_secondary_threshold} \\
                -nc {params.min_overlap} \\
                -cm {params.coverage_method} \\
                --S_algorithm {params.second_cluster_algo} \\
                -comp {params.completeness} \\
                --primary_chunksize $NUM_CHUNK \\
                ${{drep_work}}
                touch ${{drep_work}}/dRep.SMALL_BINSET_WARNING
        fi

        # Only the four items downstream rules actually read ever leave lscratch;
        # dRep's heavy intermediate churn (data/checkM, data/prodigal, etc.) is left
        # behind on node-local scratch and cleaned up by the trap above.
        mkdir -p {params.outto}/data_tables {params.outto}/figures {params.outto}/log
        cp -r ${{drep_work}}/data_tables/. {params.outto}/data_tables/
        cp -r ${{drep_work}}/figures/. {params.outto}/figures/
        cp -r ${{drep_work}}/log/. {params.outto}/log/
        cp -r ${{drep_work}}/dereplicated_genomes {params.outto}/dereplicated_genomes
        if [ -f ${{drep_work}}/dRep.SMALL_BINSET_WARNING ]; then
            cp ${{drep_work}}/dRep.SMALL_BINSET_WARNING {params.outto}/dRep.SMALL_BINSET_WARNING
        fi
        """


rule contig_annotation:
    input:
        derep_genome_info           = join(top_refine_dir, "dRep", "data_tables", "Widb.csv"),
        derep_winning_figure        = join(top_refine_dir, "dRep", "figures", "Winning_genomes.pdf"),
        derep_args                  = join(top_refine_dir, "dRep", "log", "cluster_arguments.json"),
        dRep_dir                    = join(top_refine_dir, "dRep", "dereplicated_genomes"),
    output:
        cat_bin2cls_filename        = join(top_refine_dir, "contig_annotation", "out.BAT.bin2classification.txt"),
        cat_bin2cls_official        = join(top_refine_dir, "contig_annotation", "out.BAT.bin2classification.official_names.txt"),
        cat_bin2cls_summary         = join(top_refine_dir, "contig_annotation", "out.BAT.bin2classification.summary.txt"),
    params:
        cat_dir                     = join(top_refine_dir, "contig_annotation"),
        rname                       = "contig_annotation",
        cat_db                      = "/data2/CAT",         # from <root>/config/resources.json
        tax_db                      = "/data2/CAT_tax",     # from <root>/config/resourcs.json
        diamond_exec                = "/usr/bin/diamond",
    threads: int(cluster["contig_annotation"].get("threads", default_threads)),
    singularity: metawrap_container,
    shell:
        """
        if [ ! -d "{params.cat_dir}" ]; then mkdir -p "{params.cat_dir}"; fi
        cd {params.cat_dir} && CAT bins \
            -b {input.dRep_dir} \
            -d {params.cat_db} \
            -t {params.tax_db} \
            --path_to_diamond {params.diamond_exec} \
            --bin_suffix .fa \
            -f 0.5 \
            -n {threads} \
            --force
        cd {params.cat_dir} && CAT add_names \
            -i {output.cat_bin2cls_filename} \
            -o {output.cat_bin2cls_official} \
            -t {params.tax_db} \
            --only_official
        cd {params.cat_dir} && CAT summarise \
            -i {output.cat_bin2cls_official} \
            -o {output.cat_bin2cls_summary}
        """

    
rule bbtools_index_map:
    input:
        derep_genome_info           = join(top_refine_dir, "dRep", "data_tables", "Widb.csv"),
        derep_winning_figure        = join(top_refine_dir, "dRep", "figures", "Winning_genomes.pdf"),
        derep_args                  = join(top_refine_dir, "dRep", "log", "cluster_arguments.json"),
        dRep_dir                    = join(top_refine_dir, "dRep", "dereplicated_genomes"),
        cat_bin2cls_filename        = join(top_refine_dir, "contig_annotation", "out.BAT.bin2classification.txt"),
        cat_bing2cls_official       = join(top_refine_dir, "contig_annotation", "out.BAT.bin2classification.official_names.txt"),
        cat_bing2cls_summary        = join(top_refine_dir, "contig_annotation", "out.BAT.bin2classification.summary.txt"),
        R1                          = join(top_trim_dir, "{name}", "{name}_R1_dehost.fastq"),
        R2                          = join(top_trim_dir, "{name}", "{name}_R2_dehost.fastq"),   
    output:
        index                       = directory(join(top_mags_dir, "{name}", "index")),
        statsfile                   = join(top_mags_dir, "{name}", "DNA", "{name}.statsfile"),
        scafstats                   = join(top_mags_dir, "{name}", "DNA", "{name}.scafstats"),
        covstats                    = join(top_mags_dir, "{name}", "DNA", "{name}.covstat"),
        rpkm                        = join(top_mags_dir, "{name}", "DNA", "{name}.rpkm"),
        refstats                    = join(top_mags_dir, "{name}", "DNA", "{name}.refstats"),
    params:
        rname                       = "bbtools_index_map",
        sid                         = "{name}",
        this_mag_dir                = join(top_mags_dir, "{name}"),
    singularity: metawrap_container,
    threads: int(cluster["bbtools_index_map"].get("threads", default_threads)),
    shell:
        """
	    cd {params.this_mag_dir} && bbsplit.sh -Xmx100g threads={threads} ref={input.dRep_dir}
        mkdir -p {params.this_mag_dir}/DNA
        bbsplit.sh \
            -Xmx100g \
            unpigz=t \
            threads={threads} \
            minid=0.90 \
            path={params.this_mag_dir}/index \
            ref={input.dRep_dir} \
            sortscafs=f \
            nzo=f \
            statsfile={params.this_mag_dir}/DNA/{params.sid}.statsfile \
            scafstats={params.this_mag_dir}/DNA/{params.sid}.scafstats \
            covstats={params.this_mag_dir}/DNA/{params.sid}.covstat \
            rpkm={params.this_mag_dir}/DNA/{params.sid}.rpkm \
            refstats={params.this_mag_dir}/DNA/{params.sid}.refstats \
            in={input.R1} \
            in2={input.R2}
        """


rule gtdbtk_classify:
    input:
        dRep_dir                    = join(top_refine_dir, "dRep", "dereplicated_genomes"),
    output:
        gtdbtk_dir                   = directory(join(top_tax_dir, "GTDBTK_classify_wf"))
    threads: 
        int(cluster["gtdbtk_classify"].get("threads", default_threads)),
    params:
        rname                       = "gtdbtk_classify",
        tmp_safe_dir                = join(config['options']['tmp_dir'], 'gtdbtk_classify'),
        GTDBTK_DB                   = "/data2/GTDBTK_DB",
    singularity: metawrap_container,
    shell:
        """
        # tmp dir
        if [ ! -d "{params.tmp_safe_dir}" ]; then mkdir -p "{params.tmp_safe_dir}"; fi
        tmp=$(mktemp -d -p "{params.tmp_safe_dir}")
        # GTDB-Tk parallelizes Prodigal with Python's multiprocessing.Manager(),
        # which binds a Unix domain socket under --tmpdir; AF_UNIX socket paths
        # have a hard ~108-byte kernel limit. {params.tmp_safe_dir} can be
        # arbitrarily long (e.g. a site-configured --tmp-dir), so point
        # --tmpdir at a short-named symlink (under the system tmp directory)
        # that resolves to the real tmp dir instead of using the long path
        # directly. Files still land under {params.tmp_safe_dir}.
        short_tmp=$(mktemp -u -t gtdbtk.XXXXXX)
        ln -s "${{tmp}}" "${{short_tmp}}"
        trap 'rm -rf "{params.tmp_safe_dir}"; rm -f "${{short_tmp}}"' EXIT

        # activate conda env & db path
        . /opt/conda/etc/profile.d/conda.sh && conda activate checkm
        export GTDBTK_DATA_PATH={params.GTDBTK_DB}
        checkm data setRoot /data2/CHECKM_DB
        export CHECKM_DATA_PATH="/data2/CHECKM_DB"
        gtdbtk classify_wf \
            --genome_dir {input.dRep_dir} \
            --out_dir {output.gtdbtk_dir} \
            --cpus {threads} \
            --tmpdir ${{short_tmp}} \
            --skip_ani_screen \
            --full_tree \
            --force \
            --extension fa
        """


rule gunc_detection:
    input:
        derep_genome_info           = join(top_refine_dir, "dRep", "data_tables", "Widb.csv"),
        derep_winning_figure        = join(top_refine_dir, "dRep", "figures", "Winning_genomes.pdf"),
        derep_args                  = join(top_refine_dir, "dRep", "log", "cluster_arguments.json"),
        drep_genomes                = join(top_refine_dir, "dRep", "dereplicated_genomes"),
    output:
        GUNC_detect_out             = directory(join(top_tax_dir,"GUNC_detect"))
    params:
        rname                       = "gunc_detection",
        gunc_db                     = "/data2/GUNC_DB/gunc_db_progenomes2.1.dmnd", # from <root>/config/resourcs.json
        tmp_safe_dir                = join(config['options']['tmp_dir'], 'gunc_detect'),
    singularity: metawrap_container,
    shell:
        """
        # tmp dir
        if [ ! -d "{params.tmp_safe_dir}" ]; then mkdir -p "{params.tmp_safe_dir}"; fi
        tmp=$(mktemp -d -p "{params.tmp_safe_dir}")
        # gunc parallelizes with Python's multiprocessing.Manager(), which binds a
        # Unix domain socket under --temp_dir; AF_UNIX socket paths have a hard
        # ~108-byte kernel limit. {params.tmp_safe_dir} can be arbitrarily long
        # (e.g. a site-configured --tmp-dir), so point --temp_dir at a
        # short-named symlink (under the system tmp directory) that resolves to
        # the real tmp dir instead of using the long path directly. Files still
        # land under {params.tmp_safe_dir}.
        short_tmp=$(mktemp -u -t gunc.XXXXXX)
        ln -s "${{tmp}}" "${{short_tmp}}"
        trap 'rm -rf "{params.tmp_safe_dir}"; rm -f "${{short_tmp}}"' EXIT
        # activate conda env
        . /opt/conda/etc/profile.d/conda.sh && conda activate checkm
        # run gunc
        mkdir -p {output.GUNC_detect_out}
        gunc run \
            --input_dir {input.drep_genomes} \
            --threads {threads} \
            --temp_dir ${{short_tmp}} \
            --out_dir {output.GUNC_detect_out} \
            --sensitive \
            --detailed_output \
            --db_file {params.gunc_db}
        """


