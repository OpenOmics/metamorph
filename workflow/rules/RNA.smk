# ~~~~~~~~~~
# Metawrap metagenome assembly and analysis rules
# ~~~~~~~~~~
from os.path import join
from itertools import chain
from functools import partial
from scripts.common import str_bool, list_bool, get_paired_dna, get_dna_sid

# ~~~~~~~~~~
# Constants and paths
# ~~~~~~~~~~
workpath                            = config["project"]["workpath"]
top_mapping_dir                     = join(workpath, config['project']['id'], "humann3_dna")
rna_datapath                        = config["project"].get("rna_datapath", "/dev/null")
rna_included                        = list_bool(config.get("rna", 'false'))
rna_sample_stems                    = config.get("rna", [])
rna_compressed                      = True # if accepting uncompressed fastq input
get_dna                             = partial(get_paired_dna, config)
get_dna_bugs                        = partial(get_dna_sid, config)
metawrap_container                  = config["containers"]["metawrap"]
pairedness                          = list(range(1, config['project']['nends']+1))
top_readqc_dir_rna                  = join(workpath, config['project']['id'], "metawrap_read_qc_RNA")
top_trim_dir_rna                    = join(workpath, config['project']['id'], "trimmed_reads_RNA")
top_map_dir_rna                     = join(workpath, config['project']['id'], "mapping_RNA")
humann3_dir_rna                     = join(workpath, config['project']['id'], "humann3_rna")
top_centrifuger_dir_rna             = join(workpath, config['project']['id'], "centrifuger_rna")
humann_deep_mode                    = True if "deep_profile" in config["options"] and \
                                      bool(int(config["options"]["deep_profile"])) else False
# containers
star_container              = config["containers"]["star"]
centrifuger_sylph_container = config["containers"]["centrifuger_sylph"]

rule rna_read_qc_skipBmtagger:
    input:
        R1                  = join(workpath, "rna", "{rname}_R1.fastq.gz") if rna_included else [],
        R2                  = join(workpath, "rna", "{rname}_R2.fastq.gz") if rna_included else [],
    output:
        R1_pretrim_report   = join(top_readqc_dir_rna, "{rname}", "{rname}_R1_pretrim_report.html"),
        R2_pretrim_report   = join(top_readqc_dir_rna, "{rname}", "{rname}_R2_pretrim_report.html"),
        R1_postrim_report   = join(top_readqc_dir_rna, "{rname}", "{rname}_R1_postrim_report.html"),
        R2_postrim_report   = join(top_readqc_dir_rna, "{rname}", "{rname}_R2_postrim_report.html"),
        R1_trimmed          = join(top_trim_dir_rna, "{rname}", "{rname}_R1_trimmed.fastq"),
        R2_trimmed          = join(top_trim_dir_rna, "{rname}", "{rname}_R2_trimmed.fastq"),
        R1_trimmed_gz       = join(top_trim_dir_rna, "{rname}", "{rname}_R1_trimmed.fastq.gz"),
        R2_trimmed_gz       = join(top_trim_dir_rna, "{rname}", "{rname}_R2_trimmed.fastq.gz"),
    params:
        rname               = "rna_read_qc",
        sid                 = "{rname}",
        this_qc_dir         = join(top_readqc_dir_rna, "{rname}"),
        trim_out            = join(top_trim_dir_rna, "{rname}"),
        tmp_safe_dir        = join(config['options']['tmp_dir'], 'read_qc'),
        tmpr1               = lambda _, output, input: join(config['options']['tmp_dir'], 'read_qc', str(basename(str(input.R1))).replace('_R1.', '_1.').replace('.gz', '')),
        tmpr2               = lambda _, output, input: join(config['options']['tmp_dir'], 'read_qc', str(basename(str(input.R2))).replace('_R2.', '_2.').replace('.gz', '')),
    containerized: metawrap_container,
    threads: int(cluster["rna_read_qc"].get('threads', default_threads)),
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
            ln -s {params.this_qc_dir}/final_pure_reads_1.fastq {params.trim_out}/{params.sid}_R1_trimmed.fastq
            ln -s {params.this_qc_dir}/final_pure_reads_2.fastq {params.trim_out}/{params.sid}_R2_trimmed.fastq
            pigz -9 -p {threads} -c {params.this_qc_dir}/final_pure_reads_1.fastq  > {params.trim_out}/{params.sid}_R1_trimmed.fastq.gz
            pigz -9 -p {threads} -c {params.this_qc_dir}/final_pure_reads_2.fastq > {params.trim_out}/{params.sid}_R2_trimmed.fastq.gz
            
            # collate outputs to facilitate 
            ln -s {params.this_qc_dir}/post-QC_report/final_pure_reads_1_fastqc.html {params.this_qc_dir}/{params.sid}_R1_postrim_report.html
            ln -s {params.this_qc_dir}/post-QC_report/final_pure_reads_2_fastqc.html {params.this_qc_dir}/{params.sid}_R2_postrim_report.html
            ln -s {params.this_qc_dir}/pre-QC_report/{params.sid}_1_fastqc.html {params.this_qc_dir}/{params.sid}_R1_pretrim_report.html
            ln -s {params.this_qc_dir}/pre-QC_report/{params.sid}_2_fastqc.html {params.this_qc_dir}/{params.sid}_R2_pretrim_report.html
        """

rule rna_dehost:
    """
        This module utilizes STAR to map trimmed reads to hg38 transcriptome, and select the untrimmed reads from the bam file as clean read.
    """
    input:
        R1                          = join(top_trim_dir_rna, "{rname}", "{rname}_R1_trimmed.fastq.gz"),
        R2                          = join(top_trim_dir_rna, "{rname}", "{rname}_R2_trimmed.fastq.gz"),
    output:
        R1_dehost                   = join(top_trim_dir_rna, "{rname}", "{rname}_R1_dehost.fastq.gz"),
        R2_dehost                   = join(top_trim_dir_rna, "{rname}", "{rname}_R2_dehost.fastq.gz"),

    params:
        rname                       = "rna_dehost",
        sid                         = "{rname}",
        tmp_safe_dir                = join(config['options']['tmp_dir'], 'star_dehost', "{rname}"),
        tmp_safe_dir_star           = join(config['options']['tmp_dir'], 'star_dehost', "{rname}", "STAR"),
        hg38_star_idx               = "/data2/star",
        hg38_gtf                    = join("/data2/gtf","gencode.v30.annotation.gtf")
    containerized: star_container,
    threads: int(cluster["rna_dehost"].get('threads', default_threads)),
    shell: 
        """
            # safe temp directory
            if [ ! -d "{params.tmp_safe_dir}" ]; then mkdir -p "{params.tmp_safe_dir}"; fi
            tmp=$(mktemp -d -p "{params.tmp_safe_dir}")
            trap 'rm -rf "{params.tmp_safe_dir}"' EXIT

            # run star to map reads to hg38
            ## Optimal readlength for sjdbOverhang = max(ReadLength) - 1 [Default: 100]
            readlength=$(
                        zcat {input.R1} | \
                        awk -v maxlen=100 'NR%4==2 {{if (length($1) > maxlen+0) maxlen=length($1)}}; \
                        END {{print maxlen-1}}'
                    )

            STAR \
            --runThreadN {threads} \
            --sjdbOverhang ${{readlength}} \
            --twopassMode Basic \
            --genomeDir {params.hg38_star_idx} \
            --sjdbGTFfile {params.hg38_gtf} \
            --readFilesIn {input.R1} {input.R2} \
            --readFilesCommand zcat \
            --outTmpDir {params.tmp_safe_dir_star} \
            --outFileNamePrefix {params.tmp_safe_dir_star}_ \
            --outReadsUnmapped Fastx 
            
            # The unmapped reads will be named {params.tmp_safe_dir}_Unmapped.out.mate1 and {params.tmp_safe_dir}_Unmapped.out.mate2. Their hearder lines need to be changed
            # read 1
            sed 's/0:N:  /1:N:/' {params.tmp_safe_dir_star}_Unmapped.out.mate1 > {params.tmp_safe_dir_star}_Unmapped.out.correctHeaderLine.mate1
            # read 2
            sed 's/1:N:  /2:N:/' {params.tmp_safe_dir_star}_Unmapped.out.mate2 > {params.tmp_safe_dir_star}_Unmapped.out.correctHeaderLine.mate2

            # compress them
            gzip -9 -c {params.tmp_safe_dir_star}_Unmapped.out.correctHeaderLine.mate1 > {output.R1_dehost}
            gzip -9 -c {params.tmp_safe_dir_star}_Unmapped.out.correctHeaderLine.mate2 > {output.R2_dehost}

        """

rule rna_centrifuger:
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
        R1_dehost                   = join(top_trim_dir_rna, "{name}", "{name}_R1_dehost.fastq.gz"),
        R2_dehost                   = join(top_trim_dir_rna, "{name}", "{name}_R2_dehost.fastq.gz"),
    output:
        classification              = join(top_centrifuger_dir_rna, "{name}_centrifuger_classification.tsv"),
        centrifuger_quant           = join(top_centrifuger_dir_rna, "{name}_centrifuger_quantification_report.tsv"),
        metaphlan_quant             = join(top_centrifuger_dir_rna, "{name}_centrifuger_quantification_report_mpl.tsv"),
    params:
        rname                       = "rna_centrifuger",
        sid                         = "{name}",
        centrifuger_idx_prefix      = "/data2/centrifuger/gtdb_r226+refseq_hvfc/cfr_gtdb_r226+refseq_hvfc",
    container: centrifuger_sylph_container,
    threads: int(cluster["rna_centrifuger"].get('threads', default_threads)),
    shell: 
        """
        # Runs centrifuger on RNA to classify taxonomic
        # composition of host removed reads
        centrifuger \\
            -t {threads} \\
            -x {params.centrifuger_idx_prefix} \\
            -1 {input.R1_dehost} \\
            -2 {input.R2_dehost}  \\
        > {output.classification}
        # Runs centrifuger quant on RNA to quantify
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

rule rna_humann_classify:
    input:
        R1                  = join(top_trim_dir_rna, "{rname}", "{rname}_R1_dehost.fastq.gz"),
        R2                  = join(top_trim_dir_rna, "{rname}", "{rname}_R2_dehost.fastq.gz"),
        DNA_bug_list        = lambda wildcards: join(top_mapping_dir, f"{get_dna_bugs(wildcards.rname)}_bugs_list.tsv")
    output:
        hm3_gene_fam        = join(humann3_dir_rna, '{rname}_genefamilies.tsv'),
        hm3_path_abd        = join(humann3_dir_rna, '{rname}_pathabundance.tsv'),
        hm3_path_cov        = join(humann3_dir_rna, '{rname}_pathcoverage.tsv'),
        humann_log          = join(humann3_dir_rna, '{rname}_humann3.log'),
        humann_config       = join(humann3_dir_rna, '{rname}_humann3.conf'),
    params:
        rname               = "rna_humann_classify",
        sid                 = "{rname}",
        tmpread             = join(config['options']['tmp_dir'], 'rna_map', "{rname}_concat.fastq.gz"),
        tmp_safe_dir        = join(config['options']['tmp_dir'], 'rna_map'),
        hm3_map_dir         = humann3_dir_rna,
        uniref_db           = "/data2/uniref",      # from <root>/config/resources.json
        chocophlan_db       = "/data2/chocophlan",  # from <root>/config/resources.json
        util_map_db         = "/data2/um",          # from <root>/config/resources.json
        metaphlan_db        = "/data2/metaphlan",   # from <root>/config/resources.json
        deep_mode           = "\\\n        --bypass-translated-search" if not config["options"]["shallow_profile"] else "",
    threads: int(cluster["rna_humann_classify"].get('threads', default_threads)),
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
        --remove-temp-output \
        --input-format fastq.gz {params.deep_mode} \
        --taxonomic-profile {input.DNA_bug_list} \
        --output-basename {params.sid} \
        --log-level DEBUG \
        --o-log {output.humann_log} \
        --output {params.hm3_map_dir}
        """

rule rna_humann_summarize:
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
        path_abd            = expand(join(humann3_dir_rna, '{rname}_pathabundance.tsv'), rname=rna_sample_stems),
  
    output:
        mer_path_abd_str    = join(humann3_dir_rna, 'merged_pathabundance.cpm.stratified.rna.tsv'),
        mer_path_abd_unstr  = join(humann3_dir_rna, 'merged_pathabundance.cpm.unstratified.rna.tsv'),
    params:
        rname               = "rna_humann_summarize",
        hm3_map_dir         = humann3_dir_rna,
        stratified_dir      = join(humann3_dir_rna,'stratified'),
        unstratified_dir    = join(humann3_dir_rna,'unstratified'),
    containerized: metawrap_container,
    threads: int(cluster["rna_humann_summarize"].get('threads', default_threads)),
    shell:
        """
        . /opt/conda/etc/profile.d/conda.sh && conda activate bb3

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

rule map_to_rna_to_mag:
    input:
        dna_input                   = lambda wc: get_dna(wc.rname) if rna_included else [],
        R1                          = join(top_trim_dir_rna, "{rname}", "{rname}_R1_trimmed.fastq.gz") if rna_included else [],
        R2                          = join(top_trim_dir_rna, "{rname}", "{rname}_R2_trimmed.fastq.gz") if rna_included else [],
    output:
        aligned_rna                 = join(top_map_dir_rna, "{rname}", "{rname}.RNA.aligned.sam"),
        statsfile                   = join(top_map_dir_rna, "{rname}", "{rname}.RNA.statsfile"),
        scafstats                   = join(top_map_dir_rna, "{rname}", "{rname}.RNA.scafstats"),
        covstats                    = join(top_map_dir_rna, "{rname}", "{rname}.RNA.covstat"),
        rpkm                        = join(top_map_dir_rna, "{rname}", "{rname}.RNA.rpkm"),
        refstats                    = join(top_map_dir_rna, "{rname}", "{rname}.RNA.refstats"),
    params:
        rname                       = "map_to_rna_to_mag",
        sid                         = "{rname}",
        map_rna_dir                 = join(top_map_dir_rna, "{rname}"),
        minid                       = "0.90",
        mag_dir                     = lambda wc, input: input.dna_input[0] if input.dna_input else [],
        mag_idx                     = lambda wc, input: input.dna_input[1] if input.dna_input else [],
    threads: int(cluster["map_to_rna_to_mag"].get('threads', default_threads)),
    containerized: metawrap_container,
    shell:
        """
        mkdir -p {params.map_rna_dir}
		bbsplit.sh \
            -da \
            -Xmx100g \
            unpigz=t \
            threads={threads} \
            minid={params.minid} \
            outm={output.aligned_rna} \
            mappedonly=t \
            path={params.mag_idx} \
            ref={params.mag_dir} \
            sortscafs=f \
            nzo=f \
            statsfile={output.statsfile} \
            scafstats={output.scafstats} \
            covstats={output.covstats} \
            rpkm={output.rpkm} \
            refstats={output.refstats} \
            in={input.R1} \
            in2={input.R2}
        """
