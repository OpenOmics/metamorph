# ~~~~~~~~~~
# Metawrap metagenome assembly and analysis rules
# ~~~~~~~~~~
from os.path import join
from itertools import chain
from functools import partial
from scripts.common import str_bool, list_bool, get_paired_dna

# ~~~~~~~~~~
# Constants and paths
# ~~~~~~~~~~
workpath                            = config["project"]["workpath"]
rna_datapath                        = config["project"].get("rna_datapath", "/dev/null")
rna_included                        = list_bool(config.get("rna", 'false'))
rna_sample_stems                    = config.get("rna", [])
rna_compressed                      = True # if accepting uncompressed fastq input
get_dna                             = partial(get_paired_dna, config)
metawrap_container                  = config["containers"]["metawrap"]
pairedness                          = list(range(1, config['project']['nends']+1))
top_readqc_dir_rna                  = join(workpath, config['project']['id'], "metawrap_read_qc_RNA")
top_trim_dir_rna                    = join(workpath, config['project']['id'], "trimmed_reads_RNA")
top_map_dir_rna                     = join(workpath, config['project']['id'], "mapping_RNA")
humann_deep_mode                    = True if "deep_profile" in config["options"] and \
                                      bool(int(config["options"]["deep_profile"])) else False


rule rna_read_qc:
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

            # read quality control, host removal
            # TODO: add support for mouse reads (mm10 genome prefix, "-x mm10")
            mw read_qc -1 {params.tmpr1} -2 {params.tmpr2} -t {threads} -o {params.this_qc_dir}

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


rule rna_humann_classify:
    input:
        R1                  = join(top_trim_dir_rna, "{rname}", "{rname}_R1_trimmed.fastq.gz"),
        R2                  = join(top_trim_dir_rna, "{rname}", "{rname}_R2_trimmed.fastq.gz"),
    output:
        hm3_gene_fam        = join(top_map_dir_rna, "{rname}", 'humann3', '{rname}_genefamilies.tsv'),
        hm3_path_abd        = join(top_map_dir_rna, "{rname}", 'humann3', '{rname}_pathabundance.tsv'),
        hm3_path_cov        = join(top_map_dir_rna, "{rname}", 'humann3', '{rname}_pathcoverage.tsv'),
        humann_log          = join(top_map_dir_rna, "{rname}", 'humann3.log'),
        humann_config       = join(top_map_dir_rna, "{rname}", 'humann3.conf'),
    params:
        rname               = "rna_humann_classify",
        sid                 = "{rname}",
        tmpread             = join(config['options']['tmp_dir'], 'rna_map', "{rname}_concat.fastq.gz"),
        tmp_safe_dir        = join(config['options']['tmp_dir'], 'rna_map'),
        hm3_map_dir         = join(top_map_dir_rna, "{rname}", 'humann3'),
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
        --metaphlan-options "--bowtie2db {params.metaphlan_db} --nproc {threads}" \
        --output-basename {params.sid} \
        --log-level DEBUG \
        --o-log {output.humann_log} \
        --output {params.hm3_map_dir}
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
