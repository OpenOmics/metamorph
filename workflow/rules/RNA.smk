# ~~~~~~~~~~
# Metawrap metagenome assembly and analysis rules
# ~~~~~~~~~~
from os.path import join
from itertools import chain
from scripts.common import str_bool, list_bool

# ~~~~~~~~~~
# Constants and paths
# ~~~~~~~~~~
workpath                            = config["project"]["workpath"]
rna_datapath                        = config["project"].get("rna_datapath", "/dev/null")
rna_included                        = list_bool(config.get("rna", 'false'))
# rna_coasm                         = str_bool(config["options"].get("rnacoa", 'False'))
rna_coasm                           = False
rna_sample_stems                    = config.get("rna", [])
rna_compressed                      = True # if accepting uncompressed fastq input
top_readqc_dir_rna                  = join(workpath, config['project']['id'], "metawrap_read_qc_RNA")
top_trim_dir_rna                    = join(workpath, config['project']['id'], "trimmed_reads_RNA")
metawrap_container                  = config["containers"]["metawrap"]
pairedness                          = list(range(1, config['project']['nends']+1))


rule concat_rna_reads:
    input:
        all_r1_reads                = expand(join(workpath, "rna", "{rname}_R1.fastq.gz"), rname=rna_sample_stems if rna_coasm else []),
        all_r2_reads                = expand(join(workpath, "rna", "{rname}_R2.fastq.gz"), rname=rna_sample_stems if rna_coasm else []),
    output:
        big_compressed_read_r1      = join(workpath, "rna", "concatenated_R1.fastq.gz"),
        big_compressed_read_r2      = join(workpath, "rna", "concatenated_R2.fastq.gz"),
        big_read1_hash              = join(workpath, "rna", "concatenated_R1.md5"),
        big_read2_hash              = join(workpath, "rna", "concatenated_R2.md5"),
    params:
        rname                       = "concat_rna_reads",
        big_read_r1                 = join(workpath, "dna", "concatenated_R1.fastq"),
        big_read_r2                 = join(workpath, "dna", "concatenated_R2.fastq"),
        input_dir                   = workpath,
    threads: int(cluster["concat_rna_reads"].get('threads', default_threads)),
    shell: 
        """
        # concat r1
        for fastq in {params.input_dir}/*R1*fastq; do
            ext=$(echo "${{fastq: -2}}" | tr '[:upper:]' '[:lower:]')
            if [[ "$ext" == "gz" ]]; then
                zcat $fastq >> {params.big_read_r1}
            else
                cat $fastq >> {params.big_read_r1}
            fi;
        done

        # concat r2
        for fastq in {params.input_dir}/*R2*fastq; do 
            ext=$(echo "${{fastq: -2}}" | tr '[:upper:]' '[:lower:]')
            if [[ "$ext" == "gz" ]]; then
                zcat $fastq > {params.big_read_r2}
            else
                cat $fastq >> {params.big_read_r2}
            fi;
        done
        pigz -9 -p 28 -c {output.big_read_r1} > {output.big_compressed_read_r1}
        pigz -9 -p 28 -c {output.big_read_r2} > {output.big_compressed_read_r2}
        md5sum {output.big_compressed_read_r1} > {output.big_read1_hash}
        md5sum {output.big_compressed_read_r2} > {output.big_read2_hash}
        """


rule rna_read_qc:
    input:
        R1                  = join(workpath, "rna", "{rname}_R1.fastq.gz"),
        R2                  = join(workpath, "rna", "{rname}_R2.fastq.gz"),
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
