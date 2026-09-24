# <code>metamorph <b>run</b></code>

## 1. About

The `metamorph` executable is composed of several inter-related subcommands. Please see `metamorph -h` for all available options.

This part of the documentation describes the options and concepts for the <code>metamorph <b>run</b></code> sub-command in more detail. With minimal configuration, the **`run`** sub-command allows you to start running the metagenomics/metatranscriptomics pipeline.

Setting up the metamorph pipeline is fast and easy! In its most basic form, <code>metamorph <b>run</b></code> requires two arguments:

1. A sample sheet provided with `--samplesheet`
2. An output directory provided with `--output`

## 2. Synopsis

```text
$ metamorph run [--help] \
      [--mode {slurm,local}] [--job-name JOB_NAME] \
      [--master-job-node {norm,unlimited,quick}] \
      [--triggers TRIGGER [TRIGGER ...]] [--dry-run] [--silent] \
      [--singularity-cache SINGULARITY_CACHE] [--sif-cache SIF_CACHE] \
      [--tmp-dir TMP_DIR] [--threads THREADS] \
      [--assembler ASSEMBLER] [--shallow-profile] \
      [--workflow {pre-screen,read-based,assembly-based,combined}] \
      --samplesheet SHEETPATH \
      --output OUTPUT
```

Optional arguments are shown in square brackets. A user **must** provide a sample sheet with `--samplesheet` and an output directory with `--output`.

Use `metamorph run -h` to display the command-line help for this subcommand.

### 2.1 Required arguments

Each of the following arguments is required. Failure to provide a required argument results in a non-zero exit code.

`--samplesheet SHEETPATH`

> **Full path to the delimited plain-text sample sheet.**  
> *type: file*
>
> The sample sheet is parsed before the pipeline is initialized. It must contain paths to the input DNA FASTQ files and can optionally contain a second column containing matching RNA FASTQ files.
>
> ***Example:*** `--samplesheet /home/user/samplesheet.txt`

> [!NOTE]
> Supported delimiters are inferred from the file extension:
>
> | file extension  | delimiter |
> | --------------- | --------- |
> | `*.txt`         | tab       |
> | `*.tsv`         | tab       |
> | `*.csv`         | comma     |

**Supported sample sheet formats**

> [!NOTE]
> `Path` = full absolute path to input FASTQ.GZ files  
> `R1 sample X` / `R2 sample X` = paired-end read orientation  
> `n` = total number of samples

**DNA only - 1 column: `DNA`**

```text
                ________
                |  DNA  |
                |-------|
R1 sample 1     | path  |
R2 sample 1     | path  |
...             | path  |
R1 sample n     | path  |
R2 sample n     | path  |
                --------
```

**DNA & RNA - 2 columns: `DNA`, `RNA`**

```text
                ________________
                |  DNA  |  RNA  |
                |---------------|
R1 sample 1     | path  | path  |
R2 sample 1     | path  | path  |
...             | path  | path  |
R1 sample n     | path  | path  |
R2 sample n     | path  | path  |
                ----------------
```

---

`--output OUTPUT`

> **Path to an output directory.**  
> *type: path*
>
> This location is the pipeline working directory and is where metamorph creates its output files. If the directory does not exist, it is created automatically.
>
> ***Example:*** `--output /data/$USER/metamorph_out`

### 2.2 Analysis options

Each of the following arguments is optional, and do not need to be provided. 

`--workflow {pre-screen,read-based,assembly-based,combined}`

> **Select the analysis workflow.**  
> *type: string*  
> *default: `pre-screen`*
>
> Available workflows:
>
> - `pre-screen`: perform quality control, host-read removal, and run Centrifuge and Kraken2 on dehosted reads.
> - `read-based`: perform functional and taxonomic profiling on dehosted reads using HUMAnN3 and MetaPhlAn4. This stage follows pre-screening.
> - `assembly-based`: perform the assembly and binning workflow on dehosted reads. This stage also follows pre-screening.
> - `combined`: perform all steps included in the preceding modes.
>
> ***Example:*** `--workflow combined`

---

`--assembler ASSEMBLER`

> **Select the assembler(s) used by the assembly-based workflow.**  
> *type: comma-separated string*  
> *default: `megahit,metaspades`*
>
> Supported assemblers are `megahit` and `metaspades`. One or both can be supplied. The command validates the comma-separated values after parsing and rejects unsupported assembler names.
>
> ***Examples:***
>
> ```text
> --assembler megahit,metaspades
> --assembler megahit
> --assembler metaspades
> ```

---

`--shallow-profile`

> **Use shallow HUMAnN3 profiling.**  
> *type: boolean flag*  
> *default: Runs translated protein search*
>
> Disables the translated protein search of unaligned reads in HUMAnN3, which can significantly reduce runtime.
>
> ***Example:*** `--shallow-profile`

### 2.3 Orchestration options

Each of the following arguments is optional.

`--mode {slurm,local}`

> **Select the execution method.**  
> *type: string*  
> *default: `slurm`*
>
> Available modes:
>
> - `slurm`: submits the pipeline and child jobs to the [SLURM workload manager](https://slurm.schedmd.com/). It is recommended running metamorph in this mode as execution will be significantly faster in a distributed environment. This is the default mode of execution. 
> - `local`: run the pipeline locally on the current compute instance. This is useful for testing, debugging, or environments without a supported cluster job scheduler.
>
> ***Example:*** `--mode slurm`

---

`--job-name JOB_NAME`

> **Set the name of the pipeline master job.**  
> *type: string*  
> *default: `pl:metamorph`*
>
> Overrides the default name used for the master job when the pipeline is submitted to a cluster job scheduler.  By default, the name of the pipeline's master job is set to "pl:metamorph".
>
> ***Example:*** `--job-name pl_id-2026-09-30`

---

`--master-job-node {norm,unlimited,quick}`

> **Select the node/partition class used for the master job.**  
> *type: string*  
> *default: `unlimited`*
>
> The selected value is passed to the master-job submission wrapper as its job-node setting. Valid values are `norm`, `unlimited`, and `quick`.
>
> ***Example:*** `--master-job-node norm`

---

`--triggers TRIGGER [TRIGGER ...]`

> **Control Snakemake rerun triggers.**  
> *type: zero or more validated strings*  
> *default(in prority order): params,mtime,code,software-env,input*
>
> Valid trigger names are:
>
> - `mtime`
> - `code`
> - `software-env`
> - `input`
> - `params`
>
> Comma-separated trigger names are also accepted within a value, so both of the following forms are valid:
>
> ```text
> --triggers mtime code
> --triggers mtime,code
> ```
>
> See the [Snakemake command-line documentation](https://snakemake.readthedocs.io/en/stable/executing/cli.html) for details about rerun triggers. Triggers are tied to the version of snakemake you are using.
>
> ***Example:*** `--triggers mtime,code`

---

`--dry-run`

> **Dry run the pipeline.**  
> *type: boolean flag*  
> *default: No dryrun is performed*
>
> Displays the steps that remain or would be run without executing the workflow. 
>
> ***Example:*** `--dry-run`

---

`--silent`

> **Silence standard output from master-job submission.**  
> *type: boolean flag*  
> *default: Standard output is more verbose*
>
> Reduces the information written to standard output when the master job is submitted. In SLURM mode, only the master job ID is printed.
>
> ***Example:*** `--silent`

---

`--singularity-cache SINGULARITY_CACHE`

> **Override the Singularity cache directory.**  
> *type: path*  
> *default: `<OUTPUT>/.singularity` when no alternate cache is provided*
>
> Singularity will cache image layers pulled from remote registries. This ultimately speeds up the process of pull an image from DockerHub if an image layer already exists in the singularity cache directory. By default, the cache is set to the value provided to the `--output` argument. Please note that this cache cannot be shared across users. Singularity strictly enforces you own the cache directory and will return a non-zero exit code if you do not own the cache directory! See the `--sif-cache` option to create a shareable resource. 
>
> ***Example:*** `--singularity-cache /data/$USER/.singularity`

---

`--sif-cache SIF_CACHE`

> **Use a local cache of Singularity Image Format (SIF) files.**  
> *type: path*  
> *default: `<OUTPUT>/.snakemake` when no local sif cache is provided*
>
> Uses a local cache of SIFs on the filesystem. This SIF cache can be shared across users if permissions are set correctly. If a SIF does not exist in the SIF cache, the image will be pulled from Dockerhub and a warning message will be displayed. The `metamorph cache` subcommand can be used to create a local SIF cache. Please see `metamorph cache` for more information. This command is extremely useful for avoiding DockerHub pull rate limits. It also remove any potential errors that could occur due to network issues or DockerHub being temporarily unavailable. We recommend running metamorph with this option when ever possible. 
>
> ***Example:*** `--sif-cache /data/OpenOmics/SIFs`

---

`--tmp-dir TMP_DIR`

> **Set the base directory for temporary/intermediate files.**  
> *type: path*  
> *default: `/lscratch/$SLURM_JOB_ID/`*
>
> Path on the file system for writing temporary output files. By default, the temporary directory is set to '/lscratch/$SLURM_JOBID' for backwards compatibility with the NIH's Biowulf cluster; however, if you are running the pipeline on another cluster, this option will need to be specified. Ideally, this path should point to a dedicated location on the filesystem for writing tmp files. On many systems, this location is set to somewhere in /scratch. If you need to inject a variable into this string that should NOT be expanded, please quote this options value in single quotes.
>
> ***Example:*** `--tmp-dir '/data/scratch/$USER/'`

---

`--threads THREADS`

> **Set the number of threads used by the pipeline's main/local processes.**  
> *type: integer*  
> *default: `2`*
>
>  Max number of threads for each process. This option is more applicable when running the pipeline with `--mode local`.  It is recommended setting this vaule to the maximum number of CPUs available on the host machine.
>
> ***Example:*** `--threads 12`

### 2.4 Miscellaneous options

`-h, --help`

> **Display command help.**  
> *type: boolean flag*
>
> Shows the command synopsis, help message, and examples, then exits.
>
> ***Example:*** `metamorph run --help`

## 3. Examples

### 3.1 Dry run with default options

```bash
# Step 1.) Grab an interactive node,
# do not run on head node!
srun -N 1 -n 1 --time=1:00:00 --mem=8gb  --cpus-per-task=2 --pty bash
module purge
module load singularity snakemake

metamorph run \
  --samplesheet /path/to/samplesheet.tsv \
  --output /data/$USER/output \
  --mode slurm \
  --dry-run
```

### 3.2 Run the default pre-screen workflow

```bash
metamorph run \
  --samplesheet /path/to/samplesheet.tsv \
  --output /data/$USER/output \
  --mode slurm
```

### 3.3 Run the combined workflow with one assembler

```bash
metamorph run \
  --samplesheet /path/to/samplesheet.tsv \
  --output /data/$USER/output \
  --mode slurm \
  --workflow combined \
  --assembler megahit
```

### 3.4 Run with shallow HUMAnN3 profiling and custom rerun triggers

```bash
metamorph run \
  --samplesheet /path/to/samplesheet.tsv \
  --output /data/$USER/output \
  --workflow read-based \
  --shallow-profile \
  --triggers params,mtime,code,software-env,input
```
