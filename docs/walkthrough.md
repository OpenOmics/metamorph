# Walkthrough

## 1. About

This page walks through setting up and running metamorph end-to-end using a small, public, real paired-end sequencing dataset, so you can see what each [`--workflow`](usage/run.md#22-analysis-options) mode actually does before pointing the pipeline at your own data. If you haven't already, read the [Introduction](manual/introduction.md) first — it explains *why* the pipeline is structured the way it is; this page focuses on *doing*.

We'll use two tiny paired-end FASTQ samples from the [`nf-core/test-datasets`](https://github.com/nf-core/test-datasets/tree/mag/test_data) `mag` branch. They're real assembled-metagenome reads with a small fraction of human (hg38) reads spiked in, which is exactly what makes them useful here — they're small enough to run in minutes, but large enough to exercise metamorph's host-removal step, not just skip it.

## 2. Prerequisites

```bash
# Grab an interactive node — do not run on the head node!
srun -N 1 -n 1 --time=1:00:00 --mem=8gb --cpus-per-task=4 --pty bash
module purge
module load singularity snakemake
```

If you're running outside of a cluster that already hosts metamorph's reference bundle (e.g. off of NIH's Biowulf), also see [`metamorph install`](usage/install.md) for pulling the reference bundle, and [`metamorph cache`](usage/cache.md) for pre-fetching containers into a shared, offline `--sif-cache`.

## 3. Getting a small, public test dataset

```bash
mkdir -p /data/$USER/metamorph_test/fastqs
cd /data/$USER/metamorph_test/fastqs

BASE="https://raw.githubusercontent.com/nf-core/test-datasets/mag/test_data"
for f in test_minigut_hg38host_R1.fastq.gz  test_minigut_hg38host_R2.fastq.gz \
         test_minigut_sample2_hg38host_R1.fastq.gz test_minigut_sample2_hg38host_R2.fastq.gz; do
    curl -sL -o "$f" "${BASE}/${f}"
done
```

This gives you two DNA samples:

| Sample                        | Read pairs | Compressed size (R1 + R2) |
| ------------------------------ | ---------- | -------------------------- |
| `test_minigut_hg38host`        | 51,000     | ~6.9 MB                    |
| `test_minigut_sample2_hg38host`| 25,993     | ~3.6 MB                    |

## 4. Building a sample sheet

metamorph's sample sheet is a plain delimited text file with one column per data type (`DNA`, and optionally `RNA`), listing the R1/R2 paths for every sample — see [`metamorph run`](usage/run.md#21-required-arguments) for the full format reference. This dataset is DNA-only, so it's a single `DNA` column:

```text title="samplesheet.tsv"
DNA
/data/$USER/metamorph_test/fastqs/test_minigut_hg38host_R1.fastq.gz
/data/$USER/metamorph_test/fastqs/test_minigut_hg38host_R2.fastq.gz
/data/$USER/metamorph_test/fastqs/test_minigut_sample2_hg38host_R1.fastq.gz
/data/$USER/metamorph_test/fastqs/test_minigut_sample2_hg38host_R2.fastq.gz
```

!!! tip
    `.tsv`/`.txt` sample sheets are tab-delimited; `.csv` sample sheets are comma-delimited. Paths must be absolute.

## 5. Dry-run the pipeline

Before running anything for real, `--dry-run` shows you the exact DAG of jobs metamorph will build from your sample sheet, without executing a single one:

```bash
./metamorph run \
    --samplesheet /data/$USER/metamorph_test/samplesheet.tsv \
    --output /data/$USER/metamorph_test/output \
    --mode local \
    --sif-cache /data/OpenOmics/SIFs \
    --dry-run
```

`--mode local` runs everything serially on the current machine instead of submitting to SLURM — this is what makes it practical to iterate on a test dataset like this one from an interactive node. Swap in `--mode slurm` (the default) once you're ready to run at scale on your own data.

## 6. Walking through each mode of execution

metamorph's `--workflow` flag controls how far the DAG extends past initial QC. All four modes below were generated from the *same* sample sheet above — only `--workflow` changes.

### 6.1 `pre-screen` (default)

```bash
./metamorph run --samplesheet samplesheet.tsv --output output_prescreen \
    --mode local --workflow pre-screen --sif-cache /data/OpenOmics/SIFs --dry-run
```

`pre-screen` performs quality control, host-read removal, and taxonomic pre-screening with Centrifuger. For our 2-sample sheet this resolves to:

```text
job                              count
-----------------------------  -------
all                                  1
bowtie2_dehost                       2
dna_centrifuger                      2
metawrap_read_qc_skipBmtagger        2
total                                 7
```

Per sample, that's three real steps chained together:

1. **`metawrap_read_qc_skipBmtagger`** — trims adapters and low-quality bases (Trim-Galore/Cutadapt via metaWRAP's `read_qc` module) and runs FastQC before/after trimming.
2. **`bowtie2_dehost`** — aligns trimmed reads to the hg38 host genome and keeps only the read pairs where *neither* mate mapped, removing host contamination. This is the step that specifically needs the hg38-spiked reads in our test set to do anything meaningful — a purely microbial test file would pass through this step with (almost) nothing removed.
3. **`dna_centrifuger`** — classifies the dehosted reads taxonomically against a GTDB+RefSeq index, producing a per-sample classification and quantification report.

This is the cheapest, fastest mode, and a good default for "does my data look reasonable" checks before committing to the more expensive modes below.

### 6.2 `read-based`

```bash
./metamorph run --samplesheet samplesheet.tsv --output output_readbased \
    --mode local --workflow read-based --sif-cache /data/OpenOmics/SIFs --dry-run
```

`read-based` runs everything in `pre-screen`, then adds functional and taxonomic profiling on the same dehosted reads using HUMAnN3/MetaPhlAn4:

```text
job                                 count
--------------------------------  -------
all                                     1
bowtie2_dehost                          2
dna_centrifuger                         2
dna_humann_classify                     2
dna_humann_diversity_calculation        1
dna_humann_summarize                    1
metawrap_read_qc_skipBmtagger           2
total                                  11
```

`dna_humann_classify` runs per sample; `dna_humann_diversity_calculation` and `dna_humann_summarize` each run once across *all* samples in the sheet — they merge per-sample gene family/pathway abundance tables and compute diversity metrics (Shannon, richness, Jaccard, UniFrac), which is why those two only appear once regardless of sample count. Use `--shallow-profile` here if you want to skip HUMAnN3's translated-search step for a faster (if less sensitive) run.

### 6.3 `assembly-based`

```bash
./metamorph run --samplesheet samplesheet.tsv --output output_assembly \
    --mode local --workflow assembly-based --sif-cache /data/OpenOmics/SIFs --dry-run
```

`assembly-based` skips HUMAnN3/MetaPhlAn4 entirely and instead assembles, bins, and dereplicates each sample's dehosted reads into candidate genomes:

```text
job                              count
-----------------------------  -------
all                                  1
bbtools_index_map                    2
bin_stats                            2
bowtie2_dehost                       2
contig_annotation                    1
cumulative_bin_stats                 1
derep_bins                           1
dna_centrifuger                      2
dna_decompress_dehost_reads          2
gtdbtk_classify                      1
gunc_detection                       1
metawrap_binning                     2
metawrap_genome_assembly             2
metawrap_read_qc_skipBmtagger        2
prep_genome_info                     1
total                               23
```

This is the most expensive of the four modes, and the job count more than triples versus `pre-screen`. Notably:

- **`dna_decompress_dehost_reads`** reappears here (and *only* here) — MEGAHIT/metaSPAdes, run inside `metawrap_genome_assembly`, need uncompressed FASTQ, so metamorph decompresses the dehosted `.fastq.gz` reads specifically for this mode rather than keeping every intermediate file uncompressed by default.
- **`metawrap_genome_assembly`** runs per sample (MEGAHIT and/or metaSPAdes, selectable via `--assembler`), followed by per-sample **`metawrap_binning`** (MetaBAT2, MaxBin2, and metaWRAP's own consensus binner).
- **`derep_bins`**, **`prep_genome_info`**, **`gtdbtk_classify`**, **`gunc_detection`**, **`contig_annotation`**, and **`cumulative_bin_stats`** each run once across the whole cohort — dereplication (dRep) and taxonomic classification (GTDB-Tk) only make sense evaluated across all samples' bins together, which is also why the job count for these doesn't scale with sample count the way per-sample steps do.

### 6.4 `combined`

```bash
./metamorph run --samplesheet samplesheet.tsv --output output_combined \
    --mode local --workflow combined --sif-cache /data/OpenOmics/SIFs --dry-run
```

`combined` is the union of all three modes above — pre-screening, read-based profiling, *and* assembly-based binning, in one DAG:

```text
job                                 count
--------------------------------  -------
all                                     1
bbtools_index_map                       2
bin_stats                               2
bowtie2_dehost                          2
contig_annotation                       1
cumulative_bin_stats                    1
derep_bins                              1
dna_centrifuger                         2
dna_decompress_dehost_reads             2
dna_humann_classify                     2
dna_humann_diversity_calculation        1
dna_humann_summarize                    1
gtdbtk_classify                         1
gunc_detection                          1
metawrap_binning                        2
metawrap_genome_assembly                2
metawrap_read_qc_skipBmtagger           2
prep_genome_info                        1
total                                  27
```

Each of the shared upstream steps (`metawrap_read_qc_skipBmtagger`, `bowtie2_dehost`) is computed exactly once and reused by both the read-based and assembly-based branches — Snakemake's DAG deduplicates on output file, so choosing `combined` over running `read-based` and `assembly-based` separately saves the cost of redoing QC and dehosting twice.

!!! note "About this environment"
    The dehosting and classification steps above depend on reference indices (hg38 Bowtie2/STAR indices, a GTDB+RefSeq Centrifuger index, HUMAnN3/MetaPhlAn4/GTDB-Tk/GUNC databases) that live outside the pipeline's source code. On Biowulf these are already in place; elsewhere, run [`metamorph install`](usage/install.md) first. Everything shown above is a `--dry-run`, so it reflects the job graph metamorph *would* execute — drop `--dry-run` once your reference bundle is in place to actually run it.

## 7. Inspecting the output directory

Regardless of mode, `--output` is a self-contained Snakemake working directory. After a run you'll find:

- `config.json` — the fully-resolved configuration generated from your sample sheet and CLI flags. Worth reading if a run behaves unexpectedly.
- `metagenome_results/metawrap_read_qc/<sample>/` — pre/post-trim FastQC reports.
- `metagenome_results/trimmed_reads/<sample>/` — trimmed and dehosted FASTQ.
- `metagenome_results/centrifuger_dna/` — taxonomic classification/quantification tables.
- `metagenome_results/humann3_dna/` (`read-based`/`combined` only) — gene family, pathway abundance/coverage, and diversity tables.
- `metagenome_results/mags/` (`assembly-based`/`combined` only) — dereplicated, taxonomically classified metagenome-assembled genomes and their per-sample statistics.

## 8. Re-running and cleaning up

If a run is interrupted (killed job, cancelled allocation), Snakemake leaves the output directory locked. Use [`metamorph unlock`](usage/unlock.md) before resubmitting to the same `--output` directory:

```bash
./metamorph unlock --output /data/$USER/metamorph_test/output
```

Because rule outputs are content-addressed by path, simply re-running the same `metamorph run` command afterward resumes from wherever the pipeline left off — completed samples/steps are not redone.
