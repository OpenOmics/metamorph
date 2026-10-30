# Introduction

## 1. About

This page is a conceptual introduction to metamorph. It is meant to be read before the [command usage](../usage/run.md) pages or the [walkthrough](../walkthrough.md) — it explains *why* the pipeline is built the way it is, rather than *how* to invoke any one sub command.

metamorph orchestrates a long chain of specialized, mostly decades-old bioinformatics tools (Trim-Galore, FastQC, Bowtie2/STAR, metaWRAP, MEGAHIT/metaSPAdes, MetaBAT2/MaxBin2, dRep, GTDB-Tk, GUNC, Centrifuger, HUMAnN3, MetaPhlAn4) behind a single, consistent interface. None of those tools were designed to interoperate, and metagenomics/metatranscriptomics data has properties that make gluing them together harder than a typical single-sample, single-reference genomics pipeline. Understanding those properties up front makes the rest of the documentation, and the pipeline's design choices, much easier to follow.

## 2. Portability

metamorph is designed to run unmodified on very different systems: a shared HPC cluster with a SLURM scheduler and a site-wide reference bundle, a single workstation with no scheduler at all, and (eventually) other cluster backends. A few design decisions make that possible:

- **Containers, not modules.** Every rule that shells out to third-party software (`containerized:` / `singularity:` in the Snakemake rules) runs inside a Singularity container pulled from a pinned image tag. This means the pipeline does not depend on whatever version of `bowtie2`, `metaspades`, or `humann3` happens to be installed on a given system, and a result produced on one cluster is reproducible on another. See [`metamorph cache`](../usage/cache.md) for pre-fetching those containers into a shared, offline SIF cache — useful for avoiding DockerHub rate limits and for air-gapped systems.
- **A pluggable execution backend.** `metamorph run --mode {slurm,local}` swaps the *orchestration* layer without changing a single rule. `slurm` submits each job through `sbatch` for distributed execution; `local` runs the same DAG serially on one machine, which is what makes it practical to test the pipeline on a laptop or an interactive node with the tiny FASTQ files described in the [walkthrough](../walkthrough.md).
- **Externalized reference data.** Large reference databases (host genome indices, GTDB/RefSeq taxonomy databases, HUMAnN3/MetaPhlAn4 databases) are never bundled with the pipeline's source code — they are pulled once with [`metamorph install`](../usage/install.md) into a `--resource-bundle` / hardcoded reference path. On systems that already host a copy (e.g. NIH's Biowulf), users skip this step entirely.
- **A generated, self-describing config.** The `run` sub command resolves the sample sheet, CLI flags, and the JSON files under `config/` into a single `config.json` written into the output directory. That output directory is a complete, standalone Snakemake working directory — it can be linted, resumed, or unlocked ([`metamorph unlock`](../usage/unlock.md)) independently of the code that generated it.

Portability is what lets the *same* sample sheet and the *same* `metamorph run` invocation produce the same directed acyclic graph of jobs whether it's submitted to a thousand-node cluster or dry-run on a single core.

## 3. Why metagenomics/metatranscriptomics needs a pipeline like this

A handful of properties specific to metagenomic and metatranscriptomic data make this domain notably harder to build a portable, reproducible pipeline for than single-organism, single-reference sequencing analysis:

**There is no single reference genome.** A human or mouse pipeline can align every read to one reference. A metagenome is an unknown mixture of many, often uncharacterized, organisms. metamorph never assumes a reference — it *assembles* each sample de novo (MEGAHIT and/or metaSPAdes, selectable via `--assembler`), bins the resulting contigs into candidate genomes (MetaBAT2, MaxBin2, and metaWRAP's consensus binning), and only then classifies and dereplicates ("MAGs") — a fundamentally different, and much more compute-intensive, path than reference-based alignment.

**Host contamination has to be actively removed, not just ignored.** Stool, tissue, and most other sample types contain host DNA/RNA alongside the microbial signal of interest. If it isn't removed, host reads can dominate assemblies, waste compute, and are a re-identification/privacy risk when host is human. Every DNA and RNA sample pair in metamorph passes through an explicit dehosting step (`bowtie2_dehost` for DNA against hg38, STAR-based dehosting for RNA) *before* anything taxonomic or functional is computed downstream.

**The tools involved have inconsistent, sometimes contradictory input expectations.** metaWRAP's `read_qc` module, for example, only accepts uncompressed FASTQ, while every other stage in the pipeline standardizes on gzip-compressed FASTQ to keep the working directory's storage footprint manageable. metamorph absorbs that inconsistency internally (decompressing on the way in, recompressing the trimmed output on the way out) so it never leaks into the sample sheet or the user-facing interface. Multiply this kind of quirk across a dozen tools and the value of centralizing it in one maintained pipeline, instead of every user re-solving it themselves, becomes clear.

**Reference databases are enormous and shared, not per-sample.** Host indices, GTDB-Tk's taxonomy database, and Centrifuger/HUMAnN3/MetaPhlAn4's classification databases are commonly hundreds of gigabytes, far larger than any individual FASTQ input. These are treated as external, site-provided resources ([`metamorph install`](../usage/install.md)) rather than pipeline inputs, so that a run's cost scales with sample size, not with the size of the reference collection it's compared against.

**Multi-omic sample pairing is optional but load-bearing when present.** A project may contain DNA only, or matched DNA + RNA per subject. When RNA is present, metamorph pairs it with its corresponding DNA sample throughout the assembly-based workflow so functional profiling (HUMAnN3) and MAG-level RNA mapping are computed against the *same* subject's assembled genomes — not a generic reference — which is what makes the metatranscriptomic results interpretable at all.

**Runtime and resource needs vary by two orders of magnitude across stages.** Read QC and dehosting finish in minutes; genome assembly, binning, and dereplication can take hours to days and need substantially more memory. The `--workflow {pre-screen,read-based,assembly-based,combined}` switch exists specifically so users can stop after the cheap, fast stages (`pre-screen`) when that's all a project needs, rather than always paying for the most expensive stages.

The rest of the manual and the [walkthrough](../walkthrough.md) build on these ideas concretely — showing what actually gets produced at each stage, using a small, public, real dataset.
