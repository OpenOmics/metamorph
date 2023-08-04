<div align="center">

  <h1 style="font-size: 250%">metamorph 🔬</h1>

  <b><i>An awesome metagenomic and metatranscriptomics pipeline</i></b><br> 
  <a href="https://github.com/OpenOmics/metamorph/actions/workflows/main.yaml">
    <img alt="tests" src="https://github.com/OpenOmics/metamorph/workflows/tests/badge.svg">
  </a>
  <a href="https://github.com/OpenOmics/metamorph/actions/workflows/docs.yml">
    <img alt="docs" src="https://github.com/OpenOmics/metamorph/workflows/docs/badge.svg">
  </a>
  <a href="https://github.com/OpenOmics/metamorph/issues">
    <img alt="GitHub issues" src="https://img.shields.io/github/issues/OpenOmics/metamorph?color=brightgreen">
  </a>
  <a href="https://github.com/OpenOmics/metamorph/blob/main/LICENSE">
    <img alt="GitHub license" src="https://img.shields.io/github/license/OpenOmics/metamorph">
  </a>

  <p>
    This is the home of the pipeline, metamorph. Its long-term goals: to provide accurate quantification, taxonomic classification, and functional profiling of assembled (bacteria and archaea) metagenomes!
  </p>

</div>  


## Overview
Welcome to metamorph's documentation! This guide is the main source of documentation for users that are getting started with the [long pipeline name](https://github.com/OpenOmics/metamorph/). 

The **`./metamorph`** pipeline is composed several inter-related sub commands to setup and run the pipeline across different systems. Each of the available sub commands perform different functions: 

<section align="center" markdown="1" style="display: flex; flex-wrap: row wrap; justify-content: space-around;">

!!! inline custom-grid-button ""

    [<code style="font-size: 1em;">metamorph <b>run</b></code>](usage/run.md)   
    Run the metamorph pipeline with your input files.

!!! inline custom-grid-button ""

    [<code style="font-size: 1em;">metamorph <b>unlock</b></code>](usage/unlock.md)  
    Unlocks a previous runs output directory.

</section>

<section align="center" markdown="1" style="display: flex; flex-wrap: row wrap; justify-content: space-around;">


!!! inline custom-grid-button ""

    [<code style="font-size: 1em;">metamorph <b>install</b></code>](usage/install.md)  
    Download remote reference files locally.


!!! inline custom-grid-button ""

    [<code style="font-size: 1em;">metamorph <b>cache</b></code>](usage/cache.md)  
    Cache remote software containers locally.  

</section>

**metamorph** is a comprehensive workflow that starts off with assembly of short read DNA sequencing data (metagenome) to build contigs for each individual sample, followed by contig binning. Good quality genomic bins from all samples are futher combined and dereplicated to generate representative metagenome-assembled genomes or MAGs. These MAGs are further classified, quantified, and used to predict functional profiles. The short read DNA sequencing data (metatranscriptome) mapping to MAGs is used for functional profiling. It relies on technologies like [Singularity<sup>1</sup>](https://singularity.lbl.gov/) to maintain the highest-level of reproducibility. The pipeline consists of a series of data processing and quality-control steps orchestrated by [Snakemake<sup>2</sup>](https://snakemake.readthedocs.io/en/stable/), a flexible and scalable workflow management system, to submit jobs to a cluster.

The pipeline is compatible with data generated from Illumina short-read sequencing technologies. As input, it accepts a set of FastQ files and can be run locally on a compute instance or on-premise using a cluster. A user can define the method or mode of execution. The pipeline can submit jobs to a cluster using a job scheduler like SLURM (more coming soon!). A hybrid approach ensures the pipeline is accessible to all users.

Before getting started, we highly recommend reading through the [usage](usage/run.md) section of each available sub command.

For more information about issues or trouble-shooting a problem, please checkout our [FAQ](faq/questions.md) prior to [opening an issue on Github](https://github.com/OpenOmics/metamorph/issues).

## Contribute 

This site is a living document, created for and by members like you. metamorph is maintained by the members of NCBR and is improved by continous feedback! We encourage you to contribute new content and make improvements to existing content via pull request to our [GitHub repository :octicons-heart-fill-24:{ .heart }](https://github.com/OpenOmics/metamorph).

## Citation

If you use this software, please cite it as below:  

=== "BibTex"

    ```
    Citation coming soon!
    ```

=== "APA"

    ```
    Citation coming soon!
    ```

## References
<sup>**1.**  Kurtzer GM, Sochat V, Bauer MW (2017). Singularity: Scientific containers for mobility of compute. PLoS ONE 12(5): e0177459.</sup>  
<sup>**2.**  Koster, J. and S. Rahmann (2018). "Snakemake-a scalable bioinformatics workflow engine." Bioinformatics 34(20): 3600.</sup>  
