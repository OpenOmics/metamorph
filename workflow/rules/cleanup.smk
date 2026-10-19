# ~~~~~~~~~~
# Terminal cleanup rule
# ~~~~~~~~~~
# Removes or compresses "cruft" left behind in the output directory once
# every other rule has finished. The actual list of files to act on lives
# in config/cleanup.json (an extensible manifest) -- nothing in this file
# needs to change to add or remove cleanup targets, only that manifest.

rule cleanup_cruft:
    """
    Guaranteed to run after every other rule: its input is `all_targets`,
    the exact same list of files `rule all` requires, so Snakemake cannot
    schedule this rule until everything else the pipeline produces already
    exists. See config/cleanup.json for the delete/compress/compress_indexed
    manifest, and workflow/scripts/cleanup_cruft.py for the logic that
    applies it.
    """
    input:
        all_targets
    output:
        touch(cruft_sentinel)
    params:
        rname         = "cleanup_cruft",
        # workpath (not the metagenome_results/ subdir) so manifest patterns
        # can reach anything under --output, e.g. dna/, rna/, logfiles/.
        project_dir   = workpath,
        manifest_path = join(workpath, "config", "cleanup.json"),
    threads: int(cluster["cleanup_cruft"].get("threads", default_threads)),
    singularity: metawrap_container,
    shell:
        """
        python3 workflow/scripts/cleanup_cruft.py \\
            --manifest {params.manifest_path} \\
            --project-dir {params.project_dir}
        """
