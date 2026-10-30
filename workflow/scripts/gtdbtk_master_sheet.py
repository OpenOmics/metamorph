#!/usr/bin/env python3
# -*- coding: UTF-8 -*-

"""
ABOUT:
    Builds one assembly-wide mastersheet from gtdbtk_classify's summary
    (one row per dereplicated bin, with its GTDB lineage), derep_bins'
    genomeInfo.csv (completeness/contamination, as scored by CheckM
    before dereplication), each bin's own dereplicated FASTA (genome
    size and contig count, counted directly from the sequence dRep
    picked as that cluster's representative), and bbtools_index_map's
    per-sample refstats files (the fraction of that sample's reads
    recruited to it). A bin missing from a given sample's refstats file
    genuinely recruited zero reads from that sample -- BBMap's refstats
    only lists references with at least one assigned read -- so it's
    filled in as 0 rather than treated as missing data.

USAGE:
    gtdbtk_master_sheet.py \\
        --gtdbtk-summary GTDBTK_classify_wf/gtdbtk.bac120.summary.tsv \\
        --genome-info metawrap_bin_refine/genomeInfo.csv \\
        --dereplicated-genomes-dir metawrap_bin_refine/dRep/dereplicated_genomes \\
        --refstats sample1=mags/sample1/DNA/sample1.refstats \\
        --refstats sample2=mags/sample2/DNA/sample2.refstats \\
        --output assembly_master_sheet.tsv
"""

import argparse
import os
import sys


RANKS = ["kingdom", "phylum", "class", "order", "family", "genus", "species"]
RANK_PREFIX = {"kingdom": "d", "phylum": "p", "class": "c", "order": "o", "family": "f", "genus": "g", "species": "s"}


def err(msg):
    print(f"[gtdbtk_master_sheet] ERROR: {msg}", file=sys.stderr)
    sys.exit(1)


def load_gtdbtk_summary(path):
    """Returns (ordered list of bin IDs, {bin_id: classification string})."""
    bins = []
    lineages = {}
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        try:
            bin_idx = header.index("user_genome")
            class_idx = header.index("classification")
        except ValueError:
            err(f"{path} is missing a 'user_genome' or 'classification' column")
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            fields = line.split("\t")
            bin_id = fields[bin_idx]
            bins.append(bin_id)
            lineages[bin_id] = fields[class_idx]
    return bins, lineages


def split_lineage(classification):
    """Returns {rank: name}, keyed by RANKS, parsed out of a d__;p__;...;s__ string."""
    by_prefix = {}
    for token in classification.split(";"):
        token = token.strip()
        if "__" not in token:
            continue
        prefix, name = token.split("__", 1)
        by_prefix[prefix] = name
    return {rank: by_prefix.get(RANK_PREFIX[rank], "") for rank in RANKS}


def load_genome_info(path):
    """Returns {bin_id: (completeness, contamination)} from derep_bins' genomeInfo.csv."""
    info = {}
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split(",")
        try:
            genome_idx = header.index("genome")
            comp_idx = header.index("completeness")
            cont_idx = header.index("contamination")
        except ValueError:
            err(f"{path} is missing a 'genome', 'completeness', or 'contamination' column")
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            fields = line.split(",")
            bin_id = fields[genome_idx]
            if bin_id.endswith(".fa"):
                bin_id = bin_id[: -len(".fa")]
            info[bin_id] = (fields[comp_idx], fields[cont_idx])
    return info


def fasta_stats(path):
    """Returns (num_contigs, genome_size) counted directly from a bin's FASTA."""
    num_contigs = 0
    genome_size = 0
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                num_contigs += 1
            else:
                genome_size += len(line)
    return num_contigs, genome_size


def load_refstats(path):
    """Returns {bin_id: %unambiguousReads string} from one bbtools refstats file."""
    fractions = {}
    with open(path) as fh:
        header = fh.readline().rstrip("\n").lstrip("#").split("\t")
        try:
            name_idx = header.index("name")
            frac_idx = header.index("%unambiguousReads")
        except ValueError:
            err(f"{path} is missing a 'name' or '%unambiguousReads' column")
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            fields = line.split("\t")
            fractions[fields[name_idx]] = fields[frac_idx]
    return fractions


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--gtdbtk-summary", required=True)
    p.add_argument("--genome-info", required=True, help="derep_bins' genomeInfo.csv")
    p.add_argument("--dereplicated-genomes-dir", required=True, help="dRep's dereplicated_genomes/ directory")
    p.add_argument("--refstats", action="append", required=True, metavar="SAMPLE=PATH",
                    help="repeatable; one per sample")
    p.add_argument("--output", required=True)
    args = p.parse_args()

    bins, lineages = load_gtdbtk_summary(args.gtdbtk_summary)
    genome_info = load_genome_info(args.genome_info)

    samples = []
    sample_fractions = {}
    for entry in args.refstats:
        if "=" not in entry:
            err(f"--refstats entry '{entry}' is not of the form SAMPLE=PATH")
        sample, path = entry.split("=", 1)
        samples.append(sample)
        sample_fractions[sample] = load_refstats(path)

    columns = ["bin_id", "lineage"] + RANKS + ["completeness", "contamination", "genome_size", "num_contigs"] + samples
    with open(args.output, "w") as out:
        out.write("\t".join(columns) + "\n")
        for bin_id in bins:
            classification = lineages[bin_id]
            ranks = split_lineage(classification)
            row = [bin_id, classification] + [ranks[rank] for rank in RANKS]

            if bin_id not in genome_info:
                err(f"{args.genome_info} has no completeness/contamination row for bin '{bin_id}'")
            completeness, contamination = genome_info[bin_id]
            fasta_path = os.path.join(args.dereplicated_genomes_dir, bin_id + ".fa")
            if not os.path.isfile(fasta_path):
                err(f"dereplicated genome FASTA not found for bin '{bin_id}': {fasta_path}")
            num_contigs, genome_size = fasta_stats(fasta_path)
            row += [completeness, contamination, str(genome_size), str(num_contigs)]

            row += [sample_fractions[sample].get(bin_id, "0.00000") for sample in samples]
            out.write("\t".join(row) + "\n")


if __name__ == "__main__":
    main()
