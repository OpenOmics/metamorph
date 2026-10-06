#!/usr/bin/env python3
# -*- coding: UTF-8 -*-

"""
ABOUT:
    Builds one row of a per-sample QC/provenance summary by joining
    measurements already produced by earlier DNA.smk rules: raw/trimmed
    read counts and module statuses from FastQC's own HTML reports,
    dehosting pair counts from bowtie2_dehost's/screen_injected_host_dna's
    summary text files, and classifier input/output and domain-level
    breakdown from Centrifuger's classification + quantification report.

    Nothing here re-derives a number a prior rule already computed except
    where that number is itself the point of a reconciliation check (e.g.
    independently counting the dehosted FASTQs to confirm they match what
    Centrifuger actually consumed, rather than assuming the DAG wiring is
    correct).

USAGE:
    sample_qc_summary.py \\
        --sample-id NAME --library-type DNA --workflow combined \\
        --pipeline-version 0.1.0-beta --pipeline-git-commit <hash> \\
        --fastqc-pretrim-r1 R1_pre.html --fastqc-pretrim-r2 R2_pre.html \\
        --fastqc-posttrim-r1 R1_post.html --fastqc-posttrim-r2 R2_post.html \\
        --hg38-dehost-summary NAME_hg38_dehost_summary.txt \\
        [--injected-host-dehost-summary NAME_injected_host_dehost_summary.txt] \\
        --classification-tsv NAME_centrifuger_classification.tsv \\
        --quant-report NAME_centrifuger_quantification_report.tsv \\
        --dehosted-r1-count 16131390 --dehosted-r2-count 16131390 \\
        --max-host-pct 20 --min-overall-retained-pct 50 \\
        --output NAME_qc_summary.tsv
"""

import argparse
import os
import re
import sys


FASTQC_MODULES = [
    ("Per base sequence quality", "per_base_quality"),
    ("Adapter Content", "adapter_content"),
    ("Sequence Duplication Levels", "duplication"),
    ("Overrepresented sequences", "overrepresented_sequences"),
]

COLUMNS = [
    "sample_id", "library_type", "workflow", "pipeline_version", "pipeline_git_commit",
    "raw_read_pairs",
    "posttrim_read_pairs", "posttrim_retained_pct",
    "dehosted_read_pairs", "dehosted_retained_pct", "overall_retained_pct",
    "host_alignment_pct",
    "classifier_input_fragments", "classified_fragments", "classified_pct",
    "bacterial_fragments", "bacterial_pct",
    "archaeal_fragments", "archaeal_pct",
    "eukaryotic_fragments", "eukaryotic_pct",
    "viral_fragments", "viral_pct",
    "host_fragments", "host_pct",
    "unclassified_fragments", "unclassified_pct",
    "fastqc_pretrim_per_base_quality", "fastqc_pretrim_adapter_content",
    "fastqc_pretrim_duplication", "fastqc_pretrim_overrepresented_sequences",
    "fastqc_posttrim_per_base_quality", "fastqc_posttrim_adapter_content",
    "fastqc_posttrim_duplication", "fastqc_posttrim_overrepresented_sequences",
    "flag_r1_r2_mismatch", "flag_unexpected_count_increase",
    "flag_classifier_input_mismatch", "flag_missing_or_empty_outputs",
    "flag_high_host_content", "flag_low_retained_reads",
]

STATUS_RANK = {"PASS": 0, "WARN": 1, "WARNING": 1, "FAIL": 2}
STATUS_LABEL = {0: "PASS", 1: "WARN", 2: "FAIL"}


def err(msg):
    print(f"[sample_qc_summary] ERROR: {msg}", file=sys.stderr)
    sys.exit(1)


def parse_fastqc(html_path):
    """Returns (total_sequences:int, {module_name: status}) from one FastQC HTML report."""
    if not os.path.isfile(html_path) or os.path.getsize(html_path) == 0:
        return None, {}
    with open(html_path, "r", errors="replace") as fh:
        html = fh.read()
    m = re.search(r"<td>Total Sequences</td><td>(\d+)</td>", html)
    total = int(m.group(1)) if m else None
    statuses = {}
    for status, name in re.findall(
        r'<li><img[^>]*alt="\[(PASS|WARN|WARNING|FAIL)\]"[^>]*>\s*<a href="#\w+">([^<]+)</a></li>',
        html,
    ):
        statuses[name.strip()] = "WARN" if status == "WARNING" else status
    return total, statuses


def worst_status(a, b):
    """Worst-of-both-mates: FAIL > WARN > PASS. Missing data reports as 'NA'."""
    ranks = [STATUS_RANK[s] for s in (a, b) if s in STATUS_RANK]
    if not ranks:
        return "NA"
    return STATUS_LABEL[max(ranks)]


def parse_dehost_summary(path):
    """Returns (input_pairs, removed_pairs, retained_pairs) from a
    bowtie2_dehost/screen_injected_host_dna *_dehost_summary.txt file."""
    if not path or not os.path.isfile(path) or os.path.getsize(path) == 0:
        return None, None, None
    text = open(path).read()
    def grab(label):
        m = re.search(rf"{label}:\s*(\d+)", text)
        return int(m.group(1)) if m else None
    return grab("Input read pairs"), grab("Removed \\(mapped\\)"), grab("Retained \\(unmapped\\)")


def count_tsv_data_rows(path):
    """Line count minus header, for a plain-text TSV. Returns None if missing/empty."""
    if not path or not os.path.isfile(path) or os.path.getsize(path) == 0:
        return None
    n = -1  # subtract header
    with open(path) as fh:
        for n, _ in enumerate(fh):
            pass
    return n if n >= 0 else 0


def parse_quant_report(path):
    """Returns dict with classified_fragments (root numReads) and the
    domain-level (+ Viruses, rank 'acellular root') numReads breakdown."""
    out = {"classified": None, "Bacteria": 0, "Archaea": 0, "Eukaryota": 0, "Viruses": 0}
    if not path or not os.path.isfile(path) or os.path.getsize(path) == 0:
        return out
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        idx = {col: i for i, col in enumerate(header)}
        for line in fh:
            fields = line.rstrip("\n").split("\t")
            if len(fields) <= max(idx.values()):
                continue
            name = fields[idx["name"]]
            tax_id = fields[idx["taxID"]]
            tax_rank = fields[idx["taxRank"]]
            num_reads = int(fields[idx["numReads"]])
            if tax_id == "1" and tax_rank == "no rank" and name == "root":
                out["classified"] = num_reads
            elif tax_rank == "domain" and name in ("Bacteria", "Archaea", "Eukaryota"):
                out[name] = num_reads
            elif tax_rank == "acellular root" and name == "Viruses":
                out["Viruses"] = num_reads
    return out


def pct(numerator, denominator):
    if numerator is None or denominator is None or denominator == 0:
        return None
    return round(numerator / denominator * 100.0, 4)


def fmt(value):
    return "NA" if value is None else str(value)


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--sample-id", required=True)
    p.add_argument("--library-type", required=True)
    p.add_argument("--workflow", required=True)
    p.add_argument("--pipeline-version", required=True)
    p.add_argument("--pipeline-git-commit", required=True)
    p.add_argument("--fastqc-pretrim-r1", required=True)
    p.add_argument("--fastqc-pretrim-r2", required=True)
    p.add_argument("--fastqc-posttrim-r1", required=True)
    p.add_argument("--fastqc-posttrim-r2", required=True)
    p.add_argument("--hg38-dehost-summary", required=True)
    p.add_argument("--injected-host-dehost-summary", default=None)
    p.add_argument("--classification-tsv", required=True)
    p.add_argument("--quant-report", required=True)
    p.add_argument("--dehosted-r1-count", type=int, required=True)
    p.add_argument("--dehosted-r2-count", type=int, required=True)
    p.add_argument("--max-host-pct", type=float, default=20.0)
    p.add_argument("--min-overall-retained-pct", type=float, default=50.0)
    p.add_argument("--output", required=True)
    args = p.parse_args()

    missing = []

    def check(label, path):
        if not path or not os.path.isfile(path) or os.path.getsize(path) == 0:
            missing.append(label)

    for label, path in [
        ("fastqc_pretrim_r1", args.fastqc_pretrim_r1), ("fastqc_pretrim_r2", args.fastqc_pretrim_r2),
        ("fastqc_posttrim_r1", args.fastqc_posttrim_r1), ("fastqc_posttrim_r2", args.fastqc_posttrim_r2),
        ("hg38_dehost_summary", args.hg38_dehost_summary),
        ("classification_tsv", args.classification_tsv), ("quant_report", args.quant_report),
    ]:
        check(label, path)
    if args.injected_host_dehost_summary is not None:
        check("injected_host_dehost_summary", args.injected_host_dehost_summary)

    raw_r1, pre_status = parse_fastqc(args.fastqc_pretrim_r1)
    raw_r2, pre_status_r2 = parse_fastqc(args.fastqc_pretrim_r2)
    post_r1, post_status = parse_fastqc(args.fastqc_posttrim_r1)
    post_r2, post_status_r2 = parse_fastqc(args.fastqc_posttrim_r2)

    raw_pairs = raw_r1
    posttrim_pairs = post_r1

    hg38_in, hg38_removed, hg38_retained = parse_dehost_summary(args.hg38_dehost_summary)
    if args.injected_host_dehost_summary:
        _, inj_removed, inj_retained = parse_dehost_summary(args.injected_host_dehost_summary)
        dehosted_pairs = inj_retained
        total_host_removed = (hg38_removed or 0) + (inj_removed or 0)
    else:
        dehosted_pairs = hg38_retained
        total_host_removed = hg38_removed

    classifier_input = count_tsv_data_rows(args.classification_tsv)
    quant = parse_quant_report(args.quant_report)
    classified = quant["classified"]
    unclassified = (classifier_input - classified) if (classifier_input is not None and classified is not None) else None

    posttrim_retained_pct = pct(posttrim_pairs, raw_pairs)
    dehosted_retained_pct = pct(dehosted_pairs, posttrim_pairs)
    overall_retained_pct = pct(dehosted_pairs, raw_pairs)
    host_alignment_pct = None if posttrim_pairs in (None, 0) else round(100.0 - (dehosted_retained_pct or 0), 4)
    classified_pct = pct(classified, classifier_input)
    host_pct = pct(total_host_removed, raw_pairs)
    bacterial_pct = pct(quant["Bacteria"], raw_pairs)
    archaeal_pct = pct(quant["Archaea"], raw_pairs)
    eukaryotic_pct = pct(quant["Eukaryota"], raw_pairs)
    viral_pct = pct(quant["Viruses"], raw_pairs)
    unclassified_pct = pct(unclassified, raw_pairs)

    # --- reconciliation flags ---
    r1_r2_mismatches = []
    if raw_r1 is not None and raw_r2 is not None and raw_r1 != raw_r2:
        r1_r2_mismatches.append("raw")
    if post_r1 is not None and post_r2 is not None and post_r1 != post_r2:
        r1_r2_mismatches.append("posttrim")
    if args.dehosted_r1_count != args.dehosted_r2_count:
        r1_r2_mismatches.append("dehosted")
    flag_r1_r2_mismatch = ",".join(r1_r2_mismatches) if r1_r2_mismatches else "none"

    increases = []
    stage_counts = [("raw", raw_pairs), ("posttrim", posttrim_pairs), ("dehosted", dehosted_pairs),
                     ("classifier_input", classifier_input)]
    for (prev_name, prev_val), (cur_name, cur_val) in zip(stage_counts, stage_counts[1:]):
        if prev_val is not None and cur_val is not None and cur_val > prev_val:
            increases.append(f"{cur_name}>{prev_name}")
    flag_unexpected_count_increase = ",".join(increases) if increases else "none"

    flag_classifier_input_mismatch = "true" if (
        classifier_input is not None and dehosted_pairs is not None and classifier_input != dehosted_pairs
    ) or (args.dehosted_r1_count != args.dehosted_r2_count) else "false"
    # also cross-check against the independently-counted dehosted FASTQ pairs
    if classifier_input is not None and args.dehosted_r1_count and classifier_input != args.dehosted_r1_count:
        flag_classifier_input_mismatch = "true"

    flag_missing_or_empty_outputs = ",".join(missing) if missing else "none"
    flag_high_host_content = "true" if (host_pct is not None and host_pct > args.max_host_pct) else "false"
    flag_low_retained_reads = "true" if (
        overall_retained_pct is not None and overall_retained_pct < args.min_overall_retained_pct
    ) else "false"

    row = {
        "sample_id": args.sample_id,
        "library_type": args.library_type,
        "workflow": args.workflow,
        "pipeline_version": args.pipeline_version,
        "pipeline_git_commit": args.pipeline_git_commit,
        "raw_read_pairs": raw_pairs,
        "posttrim_read_pairs": posttrim_pairs,
        "posttrim_retained_pct": posttrim_retained_pct,
        "dehosted_read_pairs": dehosted_pairs,
        "dehosted_retained_pct": dehosted_retained_pct,
        "overall_retained_pct": overall_retained_pct,
        "host_alignment_pct": host_alignment_pct,
        "classifier_input_fragments": classifier_input,
        "classified_fragments": classified,
        "classified_pct": classified_pct,
        "bacterial_fragments": quant["Bacteria"],
        "bacterial_pct": bacterial_pct,
        "archaeal_fragments": quant["Archaea"],
        "archaeal_pct": archaeal_pct,
        "eukaryotic_fragments": quant["Eukaryota"],
        "eukaryotic_pct": eukaryotic_pct,
        "viral_fragments": quant["Viruses"],
        "viral_pct": viral_pct,
        "host_fragments": total_host_removed,
        "host_pct": host_pct,
        "unclassified_fragments": unclassified,
        "unclassified_pct": unclassified_pct,
        "fastqc_pretrim_per_base_quality": worst_status(pre_status.get("Per base sequence quality"), pre_status_r2.get("Per base sequence quality")),
        "fastqc_pretrim_adapter_content": worst_status(pre_status.get("Adapter Content"), pre_status_r2.get("Adapter Content")),
        "fastqc_pretrim_duplication": worst_status(pre_status.get("Sequence Duplication Levels"), pre_status_r2.get("Sequence Duplication Levels")),
        "fastqc_pretrim_overrepresented_sequences": worst_status(pre_status.get("Overrepresented sequences"), pre_status_r2.get("Overrepresented sequences")),
        "fastqc_posttrim_per_base_quality": worst_status(post_status.get("Per base sequence quality"), post_status_r2.get("Per base sequence quality")),
        "fastqc_posttrim_adapter_content": worst_status(post_status.get("Adapter Content"), post_status_r2.get("Adapter Content")),
        "fastqc_posttrim_duplication": worst_status(post_status.get("Sequence Duplication Levels"), post_status_r2.get("Sequence Duplication Levels")),
        "fastqc_posttrim_overrepresented_sequences": worst_status(post_status.get("Overrepresented sequences"), post_status_r2.get("Overrepresented sequences")),
        "flag_r1_r2_mismatch": flag_r1_r2_mismatch,
        "flag_unexpected_count_increase": flag_unexpected_count_increase,
        "flag_classifier_input_mismatch": flag_classifier_input_mismatch,
        "flag_missing_or_empty_outputs": flag_missing_or_empty_outputs,
        "flag_high_host_content": flag_high_host_content,
        "flag_low_retained_reads": flag_low_retained_reads,
    }

    with open(args.output, "w") as out:
        out.write("\t".join(COLUMNS) + "\n")
        out.write("\t".join(fmt(row[c]) for c in COLUMNS) + "\n")


if __name__ == "__main__":
    main()
