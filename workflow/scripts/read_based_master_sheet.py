#!/usr/bin/env python3
# -*- coding: UTF-8 -*-

"""
ABOUT:
    Builds one cohort-level mastersheet for the read-based workflow as a
    two-sheet Excel workbook:
      - "centrifuger": one row per taxID from every sample's Centrifuger
        quantification report, its GTDB-style lineage (d__;p__;c__;o__;
        f__;g__;s__) split out into per-rank columns, and that taxon's
        relative abundance in each sample. A taxID missing from a given
        sample's report genuinely wasn't detected there, so its
        abundance is filled in as 0 rather than treated as missing data.
      - "metaphlan": humann3's merged_bugs_list.tsv (MetaPhlAn4
        clade-level relative abundance per sample), passed through as-is.
    On both sheets, the per-sample abundance columns get a shared
    red(low)-to-green(high) color scale, same as manually selecting that
    block in Excel and applying a Color Scale conditional format.

    These two are kept on separate sheets rather than joined into one
    table: Centrifuger (built against a GTDB r226 taxonomy) and
    MetaPhlAn4 (bundled with its own, differently-versioned taxonomy)
    are different databases that disagree on naming at every rank above
    species -- e.g. GTDB's "Pseudomonadota" vs. MetaPhlAn's
    "Proteobacteria" -- so there's no sound way to merge them into one
    row-per-taxon table without silently implying a correspondence that
    isn't really there.

    The Centrifuger lineage is resolved by walking each taxID's ancestor
    chain in the index's own taxonomy tree. That tree (nodes.dmp) and
    scientific names (names.dmp) aren't shipped alongside the prebuilt
    gtdb_r226+refseq_hvfc index -- they're extracted from the index
    itself via `centrifuger-inspect --taxonomy-tree`/`--name-table`,
    which is how dna_centrifuger_taxonomy produces the two files this
    script reads.

    A rank slot is left empty (e.g. "f__") when no ancestor at that rank
    exists in the chain -- this happens legitimately for ranks GTDB itself
    leaves unassigned, and for the non-GTDB branches of this index (e.g.
    viruses, which root under "acellular root" rather than any domain).

USAGE:
    read_based_master_sheet.py \\
        --nodes nodes.dmp --names names.dmp \\
        --quant-report sample1=sample1_centrifuger_quantification_report.tsv \\
        --quant-report sample2=sample2_centrifuger_quantification_report.tsv \\
        --metaphlan-bugs-list merged_bugs_list.tsv \\
        --output read_based_master_sheet.xlsx
"""

import argparse
import math
import sys

from openpyxl import Workbook
from openpyxl.formatting.rule import ColorScaleRule
from openpyxl.utils import get_column_letter

BRIGHT_GREEN = "FF00FF00"
BRIGHT_YELLOW = "FFFFFF00"
BRIGHT_RED = "FFFF0000"


def _patch_openpyxl_float_formatting():
    """
    openpyxl formats every numeric cell with its own safe_string(), which
    does "%.16g" % value -- one significant digit past IEEE-754 double
    precision. For most non-power-of-two decimals (i.e. almost all real
    abundance values) that leaks binary round-off noise straight into the
    saved XML, e.g. 0.0000794 comes back out as "7.940000000000001e-05".
    Across tens of thousands of cells, real Excel's stricter parser chokes
    on this and throws up a "repair this workbook?" prompt on open, even
    though python's own float repr of the exact same value is clean (it
    uses a shortest-round-trip algorithm instead of a fixed 16 digits).
    Swapping in repr() here for numeric values -- after this module is
    imported, before any save() -- fixes that at the source.
    """
    import openpyxl.cell._writer as _writer

    def clean_safe_string(value):
        if isinstance(value, float):
            if math.isnan(value) or math.isinf(value):
                return ""
            return repr(value)
        return _original_safe_string(value)

    _original_safe_string = _writer.safe_string
    _writer.safe_string = clean_safe_string


_patch_openpyxl_float_formatting()


RANKS = ["domain", "phylum", "class", "order", "family", "genus", "species"]
PREFIX = {rank: rank[0] for rank in RANKS}
METAPHLAN_SAMPLE_SUFFIX = "_bugs_list"


def err(msg):
    print(f"[read_based_master_sheet] ERROR: {msg}", file=sys.stderr)
    sys.exit(1)


def load_nodes(path):
    """Returns {taxid: (parent_taxid, rank)} from a centrifuger-inspect --taxonomy-tree dump."""
    nodes = {}
    with open(path) as fh:
        for line in fh:
            fields = [f.strip() for f in line.rstrip("\n").split("|")]
            if len(fields) < 3 or not fields[0]:
                continue
            taxid, parent, rank = fields[0], fields[1], fields[2]
            nodes[taxid] = (parent, rank)
    return nodes


def load_names(path):
    """Returns {taxid: scientific_name} from a centrifuger-inspect --name-table dump."""
    names = {}
    with open(path) as fh:
        for line in fh:
            fields = [f.strip() for f in line.rstrip("\n").split("|")]
            if len(fields) < 3 or not fields[0]:
                continue
            taxid, name, name_type = fields[0], fields[1], fields[2]
            if name_type == "scientific name":
                names[taxid] = name
    return names


def build_lineage(taxid, nodes, names, cache):
    if taxid in cache:
        return cache[taxid]
    found = {}
    current = taxid
    seen = set()
    while current in nodes and current not in seen:
        seen.add(current)
        parent, rank = nodes[current]
        if rank in PREFIX and rank not in found:
            found[rank] = names.get(current, "")
        if parent == current:
            break
        current = parent
    lineage = ";".join(f"{PREFIX[rank]}__{found.get(rank, '')}" for rank in RANKS)
    cache[taxid] = lineage
    return lineage


def split_lineage(lineage):
    """Returns {rank: name}, keyed by RANKS, parsed out of a d__;p__;...;s__ string."""
    by_prefix = {}
    for token in lineage.split(";"):
        if "__" not in token:
            continue
        prefix, name = token.split("__", 1)
        by_prefix[prefix] = name
    return {rank: by_prefix.get(PREFIX[rank], "") for rank in RANKS}


def load_quant_report(path):
    """Returns a list of {name, taxID, taxRank, genomeSize, abundance} dicts."""
    required = ["name", "taxID", "taxRank", "genomeSize", "abundance"]
    rows = []
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        idx = {col: i for i, col in enumerate(header)}
        for col in required:
            if col not in idx:
                err(f"{path} is missing a '{col}' column")
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            fields = line.split("\t")
            rows.append({col: fields[idx[col]] for col in required})
    return rows


def build_centrifuger_rows(nodes, names, quant_report_args):
    cache = {}
    samples = []
    taxon_info = {}
    sample_abundance = {}
    for entry in quant_report_args:
        if "=" not in entry:
            err(f"--quant-report entry '{entry}' is not of the form SAMPLE=PATH")
        sample, path = entry.split("=", 1)
        samples.append(sample)
        abundances = {}
        for row in load_quant_report(path):
            taxid = row["taxID"]
            abundances[taxid] = row["abundance"]
            if taxid not in taxon_info:
                lineage = build_lineage(taxid, nodes, names, cache)
                taxon_info[taxid] = {
                    "name": row["name"],
                    "taxRank": row["taxRank"],
                    "genomeSize": row["genomeSize"],
                    "lineage": lineage,
                    "ranks": split_lineage(lineage),
                }
        sample_abundance[sample] = abundances

    header = ["taxID", "name", "taxRank", "genomeSize", "lineage"] + RANKS + samples
    rows = [header]
    for taxid in sorted(taxon_info, key=int):
        info = taxon_info[taxid]
        row = [int(taxid), info["name"], info["taxRank"], int(info["genomeSize"]), info["lineage"]]
        row += [rank_name or None for rank_name in (info["ranks"][rank] for rank in RANKS)]
        row += [float(sample_abundance[sample].get(taxid, 0.0)) for sample in samples]
        rows.append(row)

    if samples:
        first_sample_col = len(header) - len(samples)
        rows[1:] = sorted(rows[1:], key=lambda row: row[first_sample_col], reverse=True)
    return rows, len(samples)


def build_metaphlan_rows(path):
    """Returns merge_metaphlan_tables.py's merged_bugs_list.tsv as a list of rows
    (header first), with each sample column's trailing "_bugs_list" suffix stripped."""
    rows = []
    with open(path) as fh:
        first = fh.readline()
        if not first.startswith("#"):
            err(f"{path} doesn't start with the expected '#<db_version>' comment line")
        header = fh.readline().rstrip("\n").split("\t")
        if not header or header[0] != "clade_name":
            err(f"{path}'s second line isn't the expected 'clade_name\\t<sample>...' header")
        header = [
            col[: -len(METAPHLAN_SAMPLE_SUFFIX)] if col.endswith(METAPHLAN_SAMPLE_SUFFIX) else col
            for col in header
        ]
        rows.append(header)
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            fields = line.split("\t")
            rows.append([fields[0]] + [float(v) for v in fields[1:]])

    if len(header) > 1:
        rows[1:] = sorted(rows[1:], key=lambda row: row[1], reverse=True)
    return rows, len(header) - 1


def color_scale_abundance_columns(ws, n_cols, n_samples, n_rows):
    """Applies a red(low)-yellow(mid)-green(high) color scale across the
    trailing n_samples columns (the per-sample abundance columns, in both
    sheets) -- one shared scale across all of them, same as manually
    selecting that block in Excel and applying Home > Conditional
    Formatting > Color Scales > Red-Yellow-Green. The midpoint is set by
    percentile (rank), not raw value, since these distributions are
    heavily right-skewed (a handful of abundant taxa, a long tail of rare
    ones) -- a value-midpoint scale would leave nearly everything red."""
    if n_samples == 0 or n_rows < 2:
        return
    first_col = get_column_letter(n_cols - n_samples + 1)
    last_col = get_column_letter(n_cols)
    rule = ColorScaleRule(
        start_type="min", start_color=BRIGHT_RED,
        mid_type="percentile", mid_value=50, mid_color=BRIGHT_YELLOW,
        end_type="max", end_color=BRIGHT_GREEN,
    )
    ws.conditional_formatting.add(f"{first_col}2:{last_col}{n_rows}", rule)


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--nodes", required=True, help="centrifuger-inspect --taxonomy-tree dump (nodes.dmp format)")
    p.add_argument("--names", required=True, help="centrifuger-inspect --name-table dump (names.dmp format)")
    p.add_argument("--quant-report", action="append", required=True, metavar="SAMPLE=PATH",
                    help="repeatable; one per sample")
    p.add_argument("--metaphlan-bugs-list", required=True, help="humann3's merged_bugs_list.tsv")
    p.add_argument("--output", required=True)
    args = p.parse_args()

    nodes = load_nodes(args.nodes)
    names = load_names(args.names)
    centrifuger_rows, n_centrifuger_samples = build_centrifuger_rows(nodes, names, args.quant_report)
    metaphlan_rows, n_metaphlan_samples = build_metaphlan_rows(args.metaphlan_bugs_list)

    wb = Workbook()
    ws = wb.active
    ws.title = "centrifuger"
    for row in centrifuger_rows:
        ws.append(row)
    color_scale_abundance_columns(ws, n_cols=len(centrifuger_rows[0]), n_samples=n_centrifuger_samples, n_rows=len(centrifuger_rows))

    ws = wb.create_sheet("metaphlan")
    for row in metaphlan_rows:
        ws.append(row)
    color_scale_abundance_columns(ws, n_cols=len(metaphlan_rows[0]), n_samples=n_metaphlan_samples, n_rows=len(metaphlan_rows))

    wb.save(args.output)


if __name__ == "__main__":
    main()
