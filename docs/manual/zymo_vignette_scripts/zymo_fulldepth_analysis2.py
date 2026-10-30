#!/usr/bin/env python3
import re, json
import pandas as pd
import numpy as np
from scipy import stats

BASE = "/vf/users/OpenOmics/dev/metamorph/remove_cruft_at_end/enzymatic_community_data/output_full_depth/metagenome_results"
REPLICATES = ["ERR16012159","ERR16012160","ERR16012161","ERR16012162","ERR16012163","ERR16012164","ERR16012165"]

BACTERIA = {
    "Pseudomonas aeruginosa":  "Pseudomonas_aeruginosa",
    "Escherichia coli":       "Escherichia_coli",
    "Salmonella enterica":    "Salmonella_enterica",
    "Lactobacillus fermentum":"Limosilactobacillus_fermentum",
    "Enterococcus faecalis":  "Enterococcus_faecalis",
    "Staphylococcus aureus":  "Staphylococcus_aureus",
    "Listeria monocytogenes": "Listeria_monocytogenes",
    "Bacillus subtilis":      "Bacillus_spizizenii",
}
VENDOR_GENOME_COPY = {
    "Pseudomonas aeruginosa": 6.1, "Escherichia coli": 8.5, "Salmonella enterica": 8.7,
    "Lactobacillus fermentum": 21.6, "Enterococcus faecalis": 14.6, "Staphylococcus aureus": 15.2,
    "Listeria monocytogenes": 13.9, "Bacillus subtilis": 10.3,
}
vendor_8_sum = sum(VENDOR_GENOME_COPY.values())
VENDOR_RENORM = {s: VENDOR_GENOME_COPY[s]/vendor_8_sum*100 for s in BACTERIA}

# genome sizes per winning MAG, from dRep's Widb.csv
widb = pd.read_csv(f"{BASE}/metawrap_bin_refine/dRep/data_tables/Widb.csv")
genome_size = dict(zip(widb["genome"].str.replace(".fa",""), widb["size"]))

summary_path = f"{BASE}/metawrap_kmer/GTDBTK_classify_wf/gtdbtk.bac120.summary.tsv"
gdf = pd.read_csv(summary_path, sep="\t")
def species_from_classification(c):
    m = re.search(r"s__(.+)$", c)
    return m.group(1).strip() if m else None
gdf["species_raw"] = gdf["classification"].apply(species_from_classification)
def to_vendor_label(species_raw):
    norm = species_raw.replace(" ", "_")
    for disp, base in BACTERIA.items():
        pat = re.compile(rf"^{re.escape(base)}(_[A-Z]{{1,2}})?$")
        if pat.match(norm):
            return disp
    return None
gdf["vendor_label"] = gdf["species_raw"].apply(to_vendor_label)
mag_to_label = dict(zip(gdf["user_genome"], gdf["vendor_label"]))

asm_rows_cov = {s: [] for s in BACTERIA}
asm_rows_reads = {s: [] for s in BACTERIA}
for rep in REPLICATES:
    path = f"{BASE}/mags/{rep}/DNA/{rep}.refstats"
    df = pd.read_csv(path, sep="\t")
    df.columns = [c.lstrip("#") for c in df.columns]
    cov_raw = {s: 0.0 for s in BACTERIA}
    reads_raw = {s: 0.0 for s in BACTERIA}
    for _, row in df.iterrows():
        label = mag_to_label.get(row["name"])
        if label is None:
            continue
        gsize = genome_size[row["name"]]
        cov_raw[label] += row["unambiguousMB"] * 1e6 / gsize   # coverage = bases / genome size
        reads_raw[label] += row["%unambiguousReads"]
    cov_sum = sum(cov_raw.values())
    reads_sum = sum(reads_raw.values())
    for s in BACTERIA:
        asm_rows_cov[s].append(cov_raw[s] / cov_sum * 100.0 if cov_sum > 0 else 0.0)
        asm_rows_reads[s].append(reads_raw[s] / reads_sum * 100.0 if reads_sum > 0 else 0.0)

print(f"{'Species':28s} {'target%':>8s} {'reads-based mean%':>18s} {'coverage-based mean%':>20s}")
for s in BACTERIA:
    rvals = np.array(asm_rows_reads[s]); cvals = np.array(asm_rows_cov[s])
    print(f"{s:28s} {VENDOR_RENORM[s]:8.2f} {rvals.mean():18.2f} {cvals.mean():20.2f}")

print("\nPer-replicate coverage-based %:")
for s in BACTERIA:
    print(f"{s:28s}", [round(v,2) for v in asm_rows_cov[s]])

asm_summary = {}
for s in BACTERIA:
    vals = np.array(asm_rows_cov[s])
    t, p = stats.ttest_1samp(vals, VENDOR_RENORM[s])
    asm_summary[s] = dict(mean=float(vals.mean()), sd=float(vals.std(ddof=1)), t=float(t), p=float(p), target=VENDOR_RENORM[s])

print("\n=== ASSEMBLY-BASED (coverage-normalized, final) t-tests ===")
for s in BACTERIA:
    d = asm_summary[s]
    print(f"{s:28s} target={d['target']:6.2f} mean={d['mean']:6.2f} sd={d['sd']:5.2f} t={d['t']:8.2f} p={d['p']:10.6f}")

with open("./zymo_asm_coverage_results.json","w") as f:
    json.dump(asm_summary, f, indent=2)
