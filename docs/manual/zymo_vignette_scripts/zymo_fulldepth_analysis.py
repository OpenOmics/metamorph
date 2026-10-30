#!/usr/bin/env python3
import re, glob, json
import pandas as pd
import numpy as np
from scipy import stats

BASE = "/vf/users/OpenOmics/dev/metamorph/remove_cruft_at_end/enzymatic_community_data/output_full_depth/metagenome_results"
REPLICATES = ["ERR16012159","ERR16012160","ERR16012161","ERR16012162","ERR16012163","ERR16012164","ERR16012165"]

# vendor-label -> current-valid GTDB base name used as the regex root for the
# generic suffix-collapse rule (<base>(_[A-Z]{1,2})? -> <base>)
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
YEASTS = {
    "Saccharomyces cerevisiae": "Saccharomyces_cerevisiae",
    "Cryptococcus neoformans":  "Cryptococcus_neoformans",  # confirmed absent in centrifuger output; kept for completeness
}
TARGET_ORDER = list(BACTERIA.keys())  # vendor table order; chart stack order set separately

VENDOR_GENOME_COPY = {
    "Pseudomonas aeruginosa": 6.1, "Escherichia coli": 8.5, "Salmonella enterica": 8.7,
    "Lactobacillus fermentum": 21.6, "Enterococcus faecalis": 14.6, "Staphylococcus aureus": 15.2,
    "Listeria monocytogenes": 13.9, "Bacillus subtilis": 10.3,
    "Saccharomyces cerevisiae": 0.57, "Cryptococcus neoformans": 0.37,
}
vendor_8_sum = sum(VENDOR_GENOME_COPY[s] for s in BACTERIA)
VENDOR_RENORM = {s: VENDOR_GENOME_COPY[s]/vendor_8_sum*100 for s in BACTERIA}

def collapse_sum(df, base, name_col="name", val_col="abundance"):
    pat = re.compile(rf"^{re.escape(base)}(_[A-Z]{{1,2}})?$")
    mask = df[name_col].apply(lambda x: bool(pat.match(x)))
    return df.loc[mask, val_col].sum()

# ---------------- READ-BASED (Centrifuger, full depth) ----------------
read_rows = {s: [] for s in list(BACTERIA.keys()) + list(YEASTS.keys())}
for rep in REPLICATES:
    path = f"{BASE}/centrifuger_dna/{rep}_centrifuger_quantification_report.tsv"
    df = pd.read_csv(path, sep="\t")
    df = df[df["taxRank"] == "species"]
    raw = {}
    for disp, base in {**BACTERIA, **YEASTS}.items():
        raw[disp] = collapse_sum(df, base) * 100.0
    bac_sum = sum(raw[s] for s in BACTERIA)
    for s in BACTERIA:
        read_rows[s].append(raw[s] / bac_sum * 100.0)
    for s in YEASTS:
        read_rows[s].append(raw[s])  # yeasts reported raw (not renormalized), matching existing vignette convention

read_summary = {}
for s in BACTERIA:
    vals = np.array(read_rows[s])
    t, p = stats.ttest_1samp(vals, VENDOR_RENORM[s])
    read_summary[s] = dict(mean=vals.mean(), sd=vals.std(ddof=1), t=t, p=p, target=VENDOR_RENORM[s], raw=list(vals))
for s in YEASTS:
    vals = np.array(read_rows[s])
    read_summary[s] = dict(mean=vals.mean(), sd=vals.std(ddof=1), raw=list(vals))

# ---------------- ASSEMBLY-BASED (GTDB-Tk + MAG mapping, full depth) ----------------
summary_path = f"{BASE}/metawrap_kmer/GTDBTK_classify_wf/gtdbtk.bac120.summary.tsv"
gdf = pd.read_csv(summary_path, sep="\t")
def species_from_classification(c):
    m = re.search(r"s__(.+)$", c)
    return m.group(1).strip() if m else None
gdf["species_raw"] = gdf["classification"].apply(species_from_classification)

# map each MAG's raw GTDB species string to a vendor display label via the
# SAME generic suffix-collapse rule, applied on the underscore/space-normalized name
def to_vendor_label(species_raw):
    norm = species_raw.replace(" ", "_")
    for disp, base in BACTERIA.items():
        pat = re.compile(rf"^{re.escape(base)}(_[A-Z]{{1,2}})?$")
        if pat.match(norm):
            return disp
    return None  # off-target / unmapped MAG

gdf["vendor_label"] = gdf["species_raw"].apply(to_vendor_label)
mag_to_label = dict(zip(gdf["user_genome"], gdf["vendor_label"]))
print("MAG -> vendor label map:", json.dumps(mag_to_label, indent=2))

asm_rows = {s: [] for s in BACTERIA}
for rep in REPLICATES:
    path = f"{BASE}/mags/{rep}/DNA/{rep}.refstats"
    df = pd.read_csv(path, sep="\t")
    df.columns = [c.lstrip("#") for c in df.columns]
    raw = {s: 0.0 for s in BACTERIA}
    for _, row in df.iterrows():
        label = mag_to_label.get(row["name"])
        if label is not None:
            raw[label] += row["%unambiguousReads"]
    bac_sum = sum(raw.values())
    for s in BACTERIA:
        asm_rows[s].append(raw[s] / bac_sum * 100.0 if bac_sum > 0 else 0.0)

asm_summary = {}
for s in BACTERIA:
    vals = np.array(asm_rows[s])
    t, p = stats.ttest_1samp(vals, VENDOR_RENORM[s])
    asm_summary[s] = dict(mean=vals.mean(), sd=vals.std(ddof=1), t=t, p=p, target=VENDOR_RENORM[s], raw=list(vals))

# ---------------- Print final tables ----------------
print("\n=== READ-BASED (Centrifuger, full depth, GTDB-collapse applied) ===")
print(f"{'Species':28s} {'target%':>8s} {'mean%':>8s} {'sd':>7s} {'t':>8s} {'p':>10s}")
for s in BACTERIA:
    d = read_summary[s]
    print(f"{s:28s} {d['target']:8.2f} {d['mean']:8.2f} {d['sd']:7.2f} {d['t']:8.2f} {d['p']:10.5f}")
for s in YEASTS:
    d = read_summary[s]
    print(f"{s:28s} {'(excl)':>8s} {d['mean']:8.3f} {d['sd']:7.3f}")

print("\n=== ASSEMBLY-BASED (GTDB-Tk + MAG mapping, full depth) ===")
print(f"{'Species':28s} {'target%':>8s} {'mean%':>8s} {'sd':>7s} {'t':>8s} {'p':>10s}")
for s in BACTERIA:
    d = asm_summary[s]
    print(f"{s:28s} {d['target']:8.2f} {d['mean']:8.2f} {d['sd']:7.2f} {d['t']:8.2f} {d['p']:10.5f}")

print("\n=== Per-replicate raw read-based % (post-collapse, pre-renorm... actually post-renorm, bacteria only) ===")
for s in BACTERIA:
    print(s, [round(v,2) for v in read_rows[s]])

print("\n=== Per-replicate raw assembly-based % (renorm, bacteria only) ===")
for s in BACTERIA:
    print(s, [round(v,2) for v in asm_rows[s]])

# Dump machine-readable summary for the plotting script
out = {
    "vendor_renorm": VENDOR_RENORM,
    "read_based": {s: {"mean": read_summary[s]["mean"], "sd": read_summary[s]["sd"]} for s in BACTERIA},
    "assembly_based": {s: {"mean": asm_summary[s]["mean"], "sd": asm_summary[s]["sd"]} for s in BACTERIA},
}
with open("./zymo_fulldepth_results.json", "w") as f:
    json.dump(out, f, indent=2)
print("\nWrote ./zymo_fulldepth_results.json")
