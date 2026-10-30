#!/usr/bin/env python3
# D6305 (200 ng isolated-DNA input) / D6306 (2,000 ng isolated-DNA input) ZymoBIOMICS
# Microbial Community DNA Standard, mapped to internal samples in3634_1 (D6305) /
# in3634_2 (D6306) per the submitter's own records. Same target composition as D6300
# (identical 12%/2% gDNA mix, confirmed with the submitter directly, since the public
# datasheet (ds1706) publishes no numeric table for the DNA-only product) -- same
# VENDOR_GENOME_COPY values already used for D6300 in zymo_fulldepth_analysis.py, since
# genome sizes/strains are identical.
#
# n=1 per condition (D6305 and D6306 are two different input masses of the same
# product, not replicates of each other) -- no SD/t-test is computed; both samples'
# values are reported individually alongside the target.
import re, json
import pandas as pd

BASE = "/data/OpenOmics/dev/metamorph/remove_cruft_at_end/zymo/output_zymo_data/metagenome_results"
SAMPLES = {"in3634_1": "D6305 (200 ng)", "in3634_2": "D6306 (2,000 ng)"}

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
    "Cryptococcus neoformans":  "Cryptococcus_neoformans",
}

VENDOR_GENOME_COPY = {
    "Pseudomonas aeruginosa": 6.1, "Escherichia coli": 8.5, "Salmonella enterica": 8.7,
    "Lactobacillus fermentum": 21.6, "Enterococcus faecalis": 14.6, "Staphylococcus aureus": 15.2,
    "Listeria monocytogenes": 13.9, "Bacillus subtilis": 10.3,
    "Saccharomyces cerevisiae": 0.57, "Cryptococcus neoformans": 0.37,
}
vendor_8_sum = sum(VENDOR_GENOME_COPY[s] for s in BACTERIA)
VENDOR_RENORM = {s: VENDOR_GENOME_COPY[s] / vendor_8_sum * 100 for s in BACTERIA}

def collapse_sum(df, base, name_col="name", val_col="abundance"):
    pat = re.compile(rf"^{re.escape(base)}(_[A-Z]{{1,2}})?$")
    mask = df[name_col].apply(lambda x: bool(pat.match(x)))
    return df.loc[mask, val_col].sum()

read_results = {}
for sample in SAMPLES:
    path = f"{BASE}/centrifuger_dna/{sample}_centrifuger_quantification_report.tsv"
    df = pd.read_csv(path, sep="\t")
    df = df[df["taxRank"] == "species"]
    raw = {}
    for disp, base in {**BACTERIA, **YEASTS}.items():
        raw[disp] = collapse_sum(df, base) * 100.0
    bac_sum = sum(raw[s] for s in BACTERIA)
    sample_out = {}
    for s in BACTERIA:
        pct = raw[s] / bac_sum * 100.0
        sample_out[s] = dict(pct=pct, target=VENDOR_RENORM[s], fold=pct / VENDOR_RENORM[s])
    for s in YEASTS:
        sample_out[s] = dict(pct=raw[s])  # yeast reported raw, not renormalized (matches D6300 page convention)
    read_results[sample] = sample_out

print(f"{'Species':28s} {'target%':>8s} {'D6305 %':>9s} {'fold':>7s} {'D6306 %':>9s} {'fold':>7s}")
for s in BACTERIA:
    d5 = read_results["in3634_1"][s]
    d6 = read_results["in3634_2"][s]
    print(f"{s:28s} {d5['target']:8.2f} {d5['pct']:9.2f} {d5['fold']:6.2f}x {d6['pct']:9.2f} {d6['fold']:6.2f}x")
print()
for s in YEASTS:
    print(f"{s:28s} D6305={read_results['in3634_1'][s]['pct']:.4f}%  D6306={read_results['in3634_2'][s]['pct']:.4f}%")

with open("./zymo_d6305_d6306_results.json", "w") as f:
    json.dump({"vendor_renorm": VENDOR_RENORM, "read_based": read_results}, f, indent=2)
print("\nWrote ./zymo_d6305_d6306_results.json")
