#!/usr/bin/env python3
# Reads ./refmap_refstats/<replicate>.refstats, produced outside metamorph by:
#   module load bbtools
#   bbsplit.sh path=<idx> build=1 ref=<the 10 ZymoBIOMICS.STD.refseq.v3/Genomes/*_complete_genome.fna/fasta files>
#   # then per replicate, reusing the same path=<idx>:
#   bbsplit.sh path=<idx> minid=0.90 unpigz=t sortscafs=f nzo=f \
#       refstats=<replicate>.refstats statsfile=... scafstats=... covstats=... rpkm=... \
#       in=trimmed_reads/<replicate>/<replicate>_R1_dehost.fastq.gz \
#       in2=trimmed_reads/<replicate>/<replicate>_R2_dehost.fastq.gz
# Same tool/flags/criterion as metamorph's own bbtools_index_map rule (Section 4.1 of the vignette),
# just pointed at the vendor's reference genomes instead of metamorph's recovered MAGs.
import json
import pandas as pd
import numpy as np
from scipy import stats

REFMAP_OUT = "./refmap_refstats"
REPLICATES = ["ERR16012159","ERR16012160","ERR16012161","ERR16012162","ERR16012163","ERR16012164","ERR16012165"]

# vendor display label -> reference-genome-file basename (bbsplit set name)
BACTERIA = {
    "Pseudomonas aeruginosa":  "Pseudomonas_aeruginosa_complete_genome",
    "Escherichia coli":       "Escherichia_coli_complete_genome",
    "Salmonella enterica":    "Salmonella_enterica_complete_genome",
    "Lactobacillus fermentum":"Lactobacillus_fermentum_complete_genome",
    "Enterococcus faecalis":  "Enterococcus_faecalis_complete_genome",
    "Staphylococcus aureus":  "Staphylococcus_aureus_complete_genome",
    "Listeria monocytogenes": "Listeria_monocytogenes_complete_genome",
    "Bacillus subtilis":      "Bacillus_subtilis_complete_genome",
}
YEASTS = {
    "Saccharomyces cerevisiae": "Saccharomyces_cerevisiae_complete_genome",
    "Cryptococcus neoformans":  "Cryptococcus_neoformans_complete_genome",
}
ALL10 = {**BACTERIA, **YEASTS}

GENOME_SIZE = {
    "Bacillus_subtilis_complete_genome": 4222511,
    "Enterococcus_faecalis_complete_genome": 2895646,
    "Escherichia_coli_complete_genome": 4925141,
    "Lactobacillus_fermentum_complete_genome": 2227379,
    "Listeria_monocytogenes_complete_genome": 3032011,
    "Pseudomonas_aeruginosa_complete_genome": 6846579,
    "Salmonella_enterica_complete_genome": 4824135,
    "Staphylococcus_aureus_complete_genome": 2729384,
    "Saccharomyces_cerevisiae_complete_genome": 18769914,
    "Cryptococcus_neoformans_complete_genome": 26305994,
}

VENDOR_GENOME_COPY = {
    "Pseudomonas aeruginosa": 6.1, "Escherichia coli": 8.5, "Salmonella enterica": 8.7,
    "Lactobacillus fermentum": 21.6, "Enterococcus faecalis": 14.6, "Staphylococcus aureus": 15.2,
    "Listeria monocytogenes": 13.9, "Bacillus subtilis": 10.3,
}
vendor_8_sum = sum(VENDOR_GENOME_COPY.values())
VENDOR_RENORM = {s: VENDOR_GENOME_COPY[s]/vendor_8_sum*100 for s in BACTERIA}

refmap_rows_cov = {s: [] for s in BACTERIA}
yeast_rows = {s: [] for s in YEASTS}

for rep in REPLICATES:
    df = pd.read_csv(f"{REFMAP_OUT}/{rep}.refstats", sep="\t")
    df.columns = [c.lstrip("#") for c in df.columns]
    cov = {}
    for disp, setname in ALL10.items():
        row = df[df["name"] == setname]
        if len(row) == 0:
            cov[disp] = 0.0
            continue
        mb = row["unambiguousMB"].values[0]
        cov[disp] = mb * 1e6 / GENOME_SIZE[setname]
    bac_cov_sum = sum(cov[s] for s in BACTERIA)
    for s in BACTERIA:
        refmap_rows_cov[s].append(cov[s] / bac_cov_sum * 100.0 if bac_cov_sum > 0 else 0.0)
    for s in YEASTS:
        # yeasts reported as raw share of total (bacteria+yeast) coverage, not renormalized
        all_cov_sum = bac_cov_sum + sum(cov[y] for y in YEASTS)
        yeast_rows[s].append(cov[s] / all_cov_sum * 100.0 if all_cov_sum > 0 else 0.0)

refmap_summary = {}
for s in BACTERIA:
    vals = np.array(refmap_rows_cov[s])
    t, p = stats.ttest_1samp(vals, VENDOR_RENORM[s])
    refmap_summary[s] = dict(mean=float(vals.mean()), sd=float(vals.std(ddof=1)), t=float(t), p=float(p), target=VENDOR_RENORM[s], raw=list(vals))

print(f"{'Species':28s} {'target%':>8s} {'refmap mean%':>13s} {'sd':>7s} {'t':>8s} {'p':>12s}")
for s in BACTERIA:
    d = refmap_summary[s]
    print(f"{s:28s} {d['target']:8.2f} {d['mean']:13.2f} {d['sd']:7.2f} {d['t']:8.2f} {d['p']:12.6f}")

print("\nYeast raw share of total coverage (not renormalized):")
for s in YEASTS:
    vals = np.array(yeast_rows[s])
    print(f"{s:28s} mean={vals.mean():.5f}% raw={[round(v,5) for v in vals]}")

print("\nPer-replicate coverage-based %, direct reference mapping:")
for s in BACTERIA:
    print(f"{s:28s}", [round(v,2) for v in refmap_rows_cov[s]])

out = {s: {"mean": refmap_summary[s]["mean"], "sd": refmap_summary[s]["sd"],
           "t": refmap_summary[s]["t"], "p": refmap_summary[s]["p"],
           "target": refmap_summary[s]["target"], "raw": refmap_summary[s]["raw"]} for s in BACTERIA}
with open("./zymo_refmap_results.json", "w") as f:
    json.dump(out, f, indent=2)
print("\nWrote ./zymo_refmap_results.json")
