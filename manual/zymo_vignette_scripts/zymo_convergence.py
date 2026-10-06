#!/usr/bin/env python3
import re, json
import pandas as pd
import numpy as np

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
YEASTS = {"Saccharomyces cerevisiae": "Saccharomyces_cerevisiae", "Cryptococcus neoformans": "Cryptococcus_neoformans"}
ALL10 = {**BACTERIA, **YEASTS}
VENDOR_GENOME_COPY = {
    "Pseudomonas aeruginosa": 6.1, "Escherichia coli": 8.5, "Salmonella enterica": 8.7,
    "Lactobacillus fermentum": 21.6, "Enterococcus faecalis": 14.6, "Staphylococcus aureus": 15.2,
    "Listeria monocytogenes": 13.9, "Bacillus subtilis": 10.3,
    "Saccharomyces cerevisiae": 0.57, "Cryptococcus neoformans": 0.37,
}
vendor_sum = sum(VENDOR_GENOME_COPY.values())
vendor10 = {s: VENDOR_GENOME_COPY[s]/vendor_sum*100 for s in ALL10}

def collapse_sum(df, base, name_col="name", val_col="abundance"):
    pat = re.compile(rf"^{re.escape(base)}(_[A-Z]{{1,2}})?$")
    mask = df[name_col].apply(lambda x: bool(pat.match(x)))
    return df.loc[mask, val_col].sum()

# per-replicate RAW (un-renormalized) read-based abundance across all 10 taxa
read_raw = {rep: {} for rep in REPLICATES}
for rep in REPLICATES:
    df = pd.read_csv(f"{BASE}/centrifuger_dna/{rep}_centrifuger_quantification_report.tsv", sep="\t")
    df = df[df["taxRank"] == "species"]
    for disp, base in ALL10.items():
        read_raw[rep][disp] = collapse_sum(df, base) * 100.0

# per-replicate RAW assembly-based coverage (yeasts always 0, no MAGs recovered)
widb = pd.read_csv(f"{BASE}/metawrap_bin_refine/dRep/data_tables/Widb.csv")
genome_size = dict(zip(widb["genome"].str.replace(".fa",""), widb["size"]))
gdf = pd.read_csv(f"{BASE}/metawrap_kmer/GTDBTK_classify_wf/gtdbtk.bac120.summary.tsv", sep="\t")
gdf["species_raw"] = gdf["classification"].apply(lambda c: re.search(r"s__(.+)$", c).group(1).strip())
def to_vendor_label(species_raw):
    norm = species_raw.replace(" ", "_")
    for disp, base in BACTERIA.items():
        if re.match(rf"^{re.escape(base)}(_[A-Z]{{1,2}})?$", norm):
            return disp
    return None
gdf["vendor_label"] = gdf["species_raw"].apply(to_vendor_label)
mag_to_label = dict(zip(gdf["user_genome"], gdf["vendor_label"]))

asm_raw = {rep: {s: 0.0 for s in ALL10} for rep in REPLICATES}
for rep in REPLICATES:
    df = pd.read_csv(f"{BASE}/mags/{rep}/DNA/{rep}.refstats", sep="\t")
    df.columns = [c.lstrip("#") for c in df.columns]
    for _, row in df.iterrows():
        label = mag_to_label.get(row["name"])
        if label is None: continue
        gsize = genome_size[row["name"]]
        asm_raw[rep][label] += row["unambiguousMB"] * 1e6 / gsize

def bray_curtis(a, b, species):
    num = sum(abs(a[s]-b[s]) for s in species)
    den = sum(a[s]+b[s] for s in species)
    return num/den

def cumulative_bc(raw_dict, species):
    out = []
    acc = {s: [] for s in species}
    for rep in REPLICATES:
        for s in species:
            acc[s].append(raw_dict[rep][s])
        mean_comp = {s: np.mean(acc[s]) for s in species}
        tot = sum(mean_comp.values())
        norm_comp = {s: mean_comp[s]/tot*100 for s in species}
        out.append(bray_curtis(norm_comp, vendor10, species))
    return out

read_bc = cumulative_bc(read_raw, list(ALL10.keys()))
asm_bc = cumulative_bc(asm_raw, list(ALL10.keys()))
print("Read-based cumulative Bray-Curtis (rep1, 1-2, ..., 1-7):", [round(x,4) for x in read_bc])
print("Assembly-based cumulative Bray-Curtis (rep1, 1-2, ..., 1-7):", [round(x,4) for x in asm_bc])

with open("./zymo_convergence_results.json","w") as f:
    json.dump({"read_based": read_bc, "assembly_based": asm_bc}, f, indent=2)
