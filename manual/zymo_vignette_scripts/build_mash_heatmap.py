#!/usr/bin/env python3
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

SKETCH_SIZE = 500000  # largest grid point: most stable/well-resolved distance estimates

df = pd.read_csv("/data/OpenOmics/dev/metamorph/remove_cruft_at_end/enzymatic_community_data/ZymoBIOMICS.STD.refseq.v3/Genomes/mash_grid_results.tsv", sep="\t")
df = df[df["sketch_size"] == SKETCH_SIZE]

NAME_MAP = {
    "./Pseudomonas_aeruginosa_complete_genome.fna":  "P. aeruginosa",
    "./Escherichia_coli_complete_genome.fna":        "E. coli",
    "./Salmonella_enterica_complete_genome.fna":     "S. enterica",
    "./Lactobacillus_fermentum_complete_genome.fna": "L. fermentum",
    "./Enterococcus_faecalis_complete_genome.fna":   "E. faecalis",
    "./Staphylococcus_aureus_complete_genome.fna":   "S. aureus",
    "./Listeria_monocytogenes_complete_genome.fna":  "L. monocytogenes",
    "./Bacillus_subtilis_complete_genome.fna":       "B. subtilis",
    "./Saccharomyces_cerevisiae_complete_genome.fasta": "S. cerevisiae",
    "./Cryptococcus_neoformans_complete_genome.fasta":  "C. neoformans",
}
ORDER = [
    "P. aeruginosa", "E. coli", "S. enterica", "L. fermentum", "E. faecalis",
    "S. aureus", "L. monocytogenes", "B. subtilis", "S. cerevisiae", "C. neoformans",
]

df["g1"] = df["genome1"].map(NAME_MAP)
df["g2"] = df["genome2"].map(NAME_MAP)

mat = df.pivot(index="g1", columns="g2", values="mash_distance").loc[ORDER, ORDER]

fig, ax = plt.subplots(figsize=(8.5, 7.3), dpi=100)
im = ax.imshow(mat.values, cmap="Blues", vmin=0, vmax=1)

ax.set_xticks(range(len(ORDER)))
ax.set_yticks(range(len(ORDER)))
ax.set_xticklabels(ORDER, rotation=40, ha="right", fontsize=9.5)
ax.set_yticklabels(ORDER, fontsize=9.5)

for i in range(len(ORDER)):
    for j in range(len(ORDER)):
        v = mat.values[i, j]
        color = "white" if v > 0.6 else "#0b0b0b"
        ax.text(j, i, f"{v:.2f}", ha="center", va="center", fontsize=7.8, color=color)

cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.03)
cbar.set_label("Mash distance (dissimilarity)", fontsize=10.5)

ax.set_title(f"Pairwise Mash distance, 10 ZymoBIOMICS reference genomes\n(sketch size s={SKETCH_SIZE:,})", fontsize=12.5)

plt.tight_layout()
plt.savefig("./zymo_mash_heatmap.png", bbox_inches="tight", facecolor="white")
print("saved")
