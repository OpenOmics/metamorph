#!/usr/bin/env python3
import json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

with open("./zymo_fulldepth_results.json") as f:
    rb = json.load(f)
with open("./zymo_asm_coverage_results.json") as f:
    asmr = json.load(f)
with open("./zymo_refmap_results.json") as f:
    refmap = json.load(f)

SPECIES_ORDER = [
    "Lactobacillus fermentum", "Salmonella enterica", "Escherichia coli",
    "Pseudomonas aeruginosa", "Listeria monocytogenes", "Bacillus subtilis",
    "Enterococcus faecalis", "Staphylococcus aureus",
]
COLORS = {
    "target":  "#898781",
    "read":    "#2a78d6",
    "asm":     "#1baf7a",
    "refmap":  "#eda100",
}

fig, ax = plt.subplots(figsize=(12.5, 6.5), dpi=100)
x = np.arange(len(SPECIES_ORDER))
width = 0.2

target_vals = [rb["vendor_renorm"][s] for s in SPECIES_ORDER]
read_vals = [rb["read_based"][s]["mean"] for s in SPECIES_ORDER]
asm_vals = [asmr[s]["mean"] for s in SPECIES_ORDER]
refmap_vals = [refmap[s]["mean"] for s in SPECIES_ORDER]
read_err = [rb["read_based"][s]["sd"] for s in SPECIES_ORDER]
asm_err = [asmr[s]["sd"] for s in SPECIES_ORDER]
refmap_err = [refmap[s]["sd"] for s in SPECIES_ORDER]

ax.bar(x - 1.5*width, target_vals, width, color=COLORS["target"], label="Vendor target (Genome Copy %)")
ax.bar(x - 0.5*width, read_vals, width, yerr=read_err, capsize=2.5, color=COLORS["read"], label="Read-based (Centrifuger)")
ax.bar(x + 0.5*width, asm_vals, width, yerr=asm_err, capsize=2.5, color=COLORS["asm"], label="Assembly-based (recovered MAGs)")
ax.bar(x + 1.5*width, refmap_vals, width, yerr=refmap_err, capsize=2.5, color=COLORS["refmap"], label="Direct reference-genome mapping (vendor genomes)")

ax.set_xticks(x)
ax.set_xticklabels(SPECIES_ORDER, rotation=30, ha="right", fontsize=10)
ax.set_ylabel("Composition (%, renormalized across the 8 bacteria)", fontsize=11.5)
ax.spines[["top", "right"]].set_visible(False)
ax.legend(loc="upper right", fontsize=10, frameon=False)
ax.set_title("Three independent composition estimates vs. vendor target, per species", fontsize=13.5)

plt.tight_layout()
plt.savefig("./zymo_refmap_comparison.png", bbox_inches="tight", facecolor="white")
print("saved")
