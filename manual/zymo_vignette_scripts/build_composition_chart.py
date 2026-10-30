#!/usr/bin/env python3
import json
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

with open("./vendor_bars_pixel_extracted.json") as f:
    vendor_bars = json.load(f)
with open("./zymo_fulldepth_results.json") as f:
    rb = json.load(f)
with open("./zymo_asm_coverage_results.json") as f:
    asmr = json.load(f)

species_colors = {
    "Pseudomonas aeruginosa": "#00458 6".replace(" ",""),
}
# exact species order (bottom -> top) and colors, matching the original chart precisely
STACK_ORDER = [
    ("Pseudomonas aeruginosa", "#004586"),
    ("Escherichia coli",       "#FF420E"),
    ("Salmonella enterica",    "#FFD320"),
    ("Lactobacillus fermentum","#579D1C"),
    ("Enterococcus faecalis",  "#7E0021"),
    ("Staphylococcus aureus",  "#83CAFF"),
    ("Listeria monocytogenes", "#314004"),
    ("Bacillus subtilis",      "#AECF00"),
]
LEGEND_LABELS = {
    "Bacillus subtilis": "Bacillus subtilis (G+)",
    "Listeria monocytogenes": "Listeria monocytogenes (G+)",
    "Staphylococcus aureus": "Staphylococcus aureus (G+)",
    "Enterococcus faecalis": "Enterococcus faecalis (G+)",
    "Lactobacillus fermentum": "Lactobacillus fermentum (G+)",
    "Salmonella enterica": "Salmonella enterica (G-)",
    "Escherichia coli": "Escherichia coli (G-)",
    "Pseudomonas aeruginosa": "Pseudomonas aeruginosa (G-)",
}

bars = []
for name in ["Theoretical", "ZymoBIOMICS DNA Miniprep", "HMP Protocol", "Supplier M", "Supplier Q"]:
    bars.append((name.replace("ZymoBIOMICS DNA Miniprep","ZymoBIOMICS\nDNA Miniprep").replace("HMP Protocol","HMP\nProtocol"), vendor_bars[name]))
bars.append(("metamorph\n(read-based)", {s: rb["read_based"][s]["mean"] for s, _ in STACK_ORDER}))
bars.append(("metamorph\n(assembly-based)", {s: asmr[s]["mean"] for s, _ in STACK_ORDER}))

fig, ax = plt.subplots(figsize=(13.5, 7.2), dpi=100)
x = np.arange(len(bars))
bottoms = np.zeros(len(bars))
for species, color in STACK_ORDER:
    vals = np.array([b[1][species] for b in bars])
    ax.bar(x, vals, bottom=bottoms, width=0.65, color=color, label=LEGEND_LABELS[species],
           edgecolor="white", linewidth=0.4)
    bottoms += vals

ax.set_xticks(x)
ax.set_xticklabels([b[0] for b in bars], fontsize=10.5)
ax.set_ylim(0, 100)
ax.set_ylabel("Microbial Composition (%)", fontsize=12)
ax.yaxis.set_major_formatter(lambda v, pos: f"{int(v)}%")
ax.spines[["top","right"]].set_visible(False)

# dashed separator between vendor bars (0-4) and metamorph bars (5,6)
sep_x = 4.5
ax.axvline(sep_x, color="gray", linestyle="--", linewidth=1, ymin=0, ymax=1)

# legend, top-to-bottom matches stack top-to-bottom (reverse of STACK_ORDER insertion order
# since matplotlib legend lists in the order bars were added -> reverse for top-down stack order)
handles, labels = ax.get_legend_handles_labels()
ax.legend(handles[::-1], labels[::-1], title="Species", loc="upper left",
          bbox_to_anchor=(1.01, 1.0), frameon=False, fontsize=10, title_fontsize=11)

plt.tight_layout()
plt.savefig("./zymo_composition_comparison_NEW.png", bbox_inches="tight")
print("saved")
