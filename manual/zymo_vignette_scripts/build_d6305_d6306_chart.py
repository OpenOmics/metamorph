#!/usr/bin/env python3
import json
import matplotlib
matplotlib.use("Agg")
matplotlib.rcParams["hatch.linewidth"] = 3
import matplotlib.pyplot as plt
import numpy as np

with open("./zymo_d6305_d6306_results.json") as f:
    rb = json.load(f)
with open("./zymo_d6305_d6306_assembly_results.json") as f:
    asm = json.load(f)

# same species -> color mapping as build_composition_chart.py, for visual consistency
# with the D6300 vignette's charts (same organism, same color everywhere on the site)
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
COLOR = dict(STACK_ORDER)
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
COMBINED_LABEL = asm["combined_label"]

target = {s: rb["vendor_renorm"][s] for s, _ in STACK_ORDER}

bars = [("Theoretical\n(D6300/D6305/D6306\nsame mix)", "target", None)]
for sample, tag in [("in3634_1", "D6305\n(200 ng)"), ("in3634_2", "D6306\n(2,000 ng)")]:
    bars.append((f"{tag}\nread-based", "read", sample))
    bars.append((f"{tag}\nassembly-based\n(metaWRAP/dRep/GTDB-Tk)", "asm", sample))
is_asm = np.array([kind == "asm" for _, kind, _ in bars])

fig, ax = plt.subplots(figsize=(11.5, 7.2), dpi=100)
x = np.arange(len(bars))
bottoms = np.zeros(len(bars))

for species, color in STACK_ORDER:
    vals = np.zeros(len(bars))
    for i, (_, kind, sample) in enumerate(bars):
        if kind == "target":
            vals[i] = target[species]
        elif kind == "read":
            vals[i] = rb["read_based"][sample][species]["pct"]
        elif species == "Salmonella enterica":
            vals[i] = 0.0  # its mass is drawn as part of E. coli's combined/striped slot below
        elif species == "Escherichia coli":
            vals[i] = asm["assembly_based"][sample][COMBINED_LABEL]["pct"]
        else:
            vals[i] = asm["assembly_based"][sample][species]["pct"]

    # the normal (non-"combined") segments: solid fill, one legend entry per species
    normal = np.where(is_asm & (species in ("Escherichia coli", "Salmonella enterica")), 0.0, vals)
    ax.bar(x, normal, bottom=bottoms, width=0.65, color=color, label=LEGEND_LABELS[species],
           edgecolor="white", linewidth=0.4)

    if species == "Escherichia coli":
        # the combined, unseparated E. coli + S. enterica measurement on assembly-based
        # bars only -- diagonal stripes of both species' own colors, same vertical slot
        # E. coli/Salmonella occupy in every other bar (Section 5.1); deliberately
        # unlabeled so it doesn't add a separate legend entry (the plain "Escherichia
        # coli (G-)" swatch above already identifies the color)
        combined_vals = np.where(is_asm, vals, 0.0)
        ax.bar(x, combined_vals, bottom=bottoms, width=0.65, color=COLOR["Escherichia coli"],
               edgecolor=COLOR["Salmonella enterica"], hatch="//", linewidth=0)

    bottoms += vals

ax.set_xticks(x)
ax.set_xticklabels([b[0] for b in bars], fontsize=9.2)
ax.set_ylim(0, 100)
ax.set_ylabel("Microbial Composition (%, 8 bacteria renormalized)", fontsize=11.5)
ax.yaxis.set_major_formatter(lambda v, pos: f"{int(v)}%")
ax.spines[["top", "right"]].set_visible(False)

ax.axvline(0.5, color="gray", linestyle="--", linewidth=1, ymin=0, ymax=1)

handles, labels = ax.get_legend_handles_labels()
ax.legend(handles[::-1], labels[::-1], title="Species", loc="upper left",
          bbox_to_anchor=(1.01, 1.0), frameon=False, fontsize=9.5, title_fontsize=11)

plt.tight_layout()
plt.savefig("./zymo_d6305_d6306_composition.png", bbox_inches="tight")
print("saved ./zymo_d6305_d6306_composition.png")
