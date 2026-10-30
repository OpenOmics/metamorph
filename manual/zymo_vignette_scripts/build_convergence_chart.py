#!/usr/bin/env python3
import json
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

with open("./zymo_convergence_results.json") as f: cum = json.load(f)
with open("./zymo_single_rep_bc.json") as f: single = json.load(f)

n = np.arange(1, 8)
fig, ax = plt.subplots(figsize=(9.5, 6.2), dpi=100)

ax.plot(n, cum["read_based"], "-o", color="#004586", linewidth=2.5, markersize=8,
        label="metamorph (read-based): cumulative mean vs. target")
ax.scatter(n, single["read_based"], marker="x", color="#004586", s=90, linewidths=2.2,
           label="metamorph (read-based): single replicate vs. target")

ax.plot(n, cum["assembly_based"], "-o", color="#579D1C", linewidth=2.5, markersize=8,
        label="metamorph (assembly-based): cumulative mean vs. target")
ax.scatter(n, single["assembly_based"], marker="x", color="#579D1C", s=90, linewidths=2.2,
           label="metamorph (assembly-based): single replicate vs. target")

ax.set_xlabel("Number of replicates (n)", fontsize=12)
ax.set_ylabel("Bray–Curtis dissimilarity to vendor target", fontsize=12)
ax.set_title("Does metamorph's output converge\ntoward the vendor-defined composition?", fontsize=15)
ax.set_xticks(n)
ax.set_ylim(0, 0.5)
ax.spines[["top","right"]].set_visible(False)
ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.14), ncol=1, frameon=True, fontsize=9.5)

plt.tight_layout()
plt.savefig("./zymo_convergence_NEW.png", bbox_inches="tight")
print("saved")
