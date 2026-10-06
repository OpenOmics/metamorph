#!/usr/bin/env python3
import json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

with open("./dehost_read_counts.json") as f:
    counts = json.load(f)

REPLICATES = ["ERR16012159","ERR16012160","ERR16012161","ERR16012162","ERR16012163","ERR16012164","ERR16012165"]

HG38_COLOR = "#2a78d6"
FUNGAL_COLOR = "#1baf7a"

hg38_removed_pct = []
fungal_removed_pct = []
for rep in REPLICATES:
    t, h, f = counts[rep]["trimmed"], counts[rep]["post_hg38"], counts[rep]["post_fungal"]
    hg38_removed_pct.append((t - h) / t * 100.0)
    fungal_removed_pct.append((h - f) / h * 100.0)

fig, axes = plt.subplots(1, 2, figsize=(13, 5.2), dpi=100)

x = np.arange(len(REPLICATES))

ax = axes[0]
bars = ax.bar(x, hg38_removed_pct, width=0.6, color=HG38_COLOR)
for b, v in zip(bars, hg38_removed_pct):
    ax.text(b.get_x() + b.get_width() / 2, v + 0.005, f"{v:.2f}%", ha="center", va="bottom", fontsize=9)
ax.set_title("hg38 (human) screen", fontsize=13)
ax.set_ylabel("Reads removed at this stage (%)", fontsize=11)
ax.set_ylim(0, max(hg38_removed_pct) * 1.35)
ax.set_xticks(x)
ax.set_xticklabels(REPLICATES, rotation=35, ha="right", fontsize=9)
ax.spines[["top", "right"]].set_visible(False)

ax = axes[1]
bars = ax.bar(x, fungal_removed_pct, width=0.6, color=FUNGAL_COLOR)
for b, v in zip(bars, fungal_removed_pct):
    ax.text(b.get_x() + b.get_width() / 2, v + 0.04, f"{v:.2f}%", ha="center", va="bottom", fontsize=9)
ax.set_title("Injected yeast screen", fontsize=13)
ax.set_ylabel("Reads removed at this stage (%)", fontsize=11)
ax.set_ylim(0, max(fungal_removed_pct) * 1.35)
ax.set_xticks(x)
ax.set_xticklabels(REPLICATES, rotation=35, ha="right", fontsize=9)
ax.spines[["top", "right"]].set_visible(False)

fig.suptitle("Reads removed by each dehosting stage, per replicate\n(the remainder in each sample is retained and carried forward)", fontsize=13, y=1.03)
plt.tight_layout()
plt.savefig("./zymo_dehost_removed_pct.png", bbox_inches="tight", facecolor="white")
print("saved")
