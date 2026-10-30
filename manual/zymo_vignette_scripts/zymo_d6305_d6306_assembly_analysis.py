#!/usr/bin/env python3
# metamorph's own in-run assembly-based estimate for D6305/D6306 (in3634_1/in3634_2).
# Same coverage-normalized MAG-mapping method as zymo_fulldepth_analysis2.py (D6300
# page): coverage = unambiguousMB * 1e6 / genome_size, summed per vendor label,
# renormalized to 100% across all 8 target bacteria.
#
# With only 2 samples, this run's own dRep dereplication (post metaWRAP bin_refinement
# QC) recovered genomes for only 5 of the 8 target bacteria outright
# (metawrap_bin_refine/dRep/dRep.SMALL_BINSET_WARNING is present). Two more categories
# are added by hand, each for a documented, verified reason -- neither is something
# `metamorph run` itself produces (see Section 5.1 of the vignette):
#
# 1. Enterococcus faecalis (in3634_2_bin.8, metawrap_binning/in3634_2/metabat2_bins/):
#    41.19% complete, 0% contamination -- confirmed by Mash against the vendor's own
#    reference genome (distance 0.034, far closer than any other reference) to be
#    real, uncontaminated E. faecalis sequence, just incomplete. Reported as its own
#    species, same as the other 5.
#
# 2. Escherichia coli + Salmonella enterica, COMBINED, not separable (in3634_2_bin.9,
#    same binning dir): 99.27% complete but 97.71% contamination -- Mash confirms this
#    bin is almost exactly equidistant between the E. coli and Salmonella reference
#    genomes, i.e. a genuine ~50/50 mixture of both organisms' contigs, not noise.
#    There is no per-contig classification step in metaWRAP/metamorph to split this
#    bin's coverage between the two species individually, so it is reported as ONE
#    combined category -- "E. coli + S. enterica (combined)" -- compared against the
#    SUM of both organisms' individual vendor targets, not reported as either
#    species's own number. This is NOT the same claim as "recovered both species" --
#    it's "recovered their combined mass, could not tell them apart."
#
# Both rescued bins were added to a combined 7-reference BBSplit index (the 5 real
# dRep winners + these 2) and both samples' dehosted reads remapped with the exact
# same bbsplit.sh parameters metamorph's own `bbtools_index_map` rule uses
# (minid=0.90, sortscafs=f, nzo=f) -- see ./d6305_d6306_asm7_mags/{sample}.refstats.
# GTDB-Tk never saw either rescued bin (neither reached that rule), so their vendor-
# label mapping below is set by hand from the Mash identification.
#
# Coverage denominator for both rescued bins is each one's own assembled size (same
# convention as every other genome here -- dRep's Widb.csv "size" column is likewise
# each winner's *assembled* size, not an independently-known true size): 957,885 bp
# for the E. faecalis bin, 9,120,732 bp for the combined E. coli/Salmonella bin (close
# to the sum of their two real genome sizes, ~4.9+4.8 Mb, consistent with this bin
# holding near-complete copies of both).
import re, json
import pandas as pd

BASE = "/data/OpenOmics/dev/metamorph/remove_cruft_at_end/zymo/output_zymo_data/metagenome_results"
RESCUE_MAGS_DIR = "./d6305_d6306_asm7_mags"  # 7-genome refstats (5 dRep winners + 2 rescued bins)
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
VENDOR_GENOME_COPY = {
    "Pseudomonas aeruginosa": 6.1, "Escherichia coli": 8.5, "Salmonella enterica": 8.7,
    "Lactobacillus fermentum": 21.6, "Enterococcus faecalis": 14.6, "Staphylococcus aureus": 15.2,
    "Listeria monocytogenes": 13.9, "Bacillus subtilis": 10.3,
}
vendor_8_sum = sum(VENDOR_GENOME_COPY.values())
TARGET_RENORM_8 = {s: VENDOR_GENOME_COPY[s] / vendor_8_sum * 100 for s in BACTERIA}

COMBINED_LABEL = "Escherichia coli + Salmonella enterica (combined, not separable)"

widb = pd.read_csv(f"{BASE}/metawrap_bin_refine/dRep/data_tables/Widb.csv")
genome_size = dict(zip(widb["genome"].str.replace(".fa", ""), widb["size"]))
EFAECALIS_BIN = "in3634_2_bin.8"
CHIMERA_BIN = "in3634_2_bin.9"
genome_size[EFAECALIS_BIN] = 957885
genome_size[CHIMERA_BIN] = 9120732

gdf = pd.read_csv(f"{BASE}/metawrap_kmer/GTDBTK_classify_wf/gtdbtk.bac120.summary.tsv", sep="\t")
def species_from_classification(c):
    m = re.search(r"s__(.+)$", c)
    return m.group(1).strip() if m else None
gdf["species_raw"] = gdf["classification"].apply(species_from_classification)
def to_vendor_label(species_raw):
    norm = species_raw.replace(" ", "_")
    for disp, base in BACTERIA.items():
        pat = re.compile(rf"^{re.escape(base)}(_[A-Z]{{1,2}})?$")
        if pat.match(norm):
            return disp
    return None
gdf["vendor_label"] = gdf["species_raw"].apply(to_vendor_label)
mag_to_label = dict(zip(gdf["user_genome"], gdf["vendor_label"]))
mag_to_label[EFAECALIS_BIN] = "Enterococcus faecalis"  # hand-identified via Mash, Section 5.1
mag_to_label[CHIMERA_BIN] = COMBINED_LABEL              # hand-identified via Mash, Section 5.1
print("MAG -> category map (GTDB-Tk output + 2 hand-rescued bins):")
print(json.dumps(mag_to_label, indent=2))

CATEGORIES = sorted(set(mag_to_label.values()))
print(f"\nCategories ({len(CATEGORIES)}), covering all 8 target bacteria: {CATEGORIES}")

TARGET_RENORM = {s: TARGET_RENORM_8[s] for s in BACTERIA if s not in ("Escherichia coli", "Salmonella enterica")}
TARGET_RENORM[COMBINED_LABEL] = TARGET_RENORM_8["Escherichia coli"] + TARGET_RENORM_8["Salmonella enterica"]

asm_results = {}
for sample in SAMPLES:
    df = pd.read_csv(f"{RESCUE_MAGS_DIR}/{sample}.refstats", sep="\t")
    df.columns = [c.lstrip("#") for c in df.columns]
    cov_raw = {s: 0.0 for s in CATEGORIES}
    for _, row in df.iterrows():
        label = mag_to_label.get(row["name"])
        if label is None:
            continue
        gsize = genome_size[row["name"]]
        cov_raw[label] += row["unambiguousMB"] * 1e6 / gsize
    cov_sum = sum(cov_raw.values())
    sample_out = {}
    for s in CATEGORIES:
        pct = cov_raw[s] / cov_sum * 100.0 if cov_sum > 0 else 0.0
        sample_out[s] = dict(pct=pct, target=TARGET_RENORM[s], fold=pct / TARGET_RENORM[s])
    asm_results[sample] = sample_out

print(f"\n{'Category':55s} {'target%':>8s} {'D6305 asm %':>12s} {'fold':>7s} {'D6306 asm %':>12s} {'fold':>7s}")
for s in CATEGORIES:
    d5 = asm_results["in3634_1"][s]
    d6 = asm_results["in3634_2"][s]
    print(f"{s:55s} {d5['target']:8.2f} {d5['pct']:12.2f} {d5['fold']:6.2f}x {d6['pct']:12.2f} {d6['fold']:6.2f}x")

with open("./zymo_d6305_d6306_assembly_results.json", "w") as f:
    json.dump({"categories": CATEGORIES, "combined_label": COMBINED_LABEL,
               "target_renorm": TARGET_RENORM, "assembly_based": asm_results}, f, indent=2)
print("\nWrote ./zymo_d6305_d6306_assembly_results.json")
