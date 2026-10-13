# PRA2003 - Week 5: WT - mutant difference as a function of momentum (Python version)
#
# Python version of momentumDifference.R. It reads the same saved counts
# (momentum_perfile_results.csv, momentum_bin_edges.csv and
# subsample_results.csv), does the same calculation and makes the same plots
# with matplotlib.
#
# Method (same as the R script):
#   1. For each file, pair and momentum bin, the difference per event:
#        D = count_WT / nEvents - count_mutant / nEvents
#   2. Per pair and bin: point = mean of the per-file values of D,
#      uncertainty = standard deviation of those values (sample SD, n - 1),
#      significance = point / uncertainty.
#   3. Per pair, a weighted fit of a constant through the 5 bins:
#        c = sum(w * D) / sum(w),  uncertainty = 1 / sqrt(sum(w)),  w = 1 / sigma^2
#        chi2 = sum(((D - c) / sigma)^2),  ndf = 5 - 1
#   Set 10 is left out of the main result because its momentum distribution
#   is very different from Sets 1-9. The same analysis with all 10 sets is
#   also made, for comparison.
#
# HOW TO RUN (in the folder with the three CSV files):
#   python3 momentumDifference.py
#
# Writes:
#   momentum_difference_python.png           - Sets 1-9
#   momentum_difference_all_sets_python.png  - all 10 sets

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import chi2 as chi2_distribution

excluded_sets = [10]   # left out of the main result; use [] to keep all

pair_names = ["E. coli", "Bacillus subtilis", "Pseudomonas aeruginosa",
              "Streptococcus pneumoniae", "Mycobacterium tuberculosis", "Salmonella"]

# ---- Read the saved counts, bin edges and number of events per file ----
counts = pd.read_csv("momentum_perfile_results.csv")
edges = pd.read_csv("momentum_bin_edges.csv")["Edge"].to_numpy()
events = pd.read_csv("subsample_results.csv")[["Set", "nEvents"]].drop_duplicates()

n_bins = len(edges) - 1

# Labels for the momentum bins, e.g. "0.65-1.32" and ">5.87"
bin_labels = []
for b in range(1, n_bins + 1):
    if b == 1:
        bin_labels.append(f"<{edges[1]:.2f}")
    elif b == n_bins:
        bin_labels.append(f">{edges[n_bins - 1]:.2f}")
    else:
        bin_labels.append(f"{edges[b - 1]:.2f}-{edges[b]:.2f}")

# Step 1: difference per event for every file, pair and bin
counts = counts.merge(events, on="Set")
counts["D"] = (counts["CountWT"] - counts["CountMutant"]) / counts["nEvents"]


def analyse_sets(sets_used, file_name, title):
    data = counts[counts["Set"].isin(sets_used)]

    fig, axes = plt.subplots(2, 3, figsize=(12.86, 7.86))
    print(f"\n==== {title} ====")
    print(f"{'Strain':27s} {'constant':>11s} {'uncertainty':>12s} {'sigma':>6s} {'chi2/ndf':>10s}")

    for k, pair in enumerate(pair_names):
        ax = axes.flatten()[k]
        pair_data = data[data["Strain"] == pair]

        # Step 2: mean and SD over the files, per bin
        per_bin = pair_data.groupby("MomentumBin")["D"]
        x = per_bin.mean().index.to_numpy()
        y = per_bin.mean().to_numpy()
        sigma = per_bin.std().to_numpy()       # sample SD (n - 1), the same as sd() in R

        # Step 3: weighted fit of a constant
        w = 1 / sigma**2
        constant = np.sum(w * y) / np.sum(w)
        constant_uncertainty = 1 / np.sqrt(np.sum(w))
        significance = constant / constant_uncertainty
        chi2 = np.sum(((y - constant) / sigma)**2)
        ndf = n_bins - 1
        p_value = chi2_distribution.sf(chi2, ndf)

        print(f"{pair:27s} {constant:11.3g} {constant_uncertainty:12.3g} "
              f"{significance:6.2f} {chi2:6.2f}/{ndf} (p = {p_value:.2f})")

        # Plot: fitted constant with its +- 1 sigma band, points +- SD
        ax.axhspan(constant - constant_uncertainty, constant + constant_uncertainty,
                   color="red", alpha=0.15, linewidth=0)
        ax.axhline(constant, color="red", linewidth=2)
        ax.axhline(0, color="grey", linestyle="--", linewidth=1)   # no difference
        ax.errorbar(x, y, yerr=sigma, fmt="o", color="darkblue", ecolor="black",
                    elinewidth=1, capsize=3)

        ax.set_title(f"{pair}\nconstant = {constant:.3g} +/- {constant_uncertainty:.3g} "
                     f"({significance:.1f} sigma)\nchi2/ndf = {chi2:.1f}/{ndf}",
                     fontsize=9)
        ax.set_xticks(x)
        ax.set_xticklabels(bin_labels, fontsize=7)
        ax.set_xlim(0.7, n_bins + 0.3)
        ax.set_xlabel(r"Momentum bin |p| ($10^{-20}$ kg m/s)", fontsize=8)
        ax.set_ylabel("Difference per event (WT - mutant)", fontsize=8)

    fig.suptitle(f"Difference per event vs momentum ({title})", fontweight="bold")
    fig.tight_layout()
    fig.savefig(file_name, dpi=140)
    plt.close(fig)
    print(f"Saved {file_name}")


all_sets = sorted(counts["Set"].unique())
main_sets = [s for s in all_sets if s not in excluded_sets]

# Main result: Sets 1-9
analyse_sets(main_sets, "momentum_difference_python.png",
             "sets " + ",".join(str(s) for s in main_sets))

# Comparison: all 10 sets
analyse_sets(all_sets, "momentum_difference_all_sets_python.png",
             "sets " + ",".join(str(s) for s in all_sets))
