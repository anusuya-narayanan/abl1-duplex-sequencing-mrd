#!/usr/bin/env python3
"""Compare duplex (DCS) and single-strand (SSCS) consensus variant tables.

Inputs are the two Variant Annotator tables exported from Galaxy
(one row per position: base counts per read direction, CVRG, MAF ...).

Outputs (written to --out and --fig):
  variant_calls.csv          sites passing the tutorial filter (MAF >= 0.01)
  threshold_sweep.csv        sites passing each MAF cut-off, with and without
                             a minimum number of supporting reads
  background_error.csv       noise reads per base, after removing real sites
  substitution_spectrum.csv  noise by substitution class (strand-collapsed)
  threshold_sweep.png, error_spectrum.png

Usage:
  python scripts/duplex_sweep.py \
      --dcs data/dcs_variant_annotator.tsv \
      --sscs data/sscs_variant_annotator.tsv
"""
import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

BASES = ["A", "C", "G", "T"]
COUNT_COLS = [f"+{b}" for b in BASES] + [f"-{b}" for b in BASES]
CUTOFFS = [0.0005, 0.001, 0.002, 0.003, 0.005, 0.01, 0.02]
# Strand-collapsed substitution classes: a change and its reverse complement
# are the same event seen from the other strand.
CLASSES = {
    "C>A": "C>A / G>T", "G>T": "C>A / G>T",
    "C>G": "C>G / G>C", "G>C": "C>G / G>C",
    "C>T": "C>T / G>A", "G>A": "C>T / G>A",
    "T>A": "T>A / A>T", "A>T": "T>A / A>T",
    "T>C": "T>C / A>G", "A>G": "T>C / A>G",
    "T>G": "T>G / A>C", "A>C": "T>G / A>C",
}
CLASS_ORDER = ["C>A / G>T", "C>G / G>C", "C>T / G>A",
               "T>A / A>T", "T>C / A>G", "T>G / A>C"]

# Chart styling (validated two-series palette; light surface)
COLORS = {"DCS": "#2a78d6", "SSCS": "#eb6834"}
SURFACE, INK, INK2, MUTED, GRID, AXIS = (
    "#fcfcfb", "#0b0b0b", "#52514e", "#898781", "#e1e0d9", "#c3c2b7")


def load(path):
    """Read a Variant Annotator table and add per-position allele counts."""
    df = pd.read_csv(path, sep="\t")
    df.columns = [c.lstrip("#") for c in df.columns]
    counts = np.column_stack([df[f"+{b}"] + df[f"-{b}"] for b in BASES])
    order = np.argsort(-counts, axis=1, kind="stable")  # most common base first
    rows = np.arange(len(df))
    df["major_base"] = np.array(BASES)[order[:, 0]]
    df["minor_base"] = np.array(BASES)[order[:, 1]]
    df["minor_reads"] = counts[rows, order[:, 1]]
    df["MAF"] = pd.to_numeric(df["MAF"], errors="coerce").fillna(0.0)
    return df


def variant_calls(df, label, cutoff=0.01):
    hits = df[df["MAF"] >= cutoff]
    return pd.DataFrame({
        "consensus": label,
        "chrom": hits["CHR"], "pos": hits["POS"],
        "ref_major": hits["major_base"], "alt_minor": hits["minor_base"],
        "alt_reads": hits["minor_reads"], "depth": hits["CVRG"],
        "maf": hits["MAF"],
    })


def sweep(df, label, min_reads):
    rows = []
    for t in CUTOFFS:
        passing = df["MAF"] >= t
        rows.append({
            "consensus": label, "maf_cutoff": t,
            "sites_maf_only": int(passing.sum()),
            f"sites_maf_and_min{min_reads}_reads":
                int((passing & (df["minor_reads"] >= min_reads)).sum()),
        })
    return pd.DataFrame(rows)


def background(df, label, real_positions, germline_maf):
    """Noise = every minor-allele read left after removing real sites.

    Real sites are the known variant positions plus anything at or above
    germline_maf (heterozygous SNPs sit near 0.5).
    """
    noise = df[(df["MAF"] < germline_maf) & (~df["POS"].isin(real_positions))]
    reads, bases = int(noise["minor_reads"].sum()), int(noise["CVRG"].sum())
    summary = {
        "consensus": label,
        "positions": len(noise),
        "positions_with_noise": int((noise["minor_reads"] > 0).sum()),
        "noise_reads": reads, "bases_read": bases,
        "noise_per_base": reads / bases if bases else float("nan"),
    }
    sub = noise[noise["minor_reads"] > 0].copy()
    sub["class"] = (sub["major_base"] + ">" + sub["minor_base"]).map(CLASSES)
    spec = (sub.groupby("class")["minor_reads"].sum()
            .reindex(CLASS_ORDER, fill_value=0).rename("noise_reads").reset_index())
    spec.insert(0, "consensus", label)
    spec["per_100k_bases"] = spec["noise_reads"] / bases * 1e5
    return summary, spec


def _style(ax):
    ax.set_facecolor(SURFACE)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(AXIS)
    ax.tick_params(colors=MUTED, labelsize=9)
    ax.yaxis.grid(True, color=GRID, linewidth=0.8)
    ax.set_axisbelow(True)


def plot_sweep(sw, min_reads, path):
    col = f"sites_maf_and_min{min_reads}_reads"
    panels = [("sites_maf_only", "MAF cut-off only")]
    if col in sw and sw[col].notna().all():
        panels.append((col, f"MAF cut-off and at least {min_reads} supporting reads"))
    fig, axes = plt.subplots(1, len(panels), figsize=(5.2 * len(panels), 4.2),
                             sharey=True, facecolor=SURFACE)
    axes = axes if len(panels) > 1 else [axes]
    for ax, (column, title) in zip(axes, panels):
        _style(ax)
        for label, color in COLORS.items():
            d = sw[sw["consensus"] == label]
            ax.plot(d["maf_cutoff"], d[column], color=color, linewidth=2,
                    marker="o", markersize=6, markeredgecolor=SURFACE,
                    markeredgewidth=1.5, label=label)
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xticks(CUTOFFS)
        ax.set_xticklabels([f"{c:g}" for c in CUTOFFS], rotation=0)
        ax.minorticks_off()
        ticks = [1, 2, 5, 10, 20, 50, 100, 200, 500, 1000, 2000]
        top = max(160, sw[[c for c, _ in panels]].max().max() * 1.4)
        ax.set_yticks([t for t in ticks if t <= top])
        ax.set_yticklabels([str(t) for t in ticks if t <= top])
        ax.set_ylim(0.8, top)
        ax.set_xlabel("MAF cut-off (log scale)", color=INK2, fontsize=10)
        ax.set_title(title, color=INK, fontsize=11, loc="left")
    axes[0].set_ylabel("Sites passing (log scale)", color=INK2, fontsize=10)
    axes[0].legend(frameon=False, fontsize=9, loc="upper right")
    fig.suptitle("Sites called at each MAF cut-off, DCS vs SSCS",
                 x=0.01, ha="left", color=INK, fontsize=12, fontweight="bold")
    fig.tight_layout()
    fig.savefig(path, dpi=200, facecolor=SURFACE)
    plt.close(fig)


def plot_spectrum(spec, path):
    fig, ax = plt.subplots(figsize=(8, 4.2), facecolor=SURFACE)
    _style(ax)
    width, gap = 0.38, 0.02
    x = range(len(CLASS_ORDER))
    for i, (label, color) in enumerate(COLORS.items()):
        d = spec[spec["consensus"] == label].set_index("class").reindex(CLASS_ORDER)
        offs = [xi + (i - 0.5) * (width + gap) for xi in x]
        ax.bar(offs, d["per_100k_bases"], width=width, color=color, label=label)
        # Label the oxidative-damage class (the headline difference) for both series
        v = d.loc[CLASS_ORDER[0], "per_100k_bases"]
        ax.annotate(f"{v:.1f}", (offs[0], v), textcoords="offset points",
                    xytext=(0, 3), ha="center", fontsize=9, color=INK2)
    ax.set_xticks(list(x))
    ax.set_xticklabels(CLASS_ORDER, fontsize=9)
    ax.set_ylabel("Noise reads per 100,000 bases", color=INK2, fontsize=10)
    ax.legend(frameon=False, fontsize=9)
    ax.set_title("Background noise by substitution type, DCS vs SSCS",
                 color=INK, fontsize=12, fontweight="bold", loc="left")
    fig.tight_layout()
    fig.savefig(path, dpi=200, facecolor=SURFACE)
    plt.close(fig)


def main():
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("--dcs", required=True, help="Variant Annotator table, DCS")
    p.add_argument("--sscs", required=True, help="Variant Annotator table, SSCS")
    p.add_argument("--out", default="results")
    p.add_argument("--fig", default="figures")
    p.add_argument("--min-reads", type=int, default=3)
    p.add_argument("--real-positions", type=int, nargs="*", default=[130872141],
                   help="known variant positions to leave out of the noise estimate")
    p.add_argument("--germline-maf", type=float, default=0.3)
    a = p.parse_args()

    out, fig = Path(a.out), Path(a.fig)
    out.mkdir(parents=True, exist_ok=True)
    fig.mkdir(parents=True, exist_ok=True)

    tables = {"DCS": load(a.dcs), "SSCS": load(a.sscs)}
    calls = pd.concat([variant_calls(t, k) for k, t in tables.items()])
    sw = pd.concat([sweep(t, k, a.min_reads) for k, t in tables.items()])
    bg, spec = zip(*[background(t, k, a.real_positions, a.germline_maf)
                     for k, t in tables.items()])
    bg, spec = pd.DataFrame(bg), pd.concat(spec)

    calls.to_csv(out / "variant_calls.csv", index=False)
    sw.to_csv(out / "threshold_sweep.csv", index=False)
    bg.to_csv(out / "background_error.csv", index=False)
    spec.to_csv(out / "substitution_spectrum.csv", index=False)
    plot_sweep(sw, a.min_reads, fig / "threshold_sweep.png")
    plot_spectrum(spec, fig / "error_spectrum.png")

    print(calls.to_string(index=False), "\n")
    print(sw.to_string(index=False), "\n")
    print(bg.to_string(index=False))


if __name__ == "__main__":
    main()
