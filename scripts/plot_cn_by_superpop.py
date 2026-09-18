#!/usr/bin/env python3
"""Copy-number composition by 1000 Genomes superpopulation, as a figure.

    python3 scripts/plot_cn_by_superpop.py results/kgp2504/lilra6_lilra3_cn.tsv \\
        --outdir docs/figures

Writes `cn_by_superpopulation.png` and a `-dark.png` twin, because a README is
read on both surfaces and an inverted light figure is not a dark-mode figure.

One decision worth stating, since it costs something. Copy number is an *ordered*
scale, and the textbook fill for it is a one-hue ramp — light for fewer copies,
dark for more — which makes the ordering legible from the colour alone. These
panels use **categorical hues** instead, which trades that cue away and buys two
things back: every copy number gets its own class rather than being folded into a
">= 4" tail (a one-hue ramp fits five steps between its ends at the >= 0.06 OKLCH L
separation the validator wants, not seven), and adjacent slices stay far apart for
a colourblind reader — worst adjacent pair dE 9.1 light / 8.4 dark, against a
target of 8. The legend is ordered, and the panels are read left to right.

LILRB3 is not plotted. 2,500 of 2,504 samples are CN 2, so its panel would be one
flat colour and would say only that the gene is stable — which the caption can say
in fewer pixels.
"""

from __future__ import annotations

import argparse
import collections
import csv
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.patches import Patch  # noqa: E402

# The data-viz skill's categorical slots, in its fixed order — never re-ordered and
# never cycled, so copy number N is the same hue in both panels and in any later
# figure. Validated on the adjacent pairlist (the one stacked bars use) in both
# modes: lightness band, chroma floor, CVD and normal-vision separation all pass.
# Light mode warns on contrast for aqua/yellow/magenta, whose documented relief is
# visible direct labels — which is why every segment wide enough carries its value.
PALETTE = {
    "light": ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300",
              "#4a3aa7"],
    "dark": ["#3987e5", "#d95926", "#199e70", "#c98500", "#d55181", "#008300",
             "#9085e9"],
}
INK = {
    "light": dict(surface="#fcfcfb", primary="#0b0b0b", secondary="#52514e",
                  muted="#898781"),
    "dark": dict(surface="#1a1a19", primary="#ffffff", secondary="#c3c2b7",
                 muted="#898781"),
}
SUPERPOPS = ["AFR", "AMR", "EAS", "EUR", "SAS"]
GENES = [("lilra6_cn", "LILRA6", [0, 1, 2, 3, 4, 5, 6],
          ["0", "1", "2", "3", "4", "5", "6"]),
         ("lilra3_cn", "LILRA3", [0, 1, 2], ["0", "1", "2"])]


def _text_on(fill: str) -> str:
    """Ink that survives on this fill — these hues straddle the readable boundary."""
    channels = [int(fill[i:i + 2], 16) / 255 for i in (1, 3, 5)]
    linear = [c / 12.92 if c <= 0.04045 else ((c + 0.055) / 1.055) ** 2.4
              for c in channels]
    lum = 0.2126 * linear[0] + 0.7152 * linear[1] + 0.0722 * linear[2]
    return "#0b0b0b" if lum > 0.2 else "#ffffff"


def load(path: Path) -> dict:
    rows = list(csv.DictReader(open(path), delimiter="\t"))
    counts = {gene: collections.defaultdict(collections.Counter)
              for gene, _, _, _ in GENES}
    n = collections.Counter()
    for row in rows:
        sup = row["superpopulation"]
        if sup not in SUPERPOPS:
            continue
        n[sup] += 1
        for gene, _, classes, _ in GENES:
            if not row[gene]:
                continue
            cn = int(row[gene])
            if cn not in classes:
                continue
            counts[gene][sup][cn] += 1
    return counts, n


def draw(counts, n, mode: str, out: Path) -> None:
    ink, ramp = INK[mode], PALETTE[mode]
    fig, axes = plt.subplots(1, 2, figsize=(9.6, 4.3), dpi=160,
                             facecolor=ink["surface"])
    fig.subplots_adjust(left=0.07, right=0.985, top=0.78, bottom=0.17, wspace=0.28)

    for ax, (gene, title, classes, labels) in zip(axes, GENES):
        ax.set_facecolor(ink["surface"])
        bottoms = [0.0] * len(SUPERPOPS)
        for idx, cn in enumerate(classes):
            share = [100 * counts[gene][sup][cn] / n[sup] for sup in SUPERPOPS]
            fill = ramp[idx]
            # A 2px surface edge is the gap between stacked segments.
            ax.bar(SUPERPOPS, share, bottom=bottoms, width=0.68, color=fill,
                   edgecolor=ink["surface"], linewidth=2, zorder=2)
            for x, (value, base) in enumerate(zip(share, bottoms)):
                if value >= 7:        # below this the label would not fit the segment
                    ax.text(x, base + value / 2, f"{value:.0f}", ha="center",
                            va="center", fontsize=8.5, color=_text_on(fill),
                            zorder=3)
            bottoms = [b + s for b, s in zip(bottoms, share)]

        ax.set_title(title, fontsize=12, color=ink["primary"], pad=8,
                     loc="left", fontweight="bold")
        ax.set_ylim(0, 100)
        ax.set_yticks([0, 25, 50, 75, 100])
        ax.set_yticklabels(["0", "25", "50", "75", "100%"], fontsize=8.5,
                           color=ink["muted"])
        ax.set_xticks(range(len(SUPERPOPS)))
        ax.set_xticklabels([f"{s}\n$n$={n[s]}" for s in SUPERPOPS], fontsize=9,
                           color=ink["secondary"])
        ax.tick_params(length=0)
        for side in ("top", "right", "left"):
            ax.spines[side].set_visible(False)
        ax.spines["bottom"].set_color(ink["muted"])
        ax.spines["bottom"].set_linewidth(0.8)

    fig.suptitle("LILR copy number by superpopulation", x=0.07, y=0.955,
                 ha="left", fontsize=14, color=ink["primary"], fontweight="bold")
    # One legend for the figure, not one per panel: copy number is a single scale
    # across both, and LILRA3 reaching no further than 2 is itself the finding.
    widest = max(GENES, key=lambda g: len(g[2]))
    legend = fig.legend(handles=[Patch(facecolor=ramp[i], label=lab)
                                 for i, lab in enumerate(widest[3])],
                        loc="upper right", bbox_to_anchor=(0.985, 0.90),
                        ncol=len(widest[3]), frameon=False, fontsize=9,
                        handlelength=0.9, handleheight=0.9, handletextpad=0.45,
                        columnspacing=1.1, labelcolor=ink["secondary"],
                        title="copies", title_fontsize=9, alignment="left")
    legend.get_title().set_color(ink["muted"])
    fig.text(0.07, 0.035,
             "Share of samples at each diploid copy number, 2,504 1000 Genomes "
             "samples. LILRB3 is omitted: 2,500 of the 2,504 carry 2 copies.",
             fontsize=8.5, color=ink["muted"], ha="left")

    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, facecolor=ink["surface"])
    plt.close(fig)
    print(f"wrote {out}")


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("calls", type=Path,
                   help="lilra6_lilra3_cn.tsv (sample, superpopulation, per-gene CN)")
    p.add_argument("--outdir", type=Path, default=Path("docs/figures"))
    p.add_argument("--stem", default="cn_by_superpopulation")
    args = p.parse_args()

    counts, n = load(args.calls)
    draw(counts, n, "light", args.outdir / f"{args.stem}.png")
    draw(counts, n, "dark", args.outdir / f"{args.stem}-dark.png")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
