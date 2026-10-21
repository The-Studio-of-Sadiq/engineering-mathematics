"""Shared Matplotlib style for every figure in the book.

Each script in figures/src/ does:

    from style import plt, np, PALETTE, new_figure, save
    fig, ax = new_figure()
    ...
    save(fig, "ch26-growth-rates")      # writes images/ch26-growth-rates.svg

Output is deterministic SVG (no timestamps, fixed hash salt), so re-running a
script only changes the file when the drawing changes.
"""

import os
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

ROOT = Path(__file__).resolve().parent.parent
IMAGES = ROOT / "images"

# Blue, red, green, purple, orange, teal, grey: readable on white and in print.
PALETTE = ["#2563eb", "#dc2626", "#16a34a", "#9333ea", "#ea580c", "#0891b2", "#4b5563"]
INK = "#1f2937"
MUTED = "#6b7280"
FILL = "#e5edff"

plt.rcParams.update({
    "figure.figsize": (7.0, 4.2),
    "figure.dpi": 100,
    "figure.facecolor": "white",
    "axes.facecolor": "white",
    "axes.edgecolor": INK,
    "axes.labelcolor": INK,
    "axes.titlesize": 12,
    "axes.titleweight": "bold",
    "axes.labelsize": 10.5,
    "axes.grid": True,
    "grid.color": "#e5e7eb",
    "grid.linewidth": 0.8,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.prop_cycle": matplotlib.cycler(color=PALETTE),
    "xtick.color": INK,
    "ytick.color": INK,
    "font.family": "DejaVu Sans",
    "font.size": 10,
    "legend.frameon": False,
    "lines.linewidth": 2.0,
    "mathtext.fontset": "dejavusans",
    "svg.fonttype": "path",
    "svg.hashsalt": "mathematics-universe",
})


def new_figure(nrows=1, ncols=1, **kw):
    """plt.subplots with the book's defaults; pass figsize=(w, h) to override."""
    return plt.subplots(nrows, ncols, constrained_layout=True, **kw)


def diagram(figsize=(7.0, 4.2)):
    """A blank canvas (no axes, no grid) for node-and-arrow diagrams, in data units 0..10 x 0..6."""
    fig, ax = plt.subplots(figsize=figsize, constrained_layout=True)
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 6)
    ax.set_aspect("equal")
    ax.axis("off")
    return fig, ax


def node(ax, x, y, text, color=PALETTE[0], r=0.42, fc="white", fontsize=10, lw=2):
    """A labelled circle node."""
    ax.add_patch(plt.Circle((x, y), r, fc=fc, ec=color, lw=lw, zorder=3))
    ax.text(x, y, text, ha="center", va="center", fontsize=fontsize, color=INK, zorder=4)


def box(ax, x, y, text, w=2.2, h=0.8, color=PALETTE[0], fc="white", fontsize=9.5, lw=1.8):
    """A labelled rounded box centred at (x, y)."""
    from matplotlib.patches import FancyBboxPatch
    ax.add_patch(FancyBboxPatch((x - w / 2, y - h / 2), w, h, boxstyle="round,pad=0.02,rounding_size=0.12",
                                fc=fc, ec=color, lw=lw, zorder=3))
    ax.text(x, y, text, ha="center", va="center", fontsize=fontsize, color=INK, zorder=4, wrap=True)


def arrow(ax, x0, y0, x1, y1, color=MUTED, text=None, lw=1.6, shrink=24, curve=0.0, fontsize=9):
    """An arrow from (x0, y0) to (x1, y1), optionally curved and labelled at its midpoint."""
    ax.annotate("", xy=(x1, y1), xytext=(x0, y0), zorder=2,
                arrowprops=dict(arrowstyle="-|>", color=color, lw=lw, shrinkA=shrink, shrinkB=shrink,
                                connectionstyle=f"arc3,rad={curve}"))
    if text:
        mx, my = (x0 + x1) / 2, (y0 + y1) / 2
        ax.text(mx, my + 0.18, text, ha="center", va="bottom", fontsize=fontsize, color=MUTED, zorder=4)


def _rasterize_dense(fig, max_paths=400, max_points=8000):
    """Rasterize heavy layers (meshes, contour fills, big scatters, very long lines) so the SVG stays
    small; text, axes and ordinary lines remain vector."""
    from matplotlib.collections import Collection, QuadMesh
    from matplotlib.lines import Line2D
    for ax in fig.axes:
        for a in ax.get_children():
            if isinstance(a, QuadMesh):
                a.set_rasterized(True)
            elif isinstance(a, Collection) and (
                    len(a.get_paths()) + len(a.get_offsets()) > max_paths
                    or sum(len(p.vertices) for p in a.get_paths()) > 5 * max_points):
                a.set_rasterized(True)
            elif isinstance(a, Line2D) and (
                    len(a.get_xdata()) > max_points
                    or (a.get_marker() not in (None, "None", "", " ") and len(a.get_xdata()) > max_paths)):
                a.set_rasterized(True)


def save(fig, name):
    """Write images/<name>.svg deterministically and close the figure."""
    IMAGES.mkdir(exist_ok=True)
    _rasterize_dense(fig)
    fig.savefig(IMAGES / f"{name}.svg", format="svg", bbox_inches="tight", pad_inches=0.08, dpi=200,
                metadata={"Date": None, "Creator": None})
    preview = os.environ.get("FIG_PREVIEW")  # set to a folder to also get PNGs for eyeballing
    if preview:
        Path(preview).mkdir(parents=True, exist_ok=True)
        fig.savefig(Path(preview) / f"{name}.png", dpi=110, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote images/{name}.svg")
