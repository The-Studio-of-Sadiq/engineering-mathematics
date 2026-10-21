"""GitHub social preview card (1280x640 PNG): images/social-preview.png."""

import json
import sys
import pathlib

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
from style import plt, np, PALETTE, ROOT, IMAGES  # noqa: E402

BOOK = json.loads((ROOT / "tools" / "book.json").read_text(encoding="utf-8"))
TRACK_COLORS = dict(zip("ABCDEFGHIJ", PALETTE + ["#be185d", "#0d9488", "#a16207"]))
TRACK_COLORS["G"] = "#9ca3af"  # the print grey is too dark on this background
BG, INK_DARK, SUB = "#0f172a", "#f8fafc", "#94a3b8"

fig = plt.figure(figsize=(12.8, 6.4), dpi=100, facecolor=BG)
ax = fig.add_axes([0, 0, 1, 1])
ax.set_xlim(0, 12.8)
ax.set_ylim(0, 6.4)
ax.axis("off")
ax.set_facecolor(BG)

# Faint background: a few curves from the book (sinusoid, decay, Gaussian).
x = np.linspace(0, 12.8, 600)
ax.plot(x, 0.55 + 0.18 * np.sin(2.2 * x), color="#1e293b", lw=2)
ax.plot(x, 5.95 + 0.25 * np.exp(-0.3 * x) * np.cos(4 * x), color="#1e293b", lw=2)

# Title block.
ax.text(0.7, 4.55, "The Mathematics", fontsize=46, weight="bold", color=INK_DARK, va="baseline")
ax.text(0.7, 3.65, "Universe", fontsize=46, weight="bold", color="#60a5fa", va="baseline")
ax.text(0.72, 3.05, BOOK["subtitle"], fontsize=17, color=SUB, va="baseline")
stats = "25 core chapters  ·  74 track chapters  ·  61 proofs  ·  1,400-page PDF"
ax.text(0.72, 2.35, stats, fontsize=13.5, color="#cbd5e1", va="baseline")

# Track chips.
for k, t in enumerate(BOOK["tracks"]):
    cx, cy = 0.95 + k * 0.56, 1.55
    ax.add_patch(plt.Circle((cx, cy), 0.22, fc=BG, ec=TRACK_COLORS[t["id"]], lw=2.4))
    ax.text(cx, cy, t["id"], ha="center", va="center", fontsize=12, weight="bold", color=TRACK_COLORS[t["id"]])
ax.text(0.72, 0.95, "Markdown + LaTeX  ·  Obsidian + Desmos  ·  CC BY 4.0", fontsize=12, color=SUB)

# Mini universe map on the right.
cx, cy = 10.2, 3.2
ax.add_patch(plt.Circle((cx, cy), 2.45, fc="none", ec=TRACK_COLORS["J"], lw=1.8, ls="--"))
ax.add_patch(plt.Circle((cx, cy), 1.95, fc="#111c33", ec="#334155", lw=1.2))
tracks = [t for t in BOOK["tracks"] if t["level"] == 2]
for k, t in enumerate(tracks):
    a = np.pi / 2 - k * 2 * np.pi / len(tracks)
    px, py = cx + 1.42 * np.cos(a), cy + 1.42 * np.sin(a)
    c = TRACK_COLORS[t["id"]]
    ax.plot([cx + 0.55 * np.cos(a), px - 0.33 * np.cos(a)], [cy + 0.55 * np.sin(a), py - 0.33 * np.sin(a)],
            color=c, lw=1.6, alpha=0.7)
    ax.add_patch(plt.Circle((px, py), 0.33, fc=BG, ec=c, lw=2.2, zorder=3))
    ax.text(px, py, t["id"], ha="center", va="center", fontsize=14, weight="bold", color=c, zorder=4)
ax.add_patch(plt.Circle((cx, cy), 0.55, fc="#2563eb", ec="white", lw=1.5, zorder=3))
ax.text(cx, cy + 0.1, "Core", ha="center", va="center", fontsize=12, weight="bold", color="white", zorder=4)
ax.text(cx, cy - 0.17, "1–25", ha="center", va="center", fontsize=9.5, color="#dbeafe", zorder=4)
ax.text(cx, cy + 2.6, "J · Frontier", ha="center", fontsize=10.5, color=TRACK_COLORS["J"], weight="bold")

IMAGES.mkdir(exist_ok=True)
fig.savefig(IMAGES / "social-preview.png", dpi=100, facecolor=BG, metadata={"Software": None})
plt.close(fig)
print("wrote images/social-preview.png")
