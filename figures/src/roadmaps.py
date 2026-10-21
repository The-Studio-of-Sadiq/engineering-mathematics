"""Universe map (cover / table of contents) and one roadmap per track, drawn from tools/book.json."""

import json
import sys
import pathlib
import textwrap

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
from style import plt, np, PALETTE, INK, MUTED, ROOT, save  # noqa: E402

BOOK = json.loads((ROOT / "tools" / "book.json").read_text(encoding="utf-8"))
CORE_TITLES = {}
for p in sorted(ROOT.glob("part-*/[0-9][0-9]-*.md")):
    CORE_TITLES[int(p.name[:2])] = p.read_text(encoding="utf-8").splitlines()[0].split("—", 1)[-1].strip()

TRACK_COLORS = dict(zip("ABCDEFGHIJ", PALETTE + ["#be185d", "#0d9488", "#a16207"]))


def short(title, width):
    head = title.split(":")[0]
    return "\n".join(textwrap.wrap(head, width))


def universe_map():
    fig, ax = plt.subplots(figsize=(8, 8), constrained_layout=True)
    ax.set_xlim(-5.2, 5.2)
    ax.set_ylim(-5.2, 5.2)
    ax.set_aspect("equal")
    ax.axis("off")
    # Level 3 ring.
    ax.add_patch(plt.Circle((0, 0), 4.95, fc="#fdf2f8", ec=TRACK_COLORS["J"], lw=2, ls="--"))
    ax.text(0, 4.62, "Level 3 · Research Frontier (Track J)", ha="center", fontsize=10.5, color=TRACK_COLORS["J"],
            weight="bold")
    j = next(t for t in BOOK["tracks"] if t["id"] == "J")
    frontier = [short(c["title"], 14) for c in j["chapters"]]
    for k, label in enumerate(frontier):
        a = np.pi / 2 + 0.55 + k * (2 * np.pi - 1.1) / (len(frontier) - 1)
        ax.text(4.35 * np.cos(a), 4.35 * np.sin(a), label.replace("\n", " ") if len(label) < 16 else label,
                ha="center", va="center", fontsize=6.6, color=MUTED, rotation=0)
    # Level 2 tracks.
    ax.add_patch(plt.Circle((0, 0), 3.75, fc="#f8fafc", ec="#cbd5e1", lw=1.2))
    tracks = [t for t in BOOK["tracks"] if t["level"] == 2]
    for k, t in enumerate(tracks):
        a = np.pi / 2 - k * 2 * np.pi / len(tracks)
        x, y = 2.65 * np.cos(a), 2.65 * np.sin(a)
        c = TRACK_COLORS[t["id"]]
        ax.plot([0.95 * np.cos(a), x - 0.72 * np.cos(a)], [0.95 * np.sin(a), y - 0.72 * np.sin(a)], color=c, lw=1.4, alpha=0.6)
        ax.add_patch(plt.Circle((x, y), 0.72, fc="white", ec=c, lw=2.2, zorder=3))
        ax.text(x, y + 0.36, t["id"], ha="center", va="center", fontsize=13, weight="bold", color=c, zorder=4)
        ax.text(x, y - 0.12, short(t["title"].replace("The Mathematics Behind ", "").replace("The Abstract Mathematics Behind Everything", "Abstract Mathematics").replace("Mathematics for ", ""), 15),
                ha="center", va="center", fontsize=6.3, color=INK, zorder=4, linespacing=1.05)
    ax.text(0, -3.55, "Level 2 · Specializations", ha="center", fontsize=9.5, color=MUTED, weight="bold")
    # Level 1 core.
    ax.add_patch(plt.Circle((0, 0), 0.95, fc=PALETTE[0], ec="white", lw=2, zorder=3))
    ax.text(0, 0.18, "Level 1", ha="center", fontsize=9, color="white", zorder=4)
    ax.text(0, -0.18, "Core · 25 ch.", ha="center", fontsize=9.5, color="white", weight="bold", zorder=4)
    save(fig, "universe-map")


def roadmap(t):
    chapters = t["chapters"]
    n = len(chapters)
    c = TRACK_COLORS[t["id"]]
    h = max(4.0, 0.62 * n + 1.4)
    fig, ax = plt.subplots(figsize=(8, h), constrained_layout=True)
    ax.set_xlim(0, 10)
    ax.set_ylim(-0.6, n + 0.2)
    ax.axis("off")
    ys = [n - 0.5 - i for i in range(n)]
    # Core feeder line.
    core = t["core"]
    cys = np.linspace(n - 0.5, max(n - 0.5 - (len(core) - 1) * 0.8, 0), len(core)) if len(core) > 1 else [n - 0.5]
    ax.plot([1.3, 1.3], [min(cys), max(cys)], color="#cbd5e1", lw=6, solid_capstyle="round", zorder=1)
    for y, num in zip(cys, core):
        ax.add_patch(plt.Circle((1.3, y), 0.2, fc="white", ec="#94a3b8", lw=2, zorder=3))
        ax.text(1.0, y, f"Ch. {num}", ha="right", va="center", fontsize=8.5, color=MUTED)
        ax.text(1.6, y, textwrap.shorten(CORE_TITLES.get(num, ""), 36, placeholder="…"), ha="left", va="center",
                fontsize=7.2, color=MUTED)
        ax.annotate("", xy=(4.15, ys[0]), xytext=(1.5, y), arrowprops=dict(arrowstyle="-", color="#e2e8f0", lw=1,
                    connectionstyle="arc3,rad=-0.15"), zorder=0)
    ax.text(1.3, n + 0.05, "core prerequisites", ha="center", fontsize=8, color=MUTED, style="italic")
    # Track line.
    ax.plot([4.4, 4.4], [ys[-1], ys[0]], color=c, lw=7, solid_capstyle="round", zorder=1)
    for y, ch in zip(ys, chapters):
        ax.add_patch(plt.Circle((4.4, y), 0.26, fc="white", ec=c, lw=2.6, zorder=3))
        ax.text(4.4, y, str(ch["n"]), ha="center", va="center", fontsize=7.5, weight="bold", color=INK, zorder=4)
        ax.text(4.85, y, ch["title"], ha="left", va="center", fontsize=9, color=INK)
    ax.text(4.4, n + 0.05, f"Track {t['id']}", ha="center", fontsize=9.5, weight="bold", color=c)
    save(fig, f"roadmap-track-{t['id'].lower()}")


universe_map()
for t in BOOK["tracks"]:
    roadmap(t)
