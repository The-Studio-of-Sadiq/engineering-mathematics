import sys, pathlib; sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
from style import plt, np, PALETTE, INK, MUTED, new_figure, diagram, node, box, arrow, save
from matplotlib.patches import FancyBboxPatch, Ellipse

# ---------------------------------------------------------------- 1. formal system vs model
fig, ax = diagram(figsize=(10.5, 5.6))
ax.set_xlim(0, 14); ax.set_ylim(0, 7.4)

# left panel: syntax
ax.add_patch(FancyBboxPatch((0.3, 0.4), 6.2, 6.4, boxstyle="round,pad=0.02,rounding_size=0.25",
                            fc="#f8fafc", ec=PALETTE[0], lw=1.6))
ax.text(3.4, 6.45, "SYNTAX: the formal system", ha="center", fontsize=11, weight="bold", color=PALETTE[0])
ax.text(3.4, 6.05, "symbols pushed around by rules; no meaning needed", ha="center", fontsize=8.8, color=MUTED)
box(ax, 1.9, 4.9, "axioms\n(finitely many\nschemes)", w=2.2, h=1.1, color=PALETTE[0])
box(ax, 4.9, 4.9, "inference rules\ne.g. modus ponens\nfrom $\\varphi$, $\\varphi\\to\\psi$ get $\\psi$", w=2.6, h=1.1,
    color=PALETTE[0], fontsize=8.6)
box(ax, 3.4, 2.6, "theorems  $T\\vdash\\varphi$\n(ends of finite proofs)", w=3.4, h=0.95, color=PALETTE[0], fc="#e5edff")
arrow(ax, 1.9, 4.35, 3.0, 3.1, shrink=3)
arrow(ax, 4.9, 4.35, 3.8, 3.1, shrink=3)
ax.text(3.4, 1.35, "a proof is a finite list of formulas,\neach an axiom or obtained by a rule",
        ha="center", fontsize=8.8, color=MUTED)

# right panel: semantics
ax.add_patch(FancyBboxPatch((7.5, 0.4), 6.2, 6.4, boxstyle="round,pad=0.02,rounding_size=0.25",
                            fc="#fefce8", ec=PALETTE[4], lw=1.6))
ax.text(10.6, 6.45, "SEMANTICS: models", ha="center", fontsize=11, weight="bold", color=PALETTE[4])
ax.text(10.6, 6.05, "structures in which sentences are true or false", ha="center", fontsize=8.8, color=MUTED)
for (x, y, lab) in [(9.0, 4.6, "$\\mathcal{M}_1$"), (10.7, 4.95, "$\\mathcal{M}_2$"), (12.4, 4.6, "$\\mathcal{M}_3$")]:
    ax.add_patch(Ellipse((x, y), 1.4, 1.0, fc="white", ec=PALETTE[4], lw=1.4))
    ax.text(x, y, lab, ha="center", va="center", fontsize=11, color=INK)
ax.text(10.85, 3.6, "every model of the axioms $T$", ha="left", fontsize=8.8, color=MUTED)
box(ax, 10.6, 2.6, "consequences  $T\\models\\varphi$\n(true in every model of $T$)", w=3.6, h=0.95,
    color=PALETTE[4], fc="#ffedd5")
arrow(ax, 10.7, 4.35, 10.6, 3.1, shrink=3, color=PALETTE[4])
ax.text(10.6, 1.35, "truth is defined by recursion\non the shape of the formula", ha="center", fontsize=8.8, color=MUTED)

# bridges
ax.annotate("", xy=(8.75, 2.75), xytext=(5.15, 2.75),
            arrowprops=dict(arrowstyle="-|>", color=PALETTE[2], lw=2.0))
ax.text(7.0, 2.95, "soundness", ha="center", fontsize=10, color=PALETTE[2], weight="bold")
ax.annotate("", xy=(5.15, 2.4), xytext=(8.75, 2.4),
            arrowprops=dict(arrowstyle="-|>", color=PALETTE[3], lw=2.0))
ax.text(7.0, 1.95, "completeness", ha="center", va="top", fontsize=10, color=PALETTE[3], weight="bold")
save(fig, "ch65-formal-system")

# ---------------------------------------------------------------- 2. map of the big theorems
fig, ax = diagram(figsize=(10.0, 5.6))
ax.set_xlim(0, 13); ax.set_ylim(0, 7.3)
box(ax, 3.0, 6.3, "Soundness\n$T\\vdash\\varphi \\Rightarrow T\\models\\varphi$", w=3.4, h=1.0, color=PALETTE[2])
box(ax, 10.0, 6.3, "Completeness (Gödel 1929)\n$T\\models\\varphi \\Rightarrow T\\vdash\\varphi$", w=4.0, h=1.0,
    color=PALETTE[3])
box(ax, 6.5, 4.2, "$\\vdash$ and $\\models$ coincide\nfor first-order logic", w=3.6, h=1.0, color=INK, fc="#f3f4f6")
arrow(ax, 3.0, 5.8, 5.6, 4.7, shrink=3)
arrow(ax, 10.0, 5.8, 7.4, 4.7, shrink=3)
box(ax, 2.2, 2.1, "Compactness\nevery finite part has a model\n$\\Rightarrow$ the whole has a model", w=3.8,
    h=1.2, color=PALETTE[0], fontsize=9)
box(ax, 6.5, 2.1, "Proofs can be\nenumerated:\nvalidity is semi-decidable", w=3.2, h=1.2, color=PALETTE[5], fontsize=9)
box(ax, 10.8, 2.1, "Löwenheim–Skolem\na countable theory with an\ninfinite model has a countable one",
    w=4.0, h=1.2, color=PALETTE[4], fontsize=9)
arrow(ax, 5.6, 3.7, 2.8, 2.7, shrink=3, text="proofs are finite", fontsize=8.5)
arrow(ax, 6.5, 3.7, 6.5, 2.7, shrink=3)
arrow(ax, 7.4, 3.7, 10.2, 2.7, shrink=3, text="Henkin's model", fontsize=8.5)
ax.text(2.2, 0.85, "non-standard models of arithmetic,\nno first-order way to say 'finite'", ha="center",
        fontsize=8.6, color=MUTED)
ax.text(10.8, 0.85, "Skolem's paradox:\na countable model of set theory", ha="center", fontsize=8.6, color=MUTED)
ax.text(6.5, 0.85, "but NOT decidable (Church, Turing 1936)", ha="center", fontsize=8.6, color=PALETTE[1])
save(fig, "ch65-logic-map")

# ---------------------------------------------------------------- 3. the MU puzzle search tree
rules = [
    lambda s: [s + "U"] if s.endswith("I") else [],
    lambda s: ["M" + s[1:] * 2],
    lambda s: [s[:i] + "U" + s[i + 3:] for i in range(1, len(s) - 2) if s[i:i + 3] == "III"],
    lambda s: [s[:i] + s[i + 2:] for i in range(1, len(s) - 1) if s[i:i + 2] == "UU"],
]
levels = [["MI"]]
seen = {"MI"}
parent = {}
for depth in range(3):
    nxt = []
    for s in levels[-1]:
        for r, rule in enumerate(rules):
            for t in rule(s):
                if t not in seen and len(t) <= 9:
                    seen.add(t); nxt.append(t); parent[t] = (s, r + 1)
    levels.append(nxt)

fig, ax = new_figure(figsize=(10.5, 4.8))
ax.axis("off")
pos = {}
for d, lev in enumerate(levels):
    n = len(lev)
    for k, s in enumerate(lev):
        x = (k + 0.5) / n
        pos[s] = (x, 1 - d * 0.3)
for t, (s, r) in parent.items():
    (x0, y0), (x1, y1) = pos[s], pos[t]
    ax.plot([x0, x1], [y0 - 0.03, y1 + 0.03], color=PALETTE[r - 1], lw=1.4, zorder=1)
    ax.text((x0 + x1) / 2, (y0 + y1) / 2, f"R{r}", fontsize=7.5, color=PALETTE[r - 1], ha="center",
            va="center", bbox=dict(fc="white", ec="none", pad=0.6), zorder=2)
for s, (x, y) in pos.items():
    i_count = s.count("I")
    ax.text(x, y, s + "\n" + f"#I = {i_count}", ha="center", va="center", fontsize=9 if len(s) < 8 else 8,
            color=INK, zorder=3, linespacing=1.3,
            bbox=dict(boxstyle="round,pad=0.3", fc="#eef2ff", ec=PALETTE[0], lw=1.0))
ax.set_xlim(-0.02, 1.02); ax.set_ylim(-0.08, 1.08)
ax.set_title("Theorems of the MIU system, three rule applications deep (strings up to 9 letters)")
ax.text(0.5, -0.06, "the number of I's is never a multiple of 3, so MU (zero I's) is never reached",
        ha="center", fontsize=9.5, color=PALETTE[1])
save(fig, "ch65-mu-tree")
