import sys, pathlib; sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
from style import plt, np, PALETTE, INK, MUTED, new_figure, diagram, node, box, arrow, save

# ---- Curry-Howard dictionary: two columns joined row by row
rows = [
    ("proposition $A$", "type $A$"),
    ("proof of $A$", "program (term) of type $A$"),
    (r"implication $A \to B$", r"function type $A \to B$"),
    (r"conjunction $A \wedge B$", r"pair type $A \times B$"),
    (r"disjunction $A \vee B$", r"tagged union $A + B$"),
    (r"truth $\top$", "unit type (one value)"),
    (r"falsehood $\bot$", "empty type (no values)"),
    (r"negation $\neg A$", r"function $A \to \bot$"),
    (r"$\forall x{:}X.\ P(x)$", r"dependent function $\Pi_{x:X} P(x)$"),
    (r"$\exists x{:}X.\ P(x)$", r"dependent pair $\Sigma_{x:X} P(x)$"),
    ("simplifying a proof", "running the program"),
]
fig, ax = diagram(figsize=(9.0, 6.6))
ax.set_xlim(0, 12); ax.set_ylim(0, 8.8)
ax.text(2.9, 8.35, "LOGIC", ha="center", fontsize=12, weight="bold", color=PALETTE[0])
ax.text(9.1, 8.35, "COMPUTATION", ha="center", fontsize=12, weight="bold", color=PALETTE[2])
for k, (l, r) in enumerate(rows):
    y = 7.7 - k * 0.7
    shade = "#fef3c7" if k in (1, 10) else "white"
    box(ax, 2.9, y, l, w=4.4, h=0.54, color=PALETTE[0], fc=shade, fontsize=9.5, lw=1.4)
    box(ax, 9.1, y, r, w=4.4, h=0.54, color=PALETTE[2], fc=shade, fontsize=9.5, lw=1.4)
    ax.annotate("", xy=(6.85, y), xytext=(5.15, y),
                arrowprops=dict(arrowstyle="<|-|>", color=MUTED, lw=1.3))
ax.text(6.0, 0.2, "each row is a precise theorem, not an analogy: the rules for one column are the rules for the other",
        ha="center", fontsize=8.5, color=MUTED)
save(fig, "ch66-curry-howard")

# ---- typing derivation tree for composition
fig, ax = diagram(figsize=(10.0, 4.6))
ax.set_xlim(0, 13); ax.set_ylim(1.3, 7.0)
TERM, TYPE = PALETTE[0], INK


def judgement(x, y, ctx, term, typ, fs=10):
    """Draw 'ctx |- term : type' with the colon at x: context grey, term blue, type black."""
    t = ax.text(x - 0.05, y, term + " :", ha="right", va="center", fontsize=fs, color=TERM)
    ax.annotate(ctx + r" $\vdash$ ", xy=(0, 0.5), xycoords=t, ha="right", va="center", fontsize=fs, color=MUTED)
    ax.text(x + 0.05, y, typ, ha="left", va="center", fontsize=fs, color=TYPE)


def rule(x0, x1, y, name):
    ax.plot([x0, x1], [y, y], color=INK, lw=1.3)
    ax.text(x1 + 0.12, y, name, ha="left", va="center", fontsize=8.5, color=PALETTE[3])


G = r"$\Gamma$"
# leaves: var rules have no premises, so a bare line sits above each leaf
judgement(1.9, 6.3, G, r"$f$", r"$A\to B$", fs=9.5); rule(0.9, 3.0, 6.62, "var")
judgement(5.3, 6.3, G, r"$x$", r"$A$", fs=9.5); rule(4.3, 5.8, 6.62, "var")
rule(0.6, 6.0, 5.95, "app")
judgement(3.5, 5.62, G, r"$f\,x$", r"$B$", fs=9.5)
judgement(9.3, 5.62, G, r"$g$", r"$B\to C$", fs=9.5); rule(8.3, 10.3, 5.95, "var")
rule(2.4, 10.4, 5.27, "app")
judgement(6.4, 4.94, G, r"$g\,(f\,x)$", r"$C$")
rule(3.6, 9.2, 4.59, r"abs (discharge $x{:}A$)")
judgement(6.6, 4.24, r"$f{:}A\to B,\ g{:}B\to C$", r"$\lambda x.\,g\,(f\,x)$", r"$A\to C$")
rule(2.6, 9.6, 3.89, r"abs (discharge $g$)")
judgement(6.6, 3.54, r"$f{:}A\to B$", r"$\lambda g.\lambda x.\,g\,(f\,x)$", r"$(B\to C)\to A\to C$")
rule(3.0, 10.6, 3.19, r"abs (discharge $f$)")
judgement(6.0, 2.84, "", r"$\lambda f.\lambda g.\lambda x.\,g\,(f\,x)$", r"$(A\to B)\to(B\to C)\to(A\to C)$")
ax.text(6.5, 2.05, r"$\Gamma = f{:}A\to B,\ g{:}B\to C,\ x{:}A$.   Grey: hypotheses.   Blue: the program.   Black: the proposition it proves.",
        ha="center", fontsize=9, color=MUTED)
ax.text(6.5, 1.6, "Erase the blue and what remains is a natural-deduction proof that implication is transitive.",
        ha="center", fontsize=9, color=MUTED)
save(fig, "ch66-typing-derivation")
