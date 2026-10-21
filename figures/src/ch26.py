import sys, pathlib; sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
from style import plt, np, PALETTE, INK, MUTED, new_figure, diagram, node, box, arrow, save

# ---------------------------------------------------------------- growth-rate crossover
costs = [
    (r"$2^n$", lambda n: 2.0 ** n, PALETTE[1]),
    (r"$50\,n^2$", lambda n: 50 * n ** 2, PALETTE[3]),
    (r"$400\,n\log_2 n$", lambda n: 400 * n * np.log2(n), PALETTE[0]),
    (r"$5000\,n$", lambda n: 5000 * n, PALETTE[2]),
]

fig, (axL, axR) = new_figure(1, 2, figsize=(9.6, 4.2))

n = np.linspace(1, 24, 400)
for lab, f, c in costs:
    axL.plot(n, f(n), color=c, label=lab)
axL.set_xlim(1, 24)
axL.set_ylim(0, 1.6e5)
axL.set_xlabel("input size $n$")
axL.set_ylabel("operations")
axL.set_title("Small inputs: the exponential looks cheap")
axL.ticklabel_format(axis="y", style="sci", scilimits=(0, 0))
# shade region where 2^n is the cheapest
mask = np.array([costs[0][1](v) <= min(f(v) for _, f, _ in costs[1:]) for v in n])
axL.fill_between(n, 0, 1.6e5, where=mask, color=PALETTE[1], alpha=0.07, lw=0)
axL.text(5.5, 1.35e5, "$2^n$ is the\nfastest here", color=PALETTE[1], ha="center", fontsize=9.5)
axL.legend(loc="center left", fontsize=9)

n2 = np.linspace(1, 120, 600)
for lab, f, c in costs:
    axR.semilogy(n2, f(n2), color=c, label=lab)
axR.set_xlim(1, 120)
axR.set_ylim(1e2, 1e12)
axR.set_xlabel("input size $n$")
axR.set_title("Larger inputs: the order of growth wins")
for x0, lab in [(13.06, "A"), (43.6, "B"), (100, "C")]:
    axR.axvline(x0, color=MUTED, lw=1, ls=":")
    axR.text(x0 + 1, 2e2, lab, color=MUTED, fontsize=10)
axR.text(60, 3e3, "A: $2^n$ overtakes $50n^2$\nB: $50n^2$ overtakes $400n\\log_2 n$\nC: $50n^2$ overtakes $5000n$",
         fontsize=8.5, color=INK)
save(fig, "ch26-growth-crossover")

# ---------------------------------------------------------------- merge-sort recursion tree
fig, ax = diagram(figsize=(9.0, 4.6))
ax.set_xlim(-0.6, 12.2)
ax.set_ylim(0.1, 6.2)
levels = 4  # sizes 8, 4, 2, 1
ys = [5.5, 4.1, 2.7, 1.3]
width = 10.0
for k in range(levels):
    cnt = 2 ** k
    size = 8 // cnt
    w = width / cnt
    for i in range(cnt):
        x = w * (i + 0.5)
        box(ax, x, ys[k], f"{size}", w=min(0.9 * w, 2.0), h=0.6, color=PALETTE[0],
            fc="#eef3ff" if k % 2 == 0 else "white", fontsize=10)
        if k + 1 < levels:
            wc = w / 2
            for c in range(2):
                xc = wc * (2 * i + c + 0.5)
                ax.plot([x, xc], [ys[k] - 0.3, ys[k + 1] + 0.3], color=MUTED, lw=1.1, zorder=1)
    ax.text(11.2, ys[k], f"{cnt} × {size}\n= 8", ha="center", va="center", fontsize=9.5, color=PALETTE[2])
ax.text(11.2, 6.1, "items per level", ha="center", fontsize=9.5, color=PALETTE[2], weight="bold")
ax.annotate("", xy=(-0.35, 1.0), xytext=(-0.35, 5.8),
            arrowprops=dict(arrowstyle="<->", color=PALETTE[1], lw=1.5))
ax.text(-0.25, 3.4, r"$\log_2 n + 1$ levels", rotation=90, ha="left", va="center", color=PALETTE[1], fontsize=9.5)
ax.text(5.0, 0.35, r"each of the $\log_2 n$ levels above the leaves merges $n$ items in total, so $T(n) = n\log_2 n$",
        ha="center", fontsize=9.5, color=INK)
save(fig, "ch26-recursion-tree")
