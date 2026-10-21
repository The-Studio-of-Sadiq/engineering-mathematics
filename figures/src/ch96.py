"""Figures for Chapter 96 — Mathematical Economics and Social Choice."""
import sys, pathlib; sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
from style import plt, np, PALETTE, INK, MUTED, new_figure, diagram, node, box, arrow, save

BLUE, RED, GREEN, PURPLE, ORANGE, TEAL, GREY = PALETTE
a, b = 0.6, 0.4          # Ann: x^0.6 y^0.4 ; Bob: x^0.4 y^0.6
W = 10.0                 # 10 apples, 10 bread in total
wx, wy = 8.0, 2.0        # Ann's endowment (Bob has the rest)
ex, ey = 6.0, 4.0        # Walrasian equilibrium allocation for Ann, price ratio 1

uA = lambda x, y: x ** a * y ** (1 - a)
uB = lambda x, y: (W - x) ** b * (W - y) ** (1 - b)
contract = lambda x: 8 * x / (18 - x)          # equal marginal rates of substitution


def ann_curve(U, xs):
    return (U / xs ** a) ** (1 / (1 - a))


def bob_curve(U, xs):
    return W - (U / (W - xs) ** b) ** (1 / (1 - b))


# ---------------------------------------------------------------- figure 1: Edgeworth box
fig, ax = plt.subplots(figsize=(6.4, 6.0), constrained_layout=True)
ax.set_xlim(0, W); ax.set_ylim(0, W); ax.set_aspect("equal"); ax.grid(False)
xs = np.linspace(0.05, 9.95, 800)
# lens of mutually better trades
A0, B0 = uA(wx, wy), uB(wx, wy)
ya, yb = ann_curve(A0, xs), bob_curve(B0, xs)
ax.fill_between(xs, ya, yb, where=(yb > ya), color=GREEN, alpha=0.13, lw=0)
for U, ls in [(A0, "-"), (uA(ex, ey), "-"), (uA(3, 3), ":")]:
    yy = ann_curve(U, xs); m = yy <= W
    ax.plot(xs[m], yy[m], color=BLUE, ls=ls, lw=1.8 if ls == "-" else 1.2)
for U, ls in [(B0, "-"), (uB(ex, ey), "-"), (uB(7.5, 7.5), ":")]:
    yy = bob_curve(U, xs); m = yy >= 0
    ax.plot(xs[m], yy[m], color=RED, ls=ls, lw=1.8 if ls == "-" else 1.2)
ax.plot(xs, contract(xs), color=PURPLE, lw=2.4)
ax.plot([0, W], [W, 0], color=INK, ls="--", lw=1.2)
ax.plot(wx, wy, "o", color=INK, ms=8, zorder=5)
ax.plot(ex, ey, "*", color=ORANGE, ms=17, mec=INK, zorder=6)
ax.annotate("endowment $\\omega$\nAnn (8, 2), Bob (2, 8)", xy=(wx, wy), xytext=(6.2, 0.55), fontsize=9.5, color=INK,
            arrowprops=dict(arrowstyle="-", color=MUTED))
ax.annotate("equilibrium: Ann (6, 4)\nprice ratio $p_x/p_y=1$", xy=(ex, ey), xytext=(6.6, 5.6), fontsize=9.5,
            color=INK, arrowprops=dict(arrowstyle="-", color=MUTED))
ax.text(1.0, 2.05, "contract curve\n(Pareto optima)", color=PURPLE, fontsize=9.5, rotation=0)
ax.text(2.0, 7.6, "budget line\nthrough $\\omega$", color=INK, fontsize=9)
ax.text(5.2, 2.9, "both\nbetter off", color=GREEN, fontsize=9, ha="center")
ax.text(0.15, -0.75, "$O_{Ann}$", color=BLUE, fontsize=11)
ax.text(9.35, 10.3, "$O_{Bob}$", color=RED, fontsize=11)
ax.set_xlabel("apples held by Ann  (Bob holds 10 − this)")
ax.set_ylabel("bread held by Ann  (Bob holds 10 − this)")
ax.plot([], [], color=BLUE, label="Ann's indifference curves")
ax.plot([], [], color=RED, label="Bob's indifference curves")
ax.legend(loc="upper left", fontsize=8.5)
save(fig, "ch96-edgeworth-box")

# ---------------------------------------------------------------- figure 2: utility possibility frontier
fig, ax = new_figure(figsize=(6.2, 4.4))
xc = np.linspace(0, W, 600); yc = contract(xc)
UA, UB = uA(xc, yc), uB(xc, yc)
ax.plot(UA, UB, color=PURPLE, lw=2.4, label="frontier: utilities along the contract curve")
ax.plot(np.sqrt(W * UA), np.sqrt(W * UB), color=TEAL, lw=1.8, ls="--",
        label=r"same frontier after relabelling utilities $u\mapsto\sqrt{10u}$")
inlens = (UA >= A0) & (UB >= B0)
ax.plot(UA[inlens], UB[inlens], color=GREEN, lw=6, alpha=0.35, label="reachable by voluntary trade")
rng = np.random.default_rng(3)
px, py = rng.uniform(0, W, 900), rng.uniform(0, W, 900)
ax.scatter(uA(px, py), uB(px, py), s=4, color=GREY, alpha=0.35, label="random allocations")
ax.plot(A0, B0, "o", color=INK, ms=8, zorder=5)
ax.annotate("endowment", xy=(A0, B0), xytext=(A0 - 2.6, B0 - 1.3), fontsize=9.5,
            arrowprops=dict(arrowstyle="-", color=MUTED))
ax.plot(uA(ex, ey), uB(ex, ey), "*", color=ORANGE, ms=16, mec=INK, zorder=6)
ax.annotate("market equilibrium", xy=(uA(ex, ey), uB(ex, ey)), xytext=(uA(ex, ey) + 0.4, uB(ex, ey) + 1.3),
            fontsize=9.5, arrowprops=dict(arrowstyle="-", color=MUTED))
ax.set_xlabel("Ann's utility $u_A$"); ax.set_ylabel("Bob's utility $u_B$")
ax.set_xlim(0, 10.5); ax.set_ylim(0, 10.5)
ax.legend(loc="lower left", fontsize=8.5, frameon=True, facecolor="white", edgecolor="white", framealpha=0.95)
ax.set_title("The frontier: nobody can gain without someone losing", fontsize=10)
save(fig, "ch96-utility-frontier")
