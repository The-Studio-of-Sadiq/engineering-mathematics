"""Chapter 46 figures: a minimax tree with alpha-beta pruning, and regret curves for online learning."""
import sys, pathlib; sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
from style import plt, np, PALETTE, INK, MUTED, new_figure, diagram, node, box, arrow, save

# ---------------------------------------------------------------- alpha-beta tree
# root MAX -> 3 MIN -> 2 MAX each -> 2 leaves each
leaves = [3, 5, 6, 9, 1, 2, 0, 4, 7, 4, 8, 1]
tree = [[[leaves[4 * i + 2 * j + k] for k in range(2)] for j in range(2)] for i in range(3)]
visited, values, pruned_edges = set(), {}, set()

def ab(t, path, maxing, alpha, beta):
    if not isinstance(t, list):
        visited.add(path); values[path] = t; return t
    best = -np.inf if maxing else np.inf
    for i, c in enumerate(t):
        v = ab(c, path + (i,), not maxing, alpha, beta)
        if maxing:
            best = max(best, v); alpha = max(alpha, best)
        else:
            best = min(best, v); beta = min(beta, best)
        if alpha >= beta:
            for k in range(i + 1, len(t)):
                pruned_edges.add(path + (k,))
            break
    visited.add(path); values[path] = best
    return best

root = ab(tree, (), True, -np.inf, np.inf)

fig, ax = plt.subplots(figsize=(10.0, 4.8), constrained_layout=True)
ax.set_xlim(-2.6, 11.6); ax.set_ylim(-0.9, 4.6); ax.axis("off")
pos = {}
def xs(path):
    # leaf index range -> x
    lo = 0; span = 12
    for d, i in enumerate(path):
        span //= [3, 2, 2][d]; lo += i * span
    return lo + (span - 1) / 2
for d in range(4):
    for path in ([()] if d == 0 else [p for p in __import__("itertools").product(*[range(n) for n in [3, 2, 2][:d]])]):
        pos[path] = (xs(path), 3.9 - 1.3 * d)

def is_pruned(path):
    return any(path[:k] in pruned_edges for k in range(1, len(path) + 1))

for path, (x, y) in pos.items():
    if path:
        px, py = pos[path[:-1]]
        cut = path in pruned_edges
        gone = is_pruned(path)
        ax.plot([px, x], [py - 0.3, y + 0.3], color=PALETTE[1] if gone else INK,
                ls="--" if gone else "-", lw=1.2 if gone else 1.6, zorder=1)
        if cut:
            mx, my = (px + x) / 2, (py + y) / 2
            ax.plot([mx - 0.18, mx + 0.18], [my - 0.12, my + 0.12], color=PALETTE[1], lw=2.6, zorder=2)
            ax.plot([mx - 0.18, mx + 0.18], [my + 0.12, my - 0.12], color=PALETTE[1], lw=2.6, zorder=2)
for path, (x, y) in pos.items():
    d = len(path); gone = is_pruned(path)
    col = MUTED if gone else [PALETTE[0], PALETTE[1], PALETTE[0], INK][d]
    txt = str(values[path]) if (path in values and not gone) else ("?" if d < 3 else str(tree[path[0]][path[1]][path[2]]))
    if d < 3 and not gone and any(path + (k,) in pruned_edges for k in range(3)):
        txt = ("≥" if d % 2 == 0 else "≤") + txt
    if d == 3:
        ax.add_patch(plt.Rectangle((x - 0.3, y - 0.3), 0.6, 0.6, fc="#f3f4f6" if gone else "white",
                                   ec=col, lw=1.6, zorder=3))
        ax.text(x, y, txt, ha="center", va="center", color=MUTED if gone else INK, fontsize=11, zorder=4)
    else:
        mark = "^" if d % 2 == 0 else "v"
        ax.scatter([x], [y], marker=mark, s=900, c="#f3f4f6" if gone else "white", edgecolors=col,
                   linewidths=2, zorder=3)
        ax.text(x, y - (0.03 if d % 2 == 0 else -0.05), txt, ha="center", va="center",
                fontsize=10.5, color=MUTED if gone else INK, zorder=4)
for d, lab in enumerate(["MAX (you)", "MIN (opponent)", "MAX (you)", "leaves (payoff\nto MAX)"]):
    ax.text(-2.55, 3.9 - 1.3 * d, lab, ha="left", va="center", fontsize=9, color=MUTED)
ax.text(11.5, -0.75, "red X = branch cut by alpha-beta; dashed = never examined", ha="right",
        fontsize=9, color=PALETTE[1])
ax.text(pos[()][0] + 0.55, pos[()][1] + 0.25, f"game value = {root}", fontsize=10, color=PALETTE[0])
save(fig, "ch46-alpha-beta")
print("visited leaves:", sum(1 for p in visited if len(p) == 3), "of", len(leaves), "value", root)

# ---------------------------------------------------------------- regret curves
T = 2000
rng = np.random.default_rng(46)
# adversarial losses for 2 experts that punish Follow-the-Leader
adv = np.zeros((T, 2)); adv[0] = [0.5, 0.0]
adv[1::2] = [0.0, 1.0]; adv[2::2] = [1.0, 0.0]
# stochastic losses for 10 experts, best expert slightly better
n = 10
means = np.full(n, 0.5); means[3] = 0.4
sto = (rng.random((T, n)) < means).astype(float)

def run(L, rule, eta=None):
    T, n = L.shape; cum = np.zeros(n); alg = 0.0; reg = np.zeros(T)
    for t in range(T):
        if rule == "ftl":
            p = (cum == cum.min()).astype(float); p /= p.sum()
        else:
            w = np.exp(-eta * (cum - cum.min())); p = w / w.sum()
        alg += p @ L[t]; cum += L[t]
        reg[t] = alg - cum.min()
    return reg

fig, (a1, a2) = new_figure(1, 2, figsize=(10.0, 3.9))
tt = np.arange(1, T + 1)
for ax, L, title in [(a1, adv, "Adversarial losses, 2 experts"), (a2, sto, "Random losses, 10 experts")]:
    m = L.shape[1]
    eta = np.sqrt(2 * np.log(m) / T)
    ax.plot(tt, run(L, "ftl"), color=PALETTE[1], label="Follow the Leader")
    ax.plot(tt, run(L, "hedge", eta), color=PALETTE[0], label=r"Hedge, $\eta=\sqrt{2\ln n/T}$")
    ax.plot(tt, np.sqrt(2 * tt * np.log(m)), color=INK, ls="--", lw=1.2, label=r"bound $\sqrt{2T\ln n}$")
    ax.set_title(title); ax.set_xlabel("round $T$"); ax.set_ylabel("regret $R_T$")
a1.legend(loc="upper left")
a1.set_ylim(-20, 1050)
save(fig, "ch46-regret")
