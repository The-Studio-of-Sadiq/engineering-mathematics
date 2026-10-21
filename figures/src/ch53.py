"""Figures for Chapter 53 — Reinforcement Learning Through Mathematics."""
import sys, pathlib; sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
from style import plt, np, PALETTE, INK, MUTED, new_figure, diagram, node, box, arrow, save

BLUE, RED, GREEN, PURPLE, ORANGE, TEAL, GREY = PALETTE

# ------------------------------------------------- figure 1: the motor-drive MDP
fig, ax = diagram(figsize=(8.6, 4.6))
ax.set_xlim(0, 11); ax.set_ylim(0.2, 6.2)
SC, SH = (2.0, 3.2), (6.6, 3.2)
node(ax, *SC, "Cool\n$C$", color=BLUE, r=0.75, fontsize=11)
node(ax, *SH, "Hot\n$H$", color=ORANGE, r=0.75, fontsize=11)
box(ax, 9.9, 3.2, "Tripped $T$\n(terminal)", w=1.8, h=1.0, color=RED, fontsize=9.5)


def action(x, y, label):
    ax.plot(x, y, "o", ms=9, color=INK, zorder=5)
    ax.text(x, y + 0.28, label, ha="center", va="bottom", fontsize=9.5, color=INK, zorder=6)


def link(p, q, text=None, color=MUTED, curve=0.0, sa=10, sb=30, tx=None, ty=None, lw=1.5):
    ax.annotate("", xy=q, xytext=p, zorder=2,
                arrowprops=dict(arrowstyle="-|>", color=color, lw=lw, shrinkA=sa, shrinkB=sb,
                                connectionstyle=f"arc3,rad={curve}"))
    if text:
        ax.text(tx, ty, text, ha="center", va="center", fontsize=8.8, color=color,
                bbox=dict(fc="white", ec="none", pad=0.6), zorder=6)


# Cool, normal (self-loop through action node above-left)
aCn = (1.0, 5.2)
link(SC, aCn, sa=30, sb=6, color=INK, lw=1.2)
action(*aCn, "normal")
link(aCn, SC, "1.0, $r=+1$", color=BLUE, curve=0.5, sa=6, sb=46, tx=0.45, ty=4.2)
# Cool, boost
aCb = (4.3, 5.2)
link(SC, aCb, sa=30, sb=6, color=INK, lw=1.2)
action(*aCb, "boost")
link(aCb, SC, "0.5, $+2$", color=BLUE, curve=0.25, sa=6, sb=46, tx=2.55, ty=5.0)
link(aCb, SH, "0.5, $+2$", color=ORANGE, curve=-0.15, sa=6, sb=46, tx=5.75, ty=4.55)
# Hot, normal
aHn = (4.3, 1.2)
link(SH, aHn, sa=30, sb=6, color=INK, lw=1.2)
action(*aHn, "")
ax.text(4.3, 0.65, "normal", ha="center", fontsize=9.5, color=INK)
link(aHn, SC, "0.5, $+1$", color=BLUE, curve=-0.2, sa=6, sb=46, tx=2.85, ty=1.75)
link(aHn, SH, "0.5, $+1$", color=ORANGE, curve=0.25, sa=6, sb=46, tx=5.9, ty=1.6)
# Hot, boost
aHb = (8.3, 5.2)
link(SH, aHb, sa=30, sb=6, color=INK, lw=1.2)
action(*aHb, "boost")
link(aHb, SH, "0.2, $+2$", color=ORANGE, curve=0.35, sa=6, sb=46, tx=7.0, ty=4.75)
link(aHb, (9.9, 3.2), "0.8, $-8$", color=RED, curve=-0.2, sa=6, sb=28, tx=9.75, ty=4.6)
ax.text(5.5, 0.25, r"black dots: actions $a$;  coloured arrows: $P(s'\mid s,a)$ and reward $r$;  $\gamma=0.9$",
        ha="center", fontsize=9, color=MUTED)
save(fig, "ch53-mdp")

# ------------------------------------------------- cliff world: Q-learning and SARSA
NR, NC = 4, 12
MOVES = [(-1, 0), (1, 0), (0, -1), (0, 1)]
START, GOAL = 36, 47


def env_step(s, a):
    r, c = divmod(s, NC)
    dr, dc = MOVES[a]
    r = min(max(r + dr, 0), NR - 1); c = min(max(c + dc, 0), NC - 1)
    if r == 3 and 1 <= c <= 10:
        return START, -100, False
    s2 = r * NC + c
    return s2, -1, s2 == GOAL


def greedy_return(Q, cap=60):
    s, G = START, 0
    for _ in range(cap):
        s, r, done = env_step(s, int(np.argmax(Q[s])))
        G += r
        if done:
            return G
    return -cap


def train(eps, alg, episodes, rng, alpha=0.5, eval_every=10):
    Q = np.zeros((48, 4))
    online, greedy = np.zeros(episodes), []

    def act(s):
        if rng.random() < eps:
            return int(rng.integers(4))
        q = Q[s]
        return int(rng.choice(np.flatnonzero(q == q.max())))

    for e in range(episodes):
        s = START; a = act(s); G = 0; done = False; n = 0
        while not done and n < 1000:
            s2, r, done = env_step(s, a); G += r; n += 1
            a2 = act(s2)
            nxt = 0.0 if done else (Q[s2].max() if alg == "q" else Q[s2, a2])
            Q[s, a] += alpha * (r + nxt - Q[s, a])
            s, a = s2, a2
        online[e] = G
        if e % eval_every == 0:
            greedy.append(greedy_return(Q))
    return online, np.array(greedy), Q


rng = np.random.default_rng(53)
EPIS, RUNS = 500, 20
curves = {}
for alg, eps in [("q", 0.01), ("q", 0.1), ("q", 0.3), ("s", 0.1)]:
    res = [train(eps, alg, EPIS, rng) for _ in range(RUNS)]
    curves[(alg, eps)] = (np.mean([r[0] for r in res], 0), np.mean([r[1] for r in res], 0), res[0][2])
    print(alg, eps, "last-100 online", curves[(alg, eps)][0][-100:].mean().round(1),
          "final greedy", curves[(alg, eps)][1][-5:].mean().round(1))


def path_of(Q):
    s, pts = START, [divmod(START, NC)]
    for _ in range(60):
        s, r, done = env_step(s, int(np.argmax(Q[s])))
        pts.append(divmod(s, NC))
        if done:
            break
    return np.array(pts)


# ------------------------------------------------- figure 2: the cliff and the two learned paths
fig, ax = new_figure(figsize=(8.4, 3.4))
ax.grid(False)
for r in range(NR):
    for c in range(NC):
        fc = "#fee2e2" if (r == 3 and 1 <= c <= 10) else "white"
        ax.add_patch(plt.Rectangle((c, NR - 1 - r), 1, 1, fc=fc, ec="#d1d5db", lw=0.8))
ax.text(6, 0.5, "the cliff: reward $-100$, back to start", ha="center", va="center", color=RED, fontsize=10)
ax.text(0.5, 0.5, "S", ha="center", va="center", fontsize=12, weight="bold", color=INK)
ax.text(11.5, 0.5, "G", ha="center", va="center", fontsize=12, weight="bold", color=INK)
for (alg, eps), col, lab, off in [(("q", 0.1), BLUE, "Q-learning (greedy policy)", -0.12),
                                  (("s", 0.1), GREEN, "SARSA (greedy policy)", 0.12)]:
    pts = path_of(curves[(alg, eps)][2])
    ax.plot(pts[:, 1] + 0.5 + off, NR - 1 - pts[:, 0] + 0.5 + off, color=col, lw=2.6, label=lab,
            marker="o", ms=3)
ax.set_xlim(0, NC); ax.set_ylim(0, NR); ax.set_aspect("equal")
ax.set_xticks([]); ax.set_yticks([])
for sp in ax.spines.values():
    sp.set_visible(False)
ax.legend(loc="upper center", bbox_to_anchor=(0.5, 1.2), ncol=2, fontsize=9.5)
save(fig, "ch53-cliff")

# ------------------------------------------------- figure 3: learning curves for different epsilon
fig, axs = new_figure(1, 2, figsize=(10, 4.0))
k = 10
smooth = lambda y: np.convolve(y, np.ones(k) / k, mode="valid")
styles = {("q", 0.01): (BLUE, "-", r"Q-learning, $\varepsilon=0.01$"),
          ("q", 0.1): (PURPLE, "-", r"Q-learning, $\varepsilon=0.1$"),
          ("q", 0.3): (ORANGE, "-", r"Q-learning, $\varepsilon=0.3$"),
          ("s", 0.1): (GREEN, "--", r"SARSA, $\varepsilon=0.1$")}
for key, (col, ls, lab) in styles.items():
    on, gr, _ = curves[key]
    axs[0].plot(np.arange(len(smooth(on))) + k, smooth(on), color=col, ls=ls, lw=1.8, label=lab)
    axs[1].plot(np.arange(len(gr)) * 10, gr, color=col, ls=ls, lw=1.8, label=lab)
axs[0].axhline(-13, color=MUTED, lw=1, ls=":")
axs[0].text(495, -10, "optimal path: $-13$", ha="right", fontsize=9, color=MUTED)
axs[0].set_ylim(-260, 0); axs[0].set_xlim(0, EPIS)
axs[0].set_xlabel("episode"); axs[0].set_ylabel("return while learning (exploring)")
axs[0].set_title("what the agent earns while it learns", fontsize=11)
axs[1].legend(fontsize=8.5, loc="lower right")
axs[1].axhline(-13, color=MUTED, lw=1, ls=":")
axs[1].set_ylim(-65, 0); axs[1].set_xlim(0, EPIS)
axs[1].set_xlabel("episode"); axs[1].set_ylabel("return of the greedy policy")
axs[1].set_title("what it would earn if it stopped exploring", fontsize=11)
save(fig, "ch53-learning")
