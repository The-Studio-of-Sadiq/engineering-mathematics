import sys, pathlib; sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
from style import plt, np, PALETTE, INK, MUTED, new_figure, diagram, node, box, arrow, save


def view(az=-55.0, el=22.0):
    """Orthographic 'drawing camera': returns a function mapping 3-D points to 2-D paper coordinates."""
    a, e = np.radians(az), np.radians(el)
    right = np.array([np.cos(a), np.sin(a), 0.0])
    fwd = np.array([-np.sin(a) * np.cos(e), np.cos(a) * np.cos(e), np.sin(e)])
    up = np.cross(right, fwd)
    up = -up if up[2] < 0 else up
    return lambda P: np.c_[np.atleast_2d(P) @ right, np.atleast_2d(P) @ up]


def rodrigues(v, n, th):
    n = n / np.linalg.norm(n)
    return v * np.cos(th) + np.cross(n, v) * np.sin(th) + n * (n @ v) * (1 - np.cos(th))


def arrow3(ax, Pr, p0, p1, color, lw=2.0, ls="-", ms=14):
    q = Pr(np.array([p0, p1]))
    ax.annotate("", xy=q[1], xytext=q[0], annotation_clip=False,
                arrowprops=dict(arrowstyle="-|>", color=color, lw=lw, ls=ls, mutation_scale=ms))
    ax.plot(q[:, 0], q[:, 1], alpha=0)
    return q[1]


# ---------------------------------------------------------------- figure 1: axis-angle rotation
Pr = view(-60, 18)
fig, ax = new_figure(figsize=(6.6, 5.0))
ax.set_aspect("equal"); ax.axis("off")
n = np.array([0.25, 0.15, 1.0]); n = n / np.linalg.norm(n)
v = np.array([1.25, -0.25, 0.55])
th = np.radians(110)
vpar = n * (n @ v)
vperp = v - vpar
# cone circle traced by the tip of v
ts = np.linspace(0, 2 * np.pi, 200)
circ = np.array([rodrigues(v, n, t) for t in ts])
q = Pr(circ); ax.plot(q[:, 0], q[:, 1], color=MUTED, lw=1, ls="--")
arc = np.array([rodrigues(v, n, t) for t in np.linspace(0, th, 60)])
q = Pr(arc); ax.plot(q[:, 0], q[:, 1], color=PALETTE[4], lw=2.4)
ax.annotate("", xy=q[-1], xytext=q[-4], arrowprops=dict(arrowstyle="-|>", color=PALETTE[4], lw=2.2, mutation_scale=16))
# axis line through origin
q = Pr(np.array([-0.5 * n, 1.75 * n])); ax.plot(q[:, 0], q[:, 1], color=INK, lw=1, ls=":")
tip = arrow3(ax, Pr, np.zeros(3), 1.6 * n, INK, lw=2.4)
ax.text(tip[0] + 0.04, tip[1], "axis $\\hat{\\mathbf{n}}$", fontsize=11, color=INK, va="center")
tv = arrow3(ax, Pr, np.zeros(3), v, PALETTE[0], lw=2.4)
ax.text(tv[0] + 0.05, tv[1] - 0.05, "$\\mathbf{v}$", fontsize=13, color=PALETTE[0], va="top")
Rv = rodrigues(v, n, th)
tr = arrow3(ax, Pr, np.zeros(3), Rv, PALETTE[1], lw=2.4)
ax.text(tr[0] - 0.06, tr[1] + 0.03, "$\\mathbf{R}\\mathbf{v}$", fontsize=13, color=PALETTE[1], ha="right", va="bottom")
# decomposition: parallel part and perpendicular parts at the circle centre
cp = arrow3(ax, Pr, vpar, v, PALETTE[2], lw=1.6, ms=11)
arrow3(ax, Pr, vpar, Rv, PALETTE[2], lw=1.6, ms=11)
qc = Pr(vpar)[0]
ax.plot(*qc, "o", color=PALETTE[2], ms=4)
ax.text(qc[0] - 0.06, qc[1] + 0.02, "$\\mathbf{v}_\\parallel$", color=PALETTE[2], fontsize=11, ha="right", va="bottom")
mid = Pr(0.5 * (vpar + v))[0]
ax.text(mid[0], mid[1] - 0.06, "$\\mathbf{v}_\\perp$", color=PALETTE[2], fontsize=11, ha="center", va="top")
am = Pr(rodrigues(v, n, th / 2))[0]
ax.text(am[0], am[1] - 0.07, "angle $\\theta$", color=PALETTE[4], fontsize=11, ha="center", va="top")
ax.plot(*Pr(np.zeros(3))[0], "o", color=INK, ms=4)
ax.set_title("Rotation by angle $\\theta$ about the axis $\\hat{\\mathbf{n}}$", fontsize=11.5)
save(fig, "ch70-axis-rotation")


# ---------------------------------------------------------------- figure 2: gimbal lock
def ring(ax, Pr, axis, radius, color, label, lab_dir, lw=3.0):
    axis = axis / np.linalg.norm(axis)
    a = np.cross(axis, [0.3, 0.5, 0.8]); a /= np.linalg.norm(a)
    b = np.cross(axis, a)
    t = np.linspace(0, 2 * np.pi, 240)
    P = radius * (np.outer(np.cos(t), a) + np.outer(np.sin(t), b))
    q = Pr(P); ax.plot(q[:, 0], q[:, 1], color=color, lw=lw, alpha=0.85)
    tip = arrow3(ax, Pr, -1.45 * radius * axis, 1.45 * radius * axis, color, lw=1.6, ms=12)
    ax.text(tip[0] + lab_dir[0], tip[1] + lab_dir[1], label, color=color, fontsize=10,
            ha="center", va="center")


Pr = view(-35, 20)
fig, axes = plt.subplots(1, 2, figsize=(9.6, 4.6), constrained_layout=True)
z, y, x = np.eye(3)[2], np.eye(3)[1], np.eye(3)[0]
cases = [("Pitch = 20°: three independent axes", np.radians(20)),
         ("Pitch = 90°: roll axis lands on yaw axis", np.radians(90))]
for ax, (title, p) in zip(axes, cases):
    ax.set_aspect("equal"); ax.axis("off")
    # yaw about world z; pitch about (yawed) y; roll about the pitched body x
    roll_axis = np.array([np.cos(p), 0.0, np.sin(p)])
    ring(ax, Pr, z, 1.25, PALETTE[0], "yaw", (0.0, 0.12))
    ring(ax, Pr, y, 1.0, PALETTE[2], "pitch", (-0.08, -0.16))
    ring(ax, Pr, roll_axis, 0.75, PALETTE[1], "roll", (0.18, 0.0) if p < 1.5 else (0.25, -0.1))
    ax.set_title(title, fontsize=11)
axes[1].text(0.5, -0.02, "yaw and roll now turn about the same line:\none degree of freedom has vanished",
             transform=axes[1].transAxes, ha="center", va="top", fontsize=9.5, color=INK)
axes[0].text(0.5, -0.02, "each ring adds a new direction of turning",
             transform=axes[0].transAxes, ha="center", va="top", fontsize=9.5, color=INK)
save(fig, "ch70-gimbal-lock")

# ---------------------------------------------------------------- figure 3: slerp vs lerp
Om = np.radians(150)
q0, q1 = np.array([1.0, 0.0]), np.array([np.cos(Om), np.sin(Om)])
ts = np.linspace(0, 1, 9)
fig, (ax, ax2) = plt.subplots(1, 2, figsize=(10.0, 4.4), constrained_layout=True,
                              gridspec_kw=dict(width_ratios=[1.15, 1]))
ax.set_aspect("equal"); ax.grid(False)
for s in ("left", "bottom"):
    ax.spines[s].set_visible(False)
ax.set_xticks([]); ax.set_yticks([])
t = np.linspace(0, 2 * np.pi, 300)
ax.plot(np.cos(t), np.sin(t), color="#d1d5db", lw=1.2)
ax.plot([q0[0], q1[0]], [q0[1], q1[1]], color=MUTED, lw=1.2, ls="--")
lerp = np.outer(1 - ts, q0) + np.outer(ts, q1)
nl = lerp / np.linalg.norm(lerp, axis=1, keepdims=True)
sl = np.array([np.cos(tt * Om) * np.array([1, 0]) + np.sin(tt * Om) * np.array([0, 1]) for tt in ts])
for a, b in zip(lerp, nl):
    ax.plot([0, b[0]], [0, b[1]], color=PALETTE[1], lw=0.6, alpha=0.5)
ax.plot(lerp[:, 0], lerp[:, 1], "o", color=MUTED, ms=4, label="lerp (on the chord)")
ax.plot(nl[:, 0] * 1.0, nl[:, 1] * 1.0, "o", color=PALETTE[1], ms=7, label="normalised lerp")
ax.plot(sl[:, 0] * 1.12, sl[:, 1] * 1.12, "s", color=PALETTE[0], ms=5.5, label="slerp (drawn just outside)")
ax.text(1.06, -0.05, "$q_0$", fontsize=12, color=INK, va="top")
ax.text(q1[0] - 0.12, q1[1] + 0.2, "$q_1$", fontsize=12, color=INK, ha="right")
ax.text(0.17, 0.05, "$\\Omega$", fontsize=12, color=INK)
aa = np.linspace(0, Om, 50); ax.plot(0.13 * np.cos(aa), 0.13 * np.sin(aa), color=INK, lw=1)
ax.set_xlim(-1.35, 1.35); ax.set_ylim(-0.45, 1.4)
ax.legend(loc="lower center", fontsize=8.5, ncol=1, bbox_to_anchor=(0.5, -0.02))
ax.set_title("Nine equally spaced values of $t$", fontsize=11)

tt = np.linspace(0, 1, 200)
lp = np.outer(1 - tt, q0) + np.outer(tt, q1)
ang = np.degrees(np.arctan2(lp[:, 1], lp[:, 0]))
ax2.plot(tt, np.degrees(Om) * tt, color=PALETTE[0], label="slerp: angle $= t\\,\\Omega$")
ax2.plot(tt, ang, color=PALETTE[1], label="normalised lerp")
ax2.set_xlabel("interpolation parameter $t$")
ax2.set_ylabel("angle travelled from $q_0$ (degrees)")
ax2.set_title("Same endpoints, different speeds", fontsize=11)
ax2.legend(loc="upper left")
save(fig, "ch70-slerp")
