import sys, pathlib; sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
from style import plt, np, PALETTE, INK, MUTED, new_figure, diagram, node, box, arrow, save
from matplotlib.patches import Polygon


def look_at(eye, target, up=(0.0, 0.0, 1.0)):
    """Orthonormal camera frame (right, up, forward) for a camera at eye looking at target."""
    eye, target, up = map(np.asarray, (eye, target, up))
    fwd = target - eye
    fwd = fwd / np.linalg.norm(fwd)
    right = np.cross(fwd, up)
    right = right / np.linalg.norm(right)
    upv = np.cross(right, fwd)
    return right, upv, fwd


def project(P, eye, frame, f):
    """Pinhole projection of world points P (n x 3) to screen coordinates (n x 2) and depth."""
    right, upv, fwd = frame
    d = np.atleast_2d(P) - eye
    xc, yc, zc = d @ right, d @ upv, d @ fwd
    return np.c_[f * xc / zc, f * yc / zc], zc


# ---------------------------------------------------------------- the scene
cube = np.array([[x, y, z] for x in (0.6, 1.8) for y in (0.6, 1.8) for z in (0.0, 1.2)])
edges = [(a, b) for a in range(8) for b in range(a + 1, 8) if np.sum(cube[a] != cube[b]) == 1]
faces = [[0, 1, 3, 2], [4, 5, 7, 6], [0, 1, 5, 4], [2, 3, 7, 6], [0, 2, 6, 4], [1, 3, 7, 5]]

cam_eye = np.array([5.2, -2.6, 2.6])
cam_frame = look_at(cam_eye, cube.mean(axis=0))
F_CAM = 1.3                          # focal length of the scene camera (world units)
HALF_W, HALF_H = 0.42, 0.28          # half-size of its image rectangle

# the observer who draws the whole scene for us (a second, wider camera)
obs_eye = np.array([3.0, -9.5, 6.0])
obs_frame = look_at(obs_eye, np.array([2.4, 0.2, 1.0]))
F_OBS = 3.0


def O(P):
    return project(P, obs_eye, obs_frame, F_OBS)[0]


def image_point(s):
    """World position of screen point s=(xs, ys) on the scene camera's image plane."""
    r, u, fw = cam_frame
    return cam_eye + F_CAM * fw + s[0] * r + s[1] * u


fig, (ax, ax2) = plt.subplots(1, 2, figsize=(10.4, 4.6), constrained_layout=True,
                              gridspec_kw=dict(width_ratios=[1.55, 1]))
ax.set_aspect("equal"); ax.axis("off")

# ground grid (z = 0)
for g in np.linspace(-1, 3, 5):
    for line in (np.array([[g, -1, 0], [g, 3, 0]]), np.array([[-1, g, 0], [3, g, 0]])):
        q = O(line); ax.plot(q[:, 0], q[:, 1], color="#e5e7eb", lw=0.8, zorder=0)

# world frame
origin = np.zeros(3)
for v, lab, c in zip(np.eye(3), ["$x_w$", "$y_w$", "$z_w$"], PALETTE[:3]):
    q = O(np.array([origin, 1.1 * v]))
    ax.annotate("", xy=q[1], xytext=q[0], arrowprops=dict(arrowstyle="-|>", color=c, lw=1.8))
    ax.text(*(q[1] + 0.012 * (q[1] - q[0]) / np.linalg.norm(q[1] - q[0])), lab, color=c, fontsize=10,
            ha="left", va="bottom")
ax.text(*O(origin)[0], "world origin  ", ha="right", va="top", fontsize=8.5, color=MUTED)

# cube (shaded faces, then edges)
for fc in faces:
    q = O(cube[fc]); ax.add_patch(Polygon(q, closed=True, fc="#dbeafe", ec="none", alpha=0.55, zorder=1))
for a, b in edges:
    q = O(cube[[a, b]]); ax.plot(q[:, 0], q[:, 1], color=PALETTE[0], lw=1.8, zorder=2)

# camera body, frame and image rectangle
s_corners = np.array([[-HALF_W, -HALF_H], [HALF_W, -HALF_H], [HALF_W, HALF_H], [-HALF_W, HALF_H]])
rect = np.array([image_point(s) for s in s_corners])
ax.add_patch(Polygon(O(rect), closed=True, fc="#fff7ed", ec=PALETTE[4], lw=1.5, alpha=0.9, zorder=3))
for c in rect:  # frustum edges
    q = O(np.array([cam_eye, cam_eye + 1.9 * (c - cam_eye)]))
    ax.plot(q[:, 0], q[:, 1], color=PALETTE[4], lw=0.8, ls="--", zorder=2)
# rays from vertices to the eye, and the projected cube drawn on the image plane
scr, depth = project(cube, cam_eye, cam_frame, F_CAM)
hits = np.array([image_point(s) for s in scr])
for P in cube:
    q = O(np.array([P, cam_eye])); ax.plot(q[:, 0], q[:, 1], color=MUTED, lw=0.5, alpha=0.6, zorder=2)
for a, b in edges:
    q = O(hits[[a, b]]); ax.plot(q[:, 0], q[:, 1], color=PALETTE[1], lw=1.3, zorder=4)
qe = O(cam_eye)[0]
ax.plot(*qe, "o", color=INK, ms=5, zorder=5)
ax.text(qe[0] + 0.01, qe[1] - 0.02, "camera\n(pinhole)", fontsize=8.5, color=INK, va="top")
for v, lab, c, ha, va in zip(cam_frame, ["$x_c$", "$y_c$", "$z_c$"], PALETTE[:3],
                             ["left", "left", "center"], ["bottom", "center", "top"]):
    q = O(np.array([cam_eye, cam_eye + 0.7 * v]))
    ax.annotate("", xy=q[1], xytext=q[0], arrowprops=dict(arrowstyle="-|>", color=c, lw=1.5), zorder=6)
    ax.plot(*q[1], alpha=0)
    ax.text(q[1][0], q[1][1] - (0.012 if lab == "$z_c$" else 0), " " + lab + " ", color=c, fontsize=9.5, ha=ha, va=va, zorder=7)
qr = O(rect)
ax.text(qr[3, 0], qr[3, 1] + 0.01, "image plane\n(screen)", color=PALETTE[4], fontsize=8.5, ha="center", va="bottom")
ax.set_title("World, camera and screen frames", fontsize=11)

# right panel: what the camera sees, in screen coordinates
ax2.set_aspect("equal")
ax2.set_xlim(-HALF_W, HALF_W); ax2.set_ylim(-HALF_H, HALF_H)
ax2.set_facecolor("#fff7ed")
order = np.argsort([-depth[fc].mean() for fc in faces])        # painter's algorithm: far faces first
for k in order:
    ax2.add_patch(Polygon(scr[faces[k]], closed=True, fc="#dbeafe", ec=PALETTE[1], lw=1.4, alpha=0.9))
for i, (p, z) in enumerate(zip(scr, depth)):
    ax2.plot(*p, "o", color=PALETTE[1], ms=3.5)
ax2.set_xlabel("$x_s = f\\,x_c / z_c$")
ax2.set_ylabel("$y_s = f\\,y_c / z_c$")
ax2.set_title("The camera's image", fontsize=11)
save(fig, "ch69-camera-frames")

# ---------------------------------------------------------------- pinhole similar triangles (side view)
fig, ax = diagram(figsize=(7.6, 3.9))
ax.set_xlim(-2.6, 10.2); ax.set_ylim(-1.9, 4.6)
ax.axhline(0, color=MUTED, lw=1, ls=":")
f, Z, Y = 2.2, 8.4, 3.2
y_img = f * Y / Z
# object (a tree-like arrow)
ax.annotate("", xy=(Z, Y), xytext=(Z, 0), arrowprops=dict(arrowstyle="-|>", color=PALETTE[2], lw=3))
ax.text(Z + 0.2, Y / 2, "object\nheight $Y$", color=PALETTE[2], fontsize=10, va="center")
# ray through pinhole, extended to back plane
ax.plot([-f, Z], [-y_img, Y], color=PALETTE[0], lw=1.6)
ax.plot(0, 0, "o", color=INK, ms=7)
ax.text(0.05, -0.25, "pinhole $C$", fontsize=9.5, ha="left", va="top", color=INK)
# virtual image plane in front, real (inverted) plane behind
for x, lab, c, yy in [(f, "virtual image plane\n(used in graphics)", PALETTE[4], y_img),
                      (-f, "film / sensor\n(image upside down)", MUTED, -y_img)]:
    ax.plot([x, x], [-1.5, 1.8], color=c, lw=2.2)
    ax.text(x, 1.95, lab, color=c, fontsize=8.5, ha="center", va="bottom")
    ax.plot(x, yy, "o", color=PALETTE[1], ms=6, zorder=5)
ax.text(f + 0.12, y_img - 0.05, "$y_s = fY/Z$", color=PALETTE[1], fontsize=10.5, va="top")
ax.text(-f - 0.12, -y_img, "$-fY/Z$", color=PALETTE[1], fontsize=10, ha="right", va="center")
# dimension lines
for x0, x1, yb, lab in [(0, f, -1.0, "$f$"), (0, Z, -1.6, "$Z$ (depth)")]:
    ax.annotate("", xy=(x1, yb), xytext=(x0, yb), arrowprops=dict(arrowstyle="<->", color=MUTED, lw=1.1))
    ax.text((x0 + x1) / 2, yb + 0.05, lab, ha="center", va="bottom", fontsize=10, color=INK)
ax.fill([0, f, f], [0, 0, y_img], color=PALETTE[4], alpha=0.18)
ax.fill([0, Z, Z], [0, 0, Y], color=PALETTE[0], alpha=0.06)
ax.text(4.8, 3.9, "similar triangles:  $y_s / f = Y / Z$", fontsize=11, color=INK, ha="center")
ax.text(-2.5, 0.1, "optical axis", fontsize=8.5, color=MUTED, ha="left", va="bottom")
save(fig, "ch69-pinhole")
