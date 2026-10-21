"""Figures for Chapter 76 — Signal Processing for Audio and Visual Computing."""
import sys, pathlib; sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
from style import plt, np, PALETTE, INK, MUTED, new_figure, diagram, node, box, arrow, save
from scipy.signal import convolve2d

BLUE, RED, GREEN, PURPLE, ORANGE, TEAL, GREY = PALETTE

# ---------------------------------------------------------------- figure 1: kernels on a synthetic image
N = 96
y, x = np.mgrid[0:N, 0:N]
img = np.zeros((N, N))
img[(x - 30) ** 2 + (y - 34) ** 2 < 18 ** 2] = 1.0           # a disc
img[55:85, 50:86] = 0.6                                      # a rectangle
img[10:40, 62:66] = 0.85                                     # a thin bar
img += 0.25 * (np.sin(2 * np.pi * x / 6) > 0) * (y > 88)     # fine stripes along the bottom
rng = np.random.default_rng(0)
img += 0.05 * rng.normal(size=img.shape)

box5 = np.ones((5, 5)) / 25
sharp = np.array([[0, -1, 0], [-1, 5, -1], [0, -1, 0]], float)
sx = np.array([[-1, 0, 1], [-2, 0, 2], [-1, 0, 1]], float)
sy = sx.T


def conv(k):
    return convolve2d(img, k, mode="same", boundary="symm")


panels = [("original", img, None),
          ("box blur 5x5", conv(box5), box5),
          ("sharpen", conv(sharp), sharp),
          ("Sobel edge magnitude", np.hypot(conv(sx), conv(sy)), sx)]

fig, axs = new_figure(2, 4, figsize=(11.0, 5.0), gridspec_kw=dict(height_ratios=[3, 1.0]))
for col, (title, im, k) in enumerate(panels):
    ax = axs[0, col]
    ax.imshow(im, cmap="gray", vmin=np.percentile(im, 1), vmax=np.percentile(im, 99.5))
    ax.set_title(title, fontsize=10.5); ax.axis("off")
    ax = axs[1, col]
    ax.axis("off")
    if k is None:
        ax.text(0.5, 0.5, "kernel:\n(none)", ha="center", va="center", fontsize=10, color=MUTED,
                transform=ax.transAxes)
        continue
    m = k.shape[0]
    ax.set_xlim(-0.5, 2.5 if m == 3 else 4.5); ax.set_ylim(m - 0.5, -0.5); ax.set_aspect("equal")
    for i in range(m):
        for j in range(m):
            v = k[i, j]
            col_ = BLUE if v > 0 else (RED if v < 0 else "white")
            ax.add_patch(plt.Rectangle((j - 0.5, i - 0.5), 1, 1, fc=col_, alpha=0.18 if v else 1, ec=INK, lw=0.8))
            s = "1" if m == 5 else f"{v:g}"
            ax.text(j, i, s, ha="center", va="center", fontsize=8 if m == 5 else 9.5, color=INK)
    if m == 5:
        ax.text(5.0, 2, r"$\times\frac{1}{25}$", ha="left", va="center", fontsize=11, color=INK)
        ax.set_xlim(-0.5, 7.0)
    if col == 3:
        ax.text(1, m - 0.1 + 0.6, "(x-kernel; y-kernel is its transpose)", ha="center", va="top", fontsize=8, color=MUTED)
save(fig, "ch76-image-kernels")

# ---------------------------------------------------------------- figure 2: spectrogram of a chirp
fs = 8000
t = np.arange(0, 2.0, 1 / fs)
f0, f1 = 200.0, 3000.0
k = (f1 - f0) / 2.0
sig = np.cos(2 * np.pi * (f0 * t + 0.5 * k * t ** 2)) + 0.6 * np.cos(2 * np.pi * 1000 * t) * (t > 1.2)


def stft(xs, L, hop):
    w = np.hanning(L)
    frames = np.array([xs[i:i + L] * w for i in range(0, len(xs) - L, hop)])
    return np.abs(np.fft.rfft(frames, axis=1)).T, (np.arange(len(frames)) * hop + L / 2) / fs


fig, axs = new_figure(1, 3, figsize=(11.0, 3.8), gridspec_kw=dict(width_ratios=[1.1, 1, 1]))
ax = axs[0]
seg = t < 0.03
ax.plot(t[seg] * 1000, sig[seg], color=BLUE, lw=1.2)
ax.set_xlabel("time (ms)"); ax.set_title("first 30 ms of the signal", fontsize=10.5)
for ax, L, name in [(axs[1], 64, "short window (8 ms)"), (axs[2], 1024, "long window (128 ms)")]:
    S, tt = stft(sig, L, L // 4)
    ff = np.fft.rfftfreq(L, 1 / fs)
    ax.imshow(20 * np.log10(S + 1e-6), origin="lower", aspect="auto", cmap="magma",
              extent=[tt[0], tt[-1], 0, fs / 2], vmin=20 * np.log10(S.max()) - 60, interpolation="nearest")
    ax.set_xlabel("time (s)"); ax.set_ylabel("frequency (Hz)")
    ax.set_title(name, fontsize=10.5); ax.grid(False)
    ax.set_ylim(0, 4000)
save(fig, "ch76-chirp-spectrogram")

# ---------------------------------------------------------------- figure 3: aliasing of a zone plate
M = 512
yy, xx = (np.mgrid[0:M, 0:M] - M / 2) / M
zone = 0.5 + 0.5 * np.cos(2 * np.pi * 110 * (xx ** 2 + yy ** 2))
D = 8
naive = zone[::D, ::D]
g = np.exp(-0.5 * (np.arange(-16, 17) / (0.6 * D)) ** 2); g /= g.sum()
blur = convolve2d(convolve2d(zone, g[None, :], mode="same", boundary="symm"), g[:, None], mode="same", boundary="symm")
pre = blur[::D, ::D]

fig, axs = new_figure(1, 3, figsize=(11.0, 4.0))
for ax, im, title in [(axs[0], zone, f"fine rings, {M}x{M}"),
                      (axs[1], naive, f"keep every {D}th pixel: false rings"),
                      (axs[2], pre, "blur first, then keep every 8th")]:
    ax.imshow(im, cmap="gray", vmin=0, vmax=1, interpolation="nearest")
    ax.set_title(title, fontsize=10.5); ax.axis("off")
save(fig, "ch76-aliasing-rings")
