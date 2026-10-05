"""Figure for Stage 3 (trajectories).

usage:  python3 tests/plot_trajectories.py [data_dir] [output.png]
Plots every known run found in data_dir:
  neutral, pusher, puller_a3, puller_a5   (production runs, --full or inputs/trajectory_*.in)
  neutral_quick, pusher_quick             (quick test runs)
default data_dir = outputs/reference/trajectories (pre-computed production runs)
"""
import os
import sys
import numpy as np
import matplotlib
import matplotlib.pyplot as plt

A_OVER_R = 0.3
R = 1 / A_OVER_R
BMAX = R - 1
U0 = 2 / 3
CASES = [
    ("neutral", "Neutral (α = 0)", "C0", "periodic wave, amplitude conserved"),
    ("pusher", "Pusher (α = −3)", "C3", "unstable → hits the wall"),
    ("puller_a3", "Puller (α = +3)", "C2", "decays onto the axis"),
    ("puller_a5", "Puller (α = +5)", "C1", "settles near the wall, tilted toward it"),
    ("neutral_quick", "Neutral (α = 0), quick test", "C0", "half-wave"),
    ("pusher_quick", "Pusher (α = −3), quick test", "C3", "hits the wall"),
]


def load(path):
    """Read a CSV written by the C++ code: '# ...' comment lines, then a header line."""
    with open(path) as fh:
        lines = [ln for ln in fh if not ln.startswith("#")]
    return np.genfromtxt(lines, delimiter=",", names=True)


def draw_path(ax, fig, d, title, col, note):
    x, z, ex, ez = d["x"], d["z"], d["ex"], d["ez"]
    zl = (-2, z.max() + 4)
    ax.axhspan(R, R + 0.5, color="0.35")
    ax.axhspan(-R - 0.5, -R, color="0.35")
    for s in (-1, 1):
        ax.axhline(s * BMAX, color="0.6", ls=":", lw=1)
    ax.axhline(0, color="0.75", lw=0.8, ls="--")
    ax.plot(z, x, color=col, lw=2.2)
    ax.set_xlim(*zl)
    ax.set_ylim(-R - 0.5, R + 0.5)
    fig.canvas.draw()
    bb = ax.get_window_extent()
    sx, sy = bb.width / (zl[1] - zl[0]), bb.height / (2 * R + 1)  # pixels per unit length
    for i in np.linspace(0, len(z) - 1, 9).astype(int):
        L = 28  # arrow length in pixels, drawn at the TRUE angle
        ax.annotate("", xy=(z[i] + L * ez[i] / sx, x[i] + L * ex[i] / sy), xytext=(z[i], x[i]),
                    arrowprops=dict(arrowstyle="-|>", color="k", lw=1.2, mutation_scale=10))
        ax.plot(z[i], x[i], "o", color=col, ms=7, mec="k", mew=0.6)
    if d["beta"][-1] > 0.9:
        ax.plot(z[-1], x[-1], "kx", ms=12, mew=2.5)
    ax.set_title(title, fontsize=11, loc="left")
    ax.text(0.99, 0.04, note, transform=ax.transAxes, ha="right", fontsize=9, style="italic",
            bbox=dict(fc="w", ec="none", alpha=0.8))
    ax.set_ylabel("x / a")
    ax.set_xlabel("z / a  (along the tube)")


def plot(data_dir="outputs/reference/trajectories", out_png=None):
    found = [c for c in CASES if os.path.exists(f"{data_dir}/{c[0]}.csv")]
    if not found:
        raise FileNotFoundError(f"no trajectory CSV files in {data_dir}")
    D = {c[0]: load(f"{data_dir}/{c[0]}.csv") for c in found}
    nrow = (len(found) + 1) // 2
    fig = plt.figure(figsize=(15, 3.6 * nrow + 3.8))
    gs = fig.add_gridspec(nrow + 1, 2, height_ratios=[1] * nrow + [1.05], hspace=0.42, wspace=0.18)
    for k, (name, title, col, note) in enumerate(found):
        draw_path(fig.add_subplot(gs[k // 2, k % 2]), fig, D[name], title + f", start β = {D[name]['beta'][0]:.1f}", col, note)
    ax1, ax2 = fig.add_subplot(gs[nrow, 0]), fig.add_subplot(gs[nrow, 1])
    for name, title, col, _ in found:
        d = D[name]
        ax1.plot(d["t"], d["x"] / BMAX, color=col, lw=1.8, label=title)
        ax2.plot(d["t"], d["Uz"] / U0, color=col, lw=1.8, label=title)
    for s in (1, -1):
        ax1.axhline(s, color="0.5", ls=":")
    ax1.set(xlabel="time  t B$_1$/a", ylabel="x / (R − a)     (±1 = touching wall)",
            title="Transverse position", ylim=(-1.1, 1.1))
    ax2.axhline(1, color="k", lw=0.8, ls="--")
    ax2.set(xlabel="time  t B$_1$/a", ylabel="axial speed  $U_z$ / U$_0$", title="Axial swimming speed")
    ax1.legend(fontsize=8, loc="lower right")
    ax2.legend(fontsize=8, loc="lower left")
    fig.suptitle(f"Squirmer trajectories in a tube, a/R = {A_OVER_R}   (transverse direction "
                 "stretched; arrows show the true swimming direction)", fontsize=12, y=0.95)
    if out_png:
        fig.savefig(out_png, dpi=120, bbox_inches="tight")
        print("saved", out_png)
    return fig


if __name__ == "__main__":
    matplotlib.use("Agg")
    d = sys.argv[1] if len(sys.argv) > 1 else "outputs/reference/trajectories"
    o = sys.argv[2] if len(sys.argv) > 2 else "outputs/tests/fig_trajectories.png"
    plot(d, o)
