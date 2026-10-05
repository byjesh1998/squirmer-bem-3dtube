"""Figure for Stage 2 (squirmer flow field and swimming kinematics).

usage:  python3 tests/plot_squirmer.py [data_dir] [output.png]
default data_dir = outputs/tests/squirmer  (written by build/test_squirmer)
"""
import sys
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, Rectangle

U0 = 2 / 3


def load(path):
    """Read a CSV written by the C++ code: '# ...' comment lines, then a header line."""
    with open(path) as fh:
        lines = [ln for ln in fh if not ln.startswith("#")]
    return np.genfromtxt(lines, delimiter=",", names=True)


def haberman_sayre(l):
    return (1 - 0.75857 * l**5) / (1 - 2.1050 * l + 2.0865 * l**3 - 1.7068 * l**5 + 0.72603 * l**6)


def field_grid(data_dir):
    f = load(f"{data_dir}/field.csv")
    xs, zs = np.unique(f["x"]), np.unique(f["z"])
    sh = (len(xs), len(zs))
    get = lambda k: f[k].reshape(sh)
    R = np.abs(xs).max() / 0.97
    return xs, zs, R, get


def draw_field(ax, xs, zs, R, ux, uz, title, cbar_fig=None):
    X, Z = np.meshgrid(xs, zs, indexing="ij")
    sp = np.hypot(ux, uz) / U0
    pc = ax.pcolormesh(Z, X, sp, shading="gouraud", cmap="viridis", vmin=0, vmax=1.2)
    ax.streamplot(zs, xs, np.nan_to_num(uz), np.nan_to_num(ux), color="w", density=1.3,
                  linewidth=0.7, arrowsize=0.8)
    ax.add_patch(Circle((0, 0), 1, fc="0.85", ec="k", zorder=3))
    ax.annotate("", xy=(0.8, 0), xytext=(-0.8, 0),
                arrowprops=dict(arrowstyle="->", lw=1.8, color="C3"), zorder=4)
    for s in (-1, 1):
        ax.add_patch(Rectangle((zs.min(), s * R), zs.max() - zs.min(), 0.25 * s, color="0.3"))
    ax.set_aspect("equal")
    ax.set_xlim(zs.min(), zs.max())
    ax.set_ylim(-R - 0.3, R + 0.3)
    ax.set_title(title, fontsize=11)
    ax.set_xlabel("z / a  (tube axis)")
    ax.set_ylabel("x / a")
    return pc


def plot(data_dir="outputs/tests/squirmer", out_png=None):
    xs, zs, R, g = field_grid(data_dir)
    fig = plt.figure(figsize=(16, 8.6))
    gs = fig.add_gridspec(2, 3, height_ratios=[0.62, 1], hspace=0.12, wspace=0.42)
    for k, (alpha, name) in enumerate([(-3, "Pusher  (α = −3)"), (0, "Neutral  (α = 0)"),
                                       (3, "Puller  (α = +3)")]):
        ax = fig.add_subplot(gs[0, k])
        ux = g("ux_B1") + alpha * g("ux_B2")
        uz = g("uz_B1") + alpha * g("uz_B2")
        pc = draw_field(ax, xs, zs, R, ux, uz, name + ",  lab frame")
        if k == 2:
            cb = fig.colorbar(pc, ax=ax, fraction=0.03, pad=0.02)
            cb.set_label("|u| / U₀")

    ax = fig.add_subplot(gs[1, 0])
    A = load(f"{data_dir}/axis.csv")
    ll = np.linspace(0.02, 0.62, 100)
    ax.semilogy(ll, haberman_sayre(ll), "k-", label="Haberman & Sayre (1958) fit")
    ax.semilogy(A["a_over_R"], A["drag_K"], "o", mfc="none", ms=8, label="BEM, towed rigid sphere")
    ax.set(xlabel="a / R", ylabel="drag / 6πμaU", title="Validation: sphere on tube axis")
    ax.legend(fontsize=9, loc="upper left")
    ax2 = ax.twinx()
    ax2.plot(A["a_over_R"], A["U_B1"] / U0, "s-", color="C3", label="squirmer U / U₀")
    rec = load(f"{data_dir}/reciprocal.csv")
    ax2.plot([0.3], [float(rec["U_reciprocal"]) / U0], "x", color="k", ms=10, mew=2,
             label="reciprocal theorem")
    ax2.set_ylabel("squirmer speed  U / U₀", color="C3")
    ax2.set_ylim(0.6, 1.02)
    ax2.legend(loc="lower left", fontsize=9)

    O = load(f"{data_dir}/offaxis.csv")
    ax = fig.add_subplot(gs[1, 1])
    for alpha, st in [(-3, "v--"), (0, "o-"), (3, "^-.")]:
        ax.plot(O["beta"], (O["UB1x"] + alpha * O["UB2x"]) / U0, st, label=f"α = {alpha:+d}")
    ax.axhline(0, color="k", lw=0.6)
    ax.set(xlabel="β = b / (R − a)", ylabel="$U_x$ / U₀   (> 0: toward nearest wall)",
           title="Transverse velocity, a/R = 0.3")
    ax.legend(fontsize=9)

    ax = fig.add_subplot(gs[1, 2])
    ax.plot(O["beta"], O["UB1z"] / U0, "o-", color="C0", label="$U_z$/U₀ (axial, any α)")
    ax.set(xlabel="β = b / (R − a)", ylabel="axial speed  $U_z$ / U₀",
           title="Axial speed and rotation, a/R = 0.3")
    ax3 = ax.twinx()
    ax3.plot(O["beta"], O["OB1y"] / U0, "s--", color="C3", label="$Ω_y$ a / U₀ (any α)")
    ax3.set_ylabel("$Ω_y$ a / U₀   (< 0: turns away from wall)", color="C3")
    h1, l1 = ax.get_legend_handles_labels()
    h2, l2 = ax3.get_legend_handles_labels()
    ax.legend(h1 + h2, l1 + l2, fontsize=9, loc="lower left")
    if out_png:
        fig.savefig(out_png, dpi=120, bbox_inches="tight")
        print("saved", out_png)
    return fig


if __name__ == "__main__":
    matplotlib.use("Agg")
    d = sys.argv[1] if len(sys.argv) > 1 else "outputs/tests/squirmer"
    o = sys.argv[2] if len(sys.argv) > 2 else "outputs/tests/fig_squirmer.png"
    plot(d, o)
