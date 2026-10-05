"""Figure for Stage 1 (pressure-driven duct flow).

usage:  python3 tests/plot_duct_flow.py [data_dir] [output.png]
default data_dir = outputs/tests/duct_flow  (written by build/test_duct_flow)
"""
import sys
import numpy as np
import matplotlib
import matplotlib.pyplot as plt

R, L, MU = 1.0, 4.0, 1.0


def load(path):
    """Read a CSV written by the C++ code: '# ...' comment lines, then a header line."""
    with open(path) as fh:
        lines = [ln for ln in fh if not ln.startswith("#")]
    return np.genfromtxt(lines, delimiter=",", names=True)


def plot(data_dir="outputs/tests/duct_flow", out_png=None):
    prof = load(f"{data_dir}/profile.csv")
    pan = load(f"{data_dir}/panels.csv")
    fq = load(f"{data_dir}/flow_rate.csv")
    conv = load(f"{data_dir}/convergence.csv")
    dP = 1.0
    G = dP / L
    exact = lambda r: G * (R**2 - r**2) / (4 * MU)

    fig, ax = plt.subplots(2, 2, figsize=(11, 8.5))
    a = ax[0, 0]
    r_ = np.linspace(0, 1, 200)
    a.plot(r_, exact(r_), "k-", label="Hagen–Poiseuille (exact)")
    a.plot(prof["r"], prof["uz"], "o", ms=5, mfc="none", label="BEM interior, z = L/2")
    out = pan["tag"] == 2
    rc = np.hypot(pan["cx"][out], pan["cy"][out])
    a.plot(rc, pan["uz"][out], "x", ms=4, alpha=0.6, label="BEM outlet cap (solved for)")
    a.set(xlabel="r / R", ylabel="$u_z$", title="Axial velocity profile (ΔP = 1)")
    a.legend()

    a = ax[0, 1]
    w = pan["tag"] == 0
    zc = pan["cz"][w]
    er = np.c_[pan["cx"][w], pan["cy"][w]]
    er /= np.linalg.norm(er, axis=1)[:, None]
    fr = pan["fx"][w] * er[:, 0] + pan["fy"][w] * er[:, 1]
    a.plot(zc, -fr, ".", ms=3, alpha=0.5, label="BEM: wall pressure  $-f_r$")
    zz = np.linspace(0, L, 50)
    a.plot(zz, dP - G * zz, "k-", label="exact p(z)")
    a.plot(zc, pan["fz"][w], ".", ms=3, alpha=0.5, label="BEM: wall shear $f_z$")
    a.axhline(-G * R / 2, color="r", lw=1, label="exact shear  −GR/2")
    a.set(xlabel="z", title="Wall traction")
    a.legend(fontsize=8)

    a = ax[1, 0]
    a.plot(fq["dP"], fq["Q"], "o", ms=8, mfc="none", label="BEM")
    d_ = np.linspace(0, fq["dP"].max() * 1.05, 10)
    a.plot(d_, np.pi * R**4 * d_ / L / (8 * MU), "k-", label=r"$Q=\pi R^4\Delta P/(8\mu L)$")
    a.set(xlabel="ΔP", ylabel="Q", title="Flow rate vs pressure drop (linearity)")
    a.legend()

    a = ax[1, 1]
    k = np.polyfit(np.log(conv["h"]), np.log(conv["rel_err_Q"]), 1)[0]
    a.loglog(conv["h"], conv["rel_err_Q"], "o-", label="BEM")
    a.loglog(conv["h"], conv["rel_err_Q"][-1] * (conv["h"] / conv["h"][-1]) ** k, "k--",
             label=f"fitted order {k:.2f}")
    a.set(xlabel="mean panel size h", ylabel="relative error in Q", title="Mesh convergence")
    a.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
    a.set_xticks([0.14, 0.2, 0.28])
    a.set_xticklabels(["0.14", "0.20", "0.28"])
    a.legend()
    fig.tight_layout()
    if out_png:
        fig.savefig(out_png, dpi=130)
        print("saved", out_png)
    return fig


if __name__ == "__main__":
    matplotlib.use("Agg")
    d = sys.argv[1] if len(sys.argv) > 1 else "outputs/tests/duct_flow"
    o = sys.argv[2] if len(sys.argv) > 2 else "outputs/tests/fig_duct_flow.png"
    plot(d, o)
