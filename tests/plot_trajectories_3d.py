"""3D trajectories (swimmer started with an out-of-plane yaw).

usage: python3 tests/plot_trajectories_3d.py [data_dir] [output.png]
expects neutral_yaw15.csv and puller_a3_yaw15.csv in data_dir (default outputs/traj3d)
"""
import sys
import numpy as np
import matplotlib
import matplotlib.pyplot as plt

R, A = 1 / 0.3, 1.0
U0 = 2 / 3
CASES = [("neutral_yaw15", "Neutral (α = 0)", "C1"), ("puller_a3_yaw15", "Puller (α = +3)", "C2")]


def load(path):
    with open(path) as fh:
        lines = [ln for ln in fh if not ln.startswith("#")]
    return np.genfromtxt(lines, delimiter=",", names=True)

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LightSource  # <-- Required for realistic lighting

# Define global variables used in your original snippet
R = 5.0  



# Define global variables used in your original snippet
R = 5.0  

def tube_wire_with_lids(ax, zmax):
    # 1. Original Cylinder Body Code
    th = np.linspace(0, 2 * np.pi, 68)
    zz = np.linspace(0, zmax, 50)
    T, Z = np.meshgrid(th, zz)
    ax.plot_surface(Z, R * np.cos(T), R * np.sin(T), color="blue", alpha=0.1, linewidth=0)
    
    # 2. ADDED: Solid Flat Circular Lids
    # Create a grid that sweeps radially from the center (0) out to the outer radius (R)
    r_mesh = np.linspace(0, R, 20)
    R_grid, T_grid = np.meshgrid(r_mesh, th)
    
    # Flat Lid at the start (Z = 0)
    ax.plot_surface(np.zeros_like(R_grid), R_grid * np.cos(T_grid), R_grid * np.sin(T_grid), 
                    color="black", alpha=0.4, linewidth=0)
    
    # Flat Lid at the end (Z = zmax)
    ax.plot_surface(np.full_like(R_grid, zmax), R_grid * np.cos(T_grid), R_grid * np.sin(T_grid), 
                    color="black", alpha=0.4, linewidth=0)

    # 3. Original Structural Rings & Axis Code
    for z0 in np.arange(0, zmax + 1, 20):
        ax.plot(np.full_like(th, z0), R * np.cos(th), R * np.sin(th), color="0.65", lw=0.6)
    ax.plot(zz, 0 * zz, 0 * zz, "--", color="0.5", lw=1.)   # tube axis



def tube_wire(ax, zmax):
    th = np.linspace(0, 2 * np.pi, 68)
    zz = np.linspace(0, zmax, 50)
    T, Z = np.meshgrid(th, zz)
    ax.plot_surface(Z, R * np.cos(T), R * np.sin(T), color="0.75", alpha=0.18, linewidth=0)
    for z0 in np.arange(0, zmax + 1, 20):
        ax.plot(np.full_like(th, z0), R * np.cos(th), R * np.sin(th), color="0.65", lw=0.6)
    ax.plot(zz, 0 * zz, 0 * zz, "--", color="0.5", lw=1.)   # tube axis


def plot(data_dir="outputs/reference/traj3d", out_png=None):
    D = {k: load(f"{data_dir}/{k}.csv") for k, _, _ in CASES}
    fig = plt.figure(figsize=(10, 10))
    gs = fig.add_gridspec(2, 2, width_ratios=[5., 2], hspace=0.015, wspace=0.015)
    th = np.linspace(0, 2 * np.pi, 200)
    for row, (k, title, col) in enumerate(CASES):
        d = D[k]
        x, y, z = d["x"], d["y"], d["z"]
        # --- 3D view (axial direction compressed so the helix is visible)
        ax = fig.add_subplot(gs[row, 0], projection="3d")
        tube_wire_with_lids(ax, z.max())
        ax.plot(z, x, y, color=col, lw=3.4)
        # ax.plot(z, 0 * x + R, y, color=col, lw=0.8, alpha=0.35)   # shadow on the back wall
        ax.scatter([z[0]], [x[0]], [y[0]], color="k", s=30)
        idx = np.linspace(0, len(z) - 1, 14).astype(int)
        ax.quiver(z[idx], x[idx], y[idx], d["ez"][idx], d["ex"][idx], d["ey"][idx],
                  length=4, normalize=True, color="k", lw=0.8)
        ax.set_box_aspect((4., 1.1, 1.1), zoom=1.)
        ax.set_xlabel("z / a  (along the tube)",labelpad=25)
        # ax.set(xlabel="z / a  (along the tube)", ylim=(-R, R), zlim=(-R, R))
        ax.set_yticks([]); ax.set_zticks([])
        ax.view_init(elev=24, azim=-60)
        ax.set_title(f"{title}: start β = {d['beta'][0]:.1f}, yaw 15°", loc="left")
        ax.grid(False) #


        ax.xaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
        ax.yaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
        ax.zaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))


        # for axis in [ax.xaxis, ax.yaxis, ax.zaxis]:
        #     axis._axinfo['axisline']['linewidth'] = 1.5  # Set thickness of the 3 lines
        #     axis._axinfo['axisline']['color'] = 'black'  # Set color of the directional lines


        # --- cross-section (looking down the tube)
        ax = fig.add_subplot(gs[row, 1])
        ax.fill(R * np.cos(th), R * np.sin(th), color="0.95")
        ax.plot(R * np.cos(th), R * np.sin(th), "k", lw=1.5)
        ax.plot((R - A) * np.cos(th), (R - A) * np.sin(th), ":", color="0.6", lw=1)
        sc = ax.scatter(x, y, c=d["t"], cmap="plasma", s=6)
        ax.plot(x[0], y[0], "ko", ms=6)
        ax.set_aspect("equal")
        ax.set(xlabel="x / a", ylabel="y / a", title="Cross-section view (colour = time)")
        fig.colorbar(sc, ax=ax, fraction=0.046, label="t")
        # --- radial distance and azimuth vs time
        # ax = fig.add_subplot(gs[row, 2])
        # r = np.hypot(x, y) / (R - A)
        # phi = np.degrees(np.unwrap(np.arctan2(y, x)))
        # phi = np.where(r > 0.03, phi, np.nan)   # azimuth undefined on the axis
        # # ax.plot(d["t"], r, color=col, lw=2, label="β = r/(R−a)")
        # ax.set(xlabel="t", ylabel="β", ylim=(0, 1), title="Distance from axis & azimuth")
        # ax2 = ax.twinx()
        # # ax2.plot(d["t"], phi - phi[0], "--", color="0.4", label="azimuth turned")
        # ax2.set_ylabel("azimuth turned (deg)")
        # h1, l1 = ax.get_legend_handles_labels(); h2, l2 = ax2.get_legend_handles_labels()
        # ax.legend(h1 + h2, l1 + l2, fontsize=8, loc="upper right")
    fig.suptitle("3D squirmer trajectories in a tube, a/R = 0.3",
                 fontsize=13, y=0.96)
    if out_png:
        fig.savefig(out_png, dpi=115, bbox_inches="tight")
        print("saved", out_png)
    return fig


if __name__ == "__main__":
    matplotlib.use("Agg")
    d = sys.argv[1] if len(sys.argv) > 1 else "outputs/reference/traj3d"
    o = sys.argv[2] if len(sys.argv) > 2 else "outputs/reference/traj3d/fig_trajectories_3d.png"
    plot(d, o)
