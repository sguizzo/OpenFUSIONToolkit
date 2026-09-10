#!/usr/bin/env python3
"""Axial electromagnetic force on the fluid, from an em_pump run with flow_on = F.

Usage:  python3 pump_force.py RUNDIR

WHAT IT READS
-------------
`oft_xdmf.0001.h5`, written by pump_sim%plot() at the end of the run.  Each frame
holds `psi` and `B*R` on the order-2 Lagrange DOFs.  `B*R` is the L2-projected
poloidal gradient -- (-dpsi/dZ, F, dpsi/dR).

EQUATION
-------------------------------
Everything is axisymmetric and the only current is toroidal, so with

    B_R = -(1/R) dpsi/dZ        B_Z = (1/R) dpsi/dR

the axial Lorentz force density is

    (J x B)_Z = -J_phi B_R = J_phi (1/R) dpsi/dZ

Integrating over dV = 2*pi*R dR dZ, the R cancels:

    F_Z = -2*pi * INT J_phi (R B_R) dR dZ        [R*B_R is B*R component 0]

With the flow frozen the fluid carries only induced current, so Ohm's law gives
J_phi straight from psi -- no need to evaluate the Grad-Shafranov operator:

    E_phi = -dA_phi/dt = -(1/R) dpsi/dt      J_phi = E_phi/eta = -(1/(eta R)) dpsi/dt

    F_Z = -(2*pi/eta) * INT (1/R) (dpsi/dt)(dpsi/dZ) dR dZ

"""

import sys
import os

import numpy as np
import h5py
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# --- Inputs ------------------------------------------------------------------
run_dir = sys.argv[1] if len(sys.argv) > 1 else "."
eta = 1.625e-3          # fluid resistivity [Ohm-m]; must match eta_fluid in oft.in
freq = 50.0             # drive frequency [Hz]; must match freq in oft.in
fluid_region = 1        # region id of the fluid
plot_frame = -1         # frame to draw the field snapshots from; -1 is the last


# --- Read Mesh --------------------------------------------------------------------
mesh = h5py.File(os.path.join(run_dir, "pump_mesh_native.h5"), "r")
nodes = mesh["mesh/R"][:]              # (npoint, 2) as (R, Z)
tris = mesh["mesh/LC"][:] - 1          # 1-based in the file
region = mesh["mesh/REG"][:]
fluid_tris = tris[region == fluid_region]
n_node = nodes.shape[0]

# Per-cell area and centroid radius, computed once.
areas, radii = [], []
for tri in fluid_tris:
    (r0, z0), (r1, z1), (r2, z2) = nodes[tri]
    edges = np.array([[r1 - r0, z1 - z0],
                      [r2 - r0, z2 - z0]])
    areas.append(0.5 * abs(np.linalg.det(edges)))
    radii.append((r0 + r1 + r2) / 3.0)
areas = np.array(areas)
radii = np.array(radii)


# --- Read psi(t) and the projected R*B_R from the plot file ------------------

plot = h5py.File(os.path.join(run_dir, "oft_xdmf.0001.h5"), "r")["mugtok_td/smesh"]
frames = sorted(k for k in plot if isinstance(plot[k], h5py.Group) and "TIME" in plot[k])
if len(frames) < 3:
    raise SystemExit(f"need at least 3 plot frames, found {len(frames)}; lower rst_freq")
times = np.array([float(plot[k]["TIME"][0]) for k in frames])
psis = np.array([plot[k]["psi"][:n_node] for k in frames])
RB_R = np.array([plot[k]["B*R"][:n_node, 0] for k in frames])


def axial_force(rb_r, dpsi_dt):
    """F_Z = (2 pi / eta) INT (1/R) (dpsi/dt) (R B_R) dR dZ, summed over cells."""
    cell_rb_r = rb_r[fluid_tris].mean(axis=1)
    cell_dpsi_dt = dpsi_dt[fluid_tris].mean(axis=1)
    return (2.0 * np.pi / eta) * np.sum(cell_dpsi_dt * cell_rb_r / radii * areas)


# Central differences in time, so the first and last samples are dropped
force = np.array([axial_force(RB_R[i], (psis[i + 1] - psis[i - 1]) / (times[i + 1] - times[i - 1]))
                  for i in range(1, len(frames) - 1)])
t = times[1:-1]

# psi oscillates at the drive frequency, so a central difference over the frame
# spacing h under-reads dpsi/dt by sinc(w*h) -- 10% at h = 2.5 ms, 50 Hz. Undo it,
# otherwise the answer depends on rst_freq rather than on the physics.
wh = 2.0 * np.pi * freq * (times[1] - times[0])
force /= np.sin(wh) / wh

bore = nodes[fluid_tris][:, :, 0].max()
head = force / (np.pi * bore ** 2)      # equivalent stall pressure rise [Pa]


# --- Report ------------------------------------------------------------------
settled = t > 0.5 * t.max()             # ignore the drive ramp
print(f"  run              {run_dir}")
print(f"  mean F_Z         {force[settled].mean():+.6e} N   (over the second half)")


# --- Plot force vs time --------------------------------------------------------------------
fig, ax = plt.subplots(figsize=(9, 4), facecolor="#fcfcfb")
ax.set_facecolor("#fcfcfb")
ax.plot(t, force, lw=2, color="#2a78d6")
ax.axhline(0, color="#52514e", lw=1, ls="--", alpha=0.6)
ax.annotate(f"mean {force[settled].mean():+.4g} N\n= {head[settled].mean():+.1f} Pa across the bore",
            xy=(t[-1], force[-1]), xytext=(-10, 12), textcoords="offset points",
            ha="right", color="#2a78d6", fontsize=9, fontweight="bold")
for side in ("top", "right"):
    ax.spines[side].set_visible(False)
for side in ("left", "bottom"):
    ax.spines[side].set_color("#e4e3df")
ax.grid(True, color="#e4e3df", lw=0.8)
ax.set_axisbelow(True)
ax.tick_params(colors="#52514e", labelsize=9)
ax.set_xlabel("time [s]", color="#52514e")
ax.set_ylabel("axial force on the fluid  F_Z [N]", color="#52514e")
ax.set_title(f"stall force,  {os.path.basename(os.path.abspath(run_dir))}",
             color="#0b0b0b", fontsize=12, loc="left", pad=8)
fig.tight_layout()
out = os.path.join(run_dir, "force.png")
fig.savefig(out, dpi=150, facecolor="#fcfcfb")
print(f"  wrote            {out}")

# --- Plot fields at the last timestep (cell centered) --------------------------------------------------------------------

import matplotlib.tri as mtri
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm

GS_EPSILON = 1.0e-12          # same guard the Fortran uses
DIVERGING = LinearSegmentedColormap.from_list(
    "bl_gy_or", ["#184f95", "#2a78d6", "#9ec5f4", "#e8e8e6", "#f6b596", "#eb6834", "#a83c15"])

# psi_dot needs a neighbour on each side, so index the frames the force uses
iframe = range(1, len(frames) - 1)[plot_frame]
psi_dot = (psis[iframe + 1] - psis[iframe - 1]) / (times[iframe + 1] - times[iframe - 1])
RB_Z = plot[frames[iframe]]["B*R"][:n_node, 2]

# Node fields -> cell values, then divide by the centroid radius
cell_psi_dot = psi_dot[fluid_tris].mean(axis=1)
cell_RB_R = RB_R[iframe][fluid_tris].mean(axis=1)
cell_RB_Z = RB_Z[fluid_tris].mean(axis=1)

J_phi = -cell_psi_dot / (eta * (radii + GS_EPSILON))       # A/m^2
B_R = cell_RB_R / (radii + GS_EPSILON)                     # T
B_Z = cell_RB_Z / (radii + GS_EPSILON)                     # T
f_Z = -J_phi * B_R                                         # N/m^3

tri_plot = mtri.Triangulation(nodes[:, 1], nodes[:, 0], triangles=fluid_tris)
panels = [(J_phi, "toroidal current density  J_phi  [A/m$^2$]"),
          (B_R,   "radial field  B_R  [T]"),
          (B_Z,   "axial field  B_Z  [T]"),
          (f_Z,   "axial force density  (JxB)_Z  [N/m$^3$]")]

fig, axes = plt.subplots(len(panels), 1, figsize=(11, 9), facecolor="#fcfcfb")
for ax, (field, label) in zip(axes, panels):
    lim = np.percentile(np.abs(field), 99.5)
    # flat shading: one value per cell, no interpolation back onto the axis
    img = ax.tripcolor(tri_plot, facecolors=field, cmap=DIVERGING,
                       norm=TwoSlopeNorm(vmin=-lim, vcenter=0.0, vmax=lim))
    bar = fig.colorbar(img, ax=ax, pad=0.01)
    bar.ax.tick_params(colors="#52514e", labelsize=8)
    bar.outline.set_visible(False)
    ax.axvline(-0.36, color="#52514e", lw=1, ls="--")     # coil stack ends
    ax.axvline(+0.36, color="#52514e", lw=1, ls="--")
    ax.set_xlim(-0.45, 0.45)
    ax.set_ylim(0.0, nodes[fluid_tris][:, :, 0].max())
    ax.set_facecolor("#fcfcfb")
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color("#e4e3df")
    ax.tick_params(colors="#52514e", labelsize=9)
    ax.set_ylabel("R [m]", color="#52514e")
    ax.set_title(label, color="#0b0b0b", fontsize=11, loc="left", pad=6)
axes[-1].set_xlabel("Z [m]   (dashed = coil stack ends; wave travels +Z)", color="#52514e")
fig.suptitle(f"fields in the fluid at t = {times[iframe]:.4f} s  (cell-centred)",
             color="#0b0b0b", fontsize=13, fontweight="bold", x=0.01, ha="left")
fig.tight_layout(rect=[0, 0, 1, 0.97])
out_fields = os.path.join(run_dir, "fields.png")
fig.savefig(out_fields, dpi=150, facecolor="#fcfcfb")
print(f"  wrote            {out_fields}")

