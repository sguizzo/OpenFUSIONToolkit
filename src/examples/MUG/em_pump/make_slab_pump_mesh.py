#!/usr/bin/env python3
"""Build a double-sided slab (flat) induction pump mesh for the MUG pump drivers.

Edit the INPUTS block below and run:  python3 make_slab_pump_mesh.py

The magnetic field boundary conditions in em_pump are conducting, not truly free-boundary.
So, scans of air_z and air_r should be performed to ensure the air box is large enough that the 
field is not affected by the boundary.
"""

import os
import sys

import numpy as np


# =============================================================================
# INPUTS
# =============================================================================

# --- Transverse positions [m], measured from the channel center -----------
chan_half     = 0.040   # channel half-width
wall_out      = 0.045   # duct wall distance, from channel center
coil_in       = 0.055   # coil inner face, from channel center
coil_out      = 0.105   # coil outer face, from channel center
yoke_out      = 0.155   # stator outer face, from channel center

# --- Axial extents [m] -------------------------------------------------------
chan_len      = 0.72    # full duct length
pump_len      = 0.72    # length of the coil stack
iron_len      = 0.72    # length of the stator iron

# --- Air ---------------------------------------------------------------------
air_x         = 0.70    # air padding beyond each back iron in the x direction
air_z         = 1.00    # air padding on either side of the pump in the z direction

# --- Winding -----------------------------------------------------------------
n_seg         = 24      # coil segments per side
seg_gap       = 0.010   # gap between adjacent segments = stator tooth width

# --- Materials ---------------------------------------------------------------
eta_wall      = 1.2e-6  # duct wall resistivity [Ohm-m]

# --- Mesh resolution [m] -----------------------------------------------------
fluid_dx      = 0.010
wall_dx       = 0.015
coil_dx       = 0.024
iron_dx       = 0.016
vac_dx        = 0.080   

# --- Output ------------------------------------------------------------------
mesh_file     = "slab_pump_mesh_native.h5"
show_plots    = True


# =============================================================================
# SETUP
# =============================================================================

tokamaker_python_path = os.getenv('OFT_ROOTPATH')
if tokamaker_python_path is not None:
    sys.path.append(os.path.join(tokamaker_python_path,'python'))

from OpenFUSIONToolkit.TokaMaker.meshing import gs_Domain
from OpenFUSIONToolkit.util import write_native_mesh

if not (0.0 < chan_half < wall_out <= coil_in < coil_out < yoke_out):
    raise SystemExit("need 0 < chan_half < wall_out <= coil_in < coil_out < yoke_out")
if air_x <= 0.0 or air_z <= 0.0:
    raise SystemExit("air_x and air_z must be positive")
x0 = yoke_out + air_x     # channel centre; the air box spans x = 0 .. 2*x0
if iron_len < pump_len:
    raise SystemExit("iron_len must be at least pump_len")

seg_pitch = pump_len / n_seg
seg_h = seg_pitch - seg_gap
if seg_h <= 0.0:
    raise SystemExit("seg_gap is larger than the segment pitch")
seg_zc = -pump_len / 2.0 + (np.arange(n_seg) + 0.5) * seg_pitch

z_half = max(chan_len, pump_len, iron_len) / 2.0
x_ext = 2.0 * x0
z_ext = z_half + air_z


# =============================================================================
# GEOMETRY
# =============================================================================

gs_mesh = gs_Domain(rextent=x_ext, zextents=[-z_ext, z_ext])

gs_mesh.define_region("air", vac_dx, "boundary")
gs_mesh.define_region("channel", fluid_dx, "plasma")
gs_mesh.define_region("wall", wall_dx, "conductor", eta=eta_wall)
for k in range(n_seg):
    gs_mesh.define_region(f"coil_{k:02d}", coil_dx, "coil")
gs_mesh.define_region("iron_yoke", iron_dx, "vacuum")

# Fluid
gs_mesh.add_rectangle(x0, 0.0, 2.0 * chan_half, chan_len, "channel", parent_name="air")

# Duct walls, one each side
for side in (-1.0, 1.0):
    gs_mesh.add_rectangle(x0 + side * 0.5 * (chan_half + wall_out), 0.0,
                          wall_out - chan_half, chan_len, "wall", parent_name="air")


def stator_contour(side):
    """Comb-shaped stator iron on one side: back iron plus the teeth between
    the coil slots, slots open toward the channel. Every coil corner is a
    vertex of this contour, so the shared edges are defined by the same points.
    side = +1 for x > x0, -1 for x < x0."""
    pts = [(coil_in, -iron_len / 2.0)]
    for zc in seg_zc:
        pts += [(coil_in,  zc - seg_h / 2.0), (coil_out, zc - seg_h / 2.0),
                (coil_out, zc + seg_h / 2.0), (coil_in,  zc + seg_h / 2.0)]
    pts += [(coil_in, iron_len / 2.0), (yoke_out, iron_len / 2.0),
            (yoke_out, -iron_len / 2.0)]
    pts = np.asarray(pts)
    pts[:, 0] = x0 + side * pts[:, 0]
    if side < 0:
        pts = pts[::-1]          # keep the same orientation as the +x side
    return pts


for side in (-1.0, 1.0):
    gs_mesh.add_polygon(stator_contour(side), "iron_yoke", parent_name="air")

# Winding: each coil region is the pair of slots at the same Z, one per side
for k, zc in enumerate(seg_zc):
    for side in (-1.0, 1.0):
        gs_mesh.add_rectangle(x0 + side * 0.5 * (coil_in + coil_out), zc,
                              coil_out - coil_in, seg_h,
                              f"coil_{k:02d}", parent_name="air")


# =============================================================================
# PLOT topology, BUILD, PLOT mesh
# =============================================================================

if show_plots:
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 1)
    for region in gs_mesh.regions:
        region.plot_segments(fig, ax)
    ax.set_aspect("equal", "box")
    ax.set_xlabel("x [m] (R in the mesher)")
    ax.set_ylabel("Z [m]")
    ax.set_title("slab pump topology")

mesh_pts, mesh_lc, mesh_reg = gs_mesh.build_mesh()

if show_plots:
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(2, 2, figsize=(8, 8), constrained_layout=True)
    gs_mesh.plot_mesh(fig, ax)
    plt.show()


# =============================================================================
# SAVE in OFT native format (build_mesh returns 0-based connectivity; the
# native reader wants 1-based)
# =============================================================================

write_native_mesh(mesh_file, mesh_pts, mesh_lc + 1, mesh_reg)


# =============================================================================
# REPORT region ids
# build_mesh re-indexes regions by TYPE (plasma, boundary, vacuum, conductor,
# coil), not by the order define_region was called in, so read them back.
# These must be consistent with the region IDs used in the Fortran driver
# (we should figure out a way to make this better in the future)
# =============================================================================

ids = {n: r["id"] for n, r in gs_mesh.region_info.items()}
print(f"\n{mesh_pts.shape[0]} points, {mesh_lc.shape[0]} cells, {max(ids.values())} regions")
for name, rid in sorted(ids.items(), key=lambda kv: kv[1]):
    if not name.startswith("COIL_"):
        print(f"  {rid:3d}  {name}")
print(f"  {ids['COIL_00']:3d} .. {ids[f'COIL_{n_seg-1:02d}']:<3d}  COIL_00 .. COIL_{n_seg-1:02d}")
print(f"\n&em_pump_options:  cyl=F  wall_region={ids['WALL']}  iron_region={ids['IRON_YOKE']}"
      f"  coil_reg0={ids['COIL_00']}  n_seg={n_seg}")
