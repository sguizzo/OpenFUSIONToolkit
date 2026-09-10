#!/usr/bin/env python3
"""Build the cylindrical-duct induction pump mesh for the MUG em_pump driver.

Edit the INPUTS block below and run:  python3 make_pump_mesh.py

"""

import os
import sys

import numpy as np


# =============================================================================
# INPUTS
# =============================================================================

# --- Radii [m], inside to out ------------------------------------------------
pipe_rout     = 0.040   # fluid radius, i.e. the duct bore (80 mm diameter)
wall_rout     = 0.045   # duct wall outer radius
coil_rin      = 0.055   # winding inner radius
coil_rout     = 0.105   # winding outer radius
yoke_rout     = 0.155   # stator outer radius; the iron runs in from here to coil_rin

# --- Axial extents [m] -------------------------------------------------------
chan_len      = 0.72    # full duct length; runs past the pump section at both ends
                        # so the outflow boundary is clear of the driven region
pump_len      = 0.72    # axial length of the coil stack, i.e. the active pump section
iron_len      = 0.72    # axial length of the stator; >= pump_len, the excess becoming end teeth

# --- Winding -----------------------------------------------------------------
n_seg         = 24      # number of coil segments
seg_gap       = 0.010   # insulating gap between adjacent segments; this IS the stator tooth width
                        # (tooth = seg_pitch - seg_h).  5 mm teeth saturate to mu_r ~ 10 at I0 = 4e4

# --- Materials ---------------------------------------------------------------
eta_wall      = 1.2e-6  # duct wall resistivity [Ohm-m]; 316 stainless-ish

# --- Mesh resolution [m] -----------------------------------------------------
fluid_dx      = 0.010   # channel
wall_dx       = 0.015   # duct wall
coil_dx       = 0.024   # coils
iron_dx       = 0.016   # iron; must resolve the tooth, so keep it below seg_gap
vac_dx        = 0.040   # vacuum

# --- Output ------------------------------------------------------------------
mesh_file     = "pump_mesh_native.h5"   # OFT native format, what the Fortran driver reads
show_plots    = True                    # draw the topology and the finished mesh


# =============================================================================
# SETUP
# =============================================================================

tokamaker_python_path = os.getenv('OFT_ROOTPATH')
if tokamaker_python_path is not None:
    sys.path.append(os.path.join(tokamaker_python_path,'python'))

from OpenFUSIONToolkit.TokaMaker.meshing import gs_Domain
from OpenFUSIONToolkit.util import write_native_mesh

if not (0.0 < pipe_rout < wall_rout <= coil_rin < coil_rout < yoke_rout):
    raise SystemExit(
        "radii must increase: 0 < pipe_rout < wall_rout <= coil_rin < coil_rout < yoke_rout"
    )
if iron_len < pump_len:
    raise SystemExit("iron_len must be at least pump_len, or the end coils fall outside the iron")

seg_pitch = pump_len / n_seg
seg_h = seg_pitch - seg_gap
if seg_h <= 0.0:
    raise SystemExit("seg_gap is larger than the segment pitch; reduce it or use fewer segments")
seg_zc = -pump_len / 2.0 + (np.arange(n_seg) + 0.5) * seg_pitch


# =============================================================================
# GEOMETRY
# =============================================================================
# Regions are declared first, then drawn.  Each piece is a rectangle in (R,Z)
# given as (centre_R, centre_Z, width, height), all as children of "air".

gs_mesh = gs_Domain()

gs_mesh.define_region("air", vac_dx, "boundary")
gs_mesh.define_region("channel", fluid_dx, "plasma")
gs_mesh.define_region("wall", wall_dx, "conductor", eta=eta_wall)
for k in range(n_seg):
    gs_mesh.define_region(f"coil_{k:02d}", coil_dx, "coil")
gs_mesh.define_region("iron_yoke", iron_dx, "vacuum")

# Fluid.  A solid cylinder on the axis: R in [0, pipe_rout], so its centre is
# pipe_rout/2 and its width is pipe_rout.
gs_mesh.add_rectangle(pipe_rout / 2.0, 0.0, pipe_rout, chan_len,
                      "channel", parent_name="air")
# Duct wall
gs_mesh.add_rectangle(0.5 * (pipe_rout + wall_rout), 0.0,
                      wall_rout - pipe_rout, chan_len,
                      "wall", parent_name="air")
# Stator iron: one comb-shaped region, back iron plus the teeth between coils,
# with the coil slots cut out and left open to the air gap. Drawn as an explicit
# contour rather than a rectangle so that every coil corner is already a vertex
# of the iron boundary -- shared edges must be defined by the same points, and
# the mesher will not insert them itself.
def stator_contour():
    pts = [(coil_rin, -iron_len / 2.0)]
    for zc in seg_zc:                       # up the tooth faces, out around each slot
        pts += [(coil_rin,  zc - seg_h / 2.0), (coil_rout, zc - seg_h / 2.0),
                (coil_rout, zc + seg_h / 2.0), (coil_rin,  zc + seg_h / 2.0)]
    pts += [(coil_rin, iron_len / 2.0), (yoke_rout, iron_len / 2.0),
            (yoke_rout, -iron_len / 2.0)]
    return np.asarray(pts)

gs_mesh.add_polygon(stator_contour(), "iron_yoke", parent_name="air")
# Winding, one region per axial segment, seated in the slots
for k, zc in enumerate(seg_zc):
    gs_mesh.add_rectangle(0.5 * (coil_rin + coil_rout), zc,
                          coil_rout - coil_rin, seg_h,
                          f"coil_{k:02d}", parent_name="air")


# =============================================================================
# PLOT: topology, then the finished mesh
# =============================================================================

if show_plots:
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 1)
    for region in gs_mesh.regions:
        region.plot_segments(fig, ax)
    ax.set_aspect("equal", "box")
    ax.set_xlabel("R [m]")
    ax.set_ylabel("Z [m]")
    ax.set_title("pump topology")

mesh_pts, mesh_lc, mesh_reg = gs_mesh.build_mesh()

if show_plots:
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(2, 2, figsize=(8, 8), constrained_layout=True)
    gs_mesh.plot_mesh(fig, ax)

if show_plots:
    import matplotlib.pyplot as plt
    plt.show()

# =============================================================================
# SAVE in OFT native format
# build_mesh returns 0-based connectivity; the native reader wants 1-based.
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
print(f"\n&em_pump_options:  wall_region={ids['WALL']}  iron_region={ids['IRON_YOKE']}"
      f"  coil_reg0={ids['COIL_00']}  n_seg={n_seg}  pump_len={pump_len}  chan_zmax={chan_len/2}")
