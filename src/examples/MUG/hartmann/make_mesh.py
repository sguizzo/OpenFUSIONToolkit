#!/usr/bin/env python3
"""Generate a graded rectangular duct mesh for the Hartmann example.

The cross-section is a rectangle, so the mesh is a tensor product and can be
built exactly here -- no external mesher needed. What matters is the *grading*:
at Hartmann number Ha the flow has two very different boundary layers,

    Hartmann layers, on the walls normal to B : delta/a = 1/Ha
    side (Shercliff) layers, on walls along B : delta/b = 1/sqrt(Ha)

At Ha=80 those differ by a factor of 9, so a uniform mesh fine enough for the
Hartmann layer wastes an enormous number of cells everywhere else. Here each
direction is graded geometrically from the walls inward, with the first cell
sized to a requested fraction of the local layer thickness.

Periodicity (optional): OFT's native reader pairs a listed node with the
boundary node having the same coordinates except in x, and requires that partner
to sit at x = 0 (see native_bset_periodic in mesh_native.F90, ref_index=1 for a
2D surface mesh). So a periodic run must (a) be periodic in x, (b) have one
periodic boundary exactly at x = 0, and (c) list the nodes on the far side.
Note the Hartmann duct itself needs no periodic direction -- all four sides are
walls -- so this is only for channel-type variants.
"""
import numpy as np, h5py


def write_native_mesh(filename, r, lc, reg, periodic_info=None):
    """Write OFT's native HDF5 mesh format.

    Same schema as OpenFUSIONToolkit.util.write_native_mesh, inlined so this
    script does not need a compiled OFT (importing the package pulls in
    liboftpy). Note LC is 1-based, and the field names are upper case.
    """
    # OFT decides a surface mesh is 2D from the SHAPE of mesh/R: native_load_smesh
    # sets is_2d = (dim_sizes(1)==2). Only then is smesh%dim=2, which in turn makes
    # native_bset_periodic use ref_index=1 (periodic in x). Written with 3 columns
    # the mesh is treated as a surface in 3D and periodicity is sought in z, which
    # can never match for a planar mesh -- it silently pairs nothing.
    r = np.asarray(r)[:, :2]
    print("Saving mesh: {0}".format(filename))
    with h5py.File(filename, 'w') as h5:
        h5.create_dataset('mesh/R', data=r, dtype='f8')
        h5.create_dataset('mesh/LC', data=lc, dtype='i4')
        h5.create_dataset('mesh/REG', data=reg, dtype='i4')
        if periodic_info is not None:
            h5.create_dataset('mesh/periodicity/nodes', data=periodic_info, dtype='i4')


def graded_half(H, n, s1):
    """n cell edges spanning [0,H], geometric, first cell of size s1 (at 0)."""
    if abs(s1*n - H) < 1e-14*H:
        return np.linspace(0.0, H, n+1), 1.0
    lo, hi = 1.0+1e-12, 4.0
    for _ in range(200):                      # bisect on the growth ratio
        r = 0.5*(lo+hi)
        first = H*(r-1.0)/(r**n - 1.0)
        if first > s1: lo = r
        else:          hi = r
    r = 0.5*(lo+hi)
    s = H*(r-1.0)/(r**n - 1.0)
    edges = np.concatenate(([0.0], np.cumsum(s*r**np.arange(n))))
    return edges*(H/edges[-1]), r


def two_sided(half, n_half, s1):
    """Symmetric graded nodes on [-half, half], clustered at both walls."""
    e, r = graded_half(half, n_half, s1)
    f = half - e[::-1]          # flip so the fine cells sit at the wall, not the centre
    return np.concatenate((-f[::-1], f[1:])), r


def build(a=0.05, b=0.05, Ha=80.0, n_a=40, n_b=20, frac=0.25,
          periodic_x=False, filename='hartmann_mesh.h5'):
    d_hart = a/Ha                 # Hartmann layer (walls normal to B, y=+-a)
    d_side = b/np.sqrt(Ha)        # side layer     (walls along  B, x=+-b)
    y, r_y = two_sided(a, n_a, frac*d_hart)
    if periodic_x:
        x = np.linspace(0.0, 2.0*b, 2*n_b+1)   # uniform; x=0 is the partner side
        r_x = 1.0
    else:
        x, r_x = two_sided(b, n_b, frac*d_side)

    nx, ny = len(x), len(y)
    X, Y = np.meshgrid(x, y, indexing='ij')
    r = np.column_stack([X.ravel(), Y.ravel(), np.zeros(nx*ny)])
    idx = lambda i, j: i*ny + j                       # 0-based

    tris = []
    for i in range(nx-1):
        for j in range(ny-1):
            p00, p10, p11, p01 = idx(i,j), idx(i+1,j), idx(i+1,j+1), idx(i,j+1)
            if (i+j) % 2 == 0:                        # alternate the diagonal
                tris += [[p00,p10,p11],[p00,p11,p01]]
            else:
                tris += [[p00,p10,p01],[p10,p11,p01]]
    lc = np.array(tris, dtype=np.int32) + 1           # 1-based for OFT
    reg = np.ones(lc.shape[0], dtype=np.int32)

    per = None
    if periodic_x:                                    # nodes on the far (x=2b) side
        per = np.array([idx(nx-1, j)+1 for j in range(ny)], dtype=np.int32)

    write_native_mesh(filename, r, lc, reg, periodic_info=per)

    # --- report what the grading actually achieves ---
    dy0, dx0 = y[1]-y[0], x[1]-x[0]
    print(f"\n  Ha = {Ha}")
    print(f"  Hartmann layer  delta_y = {d_hart:.5e}   first cell dy = {dy0:.5e}"
          f"   ({d_hart/dy0:.1f} cells, {2*d_hart/dy0:.1f} P2 nodes)")
    print(f"  side layer      delta_x = {d_side:.5e}   first cell dx = {dx0:.5e}"
          f"   ({d_side/dx0:.1f} cells, {2*d_side/dx0:.1f} P2 nodes)")
    print(f"  growth ratio    y: {r_y:.4f}   x: {r_x:.4f}"
          f"   (keep <~1.15 for accuracy)")
    print(f"  grid {nx} x {ny} nodes, {lc.shape[0]} triangles")
    uniform = int(np.ceil(2*a/dy0))
    print(f"  a uniform mesh at this wall spacing would need {uniform} cells in y"
          f" (vs {2*n_a}) -> {uniform/(2*n_a):.1f}x saving")
    if per is not None:
        print(f"  periodic in x: {len(per)} node pairs, partner side at x=0")
    return r, lc


if __name__ == '__main__':
    import argparse
    p = argparse.ArgumentParser()
    p.add_argument('--a', type=float, default=0.05)
    p.add_argument('--b', type=float, default=0.05)
    p.add_argument('--Ha', type=float, default=80.0)
    p.add_argument('--na', type=int, default=40, help='cells per half-width in y')
    p.add_argument('--nb', type=int, default=20, help='cells per half-width in x')
    p.add_argument('--frac', type=float, default=0.25,
                   help='first cell as a fraction of the layer thickness')
    p.add_argument('--periodic-x', action='store_true')
    p.add_argument('-o', '--out', default='hartmann_mesh.h5')
    A = p.parse_args()
    build(A.a, A.b, A.Ha, A.na, A.nb, A.frac, A.periodic_x, A.out)
