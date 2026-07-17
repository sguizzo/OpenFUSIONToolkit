#------------------------------------------------------------------------------
# Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
#
# SPDX-License-Identifier: LGPL-3.0-only
#------------------------------------------------------------------------------
'''! Python interface for the coupled MUG/TokaMaker (blanket_td) time-dependent solver

An extension of @ref OpenFUSIONToolkit.TokaMaker._core.TokaMaker "TokaMaker" that couples
the MUG MHD solver to an existing TokaMaker equilibrium (see @ref MUGToksim).

@authors Sophia Guizzo
@date July 2026
@ingroup doxy_oft_python
'''
import ctypes
import numpy
from .._interface import *

## @cond

# mugtok_alloc(mugtok_ptr,equil_ptr,error_str)
mugtok_alloc = ctypes_subroutine(oftpy_lib.mugtok_alloc,
    [c_void_ptr_ptr, c_void_p, c_char_p])

# mugtok_setup(mugtok_ptr,dt,lin_tol,nl_tol,mhd_flag,dens,visc,nreg,incomp,toroidal_flow,error_str)
mugtok_setup = ctypes_subroutine(oftpy_lib.mugtok_setup,
    [c_void_p, c_double, c_double, c_double,
     ctypes_numpy_array(int32,1), ctypes_numpy_array(float64,1), ctypes_numpy_array(float64,1),
     c_int, c_bool, c_bool, c_char_p])

# mugtok_step(mugtok_ptr,curr,ncoils,time,dt,nl_its,lin_its,nretry,error_str)
mugtok_step = ctypes_subroutine(oftpy_lib.mugtok_step,
    [c_void_p, ctypes_numpy_array(float64,1), c_int,
     c_double_ptr, c_double_ptr, c_int_ptr, c_int_ptr, c_int_ptr, c_char_p])

# mugtok_field_ndofs(mugtok_ptr,field_id,n,error_str)
mugtok_field_ndofs = ctypes_subroutine(oftpy_lib.mugtok_field_ndofs,
    [c_void_p, c_int, c_int_ptr, c_char_p])

# mugtok_get_field(mugtok_ptr,field_id,vals,error_str)
mugtok_get_field = ctypes_subroutine(oftpy_lib.mugtok_get_field,
    [c_void_p, c_int, ctypes_numpy_array(float64,1), c_char_p])

# mugtok_get_pmesh(mugtok_ptr,np,r_loc,nc,lc_loc,error_str)
mugtok_get_pmesh = ctypes_subroutine(oftpy_lib.mugtok_get_pmesh,
    [c_void_p, c_int_ptr, c_double_ptr_ptr, c_int_ptr, c_int_ptr_ptr, c_char_p])

# mugtok_get_gradf(mugtok_ptr,gr,gz,error_str)
mugtok_get_gradf = ctypes_subroutine(oftpy_lib.mugtok_get_gradf,
    [c_void_p, ctypes_numpy_array(float64,1), ctypes_numpy_array(float64,1), c_char_p])

# mugtok_destroy(mugtok_ptr,error_str)
mugtok_destroy = ctypes_subroutine(oftpy_lib.mugtok_destroy,
    [c_void_p, c_char_p])

## @endcond


class MUGToksim():
    '''! Coupled MUG/TokaMaker (blanket_td) time-dependent simulation.

    Built on top of an existing @ref OpenFUSIONToolkit.TokaMaker._core.TokaMaker "TokaMaker"
    instance, whose mesh, finite-element representation, and active equilibrium are reused
    (borrowed, not copied). The poloidal flux (`psi`) and coil currents live in the TokaMaker
    equilibrium; the MHD fields (velocity, pressure, F) live in this object. Only incompressible
    flow is supported, so density is constant and is not exposed as an output.
    '''

    ## MUG-owned fields -> augmented-solution-vector block index.
    ## (density is omitted: only incompressible flow is supported, so it is constant.)
    _mug_fields = {'pressure': 5, 'F': 7}

    def __init__(self, tokamaker):
        '''! Create a combined MUG/TokaMaker simulation from an existing TokaMaker instance

        @param tokamaker A TokaMaker object with an active equilibrium
        '''
        if tokamaker._tMaker_equil is None:
            raise ValueError('TokaMaker object has no active equilibrium')
        ## Parent TokaMaker object (mesh/equilibrium/coils borrowed from it)
        self._tokamaker = tokamaker
        ## OpenFUSIONToolkit environment (See @ref OpenFUSIONToolkit._core.OFT_env "OFT_env")
        self._oft_env = tokamaker._oft_env
        ## C pointer to the Fortran-side combined simulation object
        self._ptr = ctypes.c_void_p()
        ## Cached order-1 (pressure) plotting mesh (loaded on first access)
        self._pr = None
        self._plc = None
        ## Per-region MHD flag [nreg] (1 = MHD region); set by `setup()`
        self._mhd_flag = None
        error_string = self._oft_env.get_c_errorbuff()
        mugtok_alloc(ctypes.byref(self._ptr), tokamaker._tMaker_equil.c_ptr, error_string)
        if error_string.value != b'':
            raise Exception(error_string.value)

    def __del__(self):
        '''! Free the Fortran-side object (the borrowed equilibrium/mesh are left intact)'''
        if getattr(self, '_ptr', None) is not None and self._ptr:
            error_string = self._oft_env.get_c_errorbuff()
            mugtok_destroy(self._ptr, error_string)
            self._ptr = ctypes.c_void_p()

    @property
    def r(self):
        '''! Node coordinates [np,3] for the main order fields (velocity, F).

        These are the same tessellated node points as the parent `TokaMaker.r`, so the
        order-2 fields (and `psi`) can be plotted with `self.r`/`self.lc`.
        '''
        return self._tokamaker.r

    @property
    def lc(self):
        '''! Triangle list [nc,3] for the main order fields (same as `TokaMaker.lc`)'''
        return self._tokamaker.lc

    @property
    def reg(self):
        '''! Per-cell region IDs [nc] for the main order fields (same as `TokaMaker.reg`).

        Aligned with `self.lc`, so the per-cell masks returned by `get_field` for
        `velocity`, `F`, and `current` index this triangle list.
        '''
        return self._tokamaker.reg

    def _load_pmesh(self):
        '''! Load and cache the order-1 (pressure) plotting mesh'''
        if self._pr is not None:
            return
        np_loc = ctypes.c_int()
        nc_loc = ctypes.c_int()
        r_loc = c_double_ptr()
        lc_loc = c_int_ptr()
        error_string = self._oft_env.get_c_errorbuff()
        mugtok_get_pmesh(self._ptr, ctypes.byref(np_loc), ctypes.byref(r_loc),
                         ctypes.byref(nc_loc), ctypes.byref(lc_loc), error_string)
        if error_string.value != b'':
            raise Exception(error_string.value)
        self._pr = numpy.ctypeslib.as_array(r_loc, shape=(np_loc.value, 3))
        self._plc = numpy.ctypeslib.as_array(lc_loc, shape=(nc_loc.value, 3))

    @property
    def pressure_r(self):
        '''! Node coordinates [np_p,3] for the (order-1) pressure field'''
        self._load_pmesh()
        return self._pr

    @property
    def pressure_lc(self):
        '''! Triangle list [nc_p,3] for the (order-1) pressure field'''
        self._load_pmesh()
        return self._plc

    @property
    def pressure_reg(self):
        '''! Per-cell region IDs [nc_p] for the (order-1) pressure mesh.

        The pressure mesh is tessellated at order 1, so its cells map one-to-one onto the
        base mesh cells. The main-order `self.reg` expands each base region across the `k`
        order-2 sub-cells in contiguous blocks (see `tokamaker_get_mesh`), so the base
        region per pressure cell is recovered by taking every `k`-th entry.
        '''
        self._load_pmesh()
        k = self.reg.shape[0] // self._plc.shape[0]
        return self.reg[::k]

    def _region_id(self, name):
        '''! Map a region key name to its 1-based region index (via the TokaMaker region dicts)'''
        for reg_dict in (self._tokamaker._cond_dict,
                         self._tokamaker._coil_dict,
                         self._tokamaker._vac_dict):
            if name in reg_dict:
                return reg_dict[name]['reg_id']
        raise KeyError('Unknown region "{0}"'.format(name))

    def setup(self, dt, lin_tol, nl_tol, mhd_regions, density, viscosity,
              incomp=True, allow_toroidal_flow=False):
        '''! Set up the coupled solver

        @param dt Timestep [s]
        @param lin_tol Linear solver tolerance
        @param nl_tol Non-linear solver tolerance
        @param mhd_regions List of region key names where the MHD (MUG) solve is active
        @param density List of mass densities [kg/m^3], aligned to `mhd_regions`
        @param viscosity List of dynamic viscosities [Pa-s], aligned to `mhd_regions`
        @param incomp Use incompressible flow (default: True). Compressible flow
                      (`incomp=False`) is not currently supported and will raise.
        @param allow_toroidal_flow Allow toroidal (phi) flow to evolve in the MHD regions
                                   (default: False, i.e. toroidal velocity is clamped)
        '''
        if not (len(mhd_regions) == len(density) == len(viscosity)):
            raise ValueError('"mhd_regions", "density", and "viscosity" must have equal length')
        nreg = self._tokamaker.nregs
        if nreg <= 0:
            raise ValueError('TokaMaker regions are not set up (nregs <= 0)')
        mhd_flag = numpy.zeros((nreg,), dtype=numpy.int32)
        dens_reg = -numpy.ones((nreg,), dtype=numpy.float64)
        visc_reg = -numpy.ones((nreg,), dtype=numpy.float64)
        for name, dens, visc in zip(mhd_regions, density, viscosity):
            i = self._region_id(name) - 1
            mhd_flag[i] = 1
            dens_reg[i] = dens
            visc_reg[i] = visc/dens  # dynamic [Pa-s] -> kinematic [m^2/s]
        error_string = self._oft_env.get_c_errorbuff()
        mugtok_setup(self._ptr, ctypes.c_double(dt), ctypes.c_double(lin_tol), ctypes.c_double(nl_tol),
                     numpy.ascontiguousarray(mhd_flag, dtype=numpy.int32),
                     numpy.ascontiguousarray(dens_reg, dtype=numpy.float64),
                     numpy.ascontiguousarray(visc_reg, dtype=numpy.float64),
                     ctypes.c_int(nreg), ctypes.c_bool(incomp),
                     ctypes.c_bool(allow_toroidal_flow), error_string)
        if error_string.value != b'':
            raise Exception(error_string.value)
        # Remember which regions are MHD so `get_field` can build per-cell field masks
        self._mhd_flag = mhd_flag

    def step(self, time, dt, coil_currents=None):
        '''! Advance the coupled solution by one timestep

        @param time Time at the start of the step [s]
        @param dt Timestep size [s]
        @param coil_currents Coil currents at the end of the step as a dict of name->[A]
                             (if `None`, the equilibrium's present currents are held)
        @result new time, new dt, # nonlinear iterations, # linear iterations, # retries
        '''
        c_time = ctypes.c_double(time)
        c_dt = ctypes.c_double(dt)
        nl_its = ctypes.c_int()
        lin_its = ctypes.c_int()
        nretry = ctypes.c_int()
        if coil_currents is None:
            coil_currents, _ = self._tokamaker.get_coil_currents()
        currents = numpy.ascontiguousarray(self._tokamaker.coil_dict2vec(coil_currents),
                                           dtype=numpy.float64)
        error_string = self._oft_env.get_c_errorbuff()
        mugtok_step(self._ptr, currents, ctypes.c_int(self._tokamaker.ncoils),
                    ctypes.byref(c_time), ctypes.byref(c_dt),
                    ctypes.byref(nl_its), ctypes.byref(lin_its), ctypes.byref(nretry), error_string)
        if error_string.value != b'':
            raise Exception(error_string.value)
        return c_time.value, c_dt.value, nl_its.value, lin_its.value, nretry.value

    def _get_mug_field(self, field_id):
        '''! Fetch a MUG-owned field (by augmented-vector block index) at node points'''
        n = ctypes.c_int()
        error_string = self._oft_env.get_c_errorbuff()
        mugtok_field_ndofs(self._ptr, ctypes.c_int(field_id), ctypes.byref(n), error_string)
        if error_string.value != b'':
            raise Exception(error_string.value)
        vals = numpy.zeros((n.value,), dtype=numpy.float64)
        mugtok_get_field(self._ptr, ctypes.c_int(field_id), vals, error_string)
        if error_string.value != b'':
            raise Exception(error_string.value)
        return vals

    def _get_gradf(self):
        '''! Fetch the L2-projected poloidal gradient of F, (dF/dR, dF/dZ), at node points'''
        n = ctypes.c_int()
        error_string = self._oft_env.get_c_errorbuff()
        mugtok_field_ndofs(self._ptr, ctypes.c_int(7), ctypes.byref(n), error_string)
        if error_string.value != b'':
            raise Exception(error_string.value)
        gr = numpy.zeros((n.value,), dtype=numpy.float64)
        gz = numpy.zeros((n.value,), dtype=numpy.float64)
        mugtok_get_gradf(self._ptr, gr, gz, error_string)
        if error_string.value != b'':
            raise Exception(error_string.value)
        return gr, gz

    def _mhd_cell_mask(self, reg):
        '''! Per-cell boolean mask selecting the MHD regions of a per-cell region array

        @param reg Per-cell region-ID array (e.g. `self.reg` or `self.pressure_reg`)
        @result Boolean mask, `True` where the cell's region is an MHD (MUG) region
        '''
        if self._mhd_flag is None:
            raise RuntimeError('MHD regions are not defined; call `setup()` first')
        mhd_ids = numpy.flatnonzero(self._mhd_flag) + 1  # region IDs are 1-based
        return numpy.isin(reg, mhd_ids)

    def get_field(self, name):
        '''! Get the present values of a field at node points (aligned with `TokaMaker.get_mesh()`)

        @param name One of 'velocity', 'pressure', 'F' (MUG-owned),
                    'psi', 'coils' (read from the TokaMaker equilibrium), or 'current'
        @result For 'psi'/'coils', a plain NumPy array. For the MHD fields ('velocity',
                'pressure', 'F') and 'current', a tuple `(mask, field)` where `mask` is a
                per-cell boolean mask to apply at plot time (see below). 'velocity' and
                'current' fields are (3, ndof) arrays ordered (R, phi, Z).

        The per-cell masks let each field be plotted only where it is physically defined,
        exactly as for 'current' (e.g. `tripcolor(r, z, self.lc[mask], field[i, :])`):
          - 'velocity', 'pressure': `True` only in the MHD (MUG) regions.
          - 'F': `True` everywhere (F is defined over the whole mesh).
          - 'current': the conductor mask from the TokaMaker equilibrium.
        The 'velocity', 'F', and 'current' masks index `self.lc`/`self.reg`; the 'pressure'
        mask indexes the (order-1) `self.pressure_lc`/`self.pressure_reg`.

        Density is not exposed: only incompressible flow is supported, so it is constant.
        '''
        if name == 'psi':
            return self._tokamaker.get_psi(normalized=False)
        if name == 'coils':
            return self._tokamaker.get_coil_currents()[0]
        if name == 'velocity':
            vel = numpy.stack([self._get_mug_field(2),   # R
                               self._get_mug_field(3),   # phi
                               self._get_mug_field(4)],  # Z
                              axis=0)
            return self._mhd_cell_mask(self.reg), vel
        if name == 'pressure':
            # Pressure is only defined up to a constant; shift so the minimum is 0
            p = self._get_mug_field(self._mug_fields['pressure'])
            return self._mhd_cell_mask(self.pressure_reg), p - p.min()
        if name == 'F':
            # F is defined over the whole mesh, so its mask selects every cell
            F = self._get_mug_field(self._mug_fields['F'])
            return numpy.ones((self.lc.shape[0],), dtype=bool), F
        if name == 'current':
            return self._get_current()
        raise ValueError('Unknown field "{0}"; expected one of '
                         'velocity, pressure, F, psi, coils, current'.format(name))

    def _get_current(self):
        '''! Current density in the conducting structures, J [A/m^2].

        Poloidal components come from the (projected) grad(F): J_pol = (grad F x phi_hat)/(mu0 R);
        the toroidal component and the conductor mask come from the TokaMaker equilibrium via
        `calc_conductor_currents(psi)`. The mask is handled exactly as in `calc_conductor_currents`:
        a per-cell mask is returned alongside the field, to be applied at plot time (e.g.
        `tripcolor(r, z, self.lc[mask], J[i, :])`).

        @result Tuple `(mask, J)` where `mask` is the per-cell conductor mask and `J` is a
                (3, ndof) array ordered (R, phi, Z).
        '''
        # Poloidal current from projected grad(F), with 1/R guarded at the axis
        gr, gz = self._get_gradf()          # dF/dR, dF/dZ
        R = self.r[:, 0]
        JR = numpy.zeros_like(gr)
        JZ = numpy.zeros_like(gr)
        m = R > 0.0
        JR[m] = -gz[m] / (mu0 * R[m])       # -(1/(mu0 R)) dF/dZ
        JZ[m] =  gr[m] / (mu0 * R[m])       #  (1/(mu0 R)) dF/dR
        # Toroidal current + conductor mask from the TokaMaker equilibrium
        psi = self._tokamaker.get_psi(normalized=False)
        mask, Jphi = self._tokamaker._tMaker_equil.calc_conductor_currents(psi)
        J = numpy.stack([JR, Jphi, JZ], axis=0)   # (3, ndof): (R, phi, Z)
        return mask, J
