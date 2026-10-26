! sparc_vde_tq2: continue a sparc_vde_beta_kick run from one of its restart files through a second
! thermal quench and the current quench, as M3D-C1 does at 78 ms (kappa_perp 1e-6 -> 0.2 normalized,
! i.e. 1.542e20 -> 3.08e25 1/(m s), uniform, kappa_par unchanged; the thermal energy is gone in ~7 us).
! Same setup and options as sparc_vde_beta_kick (group &sparc_vde_options, read the same way), plus:
!   kappa_tq2    kappa_perp from the start of the second quench on [1/(m s)] (default 3.084e25)
!   tq2_t_start  start time [s] of the second quench (default <= 0: the restart time)
!   tq2_dur      length of the short-step quench phase [s] (default 2e-5)
!   dt_tq2       timestep during the quench phase [s] (default 5e-7)
!   dt_cq        largest timestep in the current-quench phase [s] (default 2.5e-4); dt starts at the
!                quench dt and doubles after every dt_grow_steps clean steps (default here 3) up to dt_cq
!   mass_scale_tq2  ion mass enhancement during the quench phase (default <= 0: mass_scale_post)
!   mass_scale_cq   ion mass enhancement during the current-quench phase (default <= 0: mass_scale_post)
!                The quench response is Alfvenic, so a smaller mass there follows M3D-C1 more closely;
!                visc_post stays the kinematic nu, so the viscous decay rate does not change with the mass
! Viscosity, mass scale, SOL viscosity, tolerances etc. are the post-kick values (visc_post, mass_scale_post, ...).
! Needs rst_start > 0. Phases: 3 = second thermal quench, 4 = current quench (kick_phase in the log).
#include "local.h"
PROGRAM sparc_vde_tq2
!---Runtime
USE oft_base
!---Grid
USE multigrid, ONLY: multigrid_mesh
USE multigrid_build, ONLY: multigrid_construct_surf
USE oft_mesh_native, ONLY: mesh_filename => filename
!
USE oft_la_base, ONLY: oft_vector, oft_matrix
USE oft_solver_base, ONLY: oft_solver
USE oft_solver_utils, ONLY: create_cg_solver, create_diag_pre
!
USE oft_blag_operators, ONLY: oft_blag_zerob, oft_blag_getmop, oft_blag_project
USE oft_scalar_inits, ONLY: poss_scalar_bfield
USE mhd_utils, ONLY: elec_charge, proton_mass, mu0
USE oft_io, ONLY: hdf5_field_get_sizes, hdf5_read, hdf5_field_exist, hdf5_get_group_members
USE oft_gs, ONLY: flux_func, gs_equil, gs_update_bounds, gs_test_bounds, compute_bcmat, gs_setup_walls, gs_get_qprof, gs_factory
USE oft_gs_td, ONLY: oft_tmaker_td
USE oft_gs_util, ONLY: gs_profile_load
USE oft_lag_basis, ONLY: oft_lag_setup,oft_scalar_bfem, oft_blag_eval, oft_blag_geval, oft_2D_lagrange_cast
USE fem_base, ONLY: oft_ml_fem_type
USE mugtok_td
USE oft_lu, ONLY: oft_lusolver
USE oft_native_solvers, ONLY: oft_bjprecond

IMPLICIT NONE
INTEGER(i4) :: io_unit,ierr, i, j
TYPE(multigrid_mesh) :: mg_mesh
TYPE(oft_ml_fem_type), TARGET :: ML_blagrange
CLASS(oft_scalar_bfem), POINTER :: blagrange
TYPE(oft_blag_zerob), TARGET :: blag_zerob ! setting boundary vals to zero
!---Mass matrix solver
TYPE(poss_scalar_bfield) :: field_init
CLASS(oft_solver), POINTER :: minv => NULL()
CLASS(oft_solver), POINTER :: minv_2 => NULL()
CLASS(oft_matrix), POINTER :: mop => NULL()
CLASS(oft_matrix), POINTER :: mop_2 => NULL()
CLASS(oft_vector), POINTER :: u,v, u_2, v_2
TYPE(gs_equil), TARGET :: equil
TYPE(gs_factory), TARGET :: machine
TYPE(oft_tmaker_td) :: tokamaker
REAL(r8), POINTER, DIMENSION(:) :: tmp_arr
!---Runtime options
INTEGER(i4) :: order = 3
INTEGER(i4) :: nsteps = 1
INTEGER(i4) :: rst_freq = 1
INTEGER(i4) :: ndims, nl_its, l_its, nretry
INTEGER(i4) :: npoints
integer(i4), allocatable, dimension(:) :: dim_sizes
INTEGER(i4), POINTER, DIMENSION(:) :: cell_dofs
REAL(r8) :: dt = 0.04336664911469267
REAL(r8) :: t = 0.d0
INTEGER(i4) :: kick_phase = -1, kick_phase_prev = -1 !< Current/previous kappa-kick phase (0=pre-kick, 1=kick, 2=post-kick); used to apply the phase's target dt/kappa only on entry, so a dt reduced by step() due to a solver-failure retry is not clobbered on later iterations of the same phase
REAL (r8):: ip_ratio_target = 0.205
REAL (r8):: ip_target = 7.87E6
REAL(r8), allocatable, dimension(:) :: psi_eq, psi_pert, psi_total, curr_reg, dens_reg, kappa_reg,visc_reg, D_reg, voltages
!---Coil/conductor region info read from the equilibrium file's coil_dict/cond_dict
CHARACTER(LEN=80), ALLOCATABLE, DIMENSION(:) :: coil_names, cond_names
INTEGER(i4), ALLOCATABLE, DIMENSION(:) :: coil_nturns_file
INTEGER(i4) :: reg_id_tmp, coil_id_tmp, nturns_tmp, cond_id_tmp
REAL(r8) :: eta_tmp
REAL(r8), allocatable, dimension(:,:) :: eta_reg !< Resistivity [in-plane (poloidal), out-of-plane (toroidal)] by region, in units of mu0
REAL(r8) :: eta_plate_pol = 6.93d-7 !< Poloidal (in-plane) resistivity [Ohm*m] of the PLATE_* conductors; their toroidal resistivity is left at the file value
REAL(r8), allocatable, dimension(:) :: coil_currs_fixed !< Fixed coil currents passed to every step (step() overwrites equil%coil_currs internally, so a separate reference copy is needed)
INTEGER(i4) :: nstep_total = 10000 !< Number of time steps (dt*nstep should approach the wall time ~ms)
REAL(r8) :: t_kick_start = 1.d-6 !< Time at which the kappa_perp kick begins (let the system settle first) [s]
REAL(r8) :: t_kick_dur = 6.5d-4 !< Duration of the kappa_perp kick [s] (0.3 ms)
REAL(r8) :: kappa_kick = 1.542d22 !< DUMMY: perpendicular thermal conductivity kappa=chi*n [1/(m s)] in the plasma during the kick; tune to taste
REAL(r8) :: kappa_base = 1.542d20 !< DUMMY: baseline (pre-/post-kick) perpendicular thermal conductivity kappa=chi*n [1/(m s)] in the plasma (sets kappa_reg(1) below)
REAL(r8) :: visc_kick = 500.d0 !< Kinematic viscosity [m^2/s] in the plasma during the pre-kick and kick phases (sets visc_reg(1) below)
REAL(r8) :: visc_post = 200.d0 !< Kinematic viscosity [m^2/s] in the plasma after the kick
REAL(r8) :: visc_bulk_ratio = 2.d0 !< Ratio of bulk (compressional, grad(div v)) to shear viscosity (mug%nu2_ratio; setup_mhd's default is 2, 0 turns bulk viscosity off). visc_kick/visc_post are the TOTAL viscosity, so shear = visc/(1+visc_bulk_ratio) and bulk = shear*visc_bulk_ratio
REAL(r8) :: visc_sol_mult = 1.d0 !< Post-kick only: multiply the viscosity in cold (SOL) plasma by up to visc_sol_mult, nu -> nu*(1+(visc_sol_mult-1)*(1-tanh((T-visc_sol_T0)/visc_sol_dT))/2), leaving the hot core at visc_post; 1 = off (see xmhd_2d visc_sol_mult)
REAL(r8) :: visc_sol_T0 = 50.d0 !< Middle of the visc_sol_mult transition [eV]
REAL(r8) :: visc_sol_dT = 10.d0 !< Width of the visc_sol_mult transition [eV]
LOGICAL :: visc_sol_on_psi = .FALSE. !< Key the visc_sol_mult boost on normalized flux instead of T: boost where psiN > visc_sol_psiN0 (tanh width visc_sol_dpsiN), i.e. outside the LCFS (see xmhd_2d visc_sol_on_psi)
REAL(r8) :: visc_sol_psiN0 = 1.d0 !< Middle of the psi-keyed transition (normalized flux)
REAL(r8) :: visc_sol_dpsiN = 0.002d0 !< Width of the psi-keyed transition (normalized flux)
LOGICAL :: fix_F_boundary = .TRUE. !< Hold F = R B_phi at its initial value on the outer mesh boundary (fixed TF current, as M3D-C1's default iconst_bz=1; default); .FALSE. = the natural condition that conserves the domain's toroidal flux (runs before 2026-10-26)
LOGICAL :: ohmic_heating = .FALSE. !< Ohmic heating eta*J^2 in the temperature equation (xmhd_2d option ohmic; current frozen once per step, explicit source). Default off (runs before 2026-10-26 have none)
LOGICAL :: lu_reanalyze = .FALSE. !< Redo the LU symbolic analysis at every preconditioner rebuild (mugtok_td lu_reanalyze; default: only at the first factorization). .TRUE. reproduces runs before 2026-10-26
INTEGER(i4) :: pre_nlocal = 0 !< If > 0: spatial block-Jacobi preconditioner with pre_nlocal METIS subdomains, one LU (package lu_package) each, instead of the full LU (mugtok_td pre_nlocal); 0 = full LU
LOGICAL :: pre_overlap = .FALSE. !< With pre_nlocal > 0: overlap the subdomains by one layer of unknowns (additive Schwarz)
INTEGER(i4) :: pre_part = 0 !< With pre_nlocal > 0: 0 = METIS partition; 1 = geometric 2-way split, plasma region vs the rest (all dense boundary/coil couplings in one subdomain)
REAL(r8) :: drag_rate = 0.d0 !< Post-kick linear drag rate [1/s] on the poloidal velocity (force -m_i*n*drag_rate*v_pol); 0 = off
LOGICAL :: drop_convection = .FALSE. !< If true, drop the (v.grad)v convective term from the momentum equation (setup_mhd's default); this driver keeps it by default
LOGICAL :: visc_dynamic = .FALSE. !< If true, the dynamic viscosity m_i*n*nu is uniform in the plasma (as with M3D-C1's amu/amuc): visc_kick/visc_post are then the kinematic viscosity at density visc_nref, and nu = visc*visc_nref/n pointwise
REAL(r8) :: visc_nref = 1.d20 !< Reference density [m^-3] for visc_dynamic (M3D-C1: n0_norm)
REAL(r8) :: mass_scale_kick = 5000.d0 !< Plasma (region 1) ion mass enhancement factor during the pre-kick and kick phases (passed to setup_mhd)
REAL(r8) :: mass_scale_post = 5000.d0 !< Plasma (region 1) ion mass enhancement factor after the kick
REAL(r8) :: m_i_plasma = 2.d0*proton_mass !< Physical (un-enhanced) plasma ion mass [kg]
INTEGER(i4) :: rst_start = 0 !< If >0, resume from restart file 'mugtok_NNNNN.rst' with NNNNN=rst_start (set in oft.in group sparc_vde_options)
REAL(r8) :: dt_rst = -1.d0 !< Timestep read from the restart file when resuming
REAL(r8) :: dt_min_frac = 0.125d0 !< Stop the run (then plot) once step() retries have cut dt below dt_min_frac times the current phase's starting dt; <=0 disables (set in oft.in group sparc_vde_options)
REAL(r8) :: dt_phase = 0.d0 !< Starting (nominal) dt of the current kick phase, the reference for dt_min_frac
REAL(r8) :: kappa_tq2 = 3.084d25 !< kappa_perp [1/(m s)] from the second thermal quench on (M3D-C1: 0.2 normalized)
REAL(r8) :: tq2_t_start = -1.d0 !< Start of the second thermal quench [s]; <= 0: the restart time
REAL(r8) :: tq2_dur = 2.d-5 !< Length of the short-step quench phase [s]
REAL(r8) :: dt_tq2 = 5.d-7 !< Timestep during the quench phase [s]
REAL(r8) :: dt_cq = 2.5d-4 !< Largest timestep in the current-quench phase [s]
REAL(r8) :: mass_scale_tq2 = -1.d0 !< Ion mass enhancement in the quench phase (3); <= 0: mass_scale_post
REAL(r8) :: mass_scale_cq = -1.d0 !< Ion mass enhancement in the current-quench phase (4); <= 0: mass_scale_post
REAL(r8) :: dt_stop_ref = 0.d0 !< dt at entry to the current phase: reference for the dt_min_frac stop
REAL(r8) :: dt_post = 5.d-4 !< Nominal post-kick timestep [s] (kick_phase 2); also the cap for dt_grow_steps and the reference for dt_min_frac. A restart still starts at the restart file's dt and reaches dt_post only through dt_grow_steps
INTEGER(i4) :: dt_grow_steps = 0 !< If > 0: after this many consecutive steps with no retry, double dt again (capped at the phase's nominal dt), so a dt halved by step() retries can recover; 0 = never (dt stays reduced)
LOGICAL :: newton_backtrack = .FALSE. !< Turn on backtracking (step halving while |F| grows) in the Newton solver; mugtok_td switches it off by default
LOGICAL :: resid_dump = .FALSE. !< Write the Newton residual of failed Newton attempts to the XDMF group mugtok_resid (oft_mugtok_td%resid_dump; at most 8)
INTEGER(i4) :: n_clean = 0 !< Consecutive steps without a retry (for dt_grow_steps)
LOGICAL :: plot_only = .FALSE. !< If true, set up as usual, then only convert the mugtok_*.rst files in the current folder to XDMF output (no time stepping, no new restart files)
INTEGER(i4) :: pre_freq = 1 !< Post-kick only: rebuild the approximate Jacobian/refactor the preconditioner every pre_freq steps (see oft_mugtok_td%pre_freq); pre-kick and kick phases always rebuild every step
REAL(r8) :: kappa_par = 1.542d26 !< Parallel (field-aligned) thermal conductivity kappa=chi*n [1/(m s)], held fixed for the whole run by setting chi_par_mult=kappa_par/kappa(1) in each phase (set in oft.in group sparc_vde_options; applied at the phase setup, so also on a restart)
!---Initial density/temperature: offset + scale*profile inside the plasma (profiles from
!   tokamaker_n.prof/tokamaker_t.prof, normalized to 1 on axis), offset outside.
!   n_scale is computed below so that 2*k*n*T on axis matches the equilibrium axis pressure
REAL(r8) :: n_scale !< Density scale [m^-3] (computed from axis pressure balance)
REAL(r8) :: n_offset = 6.3907395d19 !< Density offset [m^-3] (must be > 0)
REAL(r8) :: T_scale = 1.9568088d4 !< Temperature scale [eV]
REAL(r8) :: T_offset = 7.8d0 !< Temperature offset [eV]
REAL(r8) :: p_axis_eq, p_axis_nT !< Equilibrium vs. 2*k*n*T pressure on axis [Pa] (diagnostic)
REAL(r8) :: lin_tol = 1.d-7 !< Absolute tolerance of the linear (GMRES) solve (set in oft.in group sparc_vde_options)
REAL(r8) :: nl_tol = 1.d-6 !< Absolute tolerance on the Newton residual |F| (set in oft.in group sparc_vde_options). Newton stalls at a round-off floor that rises at lower mass_scale (~5e-7 at mass_scale 1000, ~1-3e-6 at 100 for dt 0.5 ms); if nl_tol is below that floor every step fails and is retried with dt halved
REAL(r8) :: mf_b0 = 1.d-3 !< Relative finite-difference step of the matrix-free Jacobian-vector product (oft_mf_matrix%b0; mugtok_td sets 1e-3, other OFT solvers 1e-5)
REAL(r8) :: T_clip = -1.d0 !< Hard temperature floor [eV] on the accepted state after every step (oft_mugtok_td%T_clip); <= 0 disables
REAL(r8) :: T_floor = -1.d0 !< Temperature floor [eV], per-cell scheme of the 3D xmhd.F90 (xmhd_2d T_floor): max(T,T_floor) in every term but dT/dt, plus a lagged per-cell heating source in cells that dropped below it; <= 0 disables
CHARACTER(LEN=7) :: lu_package = '' !< If set, LU package of the full-LU preconditioner ('mumps', 'umfpack', 'super', ...); empty = OFT's build default (UMFPACK if built, else MUMPS, else SuperLU)
REAL (r8):: coords(3), psi(1), q(1)
LOGICAL :: pm=.TRUE.
LOGICAL :: success
LOGICAL, allocatable, dimension(:) :: mhd_flag
CHARACTER(LEN=OFT_PATH_SLEN) :: filename_eq !< Equilibrium file (same as mesh file from oft.in "native_mesh_options")
CHARACTER(LEN=25) :: tmp_str
TYPE(oft_mugtok_td):: b_sim 
CLASS(flux_func), POINTER :: n_prof => NULL(), T_prof => NULL()
NAMELIST/sparc_vde_options/rst_start,visc_kick,visc_post,mass_scale_kick,mass_scale_post,pre_freq,plot_only,dt_min_frac, &
  visc_bulk_ratio,visc_dynamic,visc_nref,drop_convection,visc_sol_mult,visc_sol_T0,visc_sol_dT,lin_tol,nl_tol,mf_b0,T_clip, &
  visc_sol_on_psi,visc_sol_psiN0,visc_sol_dpsiN,drag_rate,fix_F_boundary,ohmic_heating,lu_reanalyze,pre_nlocal,pre_overlap,pre_part,lu_package,dt_grow_steps,dt_post,newton_backtrack,resid_dump, &
  T_floor,kappa_par,kappa_tq2,tq2_t_start,tq2_dur,dt_tq2,dt_cq,mass_scale_tq2,mass_scale_cq
!------------------------------------------------------------------------------
! Initialize enviroment
!---------------------------------------------------------------------------- d--
CALL oft_init
!---Read runtime options (optional group, defaults kept if absent)
OPEN(NEWUNIT=io_unit,FILE=oft_env%ifile)
READ(io_unit,sparc_vde_options,IOSTAT=ierr)
CLOSE(io_unit)
!---------------------------------------------------------------------------
! Setup grid
!---------------------------------------------------------------------------
CALL multigrid_construct_surf(mg_mesh)
filename_eq = mesh_filename
IF(TRIM(filename_eq)=='none')CALL oft_abort('No mesh/equilibrium file in "native_mesh_options"','sparc_vde_beta_kick',__FILE__)

order = 3
CALL oft_lag_setup(mg_mesh,order,ML_blag_obj=ML_blagrange,minlev=-1)
IF(.NOT.oft_2D_lagrange_cast(blagrange,ML_blagrange%current_level))CALL oft_abort("Invalid lagrange FE object","setup",__FILE__)
!---------------------------------------------------------------------------
! Read equilibrium from file
!---------------------------------------------------------------------------
CALL hdf5_field_get_sizes(TRIM(filename_eq),"tokamaker/PSI",ndims,dim_sizes)
npoints = dim_sizes(1)
ALLOCATE(psi_eq(npoints))
CALL hdf5_read(psi_eq,TRIM(filename_eq),"tokamaker/PSI",success)


psi_total = psi_eq  


CALL hdf5_get_group_members(TRIM(filename_eq),"tokamaker/coil_dict",coil_names)
CALL hdf5_get_group_members(TRIM(filename_eq),"tokamaker/cond_dict",cond_names)

CALL machine%setup(ML_blagrange)
machine%save_visit = .FALSE.

! coil_id/cond_id are 0-based indices into "tokamaker/coil_currents" and define
! the coil/conductor ordering; use them (not dict iteration order) to place
! each entry so downstream arrays (e.g. equil%coil_currs) line up with the file.
machine%ncoil_regs = SIZE(coil_names)
machine%ncoils = machine%ncoil_regs
ALLOCATE(machine%coil_regions(machine%ncoil_regs))
ALLOCATE(coil_nturns_file(machine%ncoil_regs))
DO i=1,SIZE(coil_names)
  CALL hdf5_read(coil_id_tmp,TRIM(filename_eq), &
    "tokamaker/coil_dict/"//TRIM(coil_names(i))//"/coil_id",success)
  CALL hdf5_read(reg_id_tmp,TRIM(filename_eq), &
    "tokamaker/coil_dict/"//TRIM(coil_names(i))//"/reg_id",success)
  CALL hdf5_read(nturns_tmp,TRIM(filename_eq), &
    "tokamaker/coil_dict/"//TRIM(coil_names(i))//"/nturns",success)
  machine%coil_regions(coil_id_tmp+1)%id = reg_id_tmp
  coil_nturns_file(coil_id_tmp+1) = nturns_tmp
END DO

machine%ncond_regs = SIZE(cond_names)
ALLOCATE(machine%cond_regions(machine%ncond_regs))
DO i=1,SIZE(cond_names)
  CALL hdf5_read(cond_id_tmp,TRIM(filename_eq), &
    "tokamaker/cond_dict/"//TRIM(cond_names(i))//"/cond_id",success)
  CALL hdf5_read(reg_id_tmp,TRIM(filename_eq), &
    "tokamaker/cond_dict/"//TRIM(cond_names(i))//"/reg_id",success)
  CALL hdf5_read(eta_tmp,TRIM(filename_eq), &
    "tokamaker/cond_dict/"//TRIM(cond_names(i))//"/eta",success)
  machine%cond_regions(cond_id_tmp+1)%id = reg_id_tmp
  machine%cond_regions(cond_id_tmp+1)%eta = eta_tmp/mu0
END DO

CALL gs_setup_walls(machine)

ALLOCATE(machine%coil_nturns(machine%fe_rep%mesh%nreg,machine%ncoils))
machine%coil_nturns = 0.d0
DO i=1,machine%ncoils
  machine%coil_nturns(machine%coil_regions(i)%id,i) = coil_nturns_file(i)
END DO

machine%free = .TRUE.
CALL machine%init()
CALL compute_bcmat(machine)
ALLOCATE(machine%coil_vcont(machine%ncoils))
machine%coil_vcont = 0.d0

!---------------------------------------------------------------------------
! Now, need to equilibrium object
!---------------------------------------------------------------------------

CALL equil%new(machine)
CALL equil%psi%restore_local(psi_total)

CALL gs_update_bounds(equil)


equil%Ip_target=ip_target*mu0
equil%ip_ratio_target=ip_ratio_target
! "pnorm"/"alam" are the (deprecated) TokaMaker names for p_scale/ffp_scale
CALL hdf5_read(equil%p_scale,TRIM(filename_eq),"tokamaker/pnorm",success)
CALL hdf5_read(equil%ffp_scale,TRIM(filename_eq),"tokamaker/alam",success)
tmp_str = 'tokamaker_f.prof'
CALL gs_profile_load(tmp_str,equil%I)
tmp_str = 'tokamaker_p.prof'
CALL gs_profile_load(tmp_str,equil%P)
equil%I%plasma_bounds=equil%plasma_bounds
equil%P%plasma_bounds=equil%plasma_bounds

equil%vcontrol_val = 0.d0

! coil_id ordering (used above for machine%coil_regions) matches "tokamaker/coil_currents"
ALLOCATE(equil%coil_currs(machine%ncoil_regs))
CALL hdf5_read(equil%coil_currs,TRIM(filename_eq),"tokamaker/coil_currents",success)
equil%coil_currs = equil%coil_currs*mu0
equil%mode = 0
equil%I%f_offset = 22.49d0

! Faithful (wall-timescale) setup: keep PHYSICAL plasma resistivity (Spitzer, on by
! default) so the plasma current is preserved (its L/R time ~ seconds >> the run),
! and let the resistive vessel (TokaMaker conductor regions) set the ~ms penetration
! time. dt is chosen large (quasi-static: we don't resolve every Alfven bounce, we
! damp them with viscosity) but small enough for the nonlinear solve to converge.
! Aim for dt*nstep_total ~ the wall time (~ms).
dt = 1.d-5
! Compressible MHD plasma: dens_reg is the ion mass [kg] per region (m_i);
! the plasma (region 1) value is inflated internally by mass_scale (default 100)
ALLOCATE(dens_reg(machine%mesh%nreg))
dens_reg = -1.d0
dens_reg(1) = m_i_plasma ! plasma ion mass [kg]
! visc_reg is the kinematic viscosity [m^2/s]. Elevated in the plasma to damp the
! fast (under-resolved) Alfven/sound waves so the column moves quasi-statically and
! the startup transient from an imperfect equilibrium settles quickly.
ALLOCATE(visc_reg(machine%mesh%nreg))
visc_reg = -1.d0
visc_reg(1) = visc_kick


! Density diffusivity [m^2/s] in the plasma (only used when n is evolved, see evolve_n)
ALLOCATE(D_reg(machine%mesh%nreg))
D_reg = 0.d0
D_reg(1) = 1.542d0

! Evolve the plasma (region 1) with compressible MHD; all other regions
! (conductors, vacuum) are handled by TokaMaker
ALLOCATE(mhd_flag(machine%mesh%nreg))
mhd_flag = .FALSE.
mhd_flag(1) = .TRUE.

equil%device => machine
b_sim%pm = pm               ! solver progress output (copied to oft_env%pm in setup_mhd)
b_sim%save_rst = .NOT.plot_only ! optional; default off (off in plot-only mode so setup doesn't overwrite mugtok_00000.rst)
b_sim%rst_freq = 1
! Use a complete LU of the coupled Jacobian (keeps the psi-velocity coupling that
! per-field block-Jacobi drops). Needed for good linear convergence in the stiff,
! near-ideal (physical eta) regime; costs more memory.
b_sim%use_full_lu = .TRUE.
! Region 1 as MHD: n and T seeded from the TokaMaker n/T profiles (see n_scale etc.);
! n held fixed (evolve_n) and T evolved, Spitzer eta(n,T), and plasma ion mass inflated by mass_scale.
! kappa_reg(1) below sets the baseline (pre-/post-kick) kappa_perp=kappa_base; the VDE is driven
! by transiently raising b_sim%mug%kappa(1) DURING the step loop (below), not by
! perturbing the coils.
tmp_str = 'tokamaker_n.prof'
CALL gs_profile_load(tmp_str,n_prof)
tmp_str = 'tokamaker_t.prof'
CALL gs_profile_load(tmp_str,T_prof)
! Choose n_scale so that 2*k*n*T matches the equilibrium pressure on axis (psi = plasma_bounds(2))
n_prof%plasma_bounds = equil%plasma_bounds
T_prof%plasma_bounds = equil%plasma_bounds
p_axis_eq = equil%p_scale*equil%P%f(equil%plasma_bounds(2))/mu0
n_scale = (p_axis_eq/(2.d0*elec_charge*(T_offset+T_scale*T_prof%fp(equil%plasma_bounds(2)))) - n_offset) &
  /n_prof%fp(equil%plasma_bounds(2))
WRITE(*,'(A,ES14.6)')' n_scale [m^-3] (from axis pressure balance) = ',n_scale
b_sim%n_prof => n_prof
b_sim%T_prof => T_prof
b_sim%n_scale = n_scale
b_sim%n_offset = n_offset
b_sim%T_scale = T_scale
b_sim%T_offset = T_offset
b_sim%evolve_n = .TRUE.    ! hold n fixed at its initial profile (T still evolves)
! Thermal conductivity kappa=chi*n [1/(m s)]
ALLOCATE(kappa_reg(machine%mesh%nreg))
kappa_reg = -1.d0
kappa_reg(1) = kappa_base
! Resistivity by region: default to the file's (TokaMaker) conductor values in both
! directions, then lower only the poloidal (in-plane) resistivity of the PLATE_* conductors
ALLOCATE(eta_reg(machine%mesh%nreg,2))
eta_reg = -1.d0
DO i=1,machine%ncond_regs
  eta_reg(machine%cond_regions(i)%id,:) = machine%cond_regions(i)%eta
END DO
DO i=1,SIZE(cond_names)
  IF(INDEX(TRIM(cond_names(i)),'PLATE_')/=1)CYCLE
  CALL hdf5_read(reg_id_tmp,TRIM(filename_eq), &
    "tokamaker/cond_dict/"//TRIM(cond_names(i))//"/reg_id",success)
  eta_reg(reg_id_tmp,1) = eta_plate_pol/mu0
  WRITE(*,'(2A,I4,2(A,ES12.4))')' ',TRIM(cond_names(i)),reg_id_tmp, &
    '  eta_pol [Ohm*m] = ',eta_reg(reg_id_tmp,1)*mu0,'  eta_tor [Ohm*m] = ',eta_reg(reg_id_tmp,2)*mu0
END DO
b_sim%fix_F_boundary = fix_F_boundary
WRITE(*,'(A,L2)')' Fixed F = R*B_phi on the outer boundary (fix_F_boundary): ',fix_F_boundary
b_sim%pre_nlocal = pre_nlocal
b_sim%pre_overlap = pre_overlap
b_sim%pre_part = pre_part
CALL b_sim%setup_mhd(equil, dt, lin_tol, nl_tol, mhd_flag, dens_reg, visc_reg, eta_reg=eta_reg, kappa_reg=kappa_reg, mass_scale=mass_scale_kick, D_reg=D_reg, toroidal_flow = .FALSE.)
b_sim%mug%ohmic = ohmic_heating
b_sim%lu_reanalyze = lu_reanalyze
WRITE(*,'(A,L2)')' Redo the LU symbolic analysis at every rebuild (lu_reanalyze): ',lu_reanalyze
WRITE(*,'(A,L2)')' Ohmic heating in the temperature equation (ohmic_heating): ',ohmic_heating

! Compare the seeded 2*k*n*T with the equilibrium pressure on axis (psihat=1, where the
! normalized n/T profiles are 1)
p_axis_eq = b_sim%tkmr%gs_equil%p_scale*b_sim%tkmr%gs_equil%P%f(b_sim%tkmr%gs_equil%plasma_bounds(2))/mu0
p_axis_nT = 2.d0*elec_charge*(n_offset+n_scale*n_prof%fp(b_sim%tkmr%gs_equil%plasma_bounds(2))) &
  *(T_offset+T_scale*T_prof%fp(b_sim%tkmr%gs_equil%plasma_bounds(2)))
WRITE(*,'(A,ES12.4,A,ES12.4)')' Axis pressure [Pa]: equilibrium = ',p_axis_eq,'  2kn*T = ',p_axis_nT

! TEMPORARY: disable the temperature sink (setup_mhd sets tau_sink=1.d-3, T_sink=5 eV
! for DOFs left behind as the plasma boundary recedes); tau_sink<=0 turns it off. Remove
! this line to restore the sink.
b_sim%mug%tau_sink = -1.d0

! Keep the (v.grad)v convective term in the momentum equation unless drop_convection=T
! (setup_mhd drops it by default)
b_sim%mug%drop_convection = drop_convection

! Viscosity: bulk/shear split, and optionally a uniform dynamic viscosity, nu = visc*visc_nref/n
! (nu_n holds nu*n; region 1 is updated with nu(1) at each phase change below)
b_sim%mug%nu2_ratio = visc_bulk_ratio
IF(visc_dynamic)THEN
  ALLOCATE(b_sim%mug%nu_n(SIZE(b_sim%mug%nu)))
  b_sim%mug%nu_n = b_sim%mug%nu*visc_nref
END IF
WRITE(*,'(A,L2,A,ES12.4,A,F6.2)')' Viscosity: dynamic (uniform m_i*n*nu) =',visc_dynamic, &
  '  n_ref [m^-3] =',visc_nref,'  bulk/shear ratio =',visc_bulk_ratio
b_sim%mug%visc_sol_T0 = visc_sol_T0
b_sim%mug%visc_sol_dT = visc_sol_dT
b_sim%mug%visc_sol_on_psi = visc_sol_on_psi
b_sim%mug%visc_sol_psiN0 = visc_sol_psiN0
b_sim%mug%visc_sol_dpsiN = visc_sol_dpsiN
WRITE(*,'(A,ES12.4,A)')' Linear drag on V_pol (post-kick): drag_rate =',drag_rate,' 1/s (0: off)'
IF(visc_sol_on_psi)WRITE(*,'(A,F8.4,A,F8.4,A)')' Viscosity: SOL multiplier keyed on psiN (not T): above psiN0 =',visc_sol_psiN0,' (width',visc_sol_dpsiN,')'
WRITE(*,'(A,F8.2,A,F8.2,A,F8.2,A)')' Viscosity: post-kick SOL multiplier =',visc_sol_mult, &
  '  below T0 =',visc_sol_T0,' eV (width',visc_sol_dT,' eV)'
WRITE(*,'(A,ES10.2,A,ES10.2)')' Solver tolerances: nl_tol (|F|) =',nl_tol,'  lin_tol =',lin_tol
b_sim%mfmat%b0 = mf_b0
WRITE(*,'(A,ES10.2)')' Matrix-free Jacobian: finite-difference step mf_b0 =',mf_b0
b_sim%T_clip = T_clip
WRITE(*,'(A,F8.2,A)')' Hard temperature floor T_clip =',T_clip,' eV (<= 0: off)'
b_sim%mug%T_floor = T_floor
WRITE(*,'(A,F8.2,A)')' Temperature floor T_floor (per-cell, xmhd.F90 scheme) =',T_floor,' eV (<= 0: off)'
WRITE(*,'(A,I4,A)')' dt recovery: double dt after',dt_grow_steps,' clean steps (0: off)'
WRITE(*,'(A,ES12.4,A)')' Post-kick nominal dt (dt_post) =',dt_post,' s'
WRITE(*,'(A,ES12.4,A,ES12.4)')' Parallel conductivity kappa_par =',kappa_par,' 1/(m s); post-kick kappa_par/kappa_perp =',kappa_par/kappa_base
b_sim%nksolver%backtrack = newton_backtrack
WRITE(*,'(A,L2)')' Newton backtracking:',newton_backtrack
b_sim%resid_dump = resid_dump
WRITE(*,'(A,L2)')' Residual dump on Newton failure:',resid_dump
SELECT TYPE(lu_pre => b_sim%pre)
  TYPE IS(oft_lusolver)
    IF(LEN_TRIM(lu_package)>0)lu_pre%package = lu_package
    WRITE(*,'(2A)')' Preconditioner: full LU, package ',TRIM(lu_pre%package)
  TYPE IS(oft_bjprecond)
    SELECT TYPE(sub => lu_pre%pre)
      TYPE IS(oft_lusolver)
        IF(LEN_TRIM(lu_package)>0)sub%package = lu_package
        WRITE(*,'(A,I4,A,L2,2A)')' Preconditioner: block-Jacobi with',pre_nlocal,' subdomains (pre_nlocal), overlap', &
          pre_overlap,', LU package ',TRIM(sub%package)
    END SELECT
  CLASS DEFAULT
    IF(LEN_TRIM(lu_package)>0)CALL oft_abort('lu_package needs the full-LU preconditioner (use_full_lu)','sparc_vde_beta_kick',__FILE__)
    WRITE(*,'(A)')' Preconditioner: not a single LU (block-Jacobi)'
END SELECT

! Coil currents are held fixed for the whole run; step() overwrites equil%coil_currs
! internally, so a separate fixed reference copy is needed to pass in every step.
ALLOCATE(coil_currs_fixed(machine%ncoil_regs))
coil_currs_fixed = equil%coil_currs

! Plot-only mode: write XDMF output from the existing mugtok_*.rst files and stop
! (range set by rst_start/rst_end in the mugtok_plot_options group)
IF(plot_only)THEN
  CALL b_sim%plot()
  CALL oft_finalize()
END IF

! Optionally resume a previous run from restart file 'mugtok_NNNNN.rst' (NNNNN=rst_start).
! The restart index equals the step count, so the loop continues from step rst_start+1.
IF(rst_start > 0)THEN
  CALL b_sim%restart(rst_start, t, dt_rst)
  WRITE(*,'(A,I6,A,ES14.6,A,ES12.4)')' Resuming from restart file ',rst_start,'  t = ',t,'  dt = ',dt_rst
END IF

IF(rst_start <= 0)THEN
  CALL oft_abort('needs rst_start > 0 (a post-kick restart file of a sparc_vde_beta_kick run)','sparc_vde_tq2',__FILE__)
END IF
IF(tq2_t_start <= 0.d0) tq2_t_start = t
IF(dt_grow_steps <= 0) dt_grow_steps = 3
IF(mass_scale_tq2 <= 0.d0) mass_scale_tq2 = mass_scale_post
IF(mass_scale_cq <= 0.d0) mass_scale_cq = mass_scale_post
WRITE(*,'(A,ES12.4,A,ES12.4,A,ES12.4,A,ES12.4,A,ES12.4)')' Second thermal quench: kappa_tq2 =',kappa_tq2, &
  ' from t =',tq2_t_start,' for',tq2_dur,' s at dt =',dt_tq2,'; then current quench, dt up to',dt_cq
WRITE(*,'(A,ES12.4,A,I3,A)')'   kappa_par/kappa_perp =',kappa_par/kappa_tq2,'; dt doubles after',dt_grow_steps,' clean steps'
WRITE(*,'(A,F10.2,A,F10.2,A,F10.2)')'   mass scale: before the quench',mass_scale_post,', quench',mass_scale_tq2, &
  ', current quench',mass_scale_cq
DO i=rst_start+1,nstep_total
  ! Phase 2 (post-kick, as in the parent run) until tq2_t_start, then 3 (quench: kappa_tq2, dt_tq2)
  ! for tq2_dur, then 4 (current quench: kappa_tq2 kept, dt grows towards dt_cq)
  IF(t < tq2_t_start*(1.d0-1.d-12))THEN
    kick_phase = 2
  ELSE IF(t < tq2_t_start+tq2_dur)THEN
    kick_phase = 3
  ELSE
    kick_phase = 4
  END IF
  IF(kick_phase /= kick_phase_prev)THEN
    b_sim%mug%nu(1) = visc_post
    b_sim%mug%visc_sol_mult = visc_sol_mult
    b_sim%mug%m_i(1) = m_i_plasma*mass_scale_post
    b_sim%mug%evolve_pressure = .TRUE.
    b_sim%mug%drag_rate = drag_rate
    SELECT CASE(kick_phase)
    CASE(2) ! before the quench: the parent run's post-kick settings
      b_sim%mug%kappa(1) = kappa_base
      dt_phase = dt_post
      IF(kick_phase_prev < 0) THEN
        dt = dt_rst
      ELSE
        dt = MIN(dt, dt_post)
      END IF
      b_sim%pre_freq = pre_freq
    CASE(3) ! second thermal quench: kappa_perp up, short steps, fresh preconditioner every step
      b_sim%mug%kappa(1) = kappa_tq2
      b_sim%mug%m_i(1) = m_i_plasma*mass_scale_tq2
      dt = dt_tq2
      dt_phase = dt_tq2
      b_sim%pre_freq = 1
    CASE(4) ! current quench: keep kappa_perp, let dt grow from the quench dt to dt_cq
      b_sim%mug%kappa(1) = kappa_tq2
      b_sim%mug%m_i(1) = m_i_plasma*mass_scale_cq
      dt_phase = dt_cq
      b_sim%pre_freq = 1
    END SELECT
    b_sim%mug%chi_par_mult = kappa_par/b_sim%mug%kappa(1) ! keep kappa_par fixed
    IF(visc_dynamic) b_sim%mug%nu_n(1) = b_sim%mug%nu(1)*visc_nref ! nu = nu(1)*visc_nref/n
    dt_stop_ref = dt
    kick_phase_prev = kick_phase
    n_clean = 0
    b_sim%pre_force_update = .TRUE.
    WRITE(*,'(A,I2,A,ES12.4,A,ES12.4,A,ES12.4)')' Entering phase',kick_phase,' at t =',t,': kappa_perp =', &
      b_sim%mug%kappa(1),', dt =',dt
  END IF
  write(*,*) "Step: ", i, "  t = ", t, "  kick_phase = ", kick_phase, "  kappa(1) = ", b_sim%mug%kappa(1), &
    "  kappa_par = ", b_sim%mug%kappa(1)*b_sim%mug%chi_par_mult, &
    "  nu(1) = ", b_sim%mug%nu(1), "  mass_scale = ", b_sim%mug%m_i(1)/m_i_plasma, &
    "  evolve_pressure = ", b_sim%mug%evolve_pressure
  IF(ohmic_heating) WRITE(*,'(A,ES12.4,A)') '  Ohmic power (start of step) = ', b_sim%mug%ohm_power, ' W'
  write(*,*) "opoint: ", b_sim%tkmr%gs_equil%o_point
  write(*,*) "diverted: ", b_sim%tkmr%gs_equil%diverted
  CALL b_sim%step(coil_currs_fixed, t, dt, nl_its, l_its, nretry)
  ! Stop once solver-failure retries have cut dt too far below this phase's starting dt
  ! (convergence has broken down); the restart file for this step is already written
  IF(dt_min_frac > 0.d0 .AND. dt < dt_min_frac*dt_stop_ref*(1.d0-1.d-8))THEN
    WRITE(*,'(A,ES12.4,A,ES12.4,A,I2,A,I7)')' Stopping: dt = ',dt,' is below dt_min_frac*(dt at phase entry) = ', &
      dt_min_frac*dt_stop_ref,' (phase',kick_phase,') at step ',i
    EXIT
  END IF
  ! Let a dt reduced by retries grow back toward the phase's nominal dt after dt_grow_steps clean steps
  IF(nretry == 0)THEN
    n_clean = n_clean + 1
  ELSE
    n_clean = 0
  END IF
  IF(kick_phase /= 3 .AND. dt_grow_steps > 0 .AND. n_clean >= dt_grow_steps .AND. dt < dt_phase*(1.d0-1.d-8))THEN
    dt = MIN(2.d0*dt, dt_phase)
    n_clean = 0
    WRITE(*,'(A,ES12.4,A,I4,A)')' dt grown to ',dt,' after ',dt_grow_steps,' clean steps'
  END IF
  ! IF (MOD(i,25)==0 ) THEN
  !   CALL b_sim%plot()
  ! END IF
END DO
CALL b_sim%plot()
! ALLOCATE(voltages(machine%ncoil_regs))
! voltages = 0.d0
! CALL tokamaker%setup(equil, dt, lin_tol, nl_tol, .FALSE.)
! CALL tokamaker%step(equil%coil_currs, voltages, t, dt, nl_its, l_its, nretry)



END PROGRAM sparc_vde_tq2