PROGRAM gs_driver_full
!---Runtime
USE oft_base
!---Grid
USE multigrid, ONLY: multigrid_mesh
USE multigrid_build, ONLY: multigrid_construct_surf
!
USE oft_la_base, ONLY: oft_vector, oft_matrix
USE oft_solver_base, ONLY: oft_solver
USE oft_solver_utils, ONLY: create_cg_solver, create_diag_pre
!
USE oft_blag_operators, ONLY: oft_blag_zerob, oft_blag_getmop, oft_blag_project
USE oft_scalar_inits, ONLY: poss_scalar_bfield
USE mhd_utils, ONLY: elec_charge, proton_mass, mu0
USE oft_io, ONLY: hdf5_field_get_sizes, hdf5_read, hdf5_field_exist, hdf5_get_group_members
USE oft_gs, ONLY: gs_equil, gs_update_bounds, gs_test_bounds, compute_bcmat, gs_setup_walls, gs_get_qprof, gs_factory
USE oft_gs_td, ONLY: oft_tmaker_td
USE oft_gs_util, ONLY: gs_profile_load
USE oft_lag_basis, ONLY: oft_lag_setup,oft_scalar_bfem, oft_blag_eval, oft_blag_geval, oft_2D_lagrange_cast
USE fem_base, ONLY: oft_ml_fem_type
USE mugtok_td

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
REAL(r8) :: pulse_frac !< Trapezoid ramp/flattop/ramp-down fraction in [0,1] applied to pulse_amp
INTEGER(i4) :: pulse_phase = -1, pulse_phase_prev = -1 !< Current/previous pulse-schedule phase (0=pre, 1=ramp-up, 2=flattop, 3=ramp-down, 4=post); used to apply the phase's target dt only on entry, so a dt reduced by step() due to a solver-failure retry is not clobbered on later iterations of the same phase
REAL (r8):: ip_ratio_target = 0.205
REAL (r8):: ip_target = 7.87E6
REAL(r8), allocatable, dimension(:) :: psi_eq, psi_pert, psi_total, eta_reg,curr_reg, dens_reg, kappa_reg,visc_reg, voltages
!---Coil/conductor region info read from the equilibrium file's coil_dict/cond_dict
CHARACTER(LEN=80), ALLOCATABLE, DIMENSION(:) :: coil_names, cond_names
INTEGER(i4), ALLOCATABLE, DIMENSION(:) :: coil_nturns_file
INTEGER(i4) :: reg_id_tmp, coil_id_tmp, nturns_tmp, cond_id_tmp
REAL(r8) :: eta_tmp
REAL(r8), allocatable, dimension(:) :: coil_step !< Coil currents passed to each step (baseline, then pulsed)
REAL(r8), allocatable, dimension(:) :: coil_base !< Fixed baseline coil currents (step overwrites equil%coil_currs)
INTEGER(i4) :: nstep_total = 10000 !< Number of time steps (dt*nstep should approach the wall time ~ms)
REAL(r8) :: t_pulse_start = 5.d-4 !< Time at which the coil pulse ramp-up begins (let the system settle first) [s]
REAL(r8) :: t_ramp_dur = 5.d-4 !< Duration of the ramp-up (and ramp-down) leg of the pulse [s] (0.5 ms each: 1 ms of ramp + 1 ms flattop = 2 ms pulse)
REAL(r8) :: t_flat_dur = 1.d-3 !< Duration of the flattop [s]
REAL(r8) :: pulse_amp = 1.d4*mu0 !< Peak coil-pair perturbation amplitude at the top of the flattop
LOGICAL :: do_pulse = .TRUE. !< Run WITH the coil pulse (.TRUE.) or a no-pulse baseline (.FALSE.) for common-mode subtraction
REAL(r8) :: lin_tol = 1.d-11
REAL(r8) :: nl_tol = 1.d-9
REAL (r8):: coords(3), psi(1), q(1)
LOGICAL :: pm=.TRUE.
LOGICAL :: success
LOGICAL, allocatable, dimension(:) :: mhd_flag
CHARACTER(LEN=25) :: filename_eq = 'sparc_mug_2cm.h5' !< Name of input file for mesh, fix later for variable length
CHARACTER(LEN=25) :: tmp_str
TYPE(oft_mugtok_td):: b_sim 
!------------------------------------------------------------------------------
! Initialize enviroment
!---------------------------------------------------------------------------- d--
CALL oft_init
!---------------------------------------------------------------------------
! Setup grid
!---------------------------------------------------------------------------
CALL multigrid_construct_surf(mg_mesh)

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
lin_tol = 1.d-5
nl_tol = 1.d-4
! Compressible MHD plasma: dens_reg is the ion mass [kg] per region (m_i);
! the plasma (region 1) value is inflated internally by mass_scale (default 100)
ALLOCATE(dens_reg(machine%mesh%nreg))
dens_reg = -1.d0
dens_reg(1) = 2.d0*proton_mass ! plasma ion mass [kg]
! visc_reg is the kinematic viscosity [m^2/s]. Elevated in the plasma to damp the
! fast (under-resolved) Alfven/sound waves so the column moves quasi-statically and
! the startup transient from an imperfect equilibrium settles quickly.
ALLOCATE(visc_reg(machine%mesh%nreg))
visc_reg = -1.d0
visc_reg(1) = 50.d0


! Thermal conductivity kappa=chi*n [1/(m s)]: chi=1.d-1 m^2/s at n=n_scale
ALLOCATE(kappa_reg(machine%mesh%nreg))
kappa_reg = -1.d0
kappa_reg(1) = 1.d-1*b_sim%n_scale

! Evolve the plasma (region 1) with compressible MHD; all other regions
! (conductors, vacuum) are handled by TokaMaker
ALLOCATE(mhd_flag(machine%mesh%nreg))
mhd_flag = .FALSE.
mhd_flag(1) = .TRUE.

equil%device => machine
b_sim%save_rst = .TRUE.     ! optional; default off
b_sim%rst_freq = 1
! Use a complete LU of the coupled Jacobian (keeps the psi-velocity coupling that
! per-field block-Jacobi drops). Needed for good linear convergence in the stiff,
! near-ideal (physical eta) regime; costs more memory.
b_sim%use_full_lu = .TRUE.
! Region 1 as MHD: density flat at 1e19, T = T(psi) closure, Spitzer eta(n,T),
! and plasma ion mass inflated by mass_scale (default 100).
! NOTE: set up in the UNPERTURBED equilibrium. Perturbing the coils before setup
! has no effect, because setup defines field 6 as (total equilibrium flux - coil
! vacuum flux), so any coil change is exactly cancelled in the initial total flux.
! To drive the plasma, the coils must be pulsed DURING the step loop (below).
CALL b_sim%setup_mhd(equil, dt, lin_tol, nl_tol, mhd_flag, dens_reg, visc_reg, kappa_reg=kappa_reg, mass_scale=100000.d0)

! TEMPORARY: disable the temperature sink (setup_mhd sets tau_sink=1.d-3, T_sink=5 eV
! for DOFs left behind as the plasma boundary recedes); tau_sink<=0 turns it off. Remove
! this line to restore the sink.
b_sim%mug%tau_sink = -1.d0

! Capture a FIXED baseline: step() overwrites equil%coil_currs each step with the
! prescribed values, so we cannot use it as the unperturbed reference in the loop.
ALLOCATE(coil_base(machine%ncoil_regs), coil_step(machine%ncoil_regs))
coil_base = equil%coil_currs
coil_step = coil_base
DO i=1,nstep_total
  ! Start from the fixed baseline each step, then (if do_pulse) apply a trapezoid
  ! ramp-up / flattop / ramp-down perturbation to a pair of PF coils, timed off the
  ! simulation clock `t` (not the step index), so the pulse shape is independent of
  ! dt. Because field 6 is carried over between steps (not re-derived), this changes
  ! the TOTAL flux and drives a real response. For a clean signal, run once with
  ! do_pulse=.TRUE. and once with .FALSE. from the same setup and subtract: the
  ! equilibrium-imbalance transient is common-mode and cancels, leaving only the
  ! coil-driven motion.
  pulse_frac = 0.d0
  IF(do_pulse)THEN
    IF(t < t_pulse_start)THEN
      pulse_frac = 0.d0
      pulse_phase = 0
      IF (t > 0.d0) dt = 5.d-5
    ELSE IF(t < t_pulse_start+t_ramp_dur)THEN
      pulse_frac = (t-t_pulse_start)/t_ramp_dur                    ! ramp-up
      pulse_phase = 1
    ELSE IF(t < t_pulse_start+t_ramp_dur+t_flat_dur)THEN
      pulse_frac = 1.d0                                            ! flattop
      pulse_phase = 2
    ELSE IF(t < t_pulse_start+2.d0*t_ramp_dur+t_flat_dur)THEN
      pulse_frac = 1.d0 - (t-(t_pulse_start+t_ramp_dur+t_flat_dur))/t_ramp_dur ! ramp-down
      pulse_phase = 3
    ELSE
      pulse_frac = 0.d0
      pulse_phase = 4
    END IF
    ! Only (re-)apply the phase's target dt when first entering that phase. On later
    ! iterations of the same phase, leave `dt` alone: step() writes back whatever dt it
    ! actually used (possibly halved by a solver-failure retry), and re-imposing the
    ! phase's nominal dt here every iteration would silently discard that reduction.
    IF(pulse_phase /= pulse_phase_prev)THEN
      IF(pulse_phase==1 .OR. pulse_phase==2 .OR. pulse_phase==3)THEN
        dt = 5.d-5
      ELSE IF(pulse_phase==0)THEN
        dt = 5.d-4
      ELSE
        dt = 1.d-3
      END IF
      pulse_phase_prev = pulse_phase
    END IF
  END IF
  coil_step = coil_base
  ! Perturb the VS1U/VS1L pair (coil_id 20/21 -> 1-based indices 21/22)
  coil_step(21) = coil_base(21) + pulse_frac*pulse_amp
  coil_step(22) = coil_base(22) - pulse_frac*pulse_amp
  write(*,*) "Step: ", i, "  t = ", t
  write(*,*) "Coil currents: ", coil_step(20:22)
  write(*,*) "opoint: ", b_sim%tkmr%gs_equil%o_point
  write(*,*) "diverted: ", b_sim%tkmr%gs_equil%diverted
  CALL b_sim%step(coil_step, t, dt, nl_its, l_its, nretry)
  ! IF (MOD(i,25)==0 ) THEN
  !   CALL b_sim%plot()
  ! END IF
END DO
CALL b_sim%plot()
! ALLOCATE(voltages(machine%ncoil_regs))
! voltages = 0.d0
! CALL tokamaker%setup(equil, dt, lin_tol, nl_tol, .FALSE.)
! CALL tokamaker%step(equil%coil_currs, voltages, t, dt, nl_its, l_its, nretry)



END PROGRAM gs_driver_full