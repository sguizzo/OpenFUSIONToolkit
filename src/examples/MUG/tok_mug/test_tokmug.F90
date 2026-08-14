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
USE oft_io, ONLY: hdf5_field_get_sizes, hdf5_read, hdf5_field_exist
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
REAL (r8):: ip_ratio_target = 0.205
REAL (r8):: ip_target = 7.87E6
REAL(r8), allocatable, dimension(:) :: psi_eq, psi_pert, psi_total, eta_reg,curr_reg, areas, dens_reg, visc_reg, voltages
REAL(r8), allocatable, dimension(:) :: coil_step !< Coil currents passed to each step (baseline, then pulsed)
REAL(r8), allocatable, dimension(:) :: coil_base !< Fixed baseline coil currents (step overwrites equil%coil_currs)
INTEGER(i4) :: nstep_total = 5000 !< Number of time steps (dt*nstep should approach the wall time ~ms)
INTEGER(i4) :: pulse_start = 2!< Step at which the coil pulse turns on (let the system settle first)
INTEGER(i4) :: pulse_end = 12 !< Step at which the coil pulse turns off (let the system settle first)
LOGICAL :: do_pulse = .TRUE. !< Run WITH the coil pulse (.TRUE.) or a no-pulse baseline (.FALSE.) for common-mode subtraction
REAL(r8) :: lin_tol = 1.d-11
REAL(r8) :: nl_tol = 1.d-9
REAL (r8):: coords(3), psi(1), q(1)
LOGICAL :: pm=.FALSE.
LOGICAL :: success
LOGICAL, allocatable, dimension(:) :: mhd_flag
CHARACTER(LEN=25) :: filename_eq = 'eq_o3.h5' !< Name of input file for mesh, fix later for variable length
CHARACTER(LEN=25) :: filename_pert= 'paper_pert_0115.h5' !< Name of input file for mesh, fix later for variable length
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

CALL hdf5_field_get_sizes(TRIM(filename_pert),"tokamaker/PSI",ndims,dim_sizes)
npoints = dim_sizes(1)
ALLOCATE(psi_pert(npoints))
CALL hdf5_read(psi_pert,TRIM(filename_pert),"tokamaker/PSI",success)

psi_total = psi_eq  ! - 0.1 * psi_pert


CALL machine%setup(ML_blagrange)
machine%save_visit = .FALSE.
machine%ncoil_regs = 7
machine%ncoils = 7
ALLOCATE(machine%coil_regions(machine%ncoil_regs))
machine%coil_regions(1)%id = 9
machine%coil_regions(1)%area = 1.8
machine%coil_regions(2)%id = 10
machine%coil_regions(2)%area = 0.25
machine%coil_regions(3)%id = 11
machine%coil_regions(3)%area = 0.25
machine%coil_regions(4)%id = 12
machine%coil_regions(4)%area = 0.25
machine%coil_regions(5)%id = 13
machine%coil_regions(5)%area = 0.25
machine%coil_regions(6)%id = 14
machine%coil_regions(6)%area = 0.25
machine%coil_regions(7)%id = 15
machine%coil_regions(7)%area = 0.25

machine%ncond_regs = 5
ALLOCATE(machine%cond_regions(machine%ncond_regs))
machine%cond_regions(1)%id = 4
machine%cond_regions(1)%eta = 6.9d-7*1000.d0/mu0
machine%cond_regions(2)%id = 5
machine%cond_regions(2)%eta = 1.14d-6*1000.d0/mu0
machine%cond_regions(3)%id = 6
machine%cond_regions(3)%eta = 1.14d-6*1000.d0/mu0
machine%cond_regions(4)%id = 7
machine%cond_regions(4)%eta = 6.9d-7*1000.d0/mu0
machine%cond_regions(5)%id = 8
machine%cond_regions(5)%eta = 6.9d-7*1000.d0/mu0


CALL gs_setup_walls(machine)

ALLOCATE(machine%coil_nturns(machine%fe_rep%mesh%nreg,machine%ncoils))
machine%coil_nturns = 0
DO j=1, machine%ncoils
  machine%coil_nturns(j + 8, j) = 1
END DO

machine%free = .TRUE.
CALL machine%init()
CALL compute_bcmat(machine)
ALLOCATE(machine%coil_vcont(machine%ncoils))
machine%coil_vcont = 0.d0
! write(*,*) 82
! !---------------------------------------------------------------------------
! ! Now, need to setup a tokamaker device
! !---------------------------------------------------------------------------
! CALL machine%setup(ML_blagrange)
! machine%region_info%nnonaxi = 0
! ALLOCATE(machine%region_info%reg_map(machine%fe_rep%mesh%nreg))
! machine%region_info%reg_map=0
! machine%free = .TRUE.

! IF(.NOT.ASSOCIATED(machine%ignore_rmask))THEN
!   ALLOCATE(machine%ignore_rmask(machine%fe_rep%mesh%nreg))
!   machine%ignore_rmask=.FALSE.
!   machine%ignore_rmask(5)=.TRUE.
!   machine%ignore_rmask(6)=.TRUE.
! END IF

! ! machine%ncoils = 7
! machine%ncoil_regs = 7
! ! ALLOCATE(machine%coil_nturns(machine%fe_rep%mesh%nreg,machine%ncoils))
! ! machine%coil_nturns = 0
! ! DO j=1, machine%ncoils
! !   machine%coil_nturns(j + 8, j) = 1
! ! END DO
! write(*,*) 101
! ! machine%coil_vcont = 0.d0
! write(*,*) 108
! CALL gs_setup_walls(machine)
! write(*,*) 110
! CALL machine%init()
! machine%ncond_regs = 5
! ALLOCATE(machine%cond_regions(machine%ncond_regs))
! machine%cond_regions(1)%id = 4
! machine%cond_regions(1)%eta = 6.9d-7/mu0
! machine%cond_regions(2)%id = 5
! machine%cond_regions(2)%eta = 1.14d-6/mu0
! machine%cond_regions(3)%id = 6
! machine%cond_regions(3)%eta = 1.14d-6/mu0
! machine%cond_regions(4)%id = 7
! machine%cond_regions(4)%eta = 6.9d-7/mu0
! machine%cond_regions(5)%id = 8
! machine%cond_regions(5)%eta = 6.9d-7/mu0
! ALLOCATE(machine%coil_regions(machine%ncoil_regs))
! DO j=1, machine%ncoils
!   machine%coil_regions(j)%id = 8 + j
! END DO
! CALL compute_bcmat(machine)
!---------------------------------------------------------------------------
! Now, need to equilibrium object
!---------------------------------------------------------------------------

CALL equil%new(machine)
CALL equil%psi%restore_local(psi_total)

CALL gs_update_bounds(equil)


equil%Ip_target=ip_target*mu0
equil%ip_ratio_target=ip_ratio_target
!ORDER 2
! equil%p_scale = 0.28658184156588085
! equil%ffp_scale = 2.596639717247778
! ORDER 3
equil%p_scale = 0.2866168882023207
equil%ffp_scale = 2.5969479398703386
tmp_str = 'tokamaker_f.prof'
CALL gs_profile_load(tmp_str,equil%I)
tmp_str = 'tokamaker_p.prof'
CALL gs_profile_load(tmp_str,equil%P)
equil%I%plasma_bounds=equil%plasma_bounds
equil%P%plasma_bounds=equil%plasma_bounds

equil%vcontrol_val = 0.d0

ALLOCATE(equil%coil_currs(machine%ncoil_regs))
equil%coil_currs = [-10004540.249054534, 8131096.462417697, 8130193.40495448, -2752265.1799410507, -2755148.7923723585, -1429157.477860515,  -1424955.0239338675]*mu0
equil%coil_currs = [-10004577.39633093, 8099063.653129182, 8099246.217039664, -2754682.325483237, -2754900.2110607256, -1425245.0543276411, -1424994.322781849]*mu0
ALLOCATE(areas(machine%ncoil_regs))
areas = [1.8,0.25, 0.25, 0.25, 0.25, 0.25, 0.25 ]
! equil%coil_currs = equil%coil_currs/areas
equil%mode = 0
equil%I%f_offset = 36.d0

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
visc_reg(1) = 1.d4

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
CALL b_sim%setup_mhd(equil, dt, lin_tol, nl_tol, mhd_flag, dens_reg, visc_reg)

! Capture a FIXED baseline: step() overwrites equil%coil_currs each step with the
! prescribed values, so we cannot use it as the unperturbed reference in the loop.
ALLOCATE(coil_base(machine%ncoil_regs), coil_step(machine%ncoil_regs))
coil_step = equil%coil_currs
! DO i=1,nstep_total
!   ! Start from the fixed baseline each step, then (if do_pulse) pulse a pair of PF
!   ! coils on from step `pulse_start` onward. Because field 6 is carried over between
!   ! steps (not re-derived), this changes the TOTAL flux and drives a real response.
!   ! For a clean signal, run once with do_pulse=.TRUE. and once with .FALSE. from the
!   ! same setup and subtract: the equilibrium-imbalance transient is common-mode and
!   ! cancels, leaving only the coil-driven motion.
!   ! coil_step = coil_base
!   IF(do_pulse .AND. i>=pulse_start .AND. i<=pulse_end)THEN
!     coil_step(2) = coil_step(2) + 1.d5*mu0/10.0d0
!     coil_step(3) = coil_step(3) - 1.d5*mu0/10.0d0
!     dt = 4.d-5
!   ELSE
!     dt = 4.d-5
!   END IF
!   write(*,*) "Step: ", i, "  t = ", t
!   write(*,*) "Coil currents: ", coil_step(1:3)
!   write(*,*) "opoint: ", b_sim%tkmr%gs_equil%o_point
!   CALL b_sim%step(coil_step, t, dt, nl_its, l_its, nretry)
!   ! IF (MOD(i,25)==0 ) THEN
!   !   CALL b_sim%plot()
!   ! END IF
! END DO
CALL b_sim%plot()
! ALLOCATE(voltages(machine%ncoil_regs))
! voltages = 0.d0
! CALL tokamaker%setup(equil, dt, lin_tol, nl_tol, .FALSE.)
! CALL tokamaker%step(equil%coil_currs, voltages, t, dt, nl_its, l_its, nretry)



END PROGRAM gs_driver_full