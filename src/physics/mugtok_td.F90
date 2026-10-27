MODULE mugtok_td
USE oft_base
USE oft_io, ONLY: hdf5_read, xdmf_plot_file, oft_file_exist
USE oft_quadrature, ONLY: oft_quad_type
USE oft_mesh_type, ONLY: oft_bmesh, cell_is_curved
USE multigrid, ONLY: multigrid_mesh
USE oft_gauss_quadrature, ONLY: set_quad_1d
!
USE oft_la_base, ONLY: oft_vector, oft_matrix, oft_graph, oft_graph_ptr, map_list
USE oft_deriv_matrices, ONLY: oft_noop_matrix, oft_mf_matrix
USE oft_solver_base, ONLY: oft_solver
USE oft_native_solvers, ONLY: oft_nksolver, oft_native_gmres_solver, oft_bjprecond
USE oft_solver_utils, ONLY: create_cg_solver, create_diag_pre, create_bjacobi_pre, create_native_pre
USE oft_lu, ONLY: oft_lusolver
USE oft_la_utils, ONLY: create_matrix, graph_add_dense_blocks, create_vector, create_dense_graph
USE oft_native_la, ONLY: oft_native_matrix
!
USE fem_composite, ONLY: oft_fem_comp_type, fem_graph_create
USE fem_utils, ONLY: fem_dirichlet_diag, fem_dirichlet_vec
USE oft_lag_basis, ONLY: oft_scalar_bfem, oft_blag_eval, oft_blag_geval, oft_blag_npos
USE oft_blag_operators, ONLY: oft_blag_vproject, oft_blag_getmop, oft_lag_bginterp
USE mhd_utils, ONLY: mu0
USE oft_gs, ONLY: gs_epsilon, gs_equil, gs_factory, gs_update_bounds, gs_test_bounds, flux_func, build_dels
USE oft_gs_td, ONLY: oft_tmaker_td_mfop, build_vac_op, apply_rhs
USE oft_mesh_local_util, ONLY: mesh_local_findedge
USE xmhd_2d, ONLY: oft_xmhd_2d_sim, build_approx_jacobian, xmhd_2d_freeze_eta_temp, xmhd_2d_freeze_ohmic, xmhd_2d_floor_update
USE oft_stitching, ONLY: seam_list
IMPLICIT NONE
#include "local.h"
#if !defined(TDIFF_RST_LEN)
#define TDIFF_RST_LEN 5
#endif
PRIVATE
!> L2-projected poloidal gradient of F (dF/dR, dF/dZ); helper used by the C wrapper for output
PUBLIC :: compute_gradf

!-----------------------------------------------------------------------------------------
!> Combined MUG and TokaMaker time-dependent simulation object,
!> which allows simultaneous evolution of plasma, solid, and flowing conductor regions
!-----------------------------------------------------------------------------------------
TYPE, public :: oft_mugtok_td
CLASS(oft_vector), POINTER :: u => NULL() !< Current solution vector
CLASS(oft_vector), POINTER :: rhs => NULL() !< Temporary RHS vector
CLASS(oft_vector), POINTER :: tmp => NULL() !< Temporary storage vector
CLASS(oft_vector), POINTER :: aug_vec => NULL() !<  Augmented vector template (fields + coil currents)
TYPE(oft_tmaker_td_mfop), POINTER :: tkmr => NULL() !< TokaMaker time-dependent simulation object
TYPE(oft_xmhd_2d_sim), POINTER :: mug => NULL() !< MUG time-dependent simulation object
TYPE(oft_mf_matrix), POINTER :: mfmat => NULL() !< Matrix free Jacobian operator
TYPE(oft_mugtok_td_mfop), POINTER :: nlfun => NULL() ! !< Time-advance operator 
TYPE(oft_native_gmres_solver), POINTER :: mf_solver => NULL() !< Outer linear solver
CLASS(oft_solver), POINTER :: pre => NULL() !< Preconditioner for the approximate Jacobian (complete LU, or per-field block-Jacobi/LU in MHD-plasma mode)
TYPE(oft_nksolver) :: nksolver !< Newton-Krylov solver for time-advance
TYPE(oft_fem_comp_type), POINTER :: fe_rep => NULL() !< Finite element representation for solution fields (without coil currents)
LOGICAL :: pm = .TRUE.
LOGICAL :: nl_debug = .FALSE. !< Print per-field Newton diagnostics each nonlinear iteration ! NLDEBUG
LOGICAL :: save_rst = .FALSE. !< If true, write restart files ('mugtok_NNNNN.rst') for plotting
INTEGER(i4) :: rst_freq = 1 !< Restart-file output frequency (in steps)
INTEGER(i4) :: rst_count = 0 !< Completed-step counter (used as the restart-file index)
LOGICAL, ALLOCATABLE, DIMENSION(:) :: plasma_flag !< True for DOFs in or on the border of region 1 (the plasma)
REAL(r8) :: psi_limiting = 0.d0 !< Edge (LCFS) poloidal flux captured the first time the plasma is limited (not diverted); set once and never overwritten thereafter. See mugtok_update_psi_limiting/mugtok_apply_sink
LOGICAL :: psi_limiting_set = .FALSE. !< True once `psi_limiting` has been captured
INTEGER(i4) :: F0_node = 0 !< Index of the first plasma_flag=.TRUE. DOF, used to evolve F0 (0 if none)
REAL(r8) :: f_scale_prev = 1.d0 !< FF' scale (tkmr%f_scale) from previous timestep
CLASS(flux_func), POINTER :: I_prev => NULL() !< FF' profile (gs_equil%I) from previous timestep
LOGICAL :: mhd_plasma = .FALSE. !< True when the plasma (region 1) is evolved with MHD (set by setup_mhd)
LOGICAL :: fix_F_boundary = .TRUE. !< If true (default; set before setup_mhd), hold F = R B_phi at its initial value on the outer boundary of the mesh (Dirichlet; fixed TF current, as M3D-C1's iconst_bz=1) instead of the natural no-current condition, which conserves the total toroidal flux of the domain
INTEGER(i4) :: pre_freq = 1 !< Rebuild the approximate Jacobian and refactor the preconditioner every pre_freq steps (always on a dt change, on request via pre_force_update, and before halving dt after a Newton failure with a stale preconditioner)
LOGICAL :: lu_reanalyze = .FALSE. !< Redo the LU symbolic analysis (ordering) at every preconditioner rebuild. The Jacobian's sparsity pattern is fixed, so by default only the first factorization of the preconditioner does the analysis and later rebuilds refactor numerically (MUMPS job 2 instead of 4); .TRUE. restores the old behaviour
LOGICAL :: lu_analysed = .FALSE. !< The preconditioner's LU has had its symbolic analysis (first rebuild). Needed because pre%update(new_pattern) overwrites the solver's update_graph flag even before its first factorization, so the first rebuild must request the analysis explicitly
INTEGER(i4) :: pre_age = 0 !< Number of steps since the approximate Jacobian/preconditioner was last rebuilt
LOGICAL :: pre_force_update = .TRUE. !< Force a Jacobian/preconditioner rebuild on the next step (set after changing physics coefficients, e.g. nu, m_i, kappa)
REAL(r8) :: dt_refactor_tol = 0.d0 !< If > 0: a dt change within this relative distance of the dt of the current preconditioner (mug%jac_dt) only rebuilds the vacuum operator, the residual's one dt-dependent matrix, and keeps the preconditioner; larger changes rebuild it. 0 (default): every dt change rebuilds
LOGICAL :: pre_adapt = .FALSE. !< Cost-based preconditioner refresh (replaces pre_freq): rebuild once the last step's wall time per simulated second exceeds the average since the last rebuild (which includes the rebuild's cost)
INTEGER(i4) :: pre_max_age = 0 !< With pre_adapt: also rebuild after this many steps (0: no limit)
REAL(r8) :: pre_cyc_time = 0.d0 !< With pre_adapt: wall time [s] of the steps since the last rebuild
REAL(r8) :: pre_cyc_dt = 0.d0 !< With pre_adapt: simulated time [s] of the steps since the last rebuild
REAL(r8) :: pre_last_cost = 0.d0 !< With pre_adapt: wall time per simulated second of the last step
INTEGER(i4) :: pre_cyc_n = 0 !< With pre_adapt: steps since the last rebuild
LOGICAL :: dt_adapt = .FALSE. !< Adaptive dt: after each step set dt_next = lin_its_target*sum(dt)/sum(GMRES its) over the last dt_hist steps (xmhd.F90 scheme), growth <= dt_grow_max per step, decrease >= dt_shrink_min per step, capped at dt_max. The driver applies dt_next
INTEGER(i4) :: lin_its_target = 20 !< With dt_adapt: target GMRES iterations per step
INTEGER(i4) :: dt_hist = 5 !< With dt_adapt: number of steps in the dt/GMRES history (at most 50)
REAL(r8) :: dt_max = 0.d0 !< With dt_adapt: largest dt [s] (0: no cap)
REAL(r8) :: dt_grow_max = 1.25d0 !< With dt_adapt: largest dt increase factor per step
REAL(r8) :: dt_shrink_min = 0.5d0 !< With dt_adapt: smallest dt decrease factor per step
REAL(r8) :: dt_next = 0.d0 !< With dt_adapt: dt for the next step (output of step)
INTEGER(i4) :: dt_nhist = 0 !< With dt_adapt: steps in the history (set 0 to reset it, e.g. on a phase change)
REAL(r8) :: dt_hist_dt(50) = 0.d0 !< With dt_adapt: dt of the last dt_nhist steps (most recent first)
INTEGER(i4) :: dt_hist_its(50) = 0 !< With dt_adapt: GMRES iterations of the last dt_nhist steps (most recent first)
REAL(r8) :: T_clip = -1.d0 !< If > 0, hard temperature floor [eV]: after each accepted step, T DOFs (field 5) below T_clip are set to T_clip. Not energy conserving and outside the Newton solve; off by default
LOGICAL :: resid_dump = .FALSE. !< If true, write the Newton residual F(u)-rhs of each failed Newton attempt (fields 1-7, plus the iterate's V and T) as a time step of the XDMF group "mugtok_resid" in resid_dump/oft_xdmf.*.h5 (at most 8 per run)
INTEGER(i4) :: n_resid_dumps = 0 !< Residual dumps written so far
LOGICAL :: resid_xdmf_ready = .FALSE. !< mugtok_resid XDMF output initialized
TYPE(xdmf_plot_file) :: resid_xdmf !< XDMF output for the residual dumps
LOGICAL :: use_full_lu = .FALSE. !< In MHD-plasma mode, use a complete LU of the coupled Jacobian instead of per-field block-Jacobi (better linear convergence in the stiff near-ideal regime, but much higher memory). Set before calling setup_mhd.
INTEGER(i4) :: pre_nlocal = 0 !< If > 0: spatial block-Jacobi preconditioner instead of the LU of the whole Jacobian: the Jacobian graph is split into pre_nlocal subdomains (METIS) and each is factorized with its own LU, the coupling between subdomains being dropped (takes precedence over use_full_lu; 0 = off). Set before calling setup_mhd
LOGICAL :: pre_overlap = .FALSE. !< With pre_nlocal > 0: extend each subdomain by one layer of neighbouring unknowns (additive Schwarz, oft_bjprecond loverlap)
INTEGER(i4) :: pre_part = 0 !< Partition for the spatial block-Jacobi (pre_nlocal > 0): 0 = METIS on the Jacobian graph (pre_nlocal subdomains); 1 = geometric, 2 subdomains: the unknowns of fields 1-7 belonging to plasma-region (region 1) cells vs everything else, so the dense couplings (outer-boundary psi, coil currents) all stay in one subdomain and the cut is at the plasma boundary
REAL(r8) :: lin_rtol = 1.d-14 !< Relative tolerance of the linear (GMRES) solve inside Newton: it stops when |r| < lin_tol (absolute) or |r|/|r0| < lin_rtol, whichever comes first (inexact Newton). The default leaves only the absolute tolerance active; e.g. 1e-6 lets an inexact preconditioner (spatial block-Jacobi) stop once the linear residual is a small fraction of the current |F|. Set before calling setup_mhd
LOGICAL :: pre_parts_reported = .FALSE. !< One-time report of the spatial block-Jacobi subdomains (sizes by field) done
CLASS(flux_func), POINTER :: n_prof => NULL() !< Normalized density profile nhat(psihat), evaluated via %fp (linterp values). If set before setup_mhd (together with T_prof), n and T are seeded from the profiles and n is evolved
CLASS(flux_func), POINTER :: T_prof => NULL() !< Normalized temperature profile That(psihat), evaluated via %fp (linterp values). See n_prof
REAL(r8) :: n_scale = 1.d20 !< Density scale [m^-3]: n = n_offset + n_scale*nhat inside the plasma (also used as mug%den_scale)
REAL(r8) :: n_offset = 0.d0 !< Density offset [m^-3]: n = n_offset outside the plasma
REAL(r8) :: T_scale = 1.d3 !< Temperature scale [eV]: T = T_offset + T_scale*That inside the plasma
REAL(r8) :: T_offset = 0.d0 !< Temperature offset [eV]: T = T_offset outside the plasma
LOGICAL :: evolve_n = .TRUE. !< When n_prof is set, evolve the density with the continuity equation; if false, n is held fixed at its profile-initialized values

contains
    !> Setup time-dependent simulation (Grad-Shafranov plasma + MHD conductors)
    procedure :: setup => setup_mugtok_td
    !> Setup time-dependent simulation with the plasma (region 1) evolved as a compressible MHD region
    procedure :: setup_mhd => setup_mugtok_mhd_td
    !> Delete time-dependent simulation
    procedure :: delete => delete_mugtok_td
    !> Take a timestep
    procedure :: step => step_mugtok_td
    !> Restore the state from an rst file to resume a run
    procedure :: restart => mugtok_rst_restore
    !> Plot from rst files
    procedure :: plot => mugtok_td_plot
END TYPE oft_mugtok_td

!---------------------------------------------------------
!> Combined MUG and TokaMaker time-advance operator
!---------------------------------------------------------
type, extends(oft_noop_matrix) ::oft_mugtok_td_mfop
    real(r8) :: dt = 1.E-3 !< Time step size [s]
    CLASS(oft_matrix), POINTER :: jac_op => NULL() !<Approximate Jacobian operator for preconditioning
    TYPE(oft_mugtok_td), pointer :: parent_sim => NULL() !< Pointer to parent simulation object for access to parameters
contains
    !> Apply the nonlinear function
    procedure :: apply_real => nlfun_apply
    !> Delete time-advance operator
    procedure :: delete => delete_mfop
end type oft_mugtok_td_mfop

TYPE(oft_mugtok_td), POINTER :: current_sim => NULL() !<Current mug/tokamaker simulation object
CLASS(oft_vector), POINTER :: dbg_uprev => NULL() !< Previous Newton iterate ! NLDEBUG
CLASS(oft_vector), POINTER :: dbg_wrk => NULL() !< Work vector ! NLDEBUG
INTEGER(i4) :: dbg_it = 0 !< Newton iteration counter for diagnostics ! NLDEBUG
CLASS(oft_bmesh), POINTER, PUBLIC :: mesh => NULL() !< Pointer to FE mesh
CLASS(oft_scalar_bfem), POINTER :: lag_rep => NULL()!< Pointer to single field FE representation 

CONTAINS

!---------------------------------------------------------
!> Setup combined MUG and TokaMaker simulation object
!---------------------------------------------------------
subroutine setup_mugtok_td(self, equil, dt,lin_tol,nl_tol, mhd_flag, dens_reg, visc_reg, eta_reg, incomp, toroidal_flow)
CLASS(oft_mugtok_td), INTENT(inout), TARGET :: self !< Simulation object to be setup
TYPE(gs_equil), INTENT(inout), TARGET :: equil !< Tokamaker equilibrium object
REAL(8), INTENT(in) :: dt !< Desired timestep [s]
REAL(8), INTENT(in) :: lin_tol !< Linear solver tolerance
REAL(8), INTENT(in) :: nl_tol !< Non-linear solver tolerance
LOGICAL, INTENT(in) :: mhd_flag(:) !< True for regions where MHD should be used, false for TokaMaker
REAL(8), INTENT(in) :: dens_reg(:) !< Density [kg/m^3] in each region [unused where mhd_flag is false]
REAL(8), INTENT(in) :: visc_reg(:) !< Viscosity [m^2/s] in each region [unused where mhd_flag is false]
REAL(8), INTENT(in), optional :: eta_reg(:,:)!< Electrical resistivity [in plane, out of plane], in units of mu0, in each region. Defaults to TokaMaker region resistivites
LOGICAL, INTENT(in), optional :: incomp !< If true, use incompressible MHD model (false currently not supported)
LOGICAL, INTENT(in), optional :: toroidal_flow !< Allow toroidal (phi) flow in MHD regions 
CLASS(multigrid_mesh), POINTER :: mg_mesh 
LOGICAL :: tor_flow
INTEGER(i4) :: i,j,ierr
TYPE(oft_xmhd_2d_sim), POINTER :: mhd_sim
REAL(r8), POINTER, DIMENSION(:) :: tmp_arr, vals_out
TYPE(oft_graph_ptr), ALLOCATABLE :: graphs(:,:), fe_graphs(:,:), known_graphs(:)
INTEGER(i4) :: nkgraphs
TYPE(oft_graph), TARGET :: dense_graph, dense_graph7
type(oft_1d_int), pointer, dimension(:) :: bc_nodes, f0_nodes
integer(i4), allocatable :: dense_flag(:)
LOGICAL, ALLOCATABLE :: p_dir_set(:)
INTEGER(i4), POINTER, DIMENSION(:) :: cell_dofs, cell_dofs_p
type(seam_list), pointer, dimension(:) :: stitch_tmp
type(map_list), pointer, dimension(:) :: map_tmp
CLASS(oft_vector), pointer :: tmp_vec

!--Take mesh and single field FE reps from TokaMaker equilibrium
mesh => equil%device%mesh
lag_rep=>equil%device%fe_rep
mg_mesh => equil%device%ML_fe_rep%ml_mesh

!------------------------------------------------------------------------------
! Set up TokaMaker time-dependent object
!------------------------------------------------------------------------------
ALLOCATE(self%tkmr)
!--Build the diagnostic dels_full operator over the FULL mesh, before
!  masking the MHD regions below
IF(.NOT.ASSOCIATED(equil%device%dels_full)) &
    CALL build_dels(equil%device%dels_full, equil%device, "none")
equil%device%ignore_rmask = mhd_flag !ignore regions where MHD will be used

self%tkmr%dt=dt
CALL self%tkmr%setup(equil)
CALL build_vac_op(self%tkmr,self%tkmr%vac_op)

!--Current implementation only supports I coils, throw error if V coils present
IF(self%tkmr%gs_device%ncoils > 0)THEN
    DO i=1,self%tkmr%gs_device%ncoils
        IF(self%tkmr%gs_device%Rcoils(i) > 0.d0) &
            CALL oft_abort("V coils (Rcoils>0) are not supported; all coils must be prescribed I coils", &
                           "setup_mugtok_td", __FILE__)
    END DO
END IF

!------------------------------------------------------------------------------
! Set up MUG time-dependent object
!------------------------------------------------------------------------------
ALLOCATE(mhd_sim)
ALLOCATE(mhd_sim%ignore_rmask(mesh%nreg))
mhd_sim%ignore_rmask = .NOT. mhd_flag !ignore regions where MHD will not be used
ALLOCATE(mhd_sim%eta(mesh%nreg, 2))
ALLOCATE(mhd_sim%m_i(mesh%nreg))
ALLOCATE(mhd_sim%nu(mesh%nreg))
mhd_sim%m_i = dens_reg
mhd_sim%nu = visc_reg
IF(PRESENT(eta_reg)) THEN
   mhd_sim%eta = eta_reg
ELSE
  !Default to TokaMaker region resistivity if none are provided
  mhd_sim%eta(:,1) = self%tkmr%eta_reg 
  mhd_sim%eta(:,2) = self%tkmr%eta_reg
END IF

IF (PRESENT(incomp)) THEN
    mhd_sim%incomp = incomp
ELSE
    mhd_sim%incomp = .TRUE.
END IF
!Throw error if compressible, not currently supported/validated
IF (.NOT. mhd_sim%incomp) CALL oft_abort( &
    'Compressible flow (incomp=.FALSE.) is not currently supported; use incomp=.TRUE.', &
    'setup_mugtok_td', __FILE__)

tor_flow = .FALSE.
IF (PRESENT(toroidal_flow)) tor_flow = toroidal_flow
mhd_sim%cyl_flag = .TRUE.
mhd_sim%dt = dt
mhd_sim%den_scale = 1.d0

CALL mhd_sim%setup(mg_mesh, lag_rep%order, fe_rep_in =lag_rep)
mhd_sim%sink_func => mugtok_apply_sink

!------------------------------------------------------------------------------
! Set up plasma flag and F0 node
!------------------------------------------------------------------------------
! Flag DOFs that lie in (or on the border of) region 1, the plasma
ALLOCATE(self%plasma_flag(mhd_sim%fe_rep%fields(7)%fe%ne))
self%plasma_flag = .FALSE.
ALLOCATE(cell_dofs(lag_rep%nce))
DO i = 1, mesh%nc
    IF (mesh%reg(i) == 1) THEN
        CALL lag_rep%ncdofs(i, cell_dofs) ! By DOFs (order-2)
        DO j = 1, SIZE(cell_dofs)
            self%plasma_flag(cell_dofs(j)) = .TRUE.
        END DO
    END IF
END DO
DEALLOCATE(cell_dofs)
!---Index of the first plasma DOF, which will be used to evolve F0
self%F0_node = FINDLOC(self%plasma_flag, .TRUE., DIM=1)

!------------------------------------------------------------------------------
! Set boundary conditions in MUG solve
!------------------------------------------------------------------------------
IF (ASSOCIATED(mhd_sim%n_bc)) NULLIFY(mhd_sim%n_bc)
IF (ASSOCIATED(mhd_sim%velx_bc)) NULLIFY(mhd_sim%velx_bc)
IF (ASSOCIATED(mhd_sim%vely_bc)) NULLIFY(mhd_sim%vely_bc)
IF (ASSOCIATED(mhd_sim%velz_bc)) NULLIFY(mhd_sim%velz_bc)
IF (ASSOCIATED(mhd_sim%T_bc)) NULLIFY(mhd_sim%T_bc)
IF (ASSOCIATED(mhd_sim%psi_bc)) NULLIFY(mhd_sim%psi_bc)
IF (ASSOCIATED(mhd_sim%by_bc)) NULLIFY (mhd_sim%by_bc)

ALLOCATE(mhd_sim%n_bc(mhd_sim%fe_rep%fields(1)%fe%ne))
ALLOCATE(mhd_sim%velx_bc(mhd_sim%fe_rep%fields(2)%fe%ne))
ALLOCATE(mhd_sim%vely_bc(mhd_sim%fe_rep%fields(3)%fe%ne))
ALLOCATE(mhd_sim%velz_bc(mhd_sim%fe_rep%fields(4)%fe%ne))
ALLOCATE(mhd_sim%T_bc(mhd_sim%fe_rep%fields(5)%fe%ne))
ALLOCATE(mhd_sim%psi_bc(mhd_sim%fe_rep%fields(6)%fe%ne))
ALLOCATE(mhd_sim%by_bc(mhd_sim%fe_rep%fields(7)%fe%ne))

mhd_sim%psi_bc = .FALSE. 
mhd_sim%by_bc = .FALSE. 
mhd_sim%velx_bc = .FALSE.
!Permit or restrict flow in toroidal direction
IF (tor_flow) THEN
  mhd_sim%vely_bc = .FALSE.
ELSE
  mhd_sim%vely_bc = .TRUE.
END IF
mhd_sim%velz_bc = .FALSE.
mhd_sim%n_bc = .TRUE.
mhd_sim%T_bc = .TRUE.


ALLOCATE(cell_dofs(lag_rep%nce))
ALLOCATE(cell_dofs_p(mhd_sim%fe_rep%fields(5)%fe%nce)) 
DO i = 1, mesh%nc
  call mhd_sim%fe_rep%fields(5)%fe%ncdofs(i,cell_dofs_p) ! pressure DOFs (order-1)
  call lag_rep%ncdofs(i,cell_dofs) ! velocity/psi/F DOFs
  DO j=1, SIZE(cell_dofs)
    IF (mhd_flag(mesh%reg(i))) THEN
      IF (.NOT. mhd_sim%incomp) mhd_sim%n_bc(cell_dofs(j)) = .FALSE. !Evolve density only in MHD regions (if not incompressible)
      ELSE
        mhd_sim%velx_bc(cell_dofs(j)) = .TRUE. !Pin velocity to 0 everywhere outside of MHD regions
        mhd_sim%vely_bc(cell_dofs(j)) = .TRUE.
        mhd_sim%velz_bc(cell_dofs(j)) = .TRUE.
      END IF
  END DO
  IF (mhd_flag(mesh%reg(i))) THEN
      DO j=1, SIZE(cell_dofs_p)
        mhd_sim%T_bc(cell_dofs_p(j)) = .FALSE. !Evolve pressure (or temperature) only in MHD regions
      END DO
  END IF
END DO

!If incompressible pin value of pressure at a single DOF in each MHD region to remove null space
IF (mhd_sim%incomp) THEN
    ALLOCATE(p_dir_set(mesh%nreg))
    p_dir_set = .FALSE.
    DO i=1, mesh%nc
        IF (mhd_flag(mesh%reg(i)) .AND. .NOT. p_dir_set(mesh%reg(i))) THEN
            call mhd_sim%fe_rep%fields(5)%fe%ncdofs(i,cell_dofs_p) 
            mhd_sim%T_bc(cell_dofs_p(1)) = .TRUE.
            p_dir_set(mesh%reg(i)) = .TRUE.
        END IF
    END DO
    DEALLOCATE(p_dir_set)
END IF
DEALLOCATE(cell_dofs, cell_dofs_p)

self%mug => mhd_sim

!------------------------------------------------------------------------------------
! Create solver fields, augmented with coil currents for compatibility with TokaMaker
!------------------------------------------------------------------------------------
self%fe_rep => self%mug%fe_rep

!---Create augmented vector
IF (self%tkmr%gs_device%ncoils > 0) THEN
    ALLOCATE(stitch_tmp(8),map_tmp(8))
    stitch_tmp(1)%s=>self%fe_rep%fields(1)%fe%linkage
    map_tmp(1)%m=>self%fe_rep%fields(1)%fe%map
    stitch_tmp(2)%s=>self%fe_rep%fields(2)%fe%linkage
    map_tmp(2)%m=>self%fe_rep%fields(2)%fe%map
    stitch_tmp(3)%s=>self%fe_rep%fields(3)%fe%linkage
    map_tmp(3)%m=>self%fe_rep%fields(3)%fe%map
    stitch_tmp(4)%s=>self%fe_rep%fields(4)%fe%linkage
    map_tmp(4)%m=>self%fe_rep%fields(4)%fe%map
    stitch_tmp(5)%s=>self%fe_rep%fields(5)%fe%linkage
    map_tmp(5)%m=>self%fe_rep%fields(5)%fe%map
    stitch_tmp(6)%s=>self%fe_rep%fields(6)%fe%linkage
    map_tmp(6)%m=>self%fe_rep%fields(6)%fe%map
    stitch_tmp(7)%s=>self%fe_rep%fields(7)%fe%linkage
    map_tmp(7)%m=>self%fe_rep%fields(7)%fe%map
    stitch_tmp(8)%s=>self%tkmr%gs_device%coil_stitch
    map_tmp(8)%m=>self%tkmr%gs_device%coil_map
    CALL create_vector(self%aug_vec,stitch_tmp,map_tmp)
    DEALLOCATE(stitch_tmp,map_tmp)
    CALL self%aug_vec%new(self%u)
    CALL self%aug_vec%new(self%rhs)
    CALL self%aug_vec%new(self%tmp)
ELSE
    CALL self%fe_rep%vec_create(self%u)
    call self%fe_rep%vec_create(self%rhs)
    call self%fe_rep%vec_create(self%tmp)
END IF


!------------------------------------------------------------------------------
! Set initial field values
!------------------------------------------------------------------------------
CALL self%u%set(1.d0, 1)
CALL self%u%set(0.d0, 2)
CALL self%u%set(0.d0, 3)
CALL self%u%set(0.d0, 4)
CALL self%u%set(0.d0, 5)
CALL self%u%set(self%tkmr%gs_equil%I%f_offset, 7) !Set initial F to F0 value from TokaMaker equilibrium

!Initialize field 6 (plasma poloidal flux) and field 8 (coil currents) from the
!TokaMaker equilibrium, and populate MUG's vacuum flux (mug%psi_vac) from the coils.
NULLIFY(tmp_arr, vals_out)
CALL self%tkmr%gs_device%fe_rep%vec_create(tmp_vec)
CALL tmp_vec%add(0.d0,1.d0,self%tkmr%gs_equil%psi) ! tmp_vec = total poloidal flux
IF (self%tkmr%gs_device%ncoils > 0) THEN
    CALL self%u%get_local(vals_out,8)
    DO i=1,self%tkmr%gs_device%ncoils
        !---Subtract the coil (vacuum) flux, leaving the plasma poloidal flux for field 6
        CALL tmp_vec%add(1.d0,-self%tkmr%gs_equil%coil_currs(i),self%tkmr%gs_device%psi_coil(i)%f)
        vals_out(i)=self%tkmr%gs_equil%coil_currs(i)
    END DO
    CALL self%u%restore_local(vals_out,8)
    DEALLOCATE(vals_out)
END IF
CALL tmp_vec%get_local(tmp_arr)
CALL self%u%restore_local(tmp_arr,6)

CALL build_psi_vac(self, self%tkmr%gs_equil%coil_currs)

!---Free the scratch total-flux vector used to build the initial condition
CALL tmp_vec%delete
DEALLOCATE(tmp_vec)

!------------------------------------------------------------------------------
! Build the nonlinear function, approximate Jacobian, preconditioner and solvers
!------------------------------------------------------------------------------
CALL mugtok_finalize_setup(self, dt, lin_tol, nl_tol)
end subroutine setup_mugtok_td

!---------------------------------------------------------------------------
!> Build the nonlinear function object, approximate-Jacobian structure,
!! preconditioner and solvers. Shared by both setup entry points; the F0
!! dense-coupling block is included only when F0_node>0 (Grad-Shafranov plasma).
!---------------------------------------------------------------------------
subroutine mugtok_finalize_setup(self, dt, lin_tol, nl_tol)
class(oft_mugtok_td), intent(inout), target :: self !< Simulation object
real(r8), intent(in) :: dt !< Timestep [s]
real(r8), intent(in) :: lin_tol !< Linear solver tolerance
real(r8), intent(in) :: nl_tol !< Nonlinear solver tolerance
TYPE(oft_graph_ptr), ALLOCATABLE :: graphs(:,:), fe_graphs(:,:), known_graphs(:)
INTEGER(i4) :: nkgraphs, i, j
TYPE(oft_graph), TARGET :: dense_graph, dense_graph7
type(oft_1d_int), pointer, dimension(:) :: bc_nodes, f0_nodes
integer(i4), allocatable :: dense_flag(:)
class(oft_vector), pointer :: tmp_vec

ALLOCATE(self%nlfun)
self%nlfun%dt=dt
self%nlfun%parent_sim => self

!---Scratch vector in the TokaMaker (device) space, used as a template for the
!   coil dense-coupling graphs
CALL self%tkmr%gs_device%fe_rep%vec_create(tmp_vec)

!------------------------------------------------------------------------------
! Setup structure of approximate Jacobian matrix used for preconditioning
!------------------------------------------------------------------------------
CALL fem_graph_create(self%fe_rep, fe_graphs, known_graphs, nkgraphs) !Setup graphs for base finite element fields

!Allocate overall matrix of graphs
IF (self%tkmr%gs_device%ncoils > 0) THEN
    ALLOCATE(graphs(self%fe_rep%nfields+1,self%fe_rep%nfields+1))
ELSE
    ALLOCATE(graphs(self%fe_rep%nfields,self%fe_rep%nfields))
END IF

!Populate graphs with the finite element graphs
DO i = 1, self%fe_rep%nfields
    DO j = 1, self%fe_rep%nfields
        ALLOCATE(graphs(i,j)%g)
        graphs(i,j)%g%nnz=fe_graphs(i,j)%g%nnz
        graphs(i,j)%g%nr=fe_graphs(i,j)%g%nr
        graphs(i,j)%g%nrg=fe_graphs(i,j)%g%nrg
        graphs(i,j)%g%nc=fe_graphs(i,j)%g%nc
        graphs(i,j)%g%ncg=fe_graphs(i,j)%g%ncg
        graphs(i,j)%g%kr=>fe_graphs(i,j)%g%kr
        graphs(i,j)%g%lc=>fe_graphs(i,j)%g%lc
    END DO
END DO

! Add dense regions to psi(6,6) block for free-boundary
ALLOCATE(bc_nodes(1))
bc_nodes(1)%n = lag_rep%nbe
bc_nodes(1)%v => lag_rep%lbe
ALLOCATE(dense_flag(lag_rep%ne))
dense_flag = 0
dense_flag(bc_nodes(1)%v) = 1
CALL graph_add_dense_blocks(graphs(6,6)%g,dense_graph,dense_flag,bc_nodes)
NULLIFY(graphs(6,6)%g%kr,graphs(6,6)%g%lc)
graphs(6,6)%g%nnz=dense_graph%nnz
graphs(6,6)%g%kr=>dense_graph%kr
graphs(6,6)%g%lc=>dense_graph%lc
DEALLOCATE(dense_flag, bc_nodes)

!---Add F0 coupling to the F(7,7) block: every plasma DOF is coupled to the F0 node
IF(self%F0_node > 0)THEN
  ALLOCATE(f0_nodes(1))
  f0_nodes(1)%n = 1
  ALLOCATE(f0_nodes(1)%v(1))
  f0_nodes(1)%v(1) = self%F0_node
  ALLOCATE(dense_flag(lag_rep%ne))
  dense_flag = MERGE(1, 0, self%plasma_flag)
  CALL graph_add_dense_blocks(graphs(7,7)%g,dense_graph7,dense_flag,f0_nodes)
  NULLIFY(graphs(7,7)%g%kr,graphs(7,7)%g%lc)
  graphs(7,7)%g%nnz=dense_graph7%nnz
  graphs(7,7)%g%kr=>dense_graph7%kr
  graphs(7,7)%g%lc=>dense_graph7%lc
  DEALLOCATE(dense_flag, f0_nodes(1)%v)
  DEALLOCATE(f0_nodes)
END IF

!--- Add dense blocks for coil coupling if coils are present
IF(self%tkmr%gs_device%ncoils > 0)THEN
    CALL create_dense_graph(graphs(8,6)%g,self%tkmr%gs_device%coil_vec,tmp_vec)
    CALL create_dense_graph(graphs(6,8)%g,tmp_vec,self%tkmr%gs_device%coil_vec)
    CALL create_dense_graph(graphs(8,8)%g,self%tkmr%gs_device%coil_vec,self%tkmr%gs_device%coil_vec)
END IF

!Create jacobian matrix
CALL create_matrix(self%nlfun%jac_op,graphs,self%aug_vec,self%aug_vec)

!---Geometric partition for the spatial block-Jacobi (pre_part=1): row coloring of the Jacobian
!   (as xmhd.F90), color 2 = unknowns of fields 1-7 on region-1 cells, color 1 = the rest incl. coil currents
IF((self%pre_nlocal > 0).AND.(self%pre_part == 1))THEN
    SELECT TYPE(this => self%nlfun%jac_op)
    CLASS IS(oft_native_matrix)
        ALLOCATE(this%color(this%nr))
        this%color = 1
        j = 0
        DO i=1,self%fe_rep%nfields
            !---Unknowns of this field that belong to a plasma-region (region 1) cell, on the field's own space
            !   (plasma_flag is not used: it is only allocated by setup_mugtok_td, not in MHD-plasma mode)
            BLOCK
            INTEGER(i4) :: ic, kd
            INTEGER(i4), ALLOCATABLE :: cdofs(:)
            ALLOCATE(cdofs(self%fe_rep%fields(i)%fe%nce))
            DO ic=1,mesh%nc
                IF(mesh%reg(ic) /= 1)CYCLE
                CALL self%fe_rep%fields(i)%fe%ncdofs(ic,cdofs)
                DO kd=1,SIZE(cdofs)
                    this%color(j+cdofs(kd)) = 2
                END DO
            END DO
            DEALLOCATE(cdofs)
            END BLOCK
            j = j + self%fe_rep%fields(i)%fe%ne
        END DO
        WRITE(*,'(A,I9,A,I9,A)')' Block-Jacobi geometric partition (pre_part=1):',COUNT(this%color==2), &
            ' plasma-region unknowns,',COUNT(this%color==1),' others (walls, vacuum, boundary, coils)'
    CLASS DEFAULT
        CALL oft_abort('pre_part=1 needs a native Jacobian matrix','mugtok_finalize_setup',__FILE__)
    END SELECT
END IF

DO i=1,nkgraphs
    DEALLOCATE(known_graphs(i)%g)
END DO
DEALLOCATE(graphs, known_graphs)
CALL tmp_vec%delete
DEALLOCATE(tmp_vec)
DEALLOCATE(fe_graphs)

!------------------------------------------------------------------------------
! Preconditioner for the approximate Jacobian (cf. xmhd_2d):
!  - MHD-plasma mode (default): per-field block-Jacobi with an LU factorization
!    of each field block. The monolithic LU of the coupled compressible system is
!    memory-hungry (dense free-boundary/coil blocks). Set use_full_lu=.TRUE. to
!    force a complete LU instead -- needed when the block preconditioner's linear
!    convergence is too poor in the stiff, near-ideal (low physical eta) regime.
!  - incompressible / Grad-Shafranov plasma: always a complete LU.
!  - pre_nlocal > 0 (any mode): spatial block-Jacobi, one LU per METIS subdomain
!    (optionally overlapping by one layer, pre_overlap).
!------------------------------------------------------------------------------
IF(self%pre_nlocal > 0)THEN
    CALL create_bjacobi_pre(self%pre, self%pre_nlocal) ! one LU per spatial subdomain (METIS partition of the Jacobian graph)
    SELECT TYPE(bj => self%pre)
    CLASS IS(oft_bjprecond)
        bj%loverlap = self%pre_overlap
    END SELECT
ELSE IF(self%mhd_plasma .AND. .NOT.self%use_full_lu)THEN
    CALL create_bjacobi_pre(self%pre, -1) ! one LU block per field block
ELSE
    CALL create_native_pre(self%pre, "lu")
END IF
self%pre%A=>self%nlfun%jac_op

!------------------------------------------------------------------------------
! Setup matrix free solver
!------------------------------------------------------------------------------
ALLOCATE(self%mfmat)
CALL self%mfmat%setup(self%tmp,self%nlfun)

ALLOCATE(self%mf_solver)
self%mfmat%b0=1.d-3
self%mf_solver%A=>self%mfmat
self%mf_solver%its=250
self%mf_solver%nrits=40
self%mf_solver%atol=lin_tol
self%mf_solver%rtol=self%lin_rtol
self%mf_solver%itplot=1
oft_env%pm = self%pm
self%mf_solver%pm=oft_env%pm
self%mf_solver%pre=>self%pre
!------------------------------------------------------------------------------
! Setup Newton Solver
!------------------------------------------------------------------------------
self%nksolver%A=>self%nlfun
self%nksolver%J_inv=>self%mf_solver
self%nksolver%its=20
self%nksolver%atol=nl_tol
self%nksolver%rtol=1.d-20 ! Disable relative tolerance
self%nksolver%backtrack=.FALSE.
self%nksolver%J_update=>mugtok_mfnk_update
self%nksolver%up_freq=1

!---Save initial FF' profile for F RHS
CALL snapshot_f_profile(self)

!---Write an initial restart file capturing the initial state
IF(self%save_rst) CALL mugtok_rst_save(self, 0.d0)
end subroutine mugtok_finalize_setup

!---------------------------------------------------------------------------
!> Setup a combined simulation in which the plasma (region 1) is evolved as a
!! compressible MHD region rather than by Grad-Shafranov. By default density is
!! a flat, pinned field (den_scale) and temperature is initialized from the
!! equilibrium pressure profile, T(psi)=P(psi)/(2*k*n). If `n_prof` and `T_prof`
!! are set before this call, n and T are instead initialized from those profiles
!! (n = n_offset + n_scale*nhat, T = T_offset + T_scale*That inside the plasma,
!! the offsets outside) and density is evolved with the continuity equation in
!! the MHD regions (pinned outside), unless `evolve_n=.FALSE.` (n held fixed). Temperature is evolved with the MHD
!! temperature equation in the MHD regions (pinned at their boundary and
!! outside); set `evolve_T=.FALSE.` to instead hold T to the T(psi) closure
!! evaluated pointwise inside MUG (see mugtok_T_of_psi / xmhd_2d T_func). The
!! poloidal flux (psi) and toroidal field function (F=R*B_phi) evolve with full
!! MHD; F is initialized from the equilibrium profile.
!---------------------------------------------------------------------------
subroutine setup_mugtok_mhd_td(self, equil, dt, lin_tol, nl_tol, mhd_flag, dens_reg, visc_reg, eta_reg, toroidal_flow, mass_scale, use_spitzer, kappa_reg, evolve_T, D_reg, evolve_pressure)
CLASS(oft_mugtok_td), INTENT(inout), TARGET :: self !< Simulation object to be setup
TYPE(gs_equil), INTENT(inout), TARGET :: equil !< Tokamaker equilibrium object
REAL(8), INTENT(in) :: dt !< Desired timestep [s]
REAL(8), INTENT(in) :: lin_tol !< Linear solver tolerance
REAL(8), INTENT(in) :: nl_tol !< Non-linear solver tolerance
LOGICAL, INTENT(in) :: mhd_flag(:) !< True for regions where MHD should be used (must include region 1)
REAL(8), INTENT(in) :: dens_reg(:) !< Ion mass [kg] in each region (used as m_i for the compressible model)
REAL(8), INTENT(in) :: visc_reg(:) !< Viscosity [m^2/s] in each region
REAL(8), INTENT(in), optional :: eta_reg(:,:) !< Electrical resistivity [in plane, out of plane], in units of mu0. Defaults to TokaMaker region resistivities
LOGICAL, INTENT(in), optional :: toroidal_flow !< Allow toroidal (phi) flow in MHD regions (default .TRUE.)
REAL(8), INTENT(in), optional :: mass_scale !< Factor by which the plasma (region 1) ion mass is inflated for numerical acceleration (default 100)
LOGICAL, INTENT(in), optional :: use_spitzer !< Use Spitzer resistivity eta(n,T) in the MHD regions (default .TRUE.)
REAL(8), INTENT(in), optional :: kappa_reg(:) !< Thermal conductivity kappa=chi*n [1/(m s)] in each region, used when the temperature is evolved (defaults to 1 m^2/s * n_scale in the MHD regions and 0 elsewhere)
LOGICAL, INTENT(in), optional :: evolve_T !< Evolve the temperature with the MHD temperature equation (default .TRUE.). If false, T is held to the T(psi) closure
REAL(8), INTENT(in), optional :: D_reg(:) !< Density diffusivity [m^2/s] in each region, used when the density is evolved (defaults to 1e-6 m^2/s in the MHD regions and 0 elsewhere)
LOGICAL, INTENT(in), optional :: evolve_pressure !< Add the -(T/n)*D*lap(n) term to the temperature equation so density diffusion conserves pressure (default .FALSE.; see oft_xmhd_2d_sim%evolve_pressure)
CLASS(multigrid_mesh), POINTER :: mg_mesh
LOGICAL :: tor_flow, spitzer, evolve_temp, evolve_dens
REAL(r8) :: mscale
INTEGER(i4) :: i,j
TYPE(oft_xmhd_2d_sim), POINTER :: mhd_sim
REAL(r8), POINTER, DIMENSION(:) :: tmp_arr, vals_out
INTEGER(i4), POINTER, DIMENSION(:) :: cell_dofs
type(seam_list), pointer, dimension(:) :: stitch_tmp
type(map_list), pointer, dimension(:) :: map_tmp
CLASS(oft_vector), pointer :: tmp_vec

IF(.NOT.mhd_flag(1)) CALL oft_abort('setup_mhd requires mhd_flag(1)=.TRUE. (plasma region evolved with MHD)', &
                                    'setup_mugtok_mhd_td', __FILE__)
self%mhd_plasma = .TRUE.
current_sim => self
!---n/T are seeded from profiles when both profiles are supplied, and density is
!   then evolved unless evolve_n=.FALSE. (n held at its initial profile)
IF(ASSOCIATED(self%n_prof) .NEQV. ASSOCIATED(self%T_prof)) CALL oft_abort('n_prof and T_prof must be set together', &
                                                                        'setup_mugtok_mhd_td', __FILE__)
evolve_dens = ASSOCIATED(self%n_prof) .AND. self%evolve_n
!---n = n_offset outside the plasma, and the momentum/temperature equations divide by n
IF(ASSOCIATED(self%n_prof) .AND. (self%n_offset <= 0.d0)) CALL oft_abort('n_offset must be > 0 when n_prof is set', &
                                                             'setup_mugtok_mhd_td', __FILE__)

!--Take mesh and single field FE reps from TokaMaker equilibrium
mesh => equil%device%mesh
lag_rep => equil%device%fe_rep
mg_mesh => equil%device%ML_fe_rep%ml_mesh

!------------------------------------------------------------------------------
! Set up TokaMaker time-dependent object (region 1 ignored via ignore_rmask)
!------------------------------------------------------------------------------
ALLOCATE(self%tkmr)
IF(.NOT.ASSOCIATED(equil%device%dels_full)) &
    CALL build_dels(equil%device%dels_full, equil%device, "none")
equil%device%ignore_rmask = mhd_flag !ignore regions where MHD will be used
self%tkmr%dt=dt
CALL self%tkmr%setup(equil)
CALL build_vac_op(self%tkmr,self%tkmr%vac_op)

IF(self%tkmr%gs_device%ncoils > 0)THEN
    DO i=1,self%tkmr%gs_device%ncoils
        IF(self%tkmr%gs_device%Rcoils(i) > 0.d0) &
            CALL oft_abort("V coils (Rcoils>0) are not supported; all coils must be prescribed I coils", &
                           "setup_mugtok_mhd_td", __FILE__)
    END DO
END IF

!------------------------------------------------------------------------------
! Set up MUG time-dependent object (compressible)
!------------------------------------------------------------------------------
ALLOCATE(mhd_sim)
ALLOCATE(mhd_sim%ignore_rmask(mesh%nreg))
mhd_sim%ignore_rmask = .NOT. mhd_flag !ignore regions where MHD will not be used
ALLOCATE(mhd_sim%eta(mesh%nreg, 2))
ALLOCATE(mhd_sim%m_i(mesh%nreg))
ALLOCATE(mhd_sim%nu(mesh%nreg))
ALLOCATE(mhd_sim%gamma(mesh%nreg))
ALLOCATE(mhd_sim%chi(mesh%nreg))
ALLOCATE(mhd_sim%D_diff(mesh%nreg))
mhd_sim%m_i = dens_reg
!---Inflate the plasma (region 1) ion mass to slow the fast MHD waves (larger dt)
mscale = 400.d0
IF (PRESENT(mass_scale)) mscale = mass_scale
mhd_sim%m_i(1) = mhd_sim%m_i(1)*mscale
mhd_sim%nu = visc_reg
mhd_sim%nu2_ratio = 2.d0 !< Bulk (compressive) viscosity = 2x shear viscosity, to damp compressive noise
mhd_sim%drop_convection = .TRUE. !< Drop the (v.grad)v convective term from the momentum equation
mhd_sim%gamma = 5.d0/3.d0
!---Density diffusivity: defaults to 1e-6 m^2/s in the MHD regions when density
!   is evolved (unused when density is a flat, pinned field)
mhd_sim%D_diff = 0.d0
IF (PRESENT(D_reg)) THEN
  mhd_sim%D_diff = D_reg
ELSE IF (evolve_dens) THEN
  WHERE (mhd_flag) mhd_sim%D_diff = 1.d-6
END IF
mhd_sim%evolve_pressure = .FALSE.
IF (PRESENT(evolve_pressure)) mhd_sim%evolve_pressure = evolve_pressure
!---Thermal conductivity kappa=chi*n for the evolved temperature (xmhd_2d then
!   uses chi=kappa/n pointwise). The default (chi=1 m^2/s at n=n_scale in the
!   MHD regions) gives a diffusion time L^2/chi ~ 1 s, long compared to the
!   ms-scale transients of interest, so the initial T(psi) profile is retained
!   while advection of T stays smooth on the mesh scale.
evolve_temp = .TRUE.
IF (PRESENT(evolve_T)) evolve_temp = evolve_T
mhd_sim%chi = 0.d0 !unused: replaced by kappa, or T is prescribed via the T(psi) closure
IF (evolve_temp) THEN
  ALLOCATE(mhd_sim%kappa(mesh%nreg))
  IF (PRESENT(kappa_reg)) THEN
    mhd_sim%kappa = kappa_reg
  ELSE
    mhd_sim%kappa = 0.d0
    WHERE (mhd_flag) mhd_sim%kappa = 1.d0*self%n_scale
  END IF
END IF
!---Anisotropic thermal conduction: kappa(:) above is the cross-field conductivity,
!   and field-aligned transport is faster by chi_par_mult (default 1e3)
mhd_sim%chi_aniso = .TRUE.
IF(PRESENT(eta_reg)) THEN
   mhd_sim%eta = eta_reg
ELSE
  mhd_sim%eta(:,1) = self%tkmr%eta_reg
  mhd_sim%eta(:,2) = self%tkmr%eta_reg
END IF

!---Spitzer resistivity eta(n,T) by default, using the evolved temperature (or the
!   T(psi) closure temperature when the temperature is not evolved)
spitzer = .TRUE.
IF (PRESENT(use_spitzer)) spitzer = use_spitzer
mhd_sim%use_spitzer = spitzer

tor_flow = .TRUE. !Toroidal flow evolved in MHD regions by default
IF (PRESENT(toroidal_flow)) tor_flow = toroidal_flow
mhd_sim%cyl_flag = .TRUE.
mhd_sim%incomp = .FALSE. !compressible model (so pressure = 2*k*n*T enters momentum)
mhd_sim%dt = dt
!Density field is normalized by den_scale (flat at 1.0 when not evolved, so
!den_scale is then the physical number density)
mhd_sim%den_scale = 1.d20
IF (ASSOCIATED(self%n_prof)) mhd_sim%den_scale = self%n_scale

CALL mhd_sim%setup(mg_mesh, lag_rep%order, fe_rep_in=lag_rep)
mhd_sim%sink_func => mugtok_apply_sink
!---Temperature sink term for DOFs left behind as the plasma boundary recedes
mhd_sim%tau_sink = 1.d-3 ! 1 ms relaxation timescale
mhd_sim%T_sink = 5.d0 ! [eV]

!---When the temperature is not evolved, it is a prescribed function of the local
!   fields (T = T(psi)) evaluated pointwise in MUG instead of reconstructed from
!   the field-5 DOFs. When it is evolved, field 5 carries T and the closure is
!   used only to set the initial condition (see mugtok_seed_fields).
IF (.NOT.evolve_temp) mhd_sim%T_func => mugtok_T_of_psi

!------------------------------------------------------------------------------
! Boundary conditions: psi and F evolve; velocity, (when evolved) temperature,
! and (when evolved) density evolve in the MHD regions and are pinned elsewhere
!------------------------------------------------------------------------------
IF (ASSOCIATED(mhd_sim%n_bc)) NULLIFY(mhd_sim%n_bc)
IF (ASSOCIATED(mhd_sim%velx_bc)) NULLIFY(mhd_sim%velx_bc)
IF (ASSOCIATED(mhd_sim%vely_bc)) NULLIFY(mhd_sim%vely_bc)
IF (ASSOCIATED(mhd_sim%velz_bc)) NULLIFY(mhd_sim%velz_bc)
IF (ASSOCIATED(mhd_sim%T_bc)) NULLIFY(mhd_sim%T_bc)
IF (ASSOCIATED(mhd_sim%psi_bc)) NULLIFY(mhd_sim%psi_bc)
IF (ASSOCIATED(mhd_sim%by_bc)) NULLIFY (mhd_sim%by_bc)

ALLOCATE(mhd_sim%n_bc(mhd_sim%fe_rep%fields(1)%fe%ne))
ALLOCATE(mhd_sim%velx_bc(mhd_sim%fe_rep%fields(2)%fe%ne))
ALLOCATE(mhd_sim%vely_bc(mhd_sim%fe_rep%fields(3)%fe%ne))
ALLOCATE(mhd_sim%velz_bc(mhd_sim%fe_rep%fields(4)%fe%ne))
ALLOCATE(mhd_sim%T_bc(mhd_sim%fe_rep%fields(5)%fe%ne))
ALLOCATE(mhd_sim%psi_bc(mhd_sim%fe_rep%fields(6)%fe%ne))
ALLOCATE(mhd_sim%by_bc(mhd_sim%fe_rep%fields(7)%fe%ne))

mhd_sim%psi_bc = .FALSE.
mhd_sim%by_bc = .FALSE.
!Optionally fix F (= R B_phi) on the outer mesh boundary at its initial value (fixed TF current)
IF(self%fix_F_boundary) mhd_sim%by_bc = mhd_sim%fe_rep%fields(7)%fe%be
mhd_sim%n_bc = .NOT.evolve_dens !density prescribed (flat, or fixed initial profile) unless evolved
!Temperature: evolved (free) by default when evolved at all; pinned (Dirichlet)
!below on any DOF touched by a non-MHD cell (including the MHD/non-MHD
!interface), and everywhere when not evolved (T(psi) closure)
mhd_sim%T_bc = .NOT.evolve_temp
mhd_sim%velx_bc = .FALSE.
mhd_sim%velz_bc = .FALSE.
IF (tor_flow) THEN
  mhd_sim%vely_bc = .FALSE.
ELSE
  mhd_sim%vely_bc = .TRUE.
END IF
!Pin velocity, temperature and density to their boundary values outside the MHD regions
ALLOCATE(cell_dofs(lag_rep%nce))
DO i=1,mesh%nc
  IF(.NOT.mhd_flag(mesh%reg(i)))THEN
    CALL lag_rep%ncdofs(i,cell_dofs)
    DO j=1,SIZE(cell_dofs)
      mhd_sim%velx_bc(cell_dofs(j)) = .TRUE.
      mhd_sim%vely_bc(cell_dofs(j)) = .TRUE.
      mhd_sim%velz_bc(cell_dofs(j)) = .TRUE.
      mhd_sim%T_bc(cell_dofs(j)) = .TRUE.
      mhd_sim%n_bc(cell_dofs(j)) = .TRUE.
    END DO
  END IF
END DO
DEALLOCATE(cell_dofs)

self%mug => mhd_sim

!------------------------------------------------------------------------------------
! Create solver fields, augmented with coil currents for compatibility with TokaMaker
!------------------------------------------------------------------------------------
self%fe_rep => self%mug%fe_rep
IF (self%tkmr%gs_device%ncoils > 0) THEN
    ALLOCATE(stitch_tmp(8),map_tmp(8))
    DO i=1,7
        stitch_tmp(i)%s=>self%fe_rep%fields(i)%fe%linkage
        map_tmp(i)%m=>self%fe_rep%fields(i)%fe%map
    END DO
    stitch_tmp(8)%s=>self%tkmr%gs_device%coil_stitch
    map_tmp(8)%m=>self%tkmr%gs_device%coil_map
    CALL create_vector(self%aug_vec,stitch_tmp,map_tmp)
    DEALLOCATE(stitch_tmp,map_tmp)
    CALL self%aug_vec%new(self%u)
    CALL self%aug_vec%new(self%rhs)
    CALL self%aug_vec%new(self%tmp)
ELSE
    CALL self%fe_rep%vec_create(self%u)
    CALL self%fe_rep%vec_create(self%rhs)
    CALL self%fe_rep%vec_create(self%tmp)
END IF

!------------------------------------------------------------------------------
! Set initial field values
!------------------------------------------------------------------------------
CALL self%u%set(1.d0, 1) !flat density (physical n = den_scale); overwritten below if n_prof is set
CALL self%u%set(0.d0, 2)
CALL self%u%set(0.d0, 3)
CALL self%u%set(0.d0, 4)
CALL self%u%set(0.d0, 5)
CALL self%u%set(self%tkmr%gs_equil%I%f_offset, 7) !F0 outside the plasma; overwritten inside below

!Initialize field 6 (plasma poloidal flux) and field 8 (coil currents) from the
!TokaMaker equilibrium, and populate MUG's vacuum flux (mug%psi_vac) from the coils.
NULLIFY(tmp_arr, vals_out)
CALL self%tkmr%gs_device%fe_rep%vec_create(tmp_vec)
CALL tmp_vec%add(0.d0,1.d0,self%tkmr%gs_equil%psi) ! tmp_vec = total poloidal flux
IF (self%tkmr%gs_device%ncoils > 0) THEN
    CALL self%u%get_local(vals_out,8)
    DO i=1,self%tkmr%gs_device%ncoils
        CALL tmp_vec%add(1.d0,-self%tkmr%gs_equil%coil_currs(i),self%tkmr%gs_device%psi_coil(i)%f)
        vals_out(i)=self%tkmr%gs_equil%coil_currs(i)
    END DO
    CALL self%u%restore_local(vals_out,8)
    DEALLOCATE(vals_out)
END IF
CALL tmp_vec%get_local(tmp_arr)
CALL self%u%restore_local(tmp_arr,6)
DEALLOCATE(tmp_arr)
CALL tmp_vec%delete
DEALLOCATE(tmp_vec)

CALL build_psi_vac(self, self%tkmr%gs_equil%coil_currs)

!---Update plasma bounds from the total flux, then seed F and T from the
!   equilibrium profiles inside the plasma
CALL mugtok_refresh_bounds(self, self%u)
CALL mugtok_seed_fields(self)

!---Capture the seeded temperature for the Spitzer resistivity, so the first
!   Jacobian build below sees the initial eta profile
CALL xmhd_2d_freeze_eta_temp(self%mug, self%u)
CALL xmhd_2d_freeze_ohmic(self%mug, self%u)

!------------------------------------------------------------------------------
! Build the nonlinear function, approximate Jacobian, preconditioner and solvers
!------------------------------------------------------------------------------
CALL mugtok_finalize_setup(self, dt, lin_tol, nl_tol)
end subroutine setup_mugtok_mhd_td

!---------------------------------------------------------------------------
!> Temperature closure T(psi) for the MHD-plasma model. Returns T and its
!! spatial gradient so that the momentum pressure force (2*k*n*T) reproduces
!! the equilibrium pressure P(psi)=p_scale*P(psi) inside the plasma, and zero
!! outside the LCFS. Matches the xmhd_2d_temp_func interface (uses the module
!! current_sim pointer to reach the equilibrium).
!---------------------------------------------------------------------------
subroutine mugtok_T_of_psi(coords, n, dn, psi, dpsi, by, dby, T, dT, dT_dpsi)
real(r8), intent(in) :: coords(3), n, dn(3), psi, dpsi(3), by, dby(3)
real(r8), intent(out) :: T, dT(3)
real(r8), intent(out), optional :: dT_dpsi !< Partial derivative dT/dpsi at fixed n (for the resistivity Jacobian)
type(gs_equil), pointer :: eq
real(r8) :: k_b, p_scale_i, denom, p_eq, dpdpsi
eq => current_sim%tkmr%gs_equil
k_b = current_sim%mug%k_boltz
p_scale_i = eq%p_scale
denom = 2.d0*k_b*n
IF(gs_test_bounds(eq, coords(1:2)) .AND. (psi > eq%plasma_bounds(1)))THEN
  p_eq = p_scale_i*eq%P%f(psi)/mu0
  dpdpsi = p_scale_i*eq%P%fp(psi)/mu0
  T = p_eq/denom
  !dT = d/dx [ P(psi)/(2 k n) ] via the chain rule (n may vary in space)
  dT = (dpdpsi*dpsi)/denom - p_eq*dn/(2.d0*k_b*n**2)
  !dT/dpsi at fixed n (n is prescribed, not a function of psi)
  IF(PRESENT(dT_dpsi)) dT_dpsi = dpdpsi/denom
ELSE
  p_eq = p_scale_i*eq%P%f(eq%plasma_bounds(1))/mu0
  dpdpsi = p_scale_i*eq%P%fp(eq%plasma_bounds(1))/mu0
  T = p_eq/denom
  dT = (dpdpsi*dpsi)/denom - p_eq*dn/(2.d0*k_b*n**2)
  IF(PRESENT(dT_dpsi)) dT_dpsi = 0.d0
END IF

end subroutine mugtok_T_of_psi

!---------------------------------------------------------------------------
!> Test whether a point is inside the plasma (closed flux surface) region,
!! i.e. psi is within the current plasma_bounds. Uses the module current_sim
!! pointer to reach the equilibrium. Used by mugtok_seed_fields, the T(psi)
!! closure, and mugtok_apply_sink.
!---------------------------------------------------------------------------
function mugtok_in_plasma(coords, psi) result(in_plasma)
real(r8), intent(in) :: coords(3), psi
logical :: in_plasma
type(gs_equil), pointer :: eq
eq => current_sim%tkmr%gs_equil
in_plasma = gs_test_bounds(eq, coords(1:2)) .AND. (psi > eq%plasma_bounds(1))
end function mugtok_in_plasma

!---------------------------------------------------------------------------
!> Update the equilibrium plasma bounds (and the P/I flux-function bounds)
!! from the current total poloidal flux (field 6 + vacuum). Used so the T(psi)
!! closure and the F/T seeding normalize psi consistently.
!---------------------------------------------------------------------------
subroutine mugtok_refresh_bounds(self, vec)
class(oft_mugtok_td), intent(inout) :: self !< Simulation object
class(oft_vector), intent(inout) :: vec !< Field vector whose psi (field 6) sets the bounds
type(gs_equil), pointer :: eq
real(r8), pointer, dimension(:) :: psi6, pvac, eqpsi
eq => self%tkmr%gs_equil
NULLIFY(psi6, pvac, eqpsi)
CALL vec%get_local(psi6, 6)
CALL self%mug%psi_vac%get_local(pvac)
CALL eq%psi%get_local(eqpsi)
eqpsi = psi6 + pvac ! total poloidal flux
CALL eq%psi%restore_local(eqpsi)
CALL gs_update_bounds(eq, track_opoint=.TRUE.)
eq%P%plasma_bounds = eq%plasma_bounds
eq%I%plasma_bounds = eq%plasma_bounds
DEALLOCATE(psi6, pvac, eqpsi)
end subroutine mugtok_refresh_bounds

!---------------------------------------------------------------------------
!> Seed the toroidal field function F (field 7), the temperature (field 5) and,
!! when `n_prof`/`T_prof` are set, the density (field 1) on the plasma-region
!! DOFs from the equilibrium profiles. With profiles, inside the plasma
!! n = n_offset + n_scale*nhat(psihat) and T = T_offset + T_scale*That(psihat),
!! and n = n_offset, T = T_offset everywhere else (including the DOFs pinned
!! outside the MHD regions). Without profiles, n is left flat and
!! T = P(psi)/(2*k*n) inside the plasma (0 outside). When the temperature is
!! evolved this is its initial condition (and, on the DOFs pinned outside the
!! MHD regions, its boundary value); otherwise it is carried for diagnostics
!! only. Assumes mugtok_refresh_bounds has set eq%psi to the total flux.
!---------------------------------------------------------------------------
subroutine mugtok_seed_fields(self)
class(oft_mugtok_td), intent(inout) :: self !< Simulation object
type(gs_equil), pointer :: eq
real(r8), pointer, dimension(:) :: f_vals, t_vals, n_vals, eqpsi
integer(i4), allocatable :: cell_dofs_l(:)
real(r8) :: f_pos(3), coord(3), psi_dof, F0_val, f_scale_i, p_scale_i, k_b, n_phys, p_eq
logical :: in_plasma, use_prof
integer(i4) :: ic, jd, gd
eq => self%tkmr%gs_equil
F0_val = eq%I%f_offset
f_scale_i = self%tkmr%f_scale
p_scale_i = eq%p_scale
k_b = self%mug%k_boltz
n_phys = self%mug%den_scale ! density field is flat at 1.0 (without profiles)
use_prof = ASSOCIATED(self%n_prof)
NULLIFY(f_vals, t_vals, n_vals, eqpsi)
CALL self%u%get_local(f_vals, 7)
CALL self%u%get_local(t_vals, 5)
CALL self%u%get_local(n_vals, 1)
CALL eq%psi%get_local(eqpsi) ! total poloidal flux (set by mugtok_refresh_bounds)
IF(use_prof)THEN
  !---Normalize the n/T profiles with the current plasma bounds; offsets everywhere
  !   by default (overwritten inside the plasma below)
  self%n_prof%plasma_bounds = eq%plasma_bounds
  self%T_prof%plasma_bounds = eq%plasma_bounds
  n_vals = self%n_offset/self%mug%den_scale
  t_vals = self%T_offset
END IF
ALLOCATE(cell_dofs_l(lag_rep%nce))
DO ic=1,mesh%nc
  IF(mesh%reg(ic) /= 1) CYCLE ! plasma region only
  CALL lag_rep%ncdofs(ic, cell_dofs_l)
  DO jd=1,lag_rep%nce
    gd = cell_dofs_l(jd)
    CALL oft_blag_npos(lag_rep, ic, jd, f_pos)   ! DOF logical position
    coord = mesh%log2phys(ic, f_pos)             ! -> physical (R,Z)
    psi_dof = eqpsi(gd)
    in_plasma = gs_test_bounds(eq, coord(1:2)) .AND. (psi_dof > eq%plasma_bounds(1))
    IF(in_plasma)THEN
      f_vals(gd) = SQRT(f_scale_i*eq%I%f(psi_dof) + F0_val**2)
      p_eq = p_scale_i*eq%P%f(psi_dof)/mu0
    ELSE
      f_vals(gd) = F0_val
      p_eq = 0.d0
    END IF
    IF(use_prof)THEN
      !---linterp profiles store the profile values as the derivative dofs, so
      !   %fp returns the (normalized) profile value itself at psihat
      IF(in_plasma)THEN
        n_vals(gd) = (self%n_offset + self%n_scale*self%n_prof%fp(psi_dof))/self%mug%den_scale
        t_vals(gd) = self%T_offset + self%T_scale*self%T_prof%fp(psi_dof)
      END IF
    ELSE
      t_vals(gd) = MAX(p_eq, 0.d0)/(2.d0*k_b*n_phys)
    END IF
  END DO
END DO
DEALLOCATE(cell_dofs_l)
CALL self%u%restore_local(f_vals, 7)
CALL self%u%restore_local(t_vals, 5)
CALL self%u%restore_local(n_vals, 1)
DEALLOCATE(f_vals, t_vals, n_vals, eqpsi)
end subroutine mugtok_seed_fields

!---------------------------------------------------------------------------
!> Capture the edge (LCFS) poloidal flux the first time the plasma is limited
!! (not diverted). This reference value (`psi_limiting`) is stored permanently
!! -- once captured it is never overwritten -- and is later used by
!! `mugtok_apply_sink` to identify points that were part of the (once-limited)
!! plasma but have since fallen outside the current plasma boundary. Called
!! once per accepted timestep (not per nonlinear iteration): psi_limiting is a
!! path-dependent history variable, so it must only be set from the converged
!! state of a step, not from transient intermediate Newton iterates. Assumes
!! eq%plasma_bounds/eq%diverted already reflect the current step (refreshed
!! via mugtok_refresh_bounds during the nonlinear solve in `step_mugtok_td`).
!---------------------------------------------------------------------------
subroutine mugtok_update_psi_limiting(self)
class(oft_mugtok_td), intent(inout) :: self !< Simulation object
IF(self%psi_limiting_set) RETURN
IF(self%tkmr%gs_equil%diverted) RETURN
self%psi_limiting = self%tkmr%gs_equil%plasma_bounds(1)
self%psi_limiting_set = .TRUE.
end subroutine mugtok_update_psi_limiting

!---------------------------------------------------------------------------
!> Test whether the temperature sink term should be applied at a point.
!! Matches the xmhd_2d_sink_func interface (uses the module current_sim
!! pointer to reach the equilibrium). Never applies while the plasma is
!! currently diverted. Otherwise applies where psi is within the historical
!! (once-limited) plasma boundary `psi_limiting` but the point does not
!! satisfy the current in-plasma test -- i.e. points left behind as a limited
!! plasma boundary recedes.
!---------------------------------------------------------------------------
function mugtok_apply_sink(coords, psi) result(apply_sink)
real(r8), intent(in) :: coords(3), psi
logical :: apply_sink
type(gs_equil), pointer :: eq
eq => current_sim%tkmr%gs_equil
apply_sink = .FALSE.
IF(eq%diverted) RETURN
IF(.NOT.current_sim%psi_limiting_set) RETURN
IF(psi > current_sim%psi_limiting) apply_sink = .NOT.mugtok_in_plasma(coords, psi)
end function mugtok_apply_sink


!-----------------------------------------------------------------
!> Take one timestep of the combined MUG and TokaMaker simulation
!-----------------------------------------------------------------
!---------------------------------------------------------------------------
!> Build the spatial block-Jacobi partition and factorize its subdomains now, by applying the
!! preconditioner once to a zero vector (the factorization would otherwise happen inside the first
!! linear solve, so this adds no cost), and report the subdomains: unknowns in total and per field,
!! and the time of this first factorization
!---------------------------------------------------------------------------
subroutine mugtok_report_parts(self)
CLASS(oft_mugtok_td), target, intent(inout) :: self
CLASS(oft_vector), POINTER :: v1, v2
INTEGER(i4) :: ip, jj, kf, row
INTEGER(8) :: c0, c1, crate
INTEGER(i4), ALLOCATABLE :: field_of(:), cnt(:)
NULLIFY(v1, v2)
CALL self%rhs%new(v1); CALL self%rhs%new(v2)
CALL v1%set(0.d0); CALL v2%set(0.d0)
CALL SYSTEM_CLOCK(c0, crate)
CALL self%pre%apply(v1, v2)
CALL SYSTEM_CLOCK(c1)
CALL v1%delete; CALL v2%delete
DEALLOCATE(v1, v2)
SELECT TYPE(bj => self%pre)
CLASS IS(oft_bjprecond)
    IF(ASSOCIATED(bj%parts))THEN
        ALLOCATE(field_of(bj%A%nr), cnt(bj%A%ni))
        field_of = 0
        DO kf=1,bj%A%ni
            DO jj=1,bj%A%i_map(kf)%nslice
                field_of(bj%A%i_map(kf)%offset+bj%A%i_map(kf)%slice(jj)) = kf
            END DO
        END DO
        WRITE(*,'(A,I4,A,L2,A,I9,A,F8.1,A)')' Block-Jacobi subdomains:',bj%nlocal,' (overlap',bj%loverlap,') of', &
            bj%A%nr,' unknowns; first factorization',REAL(c1-c0,8)/REAL(crate,8),' s; per subdomain: total | by field 1..n'
        DO ip=1,bj%nlocal
            cnt = 0
            DO jj=1,bj%parts(ip)%n
                row = ABS(bj%parts(ip)%v(jj))
                IF(field_of(row)>0) cnt(field_of(row)) = cnt(field_of(row)) + 1
            END DO
            WRITE(*,'(A,I3,A,I9,A,20I8)')'   subdomain',ip,':',bj%parts(ip)%n,'  |',cnt
        END DO
        DEALLOCATE(field_of, cnt)
    END IF
END SELECT
FLUSH(6)
self%pre_parts_reported = .TRUE.
end subroutine mugtok_report_parts

subroutine step_mugtok_td(self,coil_currents,time,dt,nl_its,lin_its,nretry)
CLASS(oft_mugtok_td), target, intent(inout) :: self !< Simulation object
REAL(r8), INTENT(in) :: coil_currents(:) !< Prescribed I-coil currents for this step (size ncoils)
REAL(8), INTENT(inout) :: time,dt !< Current time and desired timestep size [s]
INTEGER(4), INTENT(out) :: nl_its,lin_its,nretry !< Number of nonlinear and linear iterations, and number of retries
INTEGER(4) :: i,j,nclip,nh
INTEGER(8) :: c_step0, c_step1, c_jac0, c_rate
REAL(r8) :: t_step, t_jac, step_cost, dt_new
LOGICAL :: rebuilt
CHARACTER(LEN=12) :: pre_why
REAL(r8), pointer :: tmp_arr(:), currs_tmp(:)
CLASS(oft_vector), pointer :: tmp_vec
CLASS(oft_native_matrix), POINTER :: P => NULL()
CLASS(oft_matrix), pointer :: pre

current_sim=>self
CALL SYSTEM_CLOCK(c_step0, c_rate)
t_jac = 0.d0
rebuilt = .FALSE.
pre_why = 'reuse'

!Flux normalization for a psi-keyed SOL viscosity (xmhd_2d visc_sol_on_psi): [LCFS, axis] of the current state
self%mug%psi_norm_bounds = self%tkmr%gs_equil%plasma_bounds

!Freeze the temperature used for the Spitzer resistivity at the start of the step,
!so that eta is the same coefficient in the RHS and LHS residual passes below
CALL xmhd_2d_freeze_eta_temp(self%mug, self%u)
CALL xmhd_2d_freeze_ohmic(self%mug, self%u)

!Update plasma time advance operator
CALL self%tkmr%update()

!Update timestep if it has changed since setup. The residual depends on dt only through the
!vacuum operator (tkmr%vac_op): with dt_refactor_tol > 0 a small change rebuilds just that
!operator and keeps the preconditioner (built at mug%jac_dt); otherwise it forces a rebuild
IF(dt/=self%nlfun%dt)THEN
    dt=ABS(dt)
    self%nlfun%dt=dt
    self%mug%dt = dt
    self%tkmr%dt=dt
    IF(self%dt_refactor_tol>0.d0 .AND. self%mug%jac_dt>0.d0 .AND. &
        ABS(dt-self%mug%jac_dt)<=self%dt_refactor_tol*self%mug%jac_dt)THEN
        CALL build_vac_op(self%tkmr,self%tkmr%vac_op)
    ELSE
        self%pre_force_update=.TRUE.
        pre_why = 'dt change'
    END IF
END IF
CALL self%tkmr%update()

!Rebuild approximate jacobian and refactor preconditioner: every pre_freq steps, or with pre_adapt
!once the last step cost more wall time per simulated second than the average since the last
!rebuild (rebuild included), or after pre_max_age steps
IF(self%pre_force_update)THEN
    IF(pre_why=='reuse')pre_why = 'forced'
ELSE IF(self%pre_adapt)THEN
    IF(self%pre_max_age>0 .AND. self%pre_age>=self%pre_max_age)THEN
        pre_why = 'max age'
    ELSE IF(self%pre_cyc_n>=2 .AND. self%pre_last_cost>self%pre_cyc_time/self%pre_cyc_dt)THEN
        pre_why = 'cost'
    END IF
ELSE IF(self%pre_age>=MAX(1,self%pre_freq))THEN
    pre_why = 'pre_freq'
END IF
IF(pre_why/='reuse')THEN
    CALL SYSTEM_CLOCK(c_jac0)
    CALL build_mugtok_td_jacobian(self, self%nlfun%jac_op, self%u, update_vac = .TRUE.)
    CALL self%pre%update(self%lu_reanalyze .OR. (.NOT.self%lu_analysed))
    self%lu_analysed = .TRUE.
    IF((self%pre_nlocal > 0).AND.(.NOT.self%pre_parts_reported))CALL mugtok_report_parts(self)
    self%pre_age=0
    self%pre_force_update=.FALSE.
    rebuilt = .TRUE.
    CALL SYSTEM_CLOCK(c_step1)
    t_jac = REAL(c_step1-c_jac0,8)/REAL(c_rate,8)
ELSE
    WRITE(*,'(A,I4,A)')' Reusing preconditioner (age ',self%pre_age,' steps)'
END IF
self%pre_age=self%pre_age+1

!Per-cell temperature-floor source for this step (xmhd_2d T_floor; no-op when off)
CALL xmhd_2d_floor_update(self%mug, .FALSE.)

!Build RHS
CALL self%tmp%add(0.d0,1.d0,self%u)
CALL apply_rhs_mugtok(self%nlfun,self%u,self%rhs)

NULLIFY(currs_tmp)
DO j = 1,4
    !---Prescribe the I-coil currents using the coil current input vector
    IF(self%tkmr%gs_device%ncoils > 0)THEN
        CALL self%rhs%get_local(currs_tmp,8)
        DO i=1,self%tkmr%gs_device%ncoils
            currs_tmp(i)=coil_currents(i)
        END DO
        CALL self%rhs%restore_local(currs_tmp,8)
        CALL self%u%get_local(currs_tmp,8)
        DO i=1,self%tkmr%gs_device%ncoils
            currs_tmp(i)=coil_currents(i)
        END DO
        CALL self%u%restore_local(currs_tmp,8)
    END IF
    !Run the nonlinear solver
    dbg_it = 0 ! NLDEBUG
    CALL self%nksolver%apply(self%u,self%rhs)
    IF(self%nksolver%cits<0)THEN
        IF(self%resid_dump)CALL mugtok_dump_residual(self, time)
        CALL self%u%add(0.d0,1.d0,self%tmp)
        CALL xmhd_2d_floor_update(self%mug, .TRUE.) ! raise the floor source from this attempt's flags
        !---A stale preconditioner may be the cause: refresh it and retry at the same dt first
        IF(self%pre_age>1)THEN
            WRITE(*,*)'Newton failed with a stale preconditioner; rebuilding and retrying at same dt'
            CALL build_mugtok_td_jacobian(self, self%nlfun%jac_op, self%u, update_vac = .TRUE.)
            CALL self%pre%update(self%lu_reanalyze .OR. (.NOT.self%lu_analysed))
            self%lu_analysed = .TRUE.
            self%pre_age=1
            rebuilt = .TRUE.
            pre_why = 'retry'
            CYCLE
        END IF
        self%nlfun%dt=self%nlfun%dt/2.d0
        self%mug%dt=self%nlfun%dt
        self%tkmr%dt=self%nlfun%dt
        CALL self%tkmr%update()
        CALL build_mugtok_td_jacobian(self, self%nlfun%jac_op, self%u, update_vac = .TRUE.)
        CALL self%pre%update(self%lu_reanalyze .OR. (.NOT.self%lu_analysed))
        self%lu_analysed = .TRUE.
        self%pre_age=1
        rebuilt = .TRUE.
        pre_why = 'retry dt/2'
        CALL apply_rhs_mugtok(self%nlfun,self%u,self%rhs)
        CYCLE
    ELSE
        EXIT
    END IF
END DO

!Advance time and return iteration counts (nretry = extra Newton attempts this step: stale-preconditioner
!retries and dt halvings; j is 5 if all four attempts failed)
nretry=j-1
time=time+self%nlfun%dt
dt=self%nlfun%dt
nl_its=self%nksolver%nlits
lin_its=self%nksolver%lits

!---Cost bookkeeping for pre_adapt (a step that rebuilt the preconditioner starts a new cycle;
! with retries its time includes the failed attempts) and the dt_adapt controller
CALL SYSTEM_CLOCK(c_step1)
t_step = REAL(c_step1-c_step0,8)/REAL(c_rate,8)
step_cost = t_step/dt
IF(rebuilt)THEN
    self%pre_cyc_time = t_step
    self%pre_cyc_dt = dt
    self%pre_cyc_n = 1
ELSE
    self%pre_cyc_time = self%pre_cyc_time + t_step
    self%pre_cyc_dt = self%pre_cyc_dt + dt
    self%pre_cyc_n = self%pre_cyc_n + 1
END IF
self%pre_last_cost = step_cost
self%dt_next = dt
IF(self%dt_adapt)THEN
    IF(nretry>0)self%dt_nhist = 0 ! retries: restart the history from this step
    nh = MAX(1,MIN(self%dt_hist,SIZE(self%dt_hist_dt)))
    self%dt_hist_dt(2:nh) = self%dt_hist_dt(1:nh-1)
    self%dt_hist_its(2:nh) = self%dt_hist_its(1:nh-1)
    self%dt_hist_dt(1) = dt
    self%dt_hist_its(1) = lin_its
    self%dt_nhist = MIN(self%dt_nhist+1,nh)
    IF(SUM(self%dt_hist_its(1:self%dt_nhist))>0)THEN
        dt_new = REAL(self%lin_its_target,8)*SUM(self%dt_hist_dt(1:self%dt_nhist)) &
            /REAL(SUM(self%dt_hist_its(1:self%dt_nhist)),8)
    ELSE
        dt_new = self%dt_grow_max*dt
    END IF
    dt_new = MAX(MIN(dt_new,self%dt_grow_max*dt),self%dt_shrink_min*dt)
    IF(nretry>0)dt_new = MIN(dt_new,dt) ! no growth right after a failure
    IF(self%dt_max>0.d0)dt_new = MIN(dt_new,self%dt_max)
    self%dt_next = dt_new
END IF
IF(self%pre_adapt .OR. self%dt_adapt)WRITE(*,'(A,ES11.4,A,I3,A,I4,A,F8.1,A,F8.1,A,ES10.3,A,ES10.3,A,I4,2A,ES11.4)') &
    ' Step summary: dt =',dt,'  Newton',nl_its,'  GMRES',lin_its,'  wall',t_step,' s  jac+pre',t_jac, &
    ' s  cost',step_cost,' s/s  cycle avg',self%pre_cyc_time/self%pre_cyc_dt,' s/s  LU age',self%pre_age, &
    '  pre: ',TRIM(pre_why)//'  dt_next =',self%dt_next

!---Optional hard temperature floor on the accepted state
IF(self%T_clip>0.d0)THEN
    NULLIFY(tmp_arr)
    CALL self%u%get_local(tmp_arr,5)
    nclip=COUNT(tmp_arr<self%T_clip)
    IF(nclip>0)THEN
        WRITE(*,'(A,I8,A,ES11.3,A,F6.2,A)')' T_clip: ',nclip,' T DOFs (min',MINVAL(tmp_arr),' eV) raised to',self%T_clip,' eV'
        tmp_arr=MAX(tmp_arr,self%T_clip)
        CALL self%u%restore_local(tmp_arr,5)
    END IF
    DEALLOCATE(tmp_arr)
END IF

!Update equilibrium object to use the coil currents (field 8) 
!and psi (field 6 + coil vaccum flux) from the converged solution vector 
IF(self%tkmr%gs_device%ncoils > 0)THEN
    CALL self%u%get_local(currs_tmp,8)
    DO i=1,self%tkmr%gs_device%ncoils
        self%tkmr%gs_equil%coil_currs(i)=currs_tmp(i)
    END DO
END IF
IF(ASSOCIATED(currs_tmp))DEALLOCATE(currs_tmp)

NULLIFY(tmp_arr)
CALL self%u%get_local(tmp_arr,6)
CALL self%tkmr%gs_equil%psi%restore_local(tmp_arr)
DEALLOCATE(tmp_arr)
DO i=1,self%tkmr%gs_device%ncoils
    CALL self%tkmr%gs_equil%psi%add(1.d0,self%tkmr%gs_equil%coil_currs(i),self%tkmr%gs_device%psi_coil(i)%f)
END DO

!---Capture the historical limiting-flux reference used by the T sink term
CALL mugtok_update_psi_limiting(self)

!---Optionally write a restart file for later plotting
self%rst_count = self%rst_count + 1
IF(self%save_rst .AND. MOD(self%rst_count, self%rst_freq)==0)THEN
    CALL mugtok_rst_save(self, time)
END IF

!---Save the profile from this iteration to be used as the RHS in future timesteps
CALL snapshot_f_profile(self)
write(*,*) 'Diagnostic b: ', self%mug%nlfun%diag_vals(1)
write(*,*) 'Diagnostic v: ', self%mug%nlfun%diag_vals(2)
write(*,*) 'Diagnostic p force: ', self%mug%nlfun%diag_vals(3)
write(*,*) 'Diagnostic T: ', self%mug%nlfun%diag_vals(4)
end subroutine step_mugtok_td

!------------------------------------------------------
!Compute the RHS of the coupled MUG + TokaMaker system 
!------------------------------------------------------
subroutine apply_rhs_mugtok(self, a, b)
class(oft_mugtok_td_mfop), intent(inout) :: self !< Time-advance operator 
class(oft_vector), target, intent(inout) :: a !< Source field
class(oft_vector), intent(inout) :: b !< Result of metric function
class(oft_vector), pointer :: tmp_in, tmp_out 
REAL(r8), POINTER, DIMENSION(:) :: tmp_arr1, tmp_arr2, coil_arr

!------------------------------------------------------
!Compute MUG contribution
!------------------------------------------------------
self%parent_sim%mug%nlfun%dt = 0.d0 !Set dt = 0 for MUG RHS

!---Rebuild psi_vac from input coil currents (currents from previous timestep)
IF(self%parent_sim%tkmr%gs_device%ncoils > 0)THEN
  NULLIFY(coil_arr)
  CALL a%get_local(coil_arr, 8)
  CALL build_psi_vac(self%parent_sim, coil_arr)
  CALL a%restore_local(coil_arr, 8)
  DEALLOCATE(coil_arr)
END IF

!---MHD-plasma model: refresh plasma bounds from this iterate's flux so the
!   T(psi) closure normalizes psi consistently (psi changes every iteration)
IF(self%parent_sim%mhd_plasma) CALL mugtok_refresh_bounds(self%parent_sim, a)

!---Apply MUG RHS operator
CALL b%set(0.d0)
CALL self%parent_sim%mug%nlfun%apply_real(a,b)

!---Add F RHS in the non-MHD regions, which are currently not computed in either code
CALL add_f_terms(self, a, b, 0.d0) !use dt = 0 to get RHS 
NULLIFY(tmp_arr1)
CALL b%get_local(tmp_arr1, 6) !put psi RHS into tmp array to add to TokaMaker

!------------------------------------------------------
!Compute TokaMaker contribution
!------------------------------------------------------
IF (self%parent_sim%tkmr%gs_device%ncoils > 0) THEN
    CALL self%parent_sim%tkmr%gs_device%aug_vec%new(tmp_in)
    CALL self%parent_sim%tkmr%gs_device%aug_vec%new(tmp_out)
ELSE
    CALL lag_rep%vec_create(tmp_in)
    CALL lag_rep%vec_create(tmp_out)
END IF
CALL tmp_in%set(0.d0)
CALL tmp_out%set(0.d0)

! Extract component 6 from a and put into component 1 of tmp_in
NULLIFY(tmp_arr2)
CALL a%get_local(tmp_arr2, 6)
CALL tmp_in%restore_local(tmp_arr2, 1)
CALL a%restore_local(tmp_arr2, 6)  ! 

! Extract component 8 from a and put into component 2 of tmp_in
NULLIFY(tmp_arr2)
IF (self%parent_sim%tkmr%gs_device%ncoils > 0) THEN
    CALL a%get_local(tmp_arr2, 8)
    CALL tmp_in%restore_local(tmp_arr2, 2)
    CALL a%restore_local(tmp_arr2, 8)
    NULLIFY(tmp_arr2)
END IF

!Compute TokaMaker RHS
CALL apply_rhs(self%parent_sim%tkmr, tmp_in, tmp_out)
CALL self%parent_sim%tkmr%gs_device%zerob_bc%apply(tmp_out) !Apply BCs

!------------------------------------------------------
! Combine
!------------------------------------------------------
CALL tmp_out%get_local(tmp_arr2, 1)
tmp_arr1 = tmp_arr1 + tmp_arr2
CALL b%restore_local(tmp_arr1, 6)
NULLIFY(tmp_arr1)

! Extract component 2 from tmp_out and put in field 8
IF (self%parent_sim%tkmr%gs_device%ncoils > 0) THEN
    CALL tmp_out%get_local(tmp_arr2, 2)
    CALL b%restore_local(tmp_arr2, 8)
    NULLIFY(tmp_arr2)
END IF

! Cleanup
CALL tmp_in%delete()
CALL tmp_out%delete()
DEALLOCATE(tmp_in, tmp_out)
end subroutine apply_rhs_mugtok

!------------------------------------------------------
!Compute the LHS of the coupled MUG + TokaMaker system 
!------------------------------------------------------
subroutine nlfun_apply(self, a, b)
class(oft_mugtok_td_mfop), intent(inout) :: self !< Time-advance operator
class(oft_vector), target, intent(inout) :: a !< Source field
class(oft_vector), intent(inout) :: b !< Result of metric function
class(oft_vector), pointer :: tmp_in, tmp_out
REAL(r8), POINTER, DIMENSION(:) :: tmp_arr1, tmp_arr2, coil_arr

!----------------------------------------------------------------------------------------
!Compute TokaMaker contribution (must happen first so f_scale, plasma bounds are updated)
!----------------------------------------------------------------------------------------
NULLIFY(tmp_arr2)
CALL a%get_local(tmp_arr2, 6)

IF (self%parent_sim%tkmr%gs_device%ncoils > 0) THEN
    CALL self%parent_sim%tkmr%gs_device%aug_vec%new(tmp_in)
    CALL self%parent_sim%tkmr%gs_device%aug_vec%new(tmp_out)
ELSE
    CALL lag_rep%vec_create(tmp_in)
    CALL lag_rep%vec_create(tmp_out)
END IF

CALL tmp_in%set(0.d0)
CALL tmp_out%set(0.d0)
CALL tmp_in%restore_local(tmp_arr2, 1) !Put current psi into tokamaker vector
NULLIFY(tmp_arr2)

IF (self%parent_sim%tkmr%gs_device%ncoils > 0) THEN
    CALL a%get_local(tmp_arr2, 8)
    CALL tmp_in%restore_local(tmp_arr2, 2) !Put coil currents into tokamaker vector
    NULLIFY(tmp_arr2)
END IF

!Apply tokamaker LHS
CALL self%parent_sim%tkmr%apply_real(tmp_in, tmp_out) 
CALL tmp_out%get_local(tmp_arr2, 1) !Put psi LHS into temp vector for adding

!------------------------------------------------------
!Compute MUG contribution
!------------------------------------------------------
self%parent_sim%mug%nlfun%dt = self%dt !Use real timestep for LHS

!---Rebuild psi_vac from input coil currents (currents from previous timestep)
IF(self%parent_sim%tkmr%gs_device%ncoils > 0)THEN
  NULLIFY(coil_arr)
  CALL a%get_local(coil_arr, 8)
  CALL build_psi_vac(self%parent_sim, coil_arr)
  CALL a%restore_local(coil_arr, 8)
  DEALLOCATE(coil_arr)
END IF

!---MHD-plasma model: refresh plasma bounds from this iterate's flux so the
!   T(psi) closure normalizes psi consistently (psi changes every iteration)
IF(self%parent_sim%mhd_plasma) CALL mugtok_refresh_bounds(self%parent_sim, a)

CALL b%set(0.d0)
CALL self%parent_sim%mug%nlfun%apply_real(a,b) !Apply MUG LHS

!---Add F LHS in the non-MHD regions, which are currently not computed in either code
CALL add_f_terms(self, a, b, self%dt)
NULLIFY(tmp_arr1)
CALL b%get_local(tmp_arr1, 6)

!------------------------------------------------------
! Combine
!------------------------------------------------------
tmp_arr1 = tmp_arr1 + tmp_arr2
CALL b%restore_local(tmp_arr1, 6)
NULLIFY(tmp_arr2, tmp_arr1)

! Put coil current RHS into field 8
CALL tmp_out%get_local(tmp_arr2, 2)
CALL b%restore_local(tmp_arr2, 8)
NULLIFY(tmp_arr2)

! Cleanup
CALL tmp_in%delete()
CALL tmp_out%delete()
DEALLOCATE(tmp_in, tmp_out)
end subroutine nlfun_apply

!-------------------------------------------------------------------------------------------
!! Compute remaining contributions to F residual which are not handled in MUG or TokaMaker:
!! 1) F diffusion residual in solid conducting regions (and vacuum, with high resistivity)
!! 2) Residual for F0 node from Faraday's law (plasma toroidal flux + limiter voltage integrals)
!! 3) Residual for plasma nodes (F-F(F0_node)) to drive F to a constant in the plasma
!-------------------------------------------------------------------------------------------
subroutine add_f_terms(self, a, b, dt)
class(oft_mugtok_td_mfop), intent(inout) :: self
class(oft_vector), target, intent(inout) :: a !< Source field
class(oft_vector), intent(inout) :: b !< Result
real(r8), intent(in) :: dt !< Timestep (0 for the RHS)
real(r8), pointer, dimension(:) :: by_weights, by_res
type(oft_quad_type), pointer :: quad
integer(i4) :: i
real(r8) :: eta_fallback

!------------------------------------------------------
! F diffusion in solid conductors
!------------------------------------------------------
!---Large-but-finite resistivity used where the region eta is
! undefined/negative (vacuum)
eta_fallback = 1.d-1/mu0
quad => lag_rep%quad
NULLIFY(by_weights, by_res)
CALL a%get_local(by_weights, 7)
CALL b%get_local(by_res, 7)
BLOCK
LOGICAL :: curved
INTEGER(i4) :: m, jr
INTEGER(i4), ALLOCATABLE :: cell_dofs(:)
REAL(r8) :: by, dby(3), jac_mat(3,4), jac_det, int_factor, coords(3), eta1
REAL(r8), ALLOCATABLE :: basis_vals(:), basis_grads(:,:), by_weights_loc(:), res_loc(:)
!$omp parallel private(m,jr,curved,cell_dofs,by,dby,jac_mat,jac_det,int_factor,coords,eta1, &
!$omp   basis_vals,basis_grads,by_weights_loc,res_loc)
ALLOCATE(basis_vals(lag_rep%nce), basis_grads(3,lag_rep%nce))
ALLOCATE(by_weights_loc(lag_rep%nce), cell_dofs(lag_rep%nce), res_loc(lag_rep%nce))
!$omp do ordered
DO i=1,mesh%nc
  !---Only the non-MHD regions (those the MUG solve ignores)
  IF(.NOT.self%parent_sim%mug%ignore_rmask(mesh%reg(i)))CYCLE
  curved=cell_is_curved(mesh,i)
  CALL lag_rep%ncdofs(i,cell_dofs)
  res_loc = 0.d0
  by_weights_loc = by_weights(cell_dofs)
  eta1 = self%parent_sim%mug%eta(mesh%reg(i),1)
  IF(eta1 <= 0.d0) eta1 = eta_fallback
  DO m=1,quad%np
    IF(curved.OR.(m==1))CALL mesh%jacobian(i,quad%pts(:,m),jac_mat,jac_det)
    DO jr=1,lag_rep%nce
      CALL oft_blag_eval(lag_rep,i,jr,quad%pts(:,m),basis_vals(jr))
      CALL oft_blag_geval(lag_rep,i,jr,quad%pts(:,m),basis_grads(:,jr),jac_mat)
    END DO
    coords = mesh%log2phys(i,quad%pts(:,m))
    by = 0.d0; dby = 0.d0
    basis_grads(3,:) = basis_grads(2,:)
    basis_grads(2,:) = 0.d0
    int_factor = jac_det*quad%wts(m)
    DO jr=1,lag_rep%nce
      by = by + by_weights_loc(jr)*basis_vals(jr)
      dby = dby + by_weights_loc(jr)*basis_grads(:,jr)
    END DO
    DO jr=1,lag_rep%nce
      res_loc(jr) = res_loc(jr) &
        + basis_vals(jr)*by*int_factor/(coords(1)+gs_epsilon) &
        + dt*eta1*DOT_PRODUCT(basis_grads(:,jr),dby)*int_factor/(coords(1)+gs_epsilon)
    END DO
  END DO
  !---Add local values to the full residual vector
  !$omp ordered
  DO jr=1,lag_rep%nce
    IF(self%parent_sim%mug%by_bc(cell_dofs(jr))) CYCLE !Dirichlet F rows keep the MUG Dirichlet residual
    !$omp atomic
    by_res(cell_dofs(jr)) = by_res(cell_dofs(jr)) + res_loc(jr)
  END DO
  !$omp end ordered
END DO
DEALLOCATE(basis_vals,basis_grads,by_weights_loc,cell_dofs,res_loc)
!$omp end parallel
END BLOCK

!------------------------------------------------------
! Constrain plasma nodes and F0 toroidal flux integral
!------------------------------------------------------
IF(self%parent_sim%F0_node > 0)THEN
BLOCK
TYPE(gs_equil), POINTER :: eq
TYPE(gs_factory), POINTER :: dev
CLASS(flux_func), POINTER :: ffp 
REAL(r8) :: f_scale_use 
TYPE(oft_quad_type) :: quad_1d
REAL(r8), ALLOCATABLE :: by_res_plasma(:), psi_weights_loc(:), basis_vals_2(:), &
                         basis_vals(:), basis_grads(:,:), by_weights_loc(:), ff(:)
REAL(r8), POINTER, DIMENSION(:) :: psi_weights, pvac, eqpsi
INTEGER(i4), ALLOCATABLE :: cell_dofs_2(:), cell_b_dofs(:), elist(:,:)
INTEGER(i4) :: i, m, jr, k, je, cell, ed, nlim, F0
LOGICAL :: curved
REAL(r8) :: F0_res, psi, coords(3), jac_mat(3,4), jac_det, x1, y1, x2, y2, signed_area
REAL(r8) :: pts(2,2), dl(2), dn(3), dby(3), eta_p_loc
eq => self%parent_sim%tkmr%gs_equil
dev => self%parent_sim%tkmr%gs_device
F0 = self%parent_sim%F0_node
!---LHS (dt>0) uses the current TokaMaker profile (updated in apply_mfop)
!---RHS (dt=0) uses the previous-step snapshot
IF(dt > 0.d0)THEN
  f_scale_use = self%parent_sim%tkmr%f_scale
  ffp => eq%I
ELSE
  f_scale_use = self%parent_sim%f_scale_prev
  ffp => self%parent_sim%I_prev
END IF

!---Constrain F to a constant (= F0) over the plasma DOFs
ALLOCATE(by_res_plasma(lag_rep%ne))
by_res_plasma = by_weights - by_weights(F0)
CALL fem_dirichlet_vec(lag_rep, by_res_plasma, by_res, self%parent_sim%plasma_flag)
DEALLOCATE(by_res_plasma)

!---Get the total poloidal flux (psi + psi_vac) for the bounds test
NULLIFY(psi_weights, pvac)
CALL a%get_local(psi_weights, 6)
CALL self%parent_sim%mug%psi_vac%get_local(pvac)

NULLIFY(eqpsi)
CALL eq%psi%get_local(eqpsi)
eqpsi = psi_weights + pvac
CALL eq%psi%restore_local(eqpsi)
DEALLOCATE(eqpsi)
CALL gs_update_bounds(eq, track_opoint=.TRUE.)
ffp%plasma_bounds = eq%plasma_bounds
F0_res = 0.d0
!---Integrate toroidal flux over the plasma region (region 1)
!$omp parallel private(m,jr,curved,cell_dofs_2,psi_weights_loc,basis_vals_2, &
!$omp   jac_mat,jac_det,coords,psi) reduction(+:F0_res)
ALLOCATE(basis_vals_2(lag_rep%nce), psi_weights_loc(lag_rep%nce), cell_dofs_2(lag_rep%nce))
!$omp do
DO i=1,mesh%nc
  IF(mesh%reg(i) /= 1)CYCLE ! plasma region only
  curved = cell_is_curved(mesh,i)
  CALL lag_rep%ncdofs(i, cell_dofs_2)
  psi_weights_loc = psi_weights(cell_dofs_2) + pvac(cell_dofs_2) ! total flux
  DO m=1,quad%np
    IF(curved.OR.(m==1))CALL mesh%jacobian(i,quad%pts(:,m),jac_mat,jac_det)
    DO jr=1,lag_rep%nce
      CALL oft_blag_eval(lag_rep,i,jr,quad%pts(:,m),basis_vals_2(jr))
    END DO
    coords = mesh%log2phys(i,quad%pts(:,m))
    psi = 0.d0
    DO jr=1,lag_rep%nce
      psi = psi + psi_weights_loc(jr)*basis_vals_2(jr)
    END DO
    IF(gs_test_bounds(eq,coords(1:2)) .AND. psi > eq%plasma_bounds(1))THEN
      ! inside the plasma: F = sqrt(f_scale*F*F'(psi) + F0^2)
      F0_res = F0_res + SQRT(f_scale_use*ffp%f(psi) + by_weights(F0)**2) &
               *jac_det*quad%wts(m)/(coords(1)+gs_epsilon)
    ELSE
      ! in-region but outside the plasma:  F = F0
      F0_res = F0_res + by_weights(F0)*jac_det*quad%wts(m)/(coords(1)+gs_epsilon)
    END IF
  END DO
END DO
DEALLOCATE(basis_vals_2, psi_weights_loc, cell_dofs_2)
!$omp end parallel

!---Voltage integral around the limiter contour
nlim = dev%nlim_con
signed_area = 0.d0
ALLOCATE(elist(2,nlim))
!---For each contour segment, find the non-plasma cell and its local edge index
DO i=1,nlim
  IF(i < nlim)THEN
    je = ABS(mesh_local_findedge(mesh,[dev%lim_con(i),dev%lim_con(i+1)]))
    x1 = mesh%r(1,dev%lim_con(i));   y1 = mesh%r(2,dev%lim_con(i))
    x2 = mesh%r(1,dev%lim_con(i+1)); y2 = mesh%r(2,dev%lim_con(i+1))
    signed_area = signed_area + (x1*y2 - x2*y1) !used to figure out orientation of contour
  ELSE
    je = ABS(mesh_local_findedge(mesh,[dev%lim_con(nlim),dev%lim_con(1)]))
  END IF
  !---Pick the cell on the non-plasma (region /= 1) side of the edge
  IF(mesh%reg(mesh%lec(mesh%kec(je))) /= 1)THEN
    elist(2,i) = mesh%lec(mesh%kec(je))
  ELSE
    elist(2,i) = mesh%lec(mesh%kec(je)+1)
  END IF
  DO m=1,3
    IF(je == ABS(mesh%lce(m,elist(2,i))))THEN
      elist(1,i) = m
      EXIT
    END IF
  END DO
END DO
!---Integrate eta_p*dt*(dBy/dn)/R along the contour
CALL set_quad_1d(quad_1d, lag_rep%order+2)
!$omp parallel private(cell,ed,eta_p_loc,cell_b_dofs,by_weights_loc,pts,dl,dn,k,ff, &
!$omp   coords,jac_mat,jac_det,jr,basis_vals,basis_grads,dby) reduction(+:F0_res)
ALLOCATE(cell_b_dofs(lag_rep%nce), basis_vals(lag_rep%nce), basis_grads(3,lag_rep%nce), &
         by_weights_loc(lag_rep%nce), ff(SIZE(quad%pts,1)))
!$omp do
DO i=1,nlim
  cell = elist(2,i)
  ed = elist(1,i)
  eta_p_loc = self%parent_sim%mug%eta(mesh%reg(cell),1) ! in-plane (poloidal) resistivity
  IF(eta_p_loc <= 0.d0) eta_p_loc = eta_fallback
  CALL lag_rep%ncdofs(cell, cell_b_dofs)
  by_weights_loc = by_weights(cell_b_dofs)
  pts(:,1) = mesh%r(1:2, mesh%lc(mesh%cell_ed(1,ed),cell))
  pts(:,2) = mesh%r(1:2, mesh%lc(mesh%cell_ed(2,ed),cell))
  dl = pts(:,1) - pts(:,2)
  IF(dev%lim_con(i) == mesh%lc(mesh%cell_ed(2,ed),cell)) dl = -dl
  dn = [-dl(2), dl(1), 0.d0] 
  DO k=1,quad_1d%np
    ff = 0.d0
    ff(mesh%cell_ed(1,ed)) = quad_1d%pts(1,k)
    ff(mesh%cell_ed(2,ed)) = 1.d0 - quad_1d%pts(1,k)
    coords = mesh%log2phys(cell, ff)
    CALL mesh%jacobian(cell, ff, jac_mat, jac_det)
    DO jr=1,lag_rep%nce
      CALL oft_blag_eval(lag_rep, cell, jr, ff, basis_vals(jr))
      CALL oft_blag_geval(lag_rep, cell, jr, ff, basis_grads(:,jr), jac_mat)
    END DO
    dby = 0.d0
    DO jr=1,lag_rep%nce
      dby = dby + by_weights_loc(jr)*basis_grads(:,jr)
    END DO
    F0_res = F0_res - SIGN(1.d0,signed_area)*eta_p_loc*dt*DOT_PRODUCT(dby,dn) &
             *quad_1d%wts(k)/(coords(1)+gs_epsilon)
  END DO
END DO
DEALLOCATE(cell_b_dofs, basis_vals, basis_grads, by_weights_loc, ff)
!$omp end parallel
DEALLOCATE(elist)
!---Overwrite the F0-node residual with the computed value
by_res(F0) = F0_res
CALL a%restore_local(psi_weights, 6)
END BLOCK
END IF
CALL b%restore_local(by_res, 7)
CALL a%restore_local(by_weights, 7)
end subroutine add_f_terms

!------------------------------------------------------------------------------
!> Snapshot current FF' profile and f_scale into the sim's
!! f_scale_prev / I_prev, so the next step's RHS pass of add_f_diff uses
!! these previous-step values while the LHS uses whatever is set/computed for
!! the next step
!------------------------------------------------------------------------------
subroutine snapshot_f_profile(self)
class(oft_mugtok_td), intent(inout) :: self !< Simulation object
IF(ASSOCIATED(self%I_prev))THEN
  CALL self%I_prev%delete()
  DEALLOCATE(self%I_prev)
END IF
CALL self%tkmr%gs_equil%I%copy(self%I_prev)
self%f_scale_prev = self%tkmr%f_scale
end subroutine snapshot_f_profile

!------------------------------------------------------------------------------
!> Compute vacuum psi from a set of coil currents, for use in MUG
!------------------------------------------------------------------------------
subroutine build_psi_vac(self, coil_currents)
class(oft_mugtok_td), intent(inout) :: self !< Simulation object
real(r8), intent(in) :: coil_currents(:) !< Coil currents (size ncoils)
real(r8), pointer, dimension(:) :: pv, cw
integer(i4) :: i
IF(self%tkmr%gs_device%ncoils <= 0) RETURN
NULLIFY(pv, cw)
CALL self%mug%psi_vac%get_local(pv)
pv = 0.d0
DO i=1,self%tkmr%gs_device%ncoils
  CALL self%tkmr%gs_device%psi_coil(i)%f%get_local(cw)
  pv = pv + coil_currents(i)*cw
END DO
CALL self%mug%psi_vac%restore_local(pv)
DEALLOCATE(pv)
IF(ASSOCIATED(cw))DEALLOCATE(cw)
end subroutine build_psi_vac

!------------------------------------------------------------------------------
!> Add approximate Jacobian entries corresponding to add_f_terms
!!  1) F resistive-diffusion Jacobian in solid conductor (and vacuum) regions
!!  2) An approximate F0/F0 diagonal = integral over the plasma region of dA/R
!!     (no profile effects).
!!  3) The plasma-DOF constraint rows F(ind) - F0
!------------------------------------------------------------------------------
subroutine add_f_jac(self, mat)
class(oft_mugtok_td), intent(inout) :: self !< simulation object
class(oft_matrix), pointer, intent(inout) :: mat !< Matrix to add entries to
type(oft_quad_type), pointer :: quad
integer(i4) :: i, m, jr, jc, F0, j
integer(i4), allocatable :: cell_dofs(:)
real(r8) :: eta1, eta_fallback, dt, coords(3), jac_mat(3,4), jac_det, int_factor, F0_diag
real(r8), allocatable :: basis_vals(:), basis_grads(:,:), jac_loc(:,:)
logical :: curved
logical, allocatable :: plasma_no_F0(:)
real(r8) :: nval(1,1), pval(1,1)
!---Vacuum By diffusion Jacobian is assembled regardless of the F0 model;
!   only the F0-specific constraint block below is gated on F0_node>0.
dt = self%nlfun%dt
quad => lag_rep%quad
!------------------------------------------------------
! F diffusion in solid conductors
!------------------------------------------------------
eta_fallback = 1.d-1/mu0 !Large-but-finite resistivity in vacuum, matching add_f_diff
!$omp parallel private(m,jr,jc,curved,cell_dofs,eta1,jac_loc,basis_vals,basis_grads, &
!$omp   jac_mat,jac_det,coords,int_factor)
ALLOCATE(basis_vals(lag_rep%nce), basis_grads(3,lag_rep%nce), &
         jac_loc(lag_rep%nce,lag_rep%nce), cell_dofs(lag_rep%nce))
!$omp do
DO i=1,mesh%nc
  IF(.NOT.self%mug%ignore_rmask(mesh%reg(i))) CYCLE ! non-MHD regions only
  curved = cell_is_curved(mesh,i)
  CALL lag_rep%ncdofs(i,cell_dofs)
  eta1 = self%mug%eta(mesh%reg(i),1)
  IF(eta1 < 0.d0) eta1 = eta_fallback
  jac_loc = 0.d0
  DO m=1,quad%np
    IF(curved.OR.(m==1)) CALL mesh%jacobian(i,quad%pts(:,m),jac_mat,jac_det)
    DO jr=1,lag_rep%nce
      CALL oft_blag_eval(lag_rep,i,jr,quad%pts(:,m),basis_vals(jr))
      CALL oft_blag_geval(lag_rep,i,jr,quad%pts(:,m),basis_grads(:,jr),jac_mat)
    END DO
    coords = mesh%log2phys(i,quad%pts(:,m))
    basis_grads(3,:) = basis_grads(2,:) 
    basis_grads(2,:) = 0.d0
    int_factor = jac_det*quad%wts(m)
    DO jr=1,lag_rep%nce
      DO jc=1,lag_rep%nce
        jac_loc(jr,jc) = jac_loc(jr,jc) &
          + basis_vals(jr)*basis_vals(jc)*int_factor/(coords(1)+gs_epsilon) &
          + dt*eta1*DOT_PRODUCT(basis_grads(:,jr),basis_grads(:,jc))*int_factor/(coords(1)+gs_epsilon)
      END DO
    END DO
  END DO
  !$omp critical
  DO jr=1,lag_rep%nce
    IF(self%F0_node > 0 .AND. self%plasma_flag(cell_dofs(jr))) CYCLE !In the F0 model plasma DOFs get the F0 constraint instead; without it they get diffusion (matching add_f_res) on top of the MHD induction
    IF(self%mug%by_bc(cell_dofs(jr))) CYCLE !Dirichlet F rows keep their identity row
    CALL mat%add_values([cell_dofs(jr)], cell_dofs, RESHAPE(jac_loc(jr,:),[1,lag_rep%nce]), &
                        1, lag_rep%nce, 7, 7)
  END DO
  !$omp end critical
END DO
DEALLOCATE(basis_vals, basis_grads, jac_loc, cell_dofs)
!$omp end parallel
!--------------------------------------------------------------------
! F0 reduced-model constraint (only when the F0 model is active)
!--------------------------------------------------------------------
IF(self%F0_node > 0)THEN
  F0 = self%F0_node
  !------------------------------------------------------------------
  ! Approximate F0/F0 diagonal (integral of dA/R over the plasma region)
  !------------------------------------------------------------------
  F0_diag = 0.d0
  !$omp parallel do private(m,curved,jac_mat,jac_det,coords) reduction(+:F0_diag)
  DO i=1,mesh%nc
    IF(mesh%reg(i) /= 1) CYCLE ! plasma region only
    curved = cell_is_curved(mesh,i)
    DO m=1,quad%np
      IF(curved.OR.(m==1)) CALL mesh%jacobian(i,quad%pts(:,m),jac_mat,jac_det)
      coords = mesh%log2phys(i,quad%pts(:,m))
      F0_diag = F0_diag + jac_det*quad%wts(m)/(coords(1)+gs_epsilon)
    END DO
  END DO
  !$omp end parallel do
  pval(1,1) = F0_diag
  CALL mat%add_values([F0],[F0], pval, 1,1, 7,7)
  !------------------------------------------------------
  ! Plasma DOF constraints (F(ind) - F0)
  !------------------------------------------------------
  ALLOCATE(plasma_no_F0(SIZE(self%plasma_flag)))
  plasma_no_F0 = self%plasma_flag
  plasma_no_F0(F0) = .FALSE.
  CALL fem_dirichlet_diag(lag_rep, mat, plasma_no_F0, 7) !1 on diagonal
  nval(1,1) = -1.d0 !-1 for F0 column
  DO i=1,lag_rep%ne
    IF(lag_rep%be(i) .OR. (.NOT.plasma_no_F0(i))) CYCLE
    CALL mat%add_values([i],[F0], nval, 1,1, 7,7)
  END DO
  DO i=1,lag_rep%nbe
    IF(.NOT.lag_rep%linkage%leo(i)) CYCLE
    j = lag_rep%lbe(i)
    IF(plasma_no_F0(j)) CALL mat%add_values([j],[F0], nval, 1,1, 7,7)
  END DO
  DEALLOCATE(plasma_no_F0)
END IF
end subroutine add_f_jac

!------------------------------------------------------
! Build approximate Jacobian for preconditioning
!------------------------------------------------------
subroutine build_mugtok_td_jacobian(self, mat, a, update_vac)
class(oft_mugtok_td), intent(inout) :: self !< Simulation object
class(oft_matrix), pointer, intent(inout) :: mat !< Jacobian matrix to populate
class(oft_vector), intent(inout) :: a !< Solution for computing Jacobian
LOGICAL, INTENT(in), optional :: update_vac !< Whether to update the vacuum psi operator (default: true)
class(oft_matrix), pointer:: vac_op, mug_native
CLASS(oft_native_matrix), POINTER :: V => NULL()
CLASS(oft_native_matrix), POINTER :: M => NULL()
INTEGER(4) :: i, colcount, row_block, col_block, jp, jn
INTEGER(4), ALLOCATABLE :: cols(:)
REAL(r8), ALLOCATABLE :: vals(:)
INTEGER(4), ALLOCATABLE :: bc_rows(:)
INTEGER(4) :: nbc_rows, k
class(oft_vector), pointer :: tmp

CALL mat%zero !Zero input matrix

!Cast matrices to native types for access to entries 
select type(vac_op => self%tkmr%vac_op)
type is (oft_native_matrix)
    V => vac_op
class default
    call oft_abort("vac_op must be an oft_native_matrix", &
                   "build_mugtok_td_jacobian", __FILE__)
end select

select type(mug_jac => self%mug%jacobian)
type is (oft_native_matrix)
    M => mug_jac
class default
    call oft_abort("mug jacobian must be oft_native_matrix", &
                   "build_mugtok_td_jacobian", __FILE__)
end select

!Populate MUG and TokaMaker Jacobian matrices
self%mug%jac_dt = self%nlfun%dt
CALL build_approx_jacobian(self%mug, a)
IF (PRESENT(update_vac)) THEN
    IF (update_vac) CALL build_vac_op(self%tkmr,self%tkmr%vac_op)
ELSE
    CALL build_vac_op(self%tkmr,self%tkmr%vac_op)
END IF

!---Zero the plasma_flag rows + F0 row of the MUG Jacobian's By block before
! copying, so add_f_jac can define these rows instead
IF(self%F0_node > 0)THEN
    nbc_rows = COUNT(self%plasma_flag)
    ALLOCATE(bc_rows(nbc_rows))
    k = 0
    DO i = 1, SIZE(self%plasma_flag)
        IF(self%plasma_flag(i))THEN
            k = k + 1
            bc_rows(k) = i
        END IF
    END DO
    CALL M%zero_rows(nbc_rows, bc_rows, 7)
    DEALLOCATE(bc_rows)
END IF

!Copy MUG Jacobian for all fields 
DO row_block = 1, M%ni
    DO col_block = 1, M%nj
        DO i = 1, M%i_map(row_block)%n
            jp=M%map(row_block,col_block)%ext(1,i)
            jn=M%map(row_block,col_block)%ext(2,i)
            colcount = jn-jp+1
            ALLOCATE(cols(colcount), vals(colcount))
            cols = M%lc(jp:jn)
            vals = M%M(jp:jn)
            cols = cols - M%j_map(col_block)%offset
            CALL mat%add_values([i], cols, RESHAPE(vals, [1,colcount]), &
                1, colcount, row_block, col_block)
            DEALLOCATE(cols, vals)
        END DO
    END DO 
END DO

!Add TokaMaker Jacobian to psi(6,6) block
row_block = 1
col_block = 1
DO i = 1, V%i_map(row_block)%n
    jp=V%map(row_block,col_block)%ext(1,i)
    jn=V%map(row_block,col_block)%ext(2,i)
    colcount = jn-jp+1
    ! write(*,*) colcount
    ALLOCATE(cols(colcount), vals(colcount))
    cols = V%lc(jp:jn)
    vals = V%M(jp:jn)
    cols = cols - V%j_map(col_block)%offset
    CALL mat%add_values([i], cols, RESHAPE(vals, [1,colcount]), &
        1, colcount, 6, 6)
    DEALLOCATE(cols, vals)
END DO

!Add TokaMaker Jacobian to 6,8 block (MAYBE UNNECESSARY FOR ICOILS)
row_block = 1
col_block = 2
DO i = 1, V%i_map(row_block)%n
    jp=V%map(row_block,col_block)%ext(1,i)
    jn=V%map(row_block,col_block)%ext(2,i)
    colcount = jn-jp+1
    ALLOCATE(cols(colcount), vals(colcount))
    cols = V%lc(jp:jn)
    vals = V%M(jp:jn)
    cols = cols - V%j_map(col_block)%offset
    CALL mat%add_values([i], cols, RESHAPE(vals, [1,colcount]), &
        1, colcount, 6, 8)
    DEALLOCATE(cols, vals)
END DO

!Add TokaMaker Jacobian to 8,6 block (MAYBE UNNECESSARY FOR ICOILS)
row_block = 2
col_block = 1
DO i = 1, V%i_map(row_block)%n
    jp=V%map(row_block,col_block)%ext(1,i)
    jn=V%map(row_block,col_block)%ext(2,i)
    colcount = jn-jp+1
    ALLOCATE(cols(colcount), vals(colcount))
    cols = V%lc(jp:jn)
    vals = V%M(jp:jn)
    cols = cols - V%j_map(col_block)%offset
    CALL mat%add_values([i], cols, RESHAPE(vals, [1,colcount]), &
        1, colcount, 8, 6)
    DEALLOCATE(cols, vals)
END DO

!Add TokaMaker Jacobian to 8,8 block
row_block = 2
col_block = 2
DO i = 1, V%i_map(row_block)%n
    jp=V%map(row_block,col_block)%ext(1,i)
    jn=V%map(row_block,col_block)%ext(2,i)
    colcount = jn-jp+1
    ALLOCATE(cols(colcount), vals(colcount))
    cols = V%lc(jp:jn)
    vals = V%M(jp:jn)
    cols = cols - V%j_map(col_block)%offset
    CALL mat%add_values([i], cols, RESHAPE(vals, [1,colcount]), &
        1, colcount, 8, 8)
    DEALLOCATE(cols, vals)
END DO

!---Add missing F contributions corresponding to add_f_terms 
CALL add_f_jac(self, mat)

CALL self%aug_vec%new(tmp)
CALL mat%assemble(tmp)
CALL tmp%delete()

end subroutine build_mugtok_td_jacobian


!---------------------------
! Delete simulation object 
!---------------------------
subroutine delete_mugtok_td(self)
class(oft_mugtok_td), intent(inout) :: self !< NL operator object
DEBUG_STACK_PUSH

IF(ASSOCIATED(self%I_prev))THEN
    CALL self%I_prev%delete()
    DEALLOCATE(self%I_prev)
END IF
!
IF(ASSOCIATED(self%nlfun))THEN
    CALL self%nlfun%delete()
    DEALLOCATE(self%nlfun)
END IF
!
IF(ASSOCIATED(self%rhs))THEN
    CALL self%rhs%delete()
    CALL self%tmp%delete()
    DEALLOCATE(self%rhs,self%tmp)
    NULLIFY(self%u)

    CALL self%mfmat%delete()
    DEALLOCATE(self%mfmat)
    !
    CALL self%pre%delete()
    CALL self%mf_solver%delete()
    DEALLOCATE(self%mf_solver)
    !
    CALL self%nksolver%delete()
END IF
DEBUG_STACK_POP
end subroutine

!------------------------------------------------------------------------------
!> Save the current solution (fields 1-7 of the augmented vector) to
!> a restart file 'mugtok_NNNNN.rst', indexed by the step counter
!------------------------------------------------------------------------------
subroutine mugtok_rst_save(self, t)
class(oft_mugtok_td), intent(inout) :: self
real(r8), intent(in) :: t !< Current solution time
class(oft_vector), pointer :: mug_vec
real(r8), pointer :: tmp(:)
character(LEN=TDIFF_RST_LEN) :: rst_char
integer(i4) :: i
NULLIFY(mug_vec, tmp)
!---Copy fields (1-7) out of the augmented solution vector
CALL self%mug%fe_rep%vec_create(mug_vec)
DO i=1,7
    NULLIFY(tmp)
    CALL self%u%get_local(tmp, i)
    CALL mug_vec%restore_local(tmp, i)
    CALL self%u%restore_local(tmp, i)
    NULLIFY(tmp)
END DO
!---Write via the MUG restart functionality (handles the composite FE layout)
WRITE(rst_char,'(I5.5)') self%rst_count
CALL self%mug%rst_save(mug_vec, t, self%nlfun%dt, 'mugtok_'//rst_char//'.rst', 'U')
CALL mug_vec%delete
DEALLOCATE(mug_vec)
end subroutine mugtok_rst_save

!------------------------------------------------------------------------------
!> Restore the solution (fields 1-7), time and timestep from the restart file
!> 'mugtok_NNNNN.rst' written by mugtok_rst_save, to resume a run from that step.
!> Call after setup/setup_mhd, with the same inputs (mesh, equilibrium, coil
!> currents) as the original run. Subsequent restart files continue the numbering
!> from rst_index. Not stored in the restart file, and so reset here from the
!> restored state: psi_limiting (sink history), I_prev/f_scale_prev (FF' snapshot)
!------------------------------------------------------------------------------
subroutine mugtok_rst_restore(self, rst_index, t, dt)
class(oft_mugtok_td), intent(inout) :: self
integer(i4), intent(in) :: rst_index !< Index of the restart file to load
real(r8), intent(out) :: t !< Time of the restored solution
real(r8), intent(out) :: dt !< Timestep used for the step that produced the restored solution
class(oft_vector), pointer :: mug_vec
real(r8), pointer :: tmp(:), pvac(:)
character(LEN=TDIFF_RST_LEN) :: rst_char
character(LEN=OFT_PATH_SLEN) :: file_tmp
integer(i4) :: i
NULLIFY(mug_vec, tmp, pvac)
WRITE(rst_char,'(I5.5)') rst_index
file_tmp='mugtok_'//rst_char//'.rst'
IF(.NOT.oft_file_exist(TRIM(file_tmp))) &
  CALL oft_abort('Restart file "'//TRIM(file_tmp)//'" not found','mugtok_rst_restore',__FILE__)
CALL self%mug%fe_rep%vec_create(mug_vec)
CALL self%mug%rst_load(mug_vec, TRIM(file_tmp), 'U', t, dt)
!---The saved field 6 is the total flux; remove the vacuum (coil) flux to recover
!   the plasma flux (coil currents, and so psi_vac, are those set up from the equilibrium)
CALL mug_vec%get_local(tmp, 6)
CALL self%mug%psi_vac%get_local(pvac)
tmp = tmp - pvac
CALL mug_vec%restore_local(tmp, 6)
DEALLOCATE(tmp, pvac)
!---Copy fields (1-7) into the augmented solution vector (field 8, the coil currents, is kept)
DO i=1,7
  NULLIFY(tmp)
  CALL mug_vec%get_local(tmp, i)
  CALL self%u%restore_local(tmp, i)
  DEALLOCATE(tmp)
END DO
CALL mug_vec%delete
DEALLOCATE(mug_vec)
!---Update equilibrium psi and plasma bounds from the restored flux; search for
!   the O-point from scratch since the plasma may have moved far from its initial position
self%tkmr%gs_equil%o_point(1) = -1.d0
CALL mugtok_refresh_bounds(self, self%u)
CALL mugtok_update_psi_limiting(self)
CALL snapshot_f_profile(self)
CALL xmhd_2d_freeze_eta_temp(self%mug, self%u)
CALL xmhd_2d_freeze_ohmic(self%mug, self%u)
self%rst_count = rst_index
end subroutine mugtok_rst_restore

!------------------------------------------------------------------------------
!> Plot saved restart states to XDMF/HDF5 output
!! Runtime options are set in the main input file using the group
!! `mugtok_plot_options`.
!! **Option group:** `mugtok_plot_options`
!! |  Option    |  Description  | Type [dim] |
!! |------------|---------------|------------|
!! | `rst_start=0`      | First restart-file index to read | int |
!! | `rst_end=99999`    | Last restart-file index to read (stops earlier at the first missing file) | int |
!------------------------------------------------------------------------------
subroutine mugtok_td_plot(self)
class(oft_mugtok_td), intent(inout) :: self
class(oft_vector), pointer :: ux,uy,uz,v_lag,u
type(oft_lag_bginterp) :: grad_psi
CLASS(oft_solver), POINTER :: lminv => NULL()
class(oft_matrix), pointer :: lmop => NULL()
real(r8), pointer :: plot_vals(:), plot_vec(:,:), pvac(:)
real(r8), allocatable :: r_vals(:)
integer(i4), allocatable :: cell_dofs_plot(:)
real(r8) :: f_pos(3), coord(3)
INTEGER(i4) :: ierr, io_unit, io_stat, nplotted
INTEGER(i4) :: rst_start=0, rst_end=99999, rst_cur, rst_tmp
INTEGER(i4) :: ic, jd, gd
CHARACTER(LEN=OFT_PATH_SLEN) :: file_tmp
CHARACTER(LEN=TDIFF_RST_LEN) :: rst_char
LOGICAL :: rst_exist
real(r8) :: t
TYPE(xdmf_plot_file) :: xdmf_plot, xdmf_plot_p
namelist/mugtok_plot_options/rst_start,rst_end

! Read plotting options from input file if they exist
open(NEWUNIT=io_unit,FILE=oft_env%ifile)
read(io_unit,mugtok_plot_options,IOSTAT=ierr)
close(io_unit)
!---------------------------------------------------------------------------
! Set up grad-psi -> B projection
!---------------------------------------------------------------------------
call self%mug%fe_rep%vec_create(u)
call lag_rep%vec_create(ux)
call lag_rep%vec_create(uy)
call lag_rep%vec_create(uz)
call lag_rep%vec_create(v_lag)
call lag_rep%vec_create(grad_psi%u)
NULLIFY(lmop)
call oft_blag_getmop(lag_rep,lmop)
CALL create_cg_solver(lminv)
lminv%A=>lmop
lminv%its=-2
CALL create_diag_pre(lminv%pre)
ALLOCATE(plot_vec(3,v_lag%n))
NULLIFY(plot_vals)
CALL grad_psi%setup(lag_rep)
!---Physical R coordinate of each Lagrange DOF (used for the poloidal current)
ALLOCATE(r_vals(lag_rep%ne))
ALLOCATE(cell_dofs_plot(lag_rep%nce))
DO ic=1,mesh%nc
  CALL lag_rep%ncdofs(ic, cell_dofs_plot)
  DO jd=1,lag_rep%nce
    gd = cell_dofs_plot(jd)
    CALL oft_blag_npos(lag_rep, ic, jd, f_pos)
    coord = mesh%log2phys(ic, f_pos)
    r_vals(gd) = coord(1)
  END DO
END DO
DEALLOCATE(cell_dofs_plot)
!---------------------------------------------------------------------------
! Pass 1: main fields (n, V, psi, B, poloidal current, and T if compressible)
! on the full-order mesh
!---------------------------------------------------------------------------
CALL xdmf_plot%setup("mugtok_td")
CALL mesh%setup_io(xdmf_plot,lag_rep%order)
!---Extract vacuum flux to add to solved psi for plotting
NULLIFY(pvac)
CALL self%mug%psi_vac%get_local(pvac)
110 FORMAT (I TDIFF_RST_LEN.TDIFF_RST_LEN)
nplotted=0
rst_cur=rst_start
!Loop through restart files
DO
  IF(rst_cur > rst_end) EXIT
  WRITE(rst_char,110)rst_cur
  READ(rst_char,110,IOSTAT=io_stat)rst_tmp
  IF((io_stat/=0).OR.(rst_tmp/=rst_cur)) &
    CALL oft_abort("Step count exceeds format width","mugtok_td_plot",__FILE__)
  file_tmp='mugtok_'//rst_char//'.rst'
  rst_exist=oft_file_exist(TRIM(file_tmp))
  CALL oft_mpi_barrier(ierr)
  IF(.NOT.rst_exist) EXIT ! reached the end of the saved restart files
  CALL hdf5_read(t,TRIM(file_tmp),'t')
  CALL self%mug%rst_load(u,TRIM(file_tmp),'U')
  CALL xdmf_plot%add_timestep(t)
  nplotted=nplotted+1
  !---Density
  NULLIFY(plot_vals)
  CALL u%get_local(plot_vals,1)
  plot_vals = plot_vals*self%mug%den_scale
  CALL mesh%save_vertex_scalar(plot_vals,xdmf_plot,'n')
  !---Velocity
  CALL u%get_local(plot_vals,2)
  plot_vec(1,:)=plot_vals
  CALL u%get_local(plot_vals,3)
  plot_vec(2,:)=plot_vals
  CALL u%get_local(plot_vals,4)
  plot_vec(3,:)=plot_vals
  CALL mesh%save_vertex_vector(plot_vec,xdmf_plot,'V')
  !---Temperature (compressible only)
  IF(.NOT.self%mug%incomp)THEN
    NULLIFY(plot_vals)
    CALL u%get_local(plot_vals,5)
    CALL mesh%save_vertex_scalar(plot_vals,xdmf_plot,'T')
  END IF
  !---Total psi (psi_solved + psi_vac)
  NULLIFY(plot_vals)
  CALL u%get_local(plot_vals,6)
  CALL mesh%save_vertex_scalar(plot_vals,xdmf_plot,'psi')
  !--- B from grad(psi + psi_vac) and F
  CALL grad_psi%u%restore_local(plot_vals)
  CALL grad_psi%setup(lag_rep)
  CALL oft_blag_vproject(lag_rep,grad_psi,ux,uy,uz)
  CALL v_lag%set(0.d0)
  CALL lminv%apply(v_lag,ux)
  CALL ux%add(0.d0,1.d0,v_lag)
  CALL v_lag%set(0.d0)
  CALL lminv%apply(v_lag,uy)
  CALL uy%add(0.d0,1.d0,v_lag)
  CALL uy%get_local(plot_vals)
  plot_vec(1,:)=-plot_vals
  CALL ux%get_local(plot_vals)
  plot_vec(3,:)=plot_vals
  CALL u%get_local(plot_vals,7)
  plot_vec(2,:)=plot_vals
  IF (self%mug%cyl_flag) THEN
    CALL mesh%save_vertex_vector(plot_vec,xdmf_plot,'B*R')
  ELSE
    CALL mesh%save_vertex_vector(plot_vec,xdmf_plot,'B')
  END IF
  !--- Poloidal current: grad(F) x phi_hat / R (F = field 7 = R*B_phi)
  NULLIFY(plot_vals)
  CALL u%get_local(plot_vals,7)
  CALL grad_psi%u%restore_local(plot_vals)
  CALL grad_psi%setup(lag_rep)
  CALL oft_blag_vproject(lag_rep,grad_psi,ux,uy,uz)
  CALL v_lag%set(0.d0)
  CALL lminv%apply(v_lag,ux)
  CALL ux%add(0.d0,1.d0,v_lag)
  CALL v_lag%set(0.d0)
  CALL lminv%apply(v_lag,uy)
  CALL uy%add(0.d0,1.d0,v_lag)
  CALL uy%get_local(plot_vals)
  plot_vec(1,:)=-plot_vals/(r_vals+gs_epsilon)
  CALL ux%get_local(plot_vals)
  plot_vec(3,:)=plot_vals/(r_vals+gs_epsilon)
  plot_vec(2,:)=0.d0
  CALL mesh%save_vertex_vector(plot_vec,xdmf_plot,'Jpol')
  rst_cur=rst_cur+self%rst_freq
END DO
IF(nplotted==0.AND.oft_env%head_proc) &
  WRITE(*,'(A)')'mugtok_td_plot: no mugtok_*.rst files found, nothing to plot'
CALL self%mug%psi_vac%restore_local(pvac)
NULLIFY(pvac)
!---------------------------------------------------------------------------
! Pass 2: incompressible pressure (order-1)
!---------------------------------------------------------------------------
! IF(self%mug%incomp)THEN
!   CALL xdmf_plot_p%setup("mugtok_td_p","pressure/")
!   CALL mesh%setup_io(xdmf_plot_p,lag_rep%order-1)
!   rst_cur=rst_start
!   nplotted=0
!   DO
!     IF(rst_cur > rst_end) EXIT
!     WRITE(rst_char,110)rst_cur
!     file_tmp='mugtok_'//rst_char//'.rst'
!     rst_exist=oft_file_exist(TRIM(file_tmp))
!     CALL oft_mpi_barrier(ierr)
!     IF(.NOT.rst_exist) EXIT
!     CALL hdf5_read(t,TRIM(file_tmp),'t')
!     CALL self%mug%rst_load(u,TRIM(file_tmp),'U')
!     CALL xdmf_plot_p%add_timestep(t)
!     NULLIFY(plot_vals)
!     CALL u%get_local(plot_vals,5)
!     CALL mesh%save_vertex_scalar(plot_vals,xdmf_plot_p,'p')
!     nplotted=nplotted+1
!     rst_cur=rst_cur+self%rst_freq
!   END DO
! END IF
!---Cleanup
CALL u%delete
CALL ux%delete
CALL uy%delete
CALL uz%delete
CALL v_lag%delete
DEALLOCATE(u,ux,uy,uz,v_lag)
IF(ASSOCIATED(plot_vals))DEALLOCATE(plot_vals)
IF(ASSOCIATED(plot_vec))DEALLOCATE(plot_vec)
DEALLOCATE(r_vals)
end subroutine mugtok_td_plot

!---------------------------------------------------------------------------------
!> Compute the projected poloidal gradient of F (for poloidal current calculation)
!---------------------------------------------------------------------------------
subroutine compute_gradf(self, gradf_r, gradf_z)
class(oft_mugtok_td), intent(inout) :: self
real(r8), intent(inout) :: gradf_r(:) !< dF/dR at nodes (length lag_rep%ne)
real(r8), intent(inout) :: gradf_z(:) !< dF/dZ at nodes (length lag_rep%ne)
class(oft_vector), pointer :: ux, uy, uz, v_lag
type(oft_lag_bginterp) :: grad_psi
class(oft_matrix), pointer :: lmop => NULL()
class(oft_solver), pointer :: lminv => NULL()
real(r8), pointer, dimension(:) :: fvals

!---Set up the Lagrange mass-matrix solve 
call lag_rep%vec_create(ux)
call lag_rep%vec_create(uy)
call lag_rep%vec_create(uz)
call lag_rep%vec_create(v_lag)
call lag_rep%vec_create(grad_psi%u)
call oft_blag_getmop(lag_rep, lmop)
CALL create_cg_solver(lminv)
lminv%A => lmop
lminv%its = -2
CALL create_diag_pre(lminv%pre)
!---Project grad(F)
NULLIFY(fvals)
CALL self%u%get_local(fvals, 7)
CALL grad_psi%u%restore_local(fvals)
CALL grad_psi%setup(lag_rep)
CALL oft_blag_vproject(lag_rep, grad_psi, ux, uy, uz)
CALL v_lag%set(0.d0)
CALL lminv%apply(v_lag, ux)
CALL ux%add(0.d0, 1.d0, v_lag)
CALL v_lag%set(0.d0)
CALL lminv%apply(v_lag, uy)
CALL uy%add(0.d0, 1.d0, v_lag)
!---Copy out dF/dR and dF/dZ
CALL ux%get_local(fvals)
gradf_r = fvals
CALL uy%get_local(fvals)
gradf_z = fvals
DEALLOCATE(fvals)
!---Cleanup
CALL ux%delete; CALL uy%delete; CALL uz%delete; CALL v_lag%delete
CALL grad_psi%u%delete
DEALLOCATE(ux, uy, uz, v_lag)
NULLIFY(lminv%A)
CALL lminv%delete
DEALLOCATE(lminv)
CALL lmop%delete
DEALLOCATE(lmop)
end subroutine compute_gradf

!---------------------------
! Delete matrix-free operator 
!---------------------------
subroutine delete_mfop(self)
class(oft_mugtok_td_mfop), intent(inout) :: self !< Time-advance operator
DEBUG_STACK_PUSH

self%dt=-1.d0
IF(ASSOCIATED(self%jac_op))THEN
    CALL self%jac_op%delete()
    DEALLOCATE(self%jac_op)
END IF
DEBUG_STACK_POP
end subroutine

!---------------------------------------------------------------------------------
!> Write the Newton residual F(u)-rhs at the current (failed) iterate self%u, field by field, plus
!! the iterate's V_R, V_Z and T, as one time step of the XDMF group "mugtok_resid" in resid_dump/. Fields whose
!! DOF count differs from the plotting basis are skipped. At most 8 dumps per run.
!---------------------------------------------------------------------------------
SUBROUTINE mugtok_dump_residual(self, time)
CLASS(oft_mugtok_td), INTENT(inout) :: self
REAL(r8), INTENT(in) :: time
CLASS(oft_vector), POINTER :: wrk
REAL(r8), POINTER :: vals(:)
CHARACTER(LEN=6), PARAMETER :: rnames(7) = ['R_n   ','R_VR  ','R_Vphi','R_VZ  ','R_T   ','R_psi ','R_F   ']
CHARACTER(LEN=2), PARAMETER :: unames(3) = ['VR','VZ','T ']
INTEGER(i4), PARAMETER :: ufields(3) = [2,4,5]
INTEGER(i4) :: k
REAL(r8) :: fnorm
IF(self%n_resid_dumps >= 8)RETURN
IF(.NOT.self%resid_xdmf_ready)THEN
  CALL self%resid_xdmf%setup("mugtok_resid","resid_dump/") ! own file: the end-of-run plot truncates ./oft_xdmf.*.h5
  CALL mesh%setup_io(self%resid_xdmf,lag_rep%order)
  self%resid_xdmf_ready = .TRUE.
END IF
NULLIFY(wrk, vals)
CALL self%u%new(wrk)
CALL self%nlfun%apply(self%u,wrk)
CALL wrk%add(1.d0,-1.d0,self%rhs)
fnorm = SQRT(wrk%dot(wrk))
CALL self%resid_xdmf%add_timestep(time + 1.d-9*REAL(self%n_resid_dumps,8))
DO k=1,7
  CALL wrk%get_local(vals,k)
  IF(SIZE(vals) == lag_rep%ne)CALL mesh%save_vertex_scalar(vals,self%resid_xdmf,TRIM(rnames(k)))
  DEALLOCATE(vals)
END DO
DO k=1,3
  CALL self%u%get_local(vals,ufields(k))
  IF(SIZE(vals) == lag_rep%ne)CALL mesh%save_vertex_scalar(vals,self%resid_xdmf,TRIM(unames(k)))
  DEALLOCATE(vals)
END DO
CALL wrk%delete()
DEALLOCATE(wrk)
self%n_resid_dumps = self%n_resid_dumps + 1
WRITE(*,'(A,I3,A,ES12.4,A,ES12.4)')' Residual dump',self%n_resid_dumps,' (mugtok_resid): t =',time,'  |F| =',fnorm
END SUBROUTINE mugtok_dump_residual
!---------------------------
! Update matrix-free operatior 
!---------------------------
SUBROUTINE mugtok_mfnk_update(a)
CLASS(oft_vector), TARGET, INTENT(inout) :: a !< Current solution vector
CALL current_sim%mfmat%update(a)
IF(current_sim%nl_debug)CALL mugtok_nl_debug(current_sim,a) ! NLDEBUG
END SUBROUTINE mugtok_mfnk_update
!---------------------------------------------------------------------------------
!> NLDEBUG: Print per-field norms of the iterate, the Newton residual F(u)-rhs and
!! the update since the previous iterate, plus min(n) and min(T). Called once per
!! Newton iteration (from mugtok_mfnk_update) with the accepted iterate.
!---------------------------------------------------------------------------------
SUBROUTINE mugtok_nl_debug(self,a)
CLASS(oft_mugtok_td), INTENT(inout) :: self !< Simulation object
CLASS(oft_vector), INTENT(inout) :: a !< Current Newton iterate
CHARACTER(LEN=5), PARAMETER :: fnames(8) = ['n    ','v_R  ','v_phi','v_Z  ','T    ','psi  ','F    ','coils']
REAL(r8), POINTER, DIMENSION(:) :: vals
REAL(r8) :: nrm(3), nmin, tmin
INTEGER(i4) :: k
IF(.NOT.ASSOCIATED(dbg_wrk))THEN
  CALL a%new(dbg_wrk)
  CALL a%new(dbg_uprev)
END IF
IF(oft_env%head_proc)THEN
  WRITE(*,'(A,I4,A)')'  NLDEBUG Newton it ',dbg_it,':  field      |u|           |F(u)-rhs|    |du|'
END IF
DO k=1,a%nblocks
  nrm = -1.d0
  CALL field_norm(a,k,nrm(1))
  CALL dbg_wrk%add(0.d0,1.d0,self%mfmat%f0,-1.d0,self%rhs)
  CALL field_norm(dbg_wrk,k,nrm(2))
  IF(dbg_it>0)THEN
    CALL dbg_wrk%add(0.d0,1.d0,a,-1.d0,dbg_uprev)
    CALL field_norm(dbg_wrk,k,nrm(3))
  END IF
  IF(oft_env%head_proc)WRITE(*,'(4X,A,2X,A5,3ES14.6)')'NLDEBUG',fnames(k),nrm
END DO
NULLIFY(vals)
CALL a%get_local(vals,1)
nmin = -oft_mpi_max(-MINVAL(vals))
DEALLOCATE(vals)
CALL a%get_local(vals,5)
tmin = -oft_mpi_max(-MINVAL(vals))
DEALLOCATE(vals)
IF(oft_env%head_proc)WRITE(*,'(4X,A,ES14.6,A,ES14.6)')'NLDEBUG min(n) [x den_scale] = ',nmin,'   min(T) [eV] = ',tmin
CALL dbg_uprev%add(0.d0,1.d0,a)
dbg_it = dbg_it + 1
CONTAINS
!---L2 norm of block k of v (other blocks zeroed in a copy)
SUBROUTINE field_norm(v,k,res)
CLASS(oft_vector), INTENT(inout) :: v
INTEGER(i4), INTENT(in) :: k
REAL(r8), INTENT(out) :: res
CLASS(oft_vector), POINTER :: tmp
INTEGER(i4) :: j
CALL v%new(tmp)
CALL tmp%add(0.d0,1.d0,v)
DO j=1,v%nblocks
  IF(j/=k)CALL tmp%set(0.d0,j)
END DO
res = SQRT(tmp%dot(tmp))
CALL tmp%delete
DEALLOCATE(tmp)
END SUBROUTINE field_norm
END SUBROUTINE mugtok_nl_debug


END MODULE mugtok_td