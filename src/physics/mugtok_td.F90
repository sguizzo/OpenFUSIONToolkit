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
USE oft_native_solvers, ONLY: oft_nksolver, oft_native_gmres_solver
USE oft_solver_utils, ONLY: create_cg_solver, create_diag_pre
USE oft_lu, ONLY: oft_lusolver
USE oft_la_utils, ONLY: create_matrix, graph_add_dense_blocks, create_vector, create_dense_graph
USE oft_native_la, ONLY: oft_native_matrix
!
USE fem_composite, ONLY: oft_fem_comp_type, fem_graph_create
USE fem_utils, ONLY: fem_dirichlet_diag, fem_dirichlet_vec
USE oft_lag_basis, ONLY: oft_scalar_bfem, oft_blag_eval, oft_blag_geval
USE oft_blag_operators, ONLY: oft_blag_vproject, oft_blag_getmop, oft_lag_bginterp
USE mhd_utils, ONLY: mu0
USE oft_gs, ONLY: gs_epsilon, gs_equil, gs_factory, gs_update_bounds, gs_test_bounds, flux_func, build_dels, gs_update_mu, gs_mu_of_B
USE oft_gs_td, ONLY: oft_tmaker_td_mfop, build_vac_op, apply_rhs
USE oft_mesh_local_util, ONLY: mesh_local_findedge
USE xmhd_2d, ONLY: oft_xmhd_2d_sim, build_approx_jacobian
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
TYPE(oft_lusolver), POINTER :: pre => NULL() !< Preconditioner using approximate jacobian
TYPE(oft_nksolver) :: nksolver !< Newton-Krylov solver for time-advance
TYPE(oft_fem_comp_type), POINTER :: fe_rep => NULL() !< Finite element representation for solution fields (without coil currents)
LOGICAL :: pm = .FALSE.
LOGICAL :: save_rst = .FALSE. !< If true, write restart files ('mugtok_NNNNN.rst') for plotting
INTEGER(i4) :: rst_freq = 1 !< Restart-file output frequency (in steps)
INTEGER(i4) :: rst_count = 0 !< Completed-step counter (used as the restart-file index)
LOGICAL, ALLOCATABLE, DIMENSION(:) :: plasma_flag !< True for DOFs in or on the border of region 1 (the plasma)
INTEGER(i4) :: F0_node = 0 !< Index of the first plasma_flag=.TRUE. DOF, used to evolve F0 (0 if none)
REAL(r8) :: f_scale_prev = 1.d0 !< FF' scale (tkmr%f_scale) from previous timestep
CLASS(flux_func), POINTER :: I_prev => NULL() !< FF' profile (gs_equil%I) from previous timestep
LOGICAL :: has_iron = .FALSE. !< True if any mesh region has a magnetic susceptibility set. Detected once at setup so the permeability term costs nothing on devices without iron
LOGICAL :: mu_jac_deriv = .FALSE. !< Assemble the d(mu)/d|B| term into the approximate Jacobian (add_iron_jacobian). Preconditioner-only, so it cannot change the converged solution; set false to recover the older behaviour where jac_op's only iron content was vac_op's frozen 1/mu_cell
LOGICAL :: mu_knee_damp = .FALSE. !< Pass `knee_damp` to gs_update_mu: damp the relaxation target where the B-H curve is steep, so saturated cells fall back toward mu=1 rather than tracking mu(|B|). Only meaningful when mu_relax > 0
REAL(r8) :: mu_relax = 0.d0 !< Relaxation weight in (0,1] for refreshing the cached per-cell permeability `gs_device%mu_cell` from the current solution, once per step, before the Jacobian is built. 0 (the default) leaves mu_cell frozen at its initial value, which is what happens on a device without iron. The PHYSICS does not depend on this: add_iron_terms adds dt*(1/mu(|B|) - 1/mu_cell)*stiffness, correcting the residual for whatever is cached, so any value converges to the same solution. What it changes is how much of the true operator sits in the assembled Jacobian rather than in that correction, and hence how representative the preconditioner is. Damping matters because mu falls steeply past the B-H knee (d ln mu / d ln B is about -5 at 1.5 T and -10 at 2 T), where undamped substitution can oscillate
REAL(r8), ALLOCATABLE, DIMENSION(:) :: source_term !< Prescribed toroidal current density \f$ J_{\phi} \f$ [A/m^2] in each mesh region (defaults to zero). Set by the driver before each call to `step` to impose a known, time-varying current distribution (e.g. the drive coils of an EM pump) without giving those regions to TokaMaker as circuit coils. Uses the same normalization as `gs_coil_source`, so a region carrying current I with n turns over area A takes source_term = I*n/A

contains
    !> Setup time-dependent simulation
    procedure :: setup => setup_mugtok_td
    !> Delete time-dependent simulation
    procedure :: delete => delete_mugtok_td
    !> Take a timestep
    procedure :: step => step_mugtok_td
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
CLASS(oft_bmesh), POINTER, PUBLIC :: mesh => NULL() !< Pointer to FE mesh
CLASS(oft_scalar_bfem), POINTER :: lag_rep => NULL()!< Pointer to single field FE representation 

CONTAINS

!---------------------------------------------------------
!> Setup combined MUG and TokaMaker simulation object
!---------------------------------------------------------
subroutine setup_mugtok_td(self, equil, dt,lin_tol,nl_tol, mhd_flag, dens_reg, visc_reg, eta_reg, incomp, toroidal_flow, plasma_reg)
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
INTEGER(i4), INTENT(in), optional :: plasma_reg !< Mesh region holding the Grad-Shafranov plasma (default 1). Pass 0 for a plasma-free device (e.g. an EM pump), in which case no F0/plasma constraints are applied and region 1 is free to be an ordinary MHD region.
CLASS(multigrid_mesh), POINTER :: mg_mesh
LOGICAL :: tor_flow
INTEGER(i4) :: i,j,ierr,preg
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

!---Region holding the Grad-Shafranov plasma (0 => plasma-free device)
preg = 1
IF(PRESENT(plasma_reg)) preg = plasma_reg
IF(preg > SIZE(mhd_flag)) CALL oft_abort('"plasma_reg" exceeds the number of mesh regions', &
                                         'setup_mugtok_td', __FILE__)
IF(preg > 0)THEN
    IF(mhd_flag(preg)) CALL oft_abort('Region "plasma_reg" is flagged for MHD; use setup_mhd to evolve the '// &
                                      'Grad-Shafranov plasma with MHD, or pass plasma_reg=0 for a plasma-free device', &
                                      'setup_mugtok_td', __FILE__)
END IF

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
!---Prescribed toroidal current source, zero unless the driver sets it
IF(.NOT.ALLOCATED(self%source_term))ALLOCATE(self%source_term(mesh%nreg))
self%source_term = 0.d0
!---Detect magnetic regions once, so the per-iteration permeability hook is
!   skipped entirely on devices without iron
self%has_iron = .FALSE.
IF(ASSOCIATED(equil%device%mag_suscep)) self%has_iron = ANY(equil%device%mag_suscep > -1.d98)
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
!---TokaMaker marks the plasma region with eta=-1 and vacuum regions with a large
!   positive value, so the defaults may not always be physically meaningful, check here
DO i=1,mesh%nreg
  IF(.NOT.mhd_flag(i))CYCLE
  IF(ANY(mhd_sim%eta(i,:) <= 0.d0))THEN
    WRITE(*,'(A,I0,A,2ES11.3)')'MHD region ',i,' has non-positive resistivity: ',mhd_sim%eta(i,:)
    CALL oft_abort('Resistivity must be positive in every MHD region; pass "eta_reg" explicitly '// &
                   '(TokaMaker only assigns eta to "conductor" regions)','setup_mugtok_td',__FILE__)
  END IF
END DO

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

!------------------------------------------------------------------------------
! Set up plasma flag and F0 node
!------------------------------------------------------------------------------
! Flag DOFs that lie in (or on the border of) the plasma region. When preg==0
! (plasma-free device) nothing is flagged, F0_node stays 0, and every F0/plasma
! constraint downstream (dense Jacobian block, F0 residual, plasma-node rows)
! switches itself off.
ALLOCATE(self%plasma_flag(mhd_sim%fe_rep%fields(7)%fe%ne))
self%plasma_flag = .FALSE.
IF (preg > 0) THEN
    ALLOCATE(cell_dofs(lag_rep%nce))
    DO i = 1, mesh%nc
        IF (mesh%reg(i) == preg) THEN
            CALL lag_rep%ncdofs(i, cell_dofs) ! By DOFs (order-2)
            DO j = 1, SIZE(cell_dofs)
                self%plasma_flag(cell_dofs(j)) = .TRUE.
            END DO
        END IF
    END DO
    DEALLOCATE(cell_dofs)
END IF
!---Index of the first plasma DOF, which will be used to evolve F0 (0 if none)
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

!---Re-apply the axis regularity conditions. mhd_sim%setup already ran setup_bc,
!   but the block above nullifies and reallocates every BC flag, so that work was
!   discarded and has to be redone here.
CALL mhd_sim%apply_axis_bc()

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
    !Still create these vectors even if not used to avoid errors down the road
    CALL self%fe_rep%vec_create(self%aug_vec)
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
!---Set initial F (=R*B_phi) to the F0 value from the TokaMaker equilibrium. A
!   plasma-free device carries no toroidal field (the coils drive I_phi and the
!   induced current is also toroidal), so F starts at zero and the FF' profile
!   is not required to be present at all.
IF (self%F0_node > 0) THEN
    CALL self%u%set(self%tkmr%gs_equil%I%f_offset, 7)
ELSE
    CALL self%u%set(0.d0, 7)
END IF

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

!------------------------------------------------------------------------------
! Setup nonlinear function object
!------------------------------------------------------------------------------
ALLOCATE(self%nlfun)
self%nlfun%dt=dt
self%nlfun%parent_sim => self

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
        graphs(i,j)%g%nr=fe_graphs(i,j)%g%nr           ! <-- ADD THIS
        graphs(i,j)%g%nrg=fe_graphs(i,j)%g%nrg         ! <-- ADD THIS
        graphs(i,j)%g%nc=fe_graphs(i,j)%g%nc           ! <-- ADD THIS
        graphs(i,j)%g%ncg=fe_graphs(i,j)%g%ncg         ! <-- ADD THIS
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

DO i=1,nkgraphs
    DEALLOCATE(known_graphs(i)%g)
END DO
DEALLOCATE(graphs, known_graphs)
CALL tmp_vec%delete
DEALLOCATE(tmp_vec)
DEALLOCATE(fe_graphs)


! Preconditioner should use approximate jacobian
ALLOCATE(self%pre)
self%pre%A=>self%nlfun%jac_op 

!------------------------------------------------------------------------------
! Setup matrix free solver
!------------------------------------------------------------------------------
ALLOCATE(self%mfmat) 
CALL self%mfmat%setup(self%tmp,self%nlfun)


ALLOCATE(self%mf_solver)
self%mfmat%b0=1.d-4
self%mf_solver%A=>self%mfmat
self%mf_solver%its=1000
self%mf_solver%nrits=20
self%mf_solver%atol=lin_tol
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
end subroutine setup_mugtok_td


!-----------------------------------------------------------------
!> Take one timestep of the combined MUG and TokaMaker simulation
!-----------------------------------------------------------------
subroutine step_mugtok_td(self,coil_currents,time,dt,nl_its,lin_its,nretry)
CLASS(oft_mugtok_td), target, intent(inout) :: self !< Simulation object
REAL(r8), INTENT(in) :: coil_currents(:) !< Prescribed I-coil currents for this step (size ncoils)
REAL(8), INTENT(inout) :: time,dt !< Current time and desired timestep size [s]
INTEGER(4), INTENT(out) :: nl_its,lin_its,nretry !< Number of nonlinear and linear iterations, and number of retries
INTEGER(4) :: i,j
REAL(r8), pointer :: tmp_arr(:), currs_tmp(:)
CLASS(oft_vector), pointer :: tmp_vec
CLASS(oft_native_matrix), POINTER :: P => NULL()
CLASS(oft_matrix), pointer :: pre
LOGICAL :: mu_updated

current_sim=>self

!---Refresh the preconditioner matrix with a more up-to-date mu 
!-- The extent to which it is updated depends on mu_relax
!-- mu_cell = (1-mu_relax)*mu_cell + mu_relax*mu(|B|)
mu_updated = .FALSE.
IF(self%has_iron .AND. self%mu_relax > 0.d0)THEN
    BLOCK
    REAL(r8), POINTER, DIMENSION(:) :: psi_ind, psi_vac
    REAL(r8), ALLOCATABLE :: psi_tot(:)
    NULLIFY(psi_ind, psi_vac)
    CALL self%u%get_local(psi_ind, 6)
    CALL self%mug%psi_vac%get_local(psi_vac)
    ALLOCATE(psi_tot(SIZE(psi_ind)))
    psi_tot = psi_ind + psi_vac
    CALL self%u%restore_local(psi_ind, 6)
    CALL self%mug%psi_vac%restore_local(psi_vac)
    DEALLOCATE(psi_ind, psi_vac)
    CALL self%tkmr%gs_equil%psi%restore_local(psi_tot)
    DEALLOCATE(psi_tot)
    !---rebuild_dels=.FALSE.: `dels` is the full jacobian, which this
    !   time-dependent path never applies
    IF(self%pm)THEN
        BLOCK
        REAL(r8) :: bstats(4), muchg
        CALL gs_update_mu(self%tkmr%gs_equil, relax=self%mu_relax, rebuild_dels=.FALSE., &
                          mu_change=muchg, bpol_stats=bstats, knee_damp=self%mu_knee_damp)
        WRITE(*,'(A,5(A,ES11.4))')'  iron:','  max|B_pol|=',bstats(1),'  mean|B_pol|=',bstats(2), &
            '  min mu(B)=',bstats(3),'  mean mu(B)=',bstats(4),'  max d(mu_cell)/mu=',muchg
        END BLOCK
    ELSE
        CALL gs_update_mu(self%tkmr%gs_equil, relax=self%mu_relax, rebuild_dels=.FALSE., &
                          knee_damp=self%mu_knee_damp)
    END IF
    mu_updated = .TRUE.
    END BLOCK
END IF

!Update timestep if it has changed since setup, and build approximate jacobian for preconditioning
IF(dt/=self%nlfun%dt)THEN
    dt=ABS(dt)
    self%nlfun%dt=dt
    self%mug%dt = dt
    self%tkmr%dt = dt
    !Update plasma time advance operator
    CALL self%tkmr%update()
    CALL build_mugtok_td_jacobian(self, self%nlfun%jac_op, self%u, update_vac = .TRUE.)
ELSE
    !Update plasma time advance operator
    CALL self%tkmr%update()
    !---vac_op must be rebuilt whenever mu_cell moved
    CALL build_mugtok_td_jacobian(self, self%nlfun%jac_op, self%u, update_vac = mu_updated)
END IF

!Update preconditioner
CALL self%pre%update(.TRUE.)
CALL self%mf_solver%pre%update(.TRUE.)

!Build RHS
CALL self%tmp%add(0.d0,1.d0,self%u)
CALL apply_rhs_mugtok(self%nlfun,self%u,self%rhs)

NULLIFY(currs_tmp)
nretry=0
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
    CALL self%nksolver%apply(self%u,self%rhs)
    IF(self%nksolver%cits<0)THEN
        CALL self%u%add(0.d0,1.d0,self%tmp)
        self%nlfun%dt=self%nlfun%dt/2.d0
        !---The halved timestep has to reach every sub-operator, not just nlfun.
        !   vac_op is assembled from tkmr%dt (passed as both the dt and the
        !   global scale), so leaving it stale assembles the TokaMaker psi block
        !   for a different timestep than the MUG blocks and add_iron_terms use.
        !   The resync branch above cannot catch this: step() returns
        !   dt = nlfun%dt, so the next call arrives with dt == nlfun%dt.
        self%mug%dt = self%nlfun%dt
        self%tkmr%dt = self%nlfun%dt
        CALL self%tkmr%update()
        CALL build_mugtok_td_jacobian(self, self%nlfun%jac_op, self%u, update_vac = .TRUE.)
        CALL self%pre%update(.TRUE.)
        CALL apply_rhs_mugtok(self%nlfun,self%u,self%rhs)
        nretry=nretry+1
        CYCLE
    ELSE
        EXIT
    END IF
END DO

!Advance time and return iteration counts
time=time+self%nlfun%dt
dt=self%nlfun%dt
nl_its=self%nksolver%nlits
lin_its=self%nksolver%lits

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

!---Optionally write a restart file for later plotting 
self%rst_count = self%rst_count + 1
IF(self%save_rst .AND. MOD(self%rst_count, self%rst_freq)==0)THEN
    CALL mugtok_rst_save(self, time)
END IF

!---Save the profile from this iteration to be used as the RHS in future timesteps
CALL snapshot_f_profile(self)

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

!---Apply MUG RHS operator
CALL b%set(0.d0)
CALL self%parent_sim%mug%nlfun%apply_real(a,b)

!---Add F RHS in the non-MHD regions, which are currently not computed in either code
CALL add_f_terms(self, a, b, 0.d0) !use dt = 0 to get RHS

!---Add any prescribed toroidal current source to the psi RHS. This is a known
!   quantity rather than a function of the solution, so it appears here only and
!   is deliberately absent from nlfun_apply (the Newton residual is nlfun - rhs).
!   Uses the real dt, matching the dt scaling of the LHS operator.
CALL add_source_terms(self, b, self%dt)

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

CALL b%set(0.d0)
CALL self%parent_sim%mug%nlfun%apply_real(a,b) !Apply MUG LHS

!---Add F LHS in the non-MHD regions, which are currently not computed in either code
CALL add_f_terms(self, a, b, self%dt)

!---Correct the psi stiffness for field-dependent permeability. LHS only.
CALL add_iron_terms(self, a, b, self%dt)

NULLIFY(tmp_arr1)
CALL b%get_local(tmp_arr1, 6)

!------------------------------------------------------
! Combine
!------------------------------------------------------
tmp_arr1 = tmp_arr1 + tmp_arr2
CALL b%restore_local(tmp_arr1, 6)
NULLIFY(tmp_arr2, tmp_arr1)

! Put coil current RHS into field 8. Guarded because with no coils tmp_out is a
! single-block vector and b has no field 8 (apply_rhs_mugtok guards the same copy).
IF (self%parent_sim%tkmr%gs_device%ncoils > 0) THEN
    CALL tmp_out%get_local(tmp_arr2, 2)
    CALL b%restore_local(tmp_arr2, 8)
    NULLIFY(tmp_arr2)
END IF

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
  !---Drop any contribution to rows MUG has pinned
  DO jr=1,lag_rep%nce
    IF(self%parent_sim%mug%by_bc(cell_dofs(jr)))res_loc(jr) = 0.d0
  END DO
  !---Add local values to the full residual vector
  !$omp ordered
  DO jr=1,lag_rep%nce
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
  eta_p_loc = self%parent_sim%tkmr%eta_reg(mesh%reg(cell))
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
IF(self%F0_node <= 0)RETURN !Return if plasma free device
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
!> Add the nonlinear-permeability correction to the psi (field 6) residual
!!
!! TokaMaker assembles the Delta* stiffness into `vac_op` using the cached
!! per-cell permeability `mu_cell`, and applies it as a plain matvec. That makes
!! mu invisible to the nonlinear solver: whatever value is baked into the matrix
!! is frozen for the whole solve. This routine adds the difference between the
!! true field-dependent permeability and that cached value,
!!
!!   dt * ( 1/mu(|B|) - 1/mu_cell ) * grad(phi_i).grad(psi) / R
!!
!! so that the assembled operator plus this term equals the operator that would
!! have been built with mu evaluated at the current iterate. Two consequences:
!!
!!  - The residual is exact for any `mu_cell`. Nothing has to keep the cached
!!    value fresh for correctness; a stale one only makes this term larger and
!!    the preconditioner less representative.
!!  - Because mu is evaluated from this routine's own argument, the matrix-free
!!    Jacobian differences it along with everything else and picks up d(1/mu)/dB
!!    for free. That is what lets the existing Newton solver handle the
!!    nonlinearity directly, instead of needing an outer Picard iteration.
!!
!! Only cells in magnetic regions contribute; elsewhere mu_cell is 1 and the
!! table returns 1, so the correction vanishes identically.
!!
!! The dt factor and the sign convention match the assembled operator, which
!! build_vac_op creates with scale=dt (see add_source_terms for the same scaling).
!------------------------------------------------------------------------------
subroutine add_iron_terms(self, a, b, dt)
class(oft_mugtok_td_mfop), intent(inout) :: self !< Time-advance operator
class(oft_vector), target, intent(inout) :: a !< Source field (current iterate)
class(oft_vector), intent(inout) :: b !< Result (correction accumulated into field 6)
real(r8), intent(in) :: dt !< Timestep [s]
real(r8), pointer, dimension(:) :: psi_weights, psi_res, pvac
type(oft_quad_type), pointer :: quad
integer(i4) :: i
logical :: curved
integer(i4) :: m, jr
integer(i4), allocatable :: cell_dofs(:)
real(r8) :: jac_mat(3,4), jac_det, int_factor, coords(3), dpsi(3), Bpol, nu, nu_frozen
real(r8), allocatable :: basis_grads(:,:), psi_loc(:), res_loc(:)

IF(.NOT.self%parent_sim%has_iron)RETURN

quad => lag_rep%quad
!---Work with the total flux: with coils present field 6 holds only the induced
!   part and the coil vacuum flux lives in psi_vac, but mu depends on the total.
NULLIFY(psi_weights, psi_res, pvac)
CALL a%get_local(psi_weights, 6)
IF(self%parent_sim%tkmr%gs_device%ncoils > 0)THEN
  CALL self%parent_sim%mug%psi_vac%get_local(pvac)
  psi_weights = psi_weights + pvac
  DEALLOCATE(pvac)
END IF
CALL b%get_local(psi_res, 6)
!$omp parallel private(m,jr,curved,cell_dofs,jac_mat,jac_det,int_factor,coords, &
!$omp   dpsi,Bpol,nu,nu_frozen,basis_grads,psi_loc,res_loc)
ALLOCATE(basis_grads(3,lag_rep%nce), cell_dofs(lag_rep%nce), &
         psi_loc(lag_rep%nce), res_loc(lag_rep%nce))
!$omp do ordered
DO i=1,mesh%nc
  IF(self%parent_sim%tkmr%gs_device%mag_suscep(mesh%reg(i)) < -1.d98)CYCLE
  curved=cell_is_curved(mesh,i)
  CALL lag_rep%ncdofs(i,cell_dofs)
  res_loc = 0.d0
  psi_loc = psi_weights(cell_dofs)
  nu_frozen = 1.d0/self%parent_sim%tkmr%gs_device%mu_cell(i)
  DO m=1,quad%np
    IF(curved.OR.(m==1))CALL mesh%jacobian(i,quad%pts(:,m),jac_mat,jac_det)
    DO jr=1,lag_rep%nce
      CALL oft_blag_geval(lag_rep,i,jr,quad%pts(:,m),basis_grads(:,jr),jac_mat)
    END DO
    coords = mesh%log2phys(i,quad%pts(:,m))
    int_factor = jac_det*quad%wts(m)
    dpsi = 0.d0
    DO jr=1,lag_rep%nce
      dpsi = dpsi + psi_loc(jr)*basis_grads(:,jr)
    END DO
    !---B_pol = |grad(psi)|/R, matching the convention in gs_update_mu
    Bpol = SQRT(SUM(dpsi(1:2)**2))/(coords(1)+gs_epsilon)
    nu = 1.d0/gs_mu_of_B(Bpol)
    DO jr=1,lag_rep%nce
      res_loc(jr) = res_loc(jr) + dt*(nu-nu_frozen) &
        *DOT_PRODUCT(basis_grads(1:2,jr),dpsi(1:2))*int_factor/(coords(1)+gs_epsilon)
    END DO
  END DO
  !---Boundary rows are owned by the free-boundary BC, not the PDE
  DO jr=1,lag_rep%nce
    IF(lag_rep%be(cell_dofs(jr)))res_loc(jr) = 0.d0
  END DO
  !$omp ordered
  DO jr=1,lag_rep%nce
    !$omp atomic
    psi_res(cell_dofs(jr)) = psi_res(cell_dofs(jr)) + res_loc(jr)
  END DO
  !$omp end ordered
END DO
DEALLOCATE(basis_grads,cell_dofs,psi_loc,res_loc)
!$omp end parallel
CALL b%restore_local(psi_res, 6)
DEALLOCATE(psi_res, psi_weights)
end subroutine add_iron_terms

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!> Add a prescribed toroidal current source to the psi (field 6) right-hand side
!!
!! Imposes a known \f$ J_{\phi} \f$ in the regions flagged by `source_term`,
!! letting a driver drive current in a region without handing that region to
!! TokaMaker as a circuit coil.

!! Because the source does not depend on the solution, it belongs on the RHS
!! only (called from apply_rhs_mugtok, never from nlfun_apply) and contributes
!! nothing to the Jacobian.
!------------------------------------------------------------------------------
subroutine add_source_terms(self, b, dt)
class(oft_mugtok_td_mfop), intent(inout) :: self !< Time-advance operator
class(oft_vector), intent(inout) :: b !< Result (source is accumulated into field 6)
real(r8), intent(in) :: dt !< Timestep [s]
real(r8), pointer, dimension(:) :: psi_res
type(oft_quad_type), pointer :: quad
integer(i4) :: i
logical :: curved
integer(i4) :: m, jr
integer(i4), allocatable :: cell_dofs(:)
real(r8) :: jac_mat(3,4), jac_det, int_factor, j_phi
real(r8), allocatable :: basis_vals(:), res_loc(:)

IF(.NOT.ALLOCATED(self%parent_sim%source_term))RETURN
IF(ALL(ABS(self%parent_sim%source_term) < 1.d-30))RETURN

quad => lag_rep%quad
NULLIFY(psi_res)
CALL b%get_local(psi_res, 6)
!$omp parallel private(m,jr,curved,cell_dofs,jac_mat,jac_det,int_factor,j_phi, &
!$omp   basis_vals,res_loc)
ALLOCATE(basis_vals(lag_rep%nce), cell_dofs(lag_rep%nce), res_loc(lag_rep%nce))
!$omp do ordered
DO i=1,mesh%nc
  j_phi = self%parent_sim%source_term(mesh%reg(i))
  IF(ABS(j_phi) < 1.d-30)CYCLE
  curved=cell_is_curved(mesh,i)
  CALL lag_rep%ncdofs(i,cell_dofs)
  res_loc = 0.d0
  DO m=1,quad%np
    IF(curved.OR.(m==1))CALL mesh%jacobian(i,quad%pts(:,m),jac_mat,jac_det)
    int_factor = jac_det*quad%wts(m)
    DO jr=1,lag_rep%nce
      CALL oft_blag_eval(lag_rep,i,jr,quad%pts(:,m),basis_vals(jr))
      res_loc(jr) = res_loc(jr) + basis_vals(jr)*int_factor
    END DO
  END DO
  !---Boundary rows of the psi equation are set by the free-boundary BC rather
  !   than the PDE, so a volume source must not touch them. This is the same
  !   guard build_dels applies to its local matrix. Source regions are interior
  !   in any sensible setup, so this only bites if one is placed against the
  !   domain edge.
  DO jr=1,lag_rep%nce
    IF(lag_rep%be(cell_dofs(jr)))res_loc(jr) = 0.d0
  END DO
  !---Add local values to the full residual vector
  !$omp ordered
  DO jr=1,lag_rep%nce
    !$omp atomic
    psi_res(cell_dofs(jr)) = psi_res(cell_dofs(jr)) + dt*j_phi*res_loc(jr)
  END DO
  !$omp end ordered
END DO
DEALLOCATE(basis_vals,cell_dofs,res_loc)
!$omp end parallel
CALL b%restore_local(psi_res, 6)
DEALLOCATE(psi_res)
end subroutine add_source_terms

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
    IF(self%F0_node > 0 .AND. self%plasma_flag(cell_dofs(jr))) CYCLE !In the F0 model plasma DOFs get the F0 constraint instead; without it they get diffusion (matching add_f_terms) on top of the MHD induction
    IF(self%mug%by_bc(cell_dofs(jr))) CYCLE !Row pinned by MUG (e.g. the axis); matches the same guard in add_f_terms
    CALL mat%add_values([cell_dofs(jr)], cell_dofs, RESHAPE(jac_loc(jr,:),[1,lag_rep%nce]), &
                        1, lag_rep%nce, 7, 7)
  END DO
  !$omp end critical
END DO
DEALLOCATE(basis_vals, basis_grads, jac_loc, cell_dofs)
!$omp end parallel
!------------------------------------------------------
! F0 reduced-model constraint (only when the F0 model is active)
!------------------------------------------------------
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

!---Coil-current coupling blocks. With no coils TokaMaker's vac_op is a
!   single-block matrix (build_dels only builds the 2x2 augmented form when
!   ncoils>0) and jac_op is nfields x nfields, so both V%i_map(2) and block
!   index 8 are out of bounds here.
IF(self%tkmr%gs_device%ncoils > 0)THEN
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
END IF

!---Add missing F contributions corresponding to add_f_terms 
CALL add_f_jac(self, mat)

!---Add the Jacobian of the permeability correction that add_iron_terms puts in
!   the residual. Without this jac_op's iron block is vac_op's frozen
!   dt*(1/mu_cell)*stiffness, and mu_relax>0 amplifies the omission (see
!   add_iron_jacobian). Preconditioner-only: the residual is already exact.
IF(self%mu_jac_deriv) CALL add_iron_jacobian(self, mat, a)

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
!> Plot saved restart states to XDMF/HDF5 output
!! Runtime options are set in the main input file using the group
!! `mugtok_plot_options`.
!! **Option group:** `mugtok_plot_options`
!! |  Option    |  Description  | Type [dim] |
!! |------------|---------------|------------|
!! | `rst_start=0`      | First restart-file index to read | int |
!! | `rst_end=2000`     | Last restart-file index to read | int |
!------------------------------------------------------------------------------
subroutine mugtok_td_plot(self)
class(oft_mugtok_td), intent(inout) :: self
class(oft_vector), pointer :: ux,uy,uz,v_lag,u
type(oft_lag_bginterp) :: grad_psi
CLASS(oft_solver), POINTER :: lminv => NULL()
class(oft_matrix), pointer :: lmop => NULL()
real(r8), pointer :: plot_vals(:), plot_vec(:,:), pvac(:)
INTEGER(i4) :: ierr, io_unit, io_stat, nplotted
INTEGER(i4) :: rst_start=0, rst_end=2000, rst_cur, rst_tmp
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
!---------------------------------------------------------------------------
! Pass 1: main fields (n, V, psi, B, and T if compressible) on the full-order mesh
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
  plot_vals = plot_vals + pvac
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
  rst_cur=rst_cur+self%rst_freq
END DO
IF(nplotted==0.AND.oft_env%head_proc) &
  WRITE(*,'(A)')'mugtok_td_plot: no mugtok_*.rst files found, nothing to plot'
CALL self%mug%psi_vac%restore_local(pvac)
NULLIFY(pvac)
!---------------------------------------------------------------------------
! Pass 2: incompressible pressure (order-1)
!---------------------------------------------------------------------------
IF(self%mug%incomp)THEN
  CALL xdmf_plot_p%setup("mugtok_td_p","pressure/")
  CALL mesh%setup_io(xdmf_plot_p,lag_rep%order-1)
  rst_cur=rst_start
  nplotted=0
  DO
    IF(rst_cur > rst_end) EXIT
    WRITE(rst_char,110)rst_cur
    file_tmp='mugtok_'//rst_char//'.rst'
    rst_exist=oft_file_exist(TRIM(file_tmp))
    CALL oft_mpi_barrier(ierr)
    IF(.NOT.rst_exist) EXIT
    CALL hdf5_read(t,TRIM(file_tmp),'t')
    CALL self%mug%rst_load(u,TRIM(file_tmp),'U')
    CALL xdmf_plot_p%add_timestep(t)
    NULLIFY(plot_vals)
    CALL u%get_local(plot_vals,5)
    CALL mesh%save_vertex_scalar(plot_vals,xdmf_plot_p,'p')
    nplotted=nplotted+1
    rst_cur=rst_cur+self%rst_freq
  END DO
END IF
!---Cleanup
CALL u%delete
CALL ux%delete
CALL uy%delete
CALL uz%delete
CALL v_lag%delete
DEALLOCATE(u,ux,uy,uz,v_lag)
IF(ASSOCIATED(plot_vals))DEALLOCATE(plot_vals)
IF(ASSOCIATED(plot_vec))DEALLOCATE(plot_vec)
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

!------------------------------------------------------------------------------
!> Jacobian of the field-dependent permeability correction that add_iron_terms
!> adds to the residual
!!
!! add_iron_terms puts dt*(nu(|B|) - 1/mu_cell)*stiffness into the residual, but
!! nothing put its derivative into the approximate Jacobian: jac_op's only iron
!! content is vac_op's dt*(1/mu_cell)*stiffness. That omission is what limits
!! `mu_relax`. Because jac_op is inverted exactly (oft_lusolver) and its iron
!! block scales like 1/mu_cell, jac_op^-1 scales like mu_cell there, so raising
!! mu_relax matches the retained term and simultaneously amplifies the omitted
!! one by mu -- which is why tracking mu(|B|) converges worse, not better, once
!! cells pass the B-H knee.
!!
!! Differentiating res_j = sum_q dt*(nu - nu_frozen)*(grad(phi_j).grad(psi))*w/R
!! with respect to psi_k, using |B| = |grad(psi)|/R, gives two symmetric pieces:
!!   (A) dt*(nu - nu_frozen)*(grad(phi_j).grad(phi_k))*w/R
!!       which lifts vac_op's frozen 1/mu_cell up to the true nu(|B|), and
!!   (B) -dt*(dmu/d|B|)/mu**2 * (grad(phi_j).grad(psi))*(grad(phi_k).grad(psi))
!!                            * w/(|grad(psi)|*R**2)
!!       the field dependence itself: a rank-1 outer product per quadrature
!!       point, positive definite wherever dmu/d|B| < 0. That sign is physical --
!!       saturation stiffens the iron against further flux -- and it is the term
!!       that is 5-13x the retained one past the knee.
!!
!! Preconditioner-only: the residual already carries the exact nonlinearity, so
!! this cannot change the converged solution, only how fast it is reached.
!------------------------------------------------------------------------------
subroutine add_iron_jacobian(self, mat, a)
class(oft_mugtok_td), intent(inout) :: self !< Simulation object
class(oft_matrix), pointer, intent(inout) :: mat !< Matrix to add entries to
class(oft_vector), intent(inout) :: a !< Current iterate
real(r8), pointer, dimension(:) :: psi_weights, pvac
type(oft_quad_type), pointer :: quad
integer(i4) :: i, m, jr, jc
integer(i4), allocatable :: cell_dofs(:)
real(r8) :: jac_mat(3,4), jac_det, int_factor, coords(3), dpsi(3), Bpol, absP
real(r8) :: nu, nu_frozen, mu_loc, dmu, coefA, coefB, dt, Rloc
real(r8), allocatable :: basis_grads(:,:), psi_loc(:), jac_loc(:,:), gdotp(:)
logical :: curved
real(r8), parameter :: absP_min = 1.d-12 ! below this grad(psi) the rank-1 term vanishes faster than it blows up

IF(.NOT.self%has_iron)RETURN
dt = self%nlfun%dt
quad => lag_rep%quad

!---Total flux, matching add_iron_terms: field 6 holds only the induced part
NULLIFY(psi_weights, pvac)
CALL a%get_local(psi_weights, 6)
IF(self%tkmr%gs_device%ncoils > 0)THEN
  CALL self%mug%psi_vac%get_local(pvac)
  psi_weights = psi_weights + pvac
  DEALLOCATE(pvac)
END IF

!$omp parallel private(m,jr,jc,curved,cell_dofs,jac_mat,jac_det,int_factor,coords, &
!$omp   dpsi,Bpol,absP,nu,nu_frozen,mu_loc,dmu,coefA,coefB,Rloc,basis_grads,psi_loc, &
!$omp   jac_loc,gdotp)
ALLOCATE(basis_grads(3,lag_rep%nce), cell_dofs(lag_rep%nce), psi_loc(lag_rep%nce), &
         jac_loc(lag_rep%nce,lag_rep%nce), gdotp(lag_rep%nce))
!$omp do
DO i=1,mesh%nc
  IF(self%tkmr%gs_device%mag_suscep(mesh%reg(i)) < -1.d98)CYCLE
  curved=cell_is_curved(mesh,i)
  CALL lag_rep%ncdofs(i,cell_dofs)
  psi_loc = psi_weights(cell_dofs)
  nu_frozen = 1.d0/self%tkmr%gs_device%mu_cell(i)
  jac_loc = 0.d0
  DO m=1,quad%np
    IF(curved.OR.(m==1))CALL mesh%jacobian(i,quad%pts(:,m),jac_mat,jac_det)
    DO jr=1,lag_rep%nce
      CALL oft_blag_geval(lag_rep,i,jr,quad%pts(:,m),basis_grads(:,jr),jac_mat)
    END DO
    coords = mesh%log2phys(i,quad%pts(:,m))
    Rloc = coords(1)+gs_epsilon
    int_factor = jac_det*quad%wts(m)
    dpsi = 0.d0
    DO jr=1,lag_rep%nce
      dpsi = dpsi + psi_loc(jr)*basis_grads(:,jr)
    END DO
    absP = SQRT(SUM(dpsi(1:2)**2))
    Bpol = absP/Rloc
    mu_loc = gs_mu_of_B(Bpol, dmu_dB=dmu)
    nu = 1.d0/mu_loc
    !---(A) lift the frozen 1/mu_cell to the true nu(|B|)
    coefA = dt*(nu-nu_frozen)*int_factor/Rloc
    !---(B) field dependence of mu; rank-1 in gdotp
    IF(absP > absP_min)THEN
      coefB = -dt*dmu/(mu_loc*mu_loc)*int_factor/(absP*Rloc*Rloc)
    ELSE
      coefB = 0.d0
    END IF
    DO jr=1,lag_rep%nce
      gdotp(jr) = DOT_PRODUCT(basis_grads(1:2,jr),dpsi(1:2))
    END DO
    DO jr=1,lag_rep%nce
      DO jc=1,lag_rep%nce
        jac_loc(jr,jc) = jac_loc(jr,jc) &
          + coefA*DOT_PRODUCT(basis_grads(1:2,jr),basis_grads(1:2,jc)) &
          + coefB*gdotp(jr)*gdotp(jc)
      END DO
    END DO
  END DO
  !$omp critical
  DO jr=1,lag_rep%nce
    !---Boundary rows belong to the free-boundary BC, matching add_iron_terms
    IF(lag_rep%be(cell_dofs(jr)))CYCLE
    CALL mat%add_values([cell_dofs(jr)], cell_dofs, RESHAPE(jac_loc(jr,:),[1,lag_rep%nce]), &
                        1, lag_rep%nce, 6, 6)
  END DO
  !$omp end critical
END DO
DEALLOCATE(basis_grads, cell_dofs, psi_loc, jac_loc, gdotp)
!$omp end parallel
DEALLOCATE(psi_weights)
end subroutine add_iron_jacobian

!------------------------------------------------------------------------------
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

!---------------------------
! Update matrix-free operatior 
!---------------------------
SUBROUTINE mugtok_mfnk_update(a)
CLASS(oft_vector), TARGET, INTENT(inout) :: a !< Current solution vector
CALL current_sim%mfmat%update(a)
END SUBROUTINE mugtok_mfnk_update


END MODULE mugtok_td