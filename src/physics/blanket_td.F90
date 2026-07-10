MODULE oft_blanket_td
USE oft_base
USE oft_io, ONLY: hdf5_read, hdf5_write, oft_file_exist, &
hdf5_field_exist, oft_bin_file, xdmf_plot_file, hdf5_create_file, hdf5_create_group
USE oft_quadrature
USE oft_mesh_type, ONLY: oft_bmesh, cell_is_curved
USE multigrid, ONLY: multigrid_mesh
USE oft_gauss_quadrature, ONLY: set_quad_1d
!
USE oft_la_base, ONLY: oft_vector, oft_matrix, oft_local_mat, oft_vector_ptr, &
  vector_extrapolate, oft_graph, oft_graph_ptr, map_list
USE oft_solver_utils, ONLY: create_solver_xml, create_diag_pre
USE oft_deriv_matrices, ONLY: oft_noop_matrix, oft_mf_matrix
USE oft_solver_base, ONLY: oft_solver
USE oft_native_solvers, ONLY: oft_nksolver, oft_native_gmres_solver
USE oft_solver_utils, ONLY: create_cg_solver, create_diag_pre
USE oft_lu, ONLY: oft_lusolver
USE oft_la_utils, ONLY: create_matrix, graph_add_dense_blocks, create_identity_graph, create_vector, create_dense_graph
USE oft_native_la, ONLY: oft_native_matrix, native_matrix_cast
!
USE fem_base, ONLY: oft_ml_fem_type, fem_common_linkage
USE fem_composite, ONLY: oft_fem_comp_type, fem_graph_create
USE fem_utils, ONLY: fem_dirichlet_diag, fem_dirichlet_vec, bfem_map_flag,  bfem_interp
USE oft_lag_basis, ONLY: oft_lag_setup,oft_scalar_bfem, oft_blag_eval, oft_blag_geval, oft_2D_lagrange_cast
USE oft_blag_operators, ONLY: oft_blag_vproject,oft_blag_project, oft_blag_getmop, oft_lag_bginterp
USE oft_scalar_inits, ONLY: poss_scalar_bfield
USE mhd_utils, ONLY: mu0, elec_charge, proton_mass
USE oft_gs, ONLY: gs_epsilon, build_dels, gs_equil, gs_factory, gs_update_bounds, gs_test_bounds, set_bcmat
USE oft_gs_td, ONLY: oft_tmaker_td_mfop, tMaker_td_mfnk_update, build_vac_op, apply_rhs
USE oft_mesh_local_util, ONLY: mesh_local_findedge
USE xmhd_2d, ONLY: oft_xmhd_2d_sim, build_approx_jacobian
USE oft_stitching, ONLY: oft_seam, seam_list
USE IEEE_ARITHMETIC
IMPLICIT NONE
#include "local.h"
#if !defined(TDIFF_RST_LEN)
#define TDIFF_RST_LEN 5
#endif
PRIVATE

TYPE, public :: oft_blanket_td_sim
REAL(r8) :: lin_tol = 1.d-13 !< absolute tolerance for linear solver
REAL(r8) :: nl_tol = 1.d-11 !< Needs docs
TYPE(oft_tmaker_td_mfop), POINTER :: tkmr => NULL() !< TokaMaker time-dependent object
TYPE(oft_xmhd_2d_sim), POINTER :: mug => NULL() !< MUG time-dependent object
CLASS(oft_vector), POINTER :: u => NULL() !< current solution vector
CLASS(oft_vector), POINTER :: rhs => NULL() !< Temporary RHS vector
CLASS(oft_vector), POINTER :: tmp => NULL() !< Temporary RHS vector
CLASS(oft_vector), POINTER :: aug_vec => NULL() !<  Augmented vector template
TYPE(oft_mf_matrix), POINTER :: mfmat => NULL() !< Matrix free operator
TYPE(oft_blanket_td_mfop), POINTER :: nlfun => NULL() ! !< Time-advance operator 
TYPE(oft_native_gmres_solver), POINTER :: mf_solver => NULL() !< Outer linear solver
TYPE(oft_lusolver), POINTER :: pre => NULL() !< Preconditioner using jacobian operator
TYPE(oft_nksolver) :: nksolver !< Newton-Krylov solver for time-advance
CLASS(oft_matrix), POINTER :: jacobian => NULL() !< Needs docs
INTEGER(i4), CONTIGUOUS, POINTER, DIMENSION(:,:) :: jacobian_block_mask => NULL() !< Matrix block mask
TYPE(oft_fem_comp_type), POINTER :: fe_rep => NULL() !< Finite element representation for solution field
LOGICAL :: pm = .TRUE.
LOGICAL :: save_rst = .FALSE. !< If true, write restart files ('blanket_NNNNN.rst') for plotting
INTEGER(i4) :: rst_freq = 1 !< Restart-file output frequency (in steps)
INTEGER(i4) :: rst_count = 0 !< Completed-step counter (used as the restart-file index)
LOGICAL, ALLOCATABLE, DIMENSION(:) :: plasma_flag !< True for By-field DOFs in or on the border of region 1 (the plasma)
INTEGER(i4) :: F0_node = 0 !< Index of the first plasma_flag=.TRUE. DOF (0 if none)

contains
    !> Apply the matrix
    procedure :: setup => setup_blanket_td
    !> Needs Docs
    procedure :: delete => delete_blanket_td
    !> Needs Docs
    procedure :: step => step_blanket_td
    !> Write XDMF/HDF5 plot files from saved restart files
    procedure :: plot => blanket_td_plot
END TYPE oft_blanket_td_sim

type, extends(oft_noop_matrix) ::oft_blanket_td_mfop
    real(r8) :: dt = 1.E-3 !< Time step size [s]
    CLASS(oft_matrix), POINTER :: jac_op => NULL() !< Time-advance operator
    TYPE(oft_blanket_td_sim), pointer :: parent_sim => NULL() !< pointer to parent simulation object for access to parameters
contains
    !> Apply operator
    procedure :: apply_real => nlfun_apply
    !> Needs Docs
    procedure :: delete => delete_mfop
end type oft_blanket_td_mfop

TYPE(oft_blanket_td_sim), POINTER :: current_sim => NULL()
CLASS(oft_bmesh), POINTER, PUBLIC :: mesh => NULL()
CLASS(oft_scalar_bfem), POINTER :: lag_rep => NULL()

CONTAINS

subroutine setup_blanket_td(self, mg_mesh, equil, dt,lin_tol,nl_tol, mhd_flag, dens_reg, visc_reg, eta_reg, incomp)
CLASS(oft_blanket_td_sim), INTENT(inout), TARGET :: self
CLASS(multigrid_mesh), INTENT(in) :: mg_mesh
TYPE(gs_equil), INTENT(inout), TARGET :: equil
REAL(8), INTENT(in) :: dt !< Needs Docs
REAL(8), INTENT(in) :: lin_tol !< Needs Docs
REAL(8), INTENT(in) :: nl_tol !< Needs Docs
INTEGER(i4) :: i,j
LOGICAL, INTENT(in) :: mhd_flag(:)
REAL(8), INTENT(in) :: dens_reg(:), visc_reg(:)
REAL, INTENT(in), optional :: eta_reg(:,:)
LOGICAL, INTENT(in), optional :: incomp
TYPE(oft_tmaker_td_mfop), POINTER :: tok_sim
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
CLASS(oft_vector), pointer :: tmp_vec, psi_coil
write(*,*) 111
mesh => equil%device%mesh
lag_rep=>equil%device%fe_rep
!------------------------------------------------------------------------------
! Set up TokaMaker and MUG objects
!------------------------------------------------------------------------------
ALLOCATE(self%tkmr)
equil%device%ignore_rmask = mhd_flag

self%tkmr%dt=dt
CALL self%tkmr%setup(equil)
CALL build_vac_op(self%tkmr,self%tkmr%vac_op)

ALLOCATE(mhd_sim)
ALLOCATE(mhd_sim%ignore_rmask(mesh%nreg))
mhd_sim%ignore_rmask = .NOT. mhd_flag
ALLOCATE(mhd_sim%eta(mesh%nreg, 2))
ALLOCATE(mhd_sim%m_i(mesh%nreg))
ALLOCATE(mhd_sim%nu(mesh%nreg))
mhd_sim%m_i = dens_reg
mhd_sim%nu = visc_reg
IF(PRESENT(eta_reg)) THEN
    mhd_sim%eta = eta_reg
ELSE
    mhd_sim%eta(:,1) = self%tkmr%eta_reg
    mhd_sim%eta(:,2) = self%tkmr%eta_reg
END IF


IF (PRESENT(incomp)) THEN
    mhd_sim%incomp = incomp
ELSE
    mhd_sim%incomp = .TRUE.
END IF
mhd_sim%cyl_flag = .TRUE.
mhd_sim%dt = dt
mhd_sim%den_scale = 1.d0

CALL mhd_sim%setup(mg_mesh, lag_rep%order, fe_rep_in =lag_rep)

!------------------------------------------------------------------------------
! Set up plasma flag and F0 node
!------------------------------------------------------------------------------
! Flag By-field DOFs that lie in (or on the border of) region 1, the plasma:
! loop cells and, for any cell in region 1, mark all of its By (order-2) DOFs.
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
!---Index of the first plasma DOF (0 if region 1 has no DOFs)
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
mhd_sim%by_bc = .FALSE.  !-> deal with by later
mhd_sim%velx_bc = .FALSE. 
mhd_sim%vely_bc = .TRUE. 
mhd_sim%velz_bc = .FALSE. 

mhd_sim%n_bc = .TRUE.
mhd_sim%T_bc = .TRUE.


ALLOCATE(cell_dofs(lag_rep%nce))
ALLOCATE(cell_dofs_p(mhd_sim%fe_rep%fields(5)%fe%nce)) ! pressure lives on its own (order-1) space
DO i = 1, mesh%nc
    call mhd_sim%fe_rep%fields(5)%fe%ncdofs(i,cell_dofs_p) ! pressure DOFs (order-1)
    call lag_rep%ncdofs(i,cell_dofs) ! velocity/psi/n DOFs (order-2)
    !---Order-2 fields (psi, velocity, density): loop over the order-2 cell DOFs
    DO j=1, SIZE(cell_dofs)
        IF (mhd_flag(mesh%reg(i))) THEN
            IF (.NOT. mhd_sim%incomp) mhd_sim%n_bc(cell_dofs(j)) = .FALSE. !add contributions everywhere in MHD regions, including boundaries, if incompressible
        ELSE
            mhd_sim%velx_bc(cell_dofs(j)) = .TRUE. !pin velocity to zero everywhere in non-MHD regions, including boundaries
            mhd_sim%vely_bc(cell_dofs(j)) = .TRUE.
            mhd_sim%velz_bc(cell_dofs(j)) = .TRUE.
        END IF
    END DO
    !---Pressure (order-1) must be looped over its OWN cell-DOF count, not the
    ! order-2 count, or cell_dofs_p is indexed past the filled entries (garbage
    ! indices -> spurious unpinned pressure rows -> singular preconditioner).
    IF (mhd_flag(mesh%reg(i))) THEN
        DO j=1, SIZE(cell_dofs_p)
            mhd_sim%T_bc(cell_dofs_p(j)) = .FALSE. !free pressure everywhere in MHD regions
        END DO
    END IF
END DO

IF (mhd_sim%incomp) THEN
    ALLOCATE(p_dir_set(mesh%nreg))
    p_dir_set = .FALSE.
    DO i=1, mesh%nc
        IF (mhd_flag(mesh%reg(i)) .AND. .NOT. p_dir_set(mesh%reg(i))) THEN
            call mhd_sim%fe_rep%fields(5)%fe%ncdofs(i,cell_dofs_p) ! Get global index of local DOFs
            mhd_sim%T_bc(cell_dofs_p(1)) = .TRUE.
            p_dir_set(mesh%reg(i)) = .TRUE.
            write(*,*) "Pinning node ", cell_dofs_p(1), " in region ", mesh%reg(i)
        END IF
    END DO
    DEALLOCATE(p_dir_set)
END IF
DEALLOCATE(cell_dofs, cell_dofs_p)

self%mug => mhd_sim

!------------------------------------------------------------------------------
! Create Solver fields, augmented with coils
!------------------------------------------------------------------------------
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
CALL self%u%set(1.d3, 5)
CALL self%u%set(self%tkmr%gs_equil%I%f_offset, 7)

NULLIFY(tmp_arr, vals_out)
CALL self%tkmr%gs_device%fe_rep%vec_create(tmp_vec)
CALL self%tkmr%gs_device%fe_rep%vec_create(psi_coil)
CALL psi_coil%set(0.d0)
CALL tmp_vec%add(0.d0,1.d0,self%tkmr%gs_equil%psi)
IF (self%tkmr%gs_device%ncoils > 0) THEN
    CALL self%u%get_local(vals_out,8)
    DO i=1,self%tkmr%gs_device%ncoils
        CALL tmp_vec%add(1.d0,-self%tkmr%gs_equil%coil_currs(i),self%tkmr%gs_device%psi_coil(i)%f)
        CALL psi_coil%add(1.d0,self%tkmr%gs_equil%coil_currs(i),self%tkmr%gs_device%psi_coil(i)%f)
        vals_out(i)=self%tkmr%gs_equil%coil_currs(i)
    END DO
    CALL self%u%restore_local(vals_out,8)
    DEALLOCATE(vals_out)
END IF
CALL tmp_vec%get_local(tmp_arr)
CALL self%u%restore_local(tmp_arr,6)
CALL psi_coil%get_local(tmp_arr)
CALL self%mug%psi_vac%restore_local(tmp_arr)



!------------------------------------------------------------------------------
! Setup nl_fun object
!------------------------------------------------------------------------------
ALLOCATE(self%nlfun)
self%nlfun%dt=dt
self%nlfun%parent_sim => self

!------------------------------------------------------------------------------
! Build Jacobian matrix
!------------------------------------------------------------------------------

ALLOCATE(self%jacobian_block_mask(self%fe_rep%nfields,self%fe_rep%nfields))
self%jacobian_block_mask=1
CALL fem_graph_create(self%fe_rep, fe_graphs, known_graphs, nkgraphs, self%jacobian_block_mask)

IF (self%tkmr%gs_device%ncoils > 0) THEN
    ALLOCATE(graphs(self%fe_rep%nfields+1,self%fe_rep%nfields+1))
ELSE
    ALLOCATE(graphs(self%fe_rep%nfields,self%fe_rep%nfields))
END IF

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

! Add dense regions to psi block for boundary
ALLOCATE(bc_nodes(1))
bc_nodes(1)%n = lag_rep%nbe
bc_nodes(1)%v => lag_rep%lbe
ALLOCATE(dense_flag(lag_rep%ne))
dense_flag = 0
dense_flag(bc_nodes(1)%v) = 1
!---Add dense blocks
CALL graph_add_dense_blocks(graphs(6,6)%g,dense_graph,dense_flag,bc_nodes)
NULLIFY(graphs(6,6)%g%kr,graphs(6,6)%g%lc)
graphs(6,6)%g%nnz=dense_graph%nnz
graphs(6,6)%g%kr=>dense_graph%kr
graphs(6,6)%g%lc=>dense_graph%lc
DEALLOCATE(dense_flag, bc_nodes)

!---Add F0 coupling to the by (7,7) block: every plasma DOF is coupled to the
! F0 node (a single-node dense group {F0} appended to each plasma row).
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

IF(self%tkmr%gs_device%ncoils > 0)THEN
    CALL create_dense_graph(graphs(8,6)%g,self%tkmr%gs_device%coil_vec,tmp_vec)
    CALL create_dense_graph(graphs(6,8)%g,tmp_vec,self%tkmr%gs_device%coil_vec)
    CALL create_dense_graph(graphs(8,8)%g,self%tkmr%gs_device%coil_vec,self%tkmr%gs_device%coil_vec)
END IF 

CALL create_matrix(self%nlfun%jac_op,graphs,self%aug_vec,self%aug_vec)

DO i=1,nkgraphs
    DEALLOCATE(known_graphs(i)%g)
END DO
DEALLOCATE(graphs, known_graphs)
CALL tmp_vec%delete
CALL psi_coil%delete
DEALLOCATE(tmp_vec)
DEALLOCATE(psi_coil)
DEALLOCATE(fe_graphs)


! Preconditioner should use approximate jacobian
ALLOCATE(self%pre) !CHECKBACK
self%pre%A=>self%nlfun%jac_op 

!------------------------------------------------------------------------------
! Setup matrix free solver
!------------------------------------------------------------------------------
ALLOCATE(self%mfmat) 
CALL self%mfmat%setup(self%tmp,self%nlfun)


ALLOCATE(self%mf_solver)
self%mfmat%b0=1.d-5
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
self%nksolver%J_update=>blanket_mfnk_update
self%nksolver%up_freq=1


end subroutine setup_blanket_td



subroutine step_blanket_td(self,time,dt,nl_its,lin_its,nretry)
CLASS(oft_blanket_td_sim), target, intent(inout) :: self !< NL operator object
REAL(8), INTENT(inout) :: time,dt
INTEGER(4), INTENT(out) :: nl_its,lin_its,nretry
INTEGER(4) :: i,j,k,ierr
REAL(r8), pointer :: tmp_arr(:), tmp_arr_2(:)
REAL(r8) :: res
CLASS(oft_vector), pointer :: tmp_vec
CLASS(oft_native_matrix), POINTER :: P => NULL()
CLASS(oft_matrix), pointer :: pre

current_sim=>self

!Update plasma time advance operator
CALL self%tkmr%update()

!Update timestep if needed
IF(dt/=self%nlfun%dt)THEN
    dt=ABS(dt)
    self%nlfun%dt=dt
    self%mug%dt = dt
    CALL build_blankettd_jacobian(self, self%nlfun%jac_op, self%u, update_vac = .TRUE.)
ELSE
    CALL build_blankettd_jacobian(self, self%nlfun%jac_op, self%u, update_vac = .FALSE.)
END IF

!Update preconditioner

CALL self%pre%update(.TRUE.)
CALL self%mf_solver%pre%update(.TRUE.)

!Build RHS and apply boundary conditions
NULLIFY(tmp_arr)
NULLIFY(tmp_arr_2)
CALL self%tmp%add(0.d0,1.d0,self%u)
CALL apply_rhs_blanket(self%nlfun,self%u,self%rhs)
CALL self%rhs%get_local(tmp_arr,6)

CALL self%nlfun%apply_real(self%u,self%tmp)
CALL self%tmp%get_local(tmp_arr_2,6)

CALL self%rhs%restore_local(tmp_arr,6)
CALL self%tmp%restore_local(tmp_arr_2,6)

DO j = 1,4
    CALL self%nksolver%apply(self%u,self%rhs)
    IF(self%nksolver%cits<0)THEN
        CALL self%u%add(0.d0,1.d0,self%tmp)
        self%nlfun%dt=self%nlfun%dt/2.d0
        CALL build_blankettd_jacobian(self, self%nlfun%jac_op, self%u, update_vac = .TRUE.)
        CALL self%pre%update(.TRUE.)
        ! Call apply_rhs_blanket which now properly manages its own temporary vectors
        CALL apply_rhs_blanket(self%nlfun,self%u,self%rhs)
        CYCLE
    ELSE
        EXIT
    END IF
END DO

time=time+self%nlfun%dt
dt=self%nlfun%dt
nl_its=self%nksolver%nlits
lin_its=self%nksolver%lits

!---Optionally write a restart file for later plotting
self%rst_count = self%rst_count + 1
IF(self%save_rst .AND. MOD(self%rst_count, self%rst_freq)==0)THEN
    CALL blanket_rst_save(self, time)
END IF

end subroutine step_blanket_td

!------------------------------------------------------------------------------
!> Save the current coupled solution (MUG fields 1-7 of the augmented vector) to
!> a restart file 'blanket_NNNNN.rst', indexed by the completed-step counter.
!! Only the MUG composite fields are stored (the coil currents, field 8, are not
!! spatial fields); this is what blanket_td_plot reads back for visualization.
!------------------------------------------------------------------------------
subroutine blanket_rst_save(self, t)
class(oft_blanket_td_sim), intent(inout) :: self
real(r8), intent(in) :: t !< Current solution time
class(oft_vector), pointer :: mug_vec
real(r8), pointer :: tmp(:)
character(LEN=TDIFF_RST_LEN) :: rst_char
integer(i4) :: i
NULLIFY(mug_vec, tmp)
!---Copy MUG fields (1-7) out of the augmented solution vector
CALL self%mug%fe_rep%vec_create(mug_vec)
DO i=1,7
    NULLIFY(tmp)
    CALL self%u%get_local(tmp, i)
    CALL mug_vec%restore_local(tmp, i)
    CALL self%u%restore_local(tmp, i) ! return checkout to the source vector
    NULLIFY(tmp)
END DO
!---Write via the MUG restart machinery (handles the composite FE layout)
WRITE(rst_char,'(I5.5)') self%rst_count
CALL self%mug%rst_save(mug_vec, t, self%nlfun%dt, 'blanket_'//rst_char//'.rst', 'U')
CALL mug_vec%delete
DEALLOCATE(mug_vec)
end subroutine blanket_rst_save

!------------------------------------------------------------------------------
!> Plot saved restart states to XDMF/HDF5, modeled on xmhd_2d_plot.
!! Automatically finds and plots every 'blanket_NNNNN.rst' file in the working
!! directory (written when save_rst=.TRUE.), writing density, velocity, psi and
!! poloidal B to the root 'oft_xdmf.*' (group 'blanket_td'). For incompressible
!! runs the (order-1) pressure is written under the 'pressure/' subdirectory
!! (group 'blanket_td_p'), mirroring xmhd_2d_plot's separate pressure output;
!! run build_xdmf.py in that subdirectory to assemble it. Takes no options.
!------------------------------------------------------------------------------
subroutine blanket_td_plot(self)
class(oft_blanket_td_sim), intent(inout) :: self
class(oft_vector), pointer :: ux,uy,uz,v_lag,u
type(oft_lag_bginterp) :: grad_psi
CLASS(oft_solver), POINTER :: lminv => NULL()
class(oft_matrix), pointer :: lmop => NULL()
real(r8), pointer :: plot_vals(:), plot_vec(:,:), pvac(:)
INTEGER(i4) :: ierr, io_unit, io_stat, nfiles, k
CHARACTER(LEN=OFT_PATH_SLEN) :: line, listfile
CHARACTER(LEN=OFT_PATH_SLEN), ALLOCATABLE :: file_list(:)
real(r8) :: t
TYPE(xdmf_plot_file) :: xdmf_plot, xdmf_plot_p
!---------------------------------------------------------------------------
! Discover every blanket_NNNNN.rst file in the working directory (sorted)
!---------------------------------------------------------------------------
listfile='.blanket_rst_list.tmp'
IF(oft_env%head_proc) &
  CALL EXECUTE_COMMAND_LINE('ls blanket_?????.rst 2>/dev/null | sort > '//TRIM(listfile))
CALL oft_mpi_barrier(ierr)
nfiles=0
OPEN(NEWUNIT=io_unit,FILE=TRIM(listfile),STATUS='old',IOSTAT=io_stat)
IF(io_stat==0)THEN
  DO
    READ(io_unit,'(A)',IOSTAT=io_stat) line
    IF(io_stat/=0) EXIT
    IF(LEN_TRIM(line)>0) nfiles=nfiles+1
  END DO
  ALLOCATE(file_list(MAX(nfiles,1)))
  REWIND(io_unit)
  k=0
  DO
    READ(io_unit,'(A)',IOSTAT=io_stat) line
    IF(io_stat/=0) EXIT
    IF(LEN_TRIM(line)>0)THEN
      k=k+1
      file_list(k)=line
    END IF
  END DO
  CLOSE(io_unit,STATUS='delete')
END IF
IF(nfiles==0)THEN
  IF(oft_env%head_proc)WRITE(*,'(A)')'blanket_td_plot: no blanket_*.rst files found, nothing to plot'
  IF(ALLOCATED(file_list))DEALLOCATE(file_list)
  RETURN
END IF
IF(oft_env%head_proc)WRITE(*,'(A,I0,A)')'blanket_td_plot: plotting ',nfiles,' restart file(s)'
!---------------------------------------------------------------------------
! Create working fields (psi-gradient -> B projection uses an L2 mass solve)
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
! Pass 1: order-2 fields (n, V, psi, B, and T if compressible) on 'mesh'
!---------------------------------------------------------------------------
CALL xdmf_plot%setup("blanket_td")
CALL mesh%setup_io(xdmf_plot,lag_rep%order)
!---Vacuum poloidal flux offset (constant in time); added to psi for the total flux
NULLIFY(pvac)
CALL self%mug%psi_vac%get_local(pvac)
DO k=1,nfiles
  CALL hdf5_read(t,TRIM(file_list(k)),'t')
  CALL self%mug%rst_load(u,TRIM(file_list(k)),'U')
  CALL xdmf_plot%add_timestep(t)
  !---Density
  NULLIFY(plot_vals)
  CALL u%get_local(plot_vals,1)
  plot_vals = plot_vals*self%mug%den_scale
  CALL mesh%save_vertex_scalar(plot_vals,xdmf_plot,'n')
  !---Velocity
  CALL u%get_local(plot_vals,2)
  plot_vec(1,:)=plot_vals
  CALL u%get_local(plot_vals,3)
  plot_vec(3,:)=plot_vals
  CALL u%get_local(plot_vals,4)
  plot_vec(2,:)=plot_vals
  CALL mesh%save_vertex_vector(plot_vec,xdmf_plot,'V')
  !---Temperature (compressible only; incompressible pressure is done in pass 2)
  IF(.NOT.self%mug%incomp)THEN
    NULLIFY(plot_vals)
    CALL u%get_local(plot_vals,5)
    CALL mesh%save_vertex_scalar(plot_vals,xdmf_plot,'T')
  END IF
  !---Total poloidal flux psi + psi_vac
  NULLIFY(plot_vals)
  CALL u%get_local(plot_vals,6)
  plot_vals = plot_vals + pvac
  CALL mesh%save_vertex_scalar(plot_vals,xdmf_plot,'psi')
  !---Poloidal B from grad(psi + psi_vac) (L2-projected onto the Lagrange space)
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
  plot_vec(2,:)=plot_vals
  CALL u%get_local(plot_vals,7)
  plot_vec(3,:)=plot_vals
  IF (self%mug%cyl_flag) THEN
    CALL mesh%save_vertex_vector(plot_vec,xdmf_plot,'B*R')
  ELSE
    CALL mesh%save_vertex_vector(plot_vec,xdmf_plot,'B')
  END IF
END DO
CALL self%mug%psi_vac%restore_local(pvac)
NULLIFY(pvac)
!---------------------------------------------------------------------------
! Pass 2: the incompressible pressure lives on the order-1 space. Plot it in a
! separate pass so the order-1 mesh IO (the 'mesh_p' role) does not clash with
! the order-2 tessellation used above; written to its own file like xmhd_2d.
!---------------------------------------------------------------------------
IF(self%mug%incomp)THEN
  ! Separate output DIRECTORY: xmdf_setup truncates oft_xdmf.*.h5, so a second
  ! plot object in the same directory would wipe the main-field file above.
  CALL xdmf_plot_p%setup("blanket_td_p","pressure/")
  CALL mesh%setup_io(xdmf_plot_p,lag_rep%order-1)
  DO k=1,nfiles
    CALL hdf5_read(t,TRIM(file_list(k)),'t')
    CALL self%mug%rst_load(u,TRIM(file_list(k)),'U')
    CALL xdmf_plot_p%add_timestep(t)
    NULLIFY(plot_vals)
    CALL u%get_local(plot_vals,5)
    CALL mesh%save_vertex_scalar(plot_vals,xdmf_plot_p,'p')
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
DEALLOCATE(file_list)
end subroutine blanket_td_plot

subroutine apply_rhs_blanket(self, a, b)
class(oft_blanket_td_mfop), intent(inout) :: self
class(oft_vector), target, intent(inout) :: a !< Source field
class(oft_vector), intent(inout) :: b !< Result of metric function
class(oft_vector), pointer :: tmp_in, tmp_out !< Temporary vectors
REAL(r8), POINTER, DIMENSION(:) :: tmp_arr1, tmp_arr2

self%parent_sim%mug%nlfun%dt = 0.d0
CALL b%set(0.d0)
CALL self%parent_sim%mug%nlfun%apply_real(a,b)
!---Add By(=F) diffusion in the non-MHD regions (dt=0 -> mass/RHS term only)
CALL add_f_diff(self, a, b, 0.d0)
NULLIFY(tmp_arr1)
NULLIFY(tmp_arr2)
CALL b%get_local(tmp_arr1, 6) 
! where (self%parent_sim%mug%psi_bc)
!     tmp_arr1 = 0.0
! end where

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
CALL a%get_local(tmp_arr2, 6)
CALL tmp_in%restore_local(tmp_arr2, 1)
CALL a%restore_local(tmp_arr2, 6)  ! Restore to a before reusing tmp_arr2
NULLIFY(tmp_arr2)

! Extract component 8 from a and put into component 2 of tmp_in (if coils exist)
IF (self%parent_sim%tkmr%gs_device%ncoils > 0) THEN
    CALL a%get_local(tmp_arr2, 8)
    CALL tmp_in%restore_local(tmp_arr2, 2)
    CALL a%restore_local(tmp_arr2, 8)  ! Restore to a before proceeding
    NULLIFY(tmp_arr2)
END IF

CALL apply_rhs(self%parent_sim%tkmr, tmp_in, tmp_out)
CALL tmp_out%get_local(tmp_arr2, 1)
CALL tmp_out%restore_local(tmp_arr2, 1)
CALL self%parent_sim%tkmr%gs_device%zerob_bc%apply(tmp_out)

! Extract component 1 from tmp_out
CALL tmp_out%get_local(tmp_arr2, 1)
! tmp_arr1 = 0.d0
! write(*,*) 'mug cont to RHS: ' , tmp_arr1(21202)
! write(*,*) 'tok cont to RHS: ' , tmp_arr2(21202)
tmp_arr1 = tmp_arr1 + tmp_arr2
CALL b%restore_local(tmp_arr1, 6)
CALL tmp_out%restore_local(tmp_arr2, 1)  ! Restore to tmp_out

NULLIFY(tmp_arr2)

! Extract component 2 from tmp_out (if coils exist)
IF (self%parent_sim%tkmr%gs_device%ncoils > 0) THEN
    CALL tmp_out%get_local(tmp_arr2, 2)
    CALL b%restore_local(tmp_arr2, 8)
    CALL tmp_out%restore_local(tmp_arr2, 2)  ! Restore to tmp_out
    NULLIFY(tmp_arr2)
END IF

! Cleanup
CALL tmp_in%delete()
CALL tmp_out%delete()
DEALLOCATE(tmp_in, tmp_out)
NULLIFY(tmp_arr1)
end subroutine apply_rhs_blanket

subroutine nlfun_apply(self, a, b)
class(oft_blanket_td_mfop), intent(inout) :: self
class(oft_vector), target, intent(inout) :: a !< Source field
class(oft_vector), intent(inout) :: b !< Result of metric function
class(oft_vector), pointer :: tmp_in, tmp_out !< Temporary vectors
REAL(r8), POINTER, DIMENSION(:) :: tmp_arr1, tmp_arr2

self%parent_sim%mug%nlfun%dt = self%dt
CALL b%set(0.d0)
CALL self%parent_sim%mug%nlfun%apply_real(a,b)
!---Add By(=F) diffusion in the non-MHD regions (mass + dt*diffusion)
CALL add_f_diff(self, a, b, self%dt)
NULLIFY(tmp_arr1)
NULLIFY(tmp_arr2)
CALL b%get_local(tmp_arr1, 6)
! tmp_arr1 = tmp_arr1
! where (self%parent_sim%mug%psi_bc)
!     tmp_arr1 = 0.0
! end where
! Extract component 6 from a
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
CALL tmp_in%restore_local(tmp_arr2, 1)
CALL a%restore_local(tmp_arr2, 6)  ! Restore to a before reusing tmp_arr2
NULLIFY(tmp_arr2)

IF (self%parent_sim%tkmr%gs_device%ncoils > 0) THEN
    CALL a%get_local(tmp_arr2, 8)
    CALL tmp_in%restore_local(tmp_arr2, 2)
    CALL a%restore_local(tmp_arr2, 8)  ! Restore to a
    NULLIFY(tmp_arr2)
END IF

CALL self%parent_sim%tkmr%apply_real(tmp_in, tmp_out) 
! Extract component 1 from tmp_out
CALL tmp_out%get_local(tmp_arr2, 1)
! tmp_arr1 = 0.d0
! write(*,*) 'mug cont to LHS: ' , tmp_arr1(21202)
! write(*,*) 'tok cont to LHS: ' , tmp_arr2(21202)
tmp_arr1 = tmp_arr1 + tmp_arr2
CALL b%restore_local(tmp_arr1, 6)
CALL tmp_out%restore_local(tmp_arr2, 1)  ! Restore to tmp_out
NULLIFY(tmp_arr2)

! Extract component 2 from tmp_out
CALL tmp_out%get_local(tmp_arr2, 2)
CALL b%restore_local(tmp_arr2, 8)
CALL tmp_out%restore_local(tmp_arr2, 2)  ! Restore to tmp_out
NULLIFY(tmp_arr2)

! Cleanup
CALL tmp_in%delete()
CALL tmp_out%delete()
DEALLOCATE(tmp_in, tmp_out)
NULLIFY(tmp_arr1)
end subroutine nlfun_apply

!------------------------------------------------------------------------------
!> Add the By(=F) diffusion residual in the NON-MHD regions (mhd_flag=.FALSE.),
!! which the MUG solve skips via ignore_rmask. Cylindrical only; mirrors the
!! xmhd_2d By residual but keeps ONLY the mass and resistive-diffusion terms
!! (no velocity/fluid terms):
!!   res(jr,7) += basis(jr)*by/R + dt*eta*grad(basis(jr)).dby/R
!! Called with dt=self%dt from nlfun_apply and dt=0 from apply_rhs_blanket, so
!! the by/R mass term forms a consistent backward-Euler time derivative.
!------------------------------------------------------------------------------
subroutine add_f_diff(self, a, b, dt)
class(oft_blanket_td_mfop), intent(inout) :: self
class(oft_vector), target, intent(inout) :: a !< Source field (uses By, field 7)
class(oft_vector), intent(inout) :: b !< Result; By diffusion is ADDED to field 7
real(r8), intent(in) :: dt !< Timestep (0 for the RHS/mass-only pass)
real(r8), pointer, dimension(:) :: by_weights, by_res
type(oft_quad_type), pointer :: quad
integer(i4) :: i
real(r8) :: eta_fallback
!---Large-but-finite resistivity (0.1 ohm-m) used where the region eta is
! undefined/negative (vacuum), converted to eta_reg (magnetic diffusivity) units.
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
ALLOCATE(basis_vals(lag_rep%nce), basis_grads(3,lag_rep%nce))
ALLOCATE(by_weights_loc(lag_rep%nce), cell_dofs(lag_rep%nce), res_loc(lag_rep%nce))
DO i=1,mesh%nc
  !---Only the non-MHD regions (those the MUG solve ignores)
  IF(.NOT.self%parent_sim%mug%ignore_rmask(mesh%reg(i)))CYCLE
  curved=cell_is_curved(mesh,i)
  CALL lag_rep%ncdofs(i,cell_dofs)
  res_loc = 0.d0
  by_weights_loc = by_weights(cell_dofs)
  !---Resistivity for this region (fallback to a large finite value in vacuum)
  eta1 = self%parent_sim%mug%eta(mesh%reg(i),1)
  IF(eta1 < 0.d0) eta1 = eta_fallback
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
  DO jr=1,lag_rep%nce
    by_res(cell_dofs(jr)) = by_res(cell_dofs(jr)) + res_loc(jr)
  END DO
END DO
DEALLOCATE(basis_vals,basis_grads,by_weights_loc,cell_dofs,res_loc)
END BLOCK
!---Plasma flux-function constraint (nlfun only, dt>0): overwrite the By residual
! at plasma DOFs with (By - By(F0_node)) so that By is driven to a constant equal
! to its value at the F0 node, then set the F0-node residual to the physical F0
! residual (plasma toroidal-flux + limiter-contour voltage integrals).
IF(dt > 0.d0 .AND. self%parent_sim%F0_node > 0)THEN
BLOCK
TYPE(gs_equil), POINTER :: eq
TYPE(gs_factory), POINTER :: dev
TYPE(oft_quad_type) :: quad_1d
REAL(r8), ALLOCATABLE :: by_res_plasma(:), psi_weights_loc(:), basis_vals_2(:), &
                         basis_vals(:), basis_grads(:,:), by_weights_loc(:), ff(:)
REAL(r8), POINTER, DIMENSION(:) :: psi_weights, pvac
INTEGER(i4), ALLOCATABLE :: cell_dofs_2(:), cell_b_dofs(:), elist(:,:)
INTEGER(i4) :: i, m, jr, k, je, cell, ed, nlim, F0
LOGICAL :: curved
REAL(r8) :: F0_res, psi, coords(3), jac_mat(3,4), jac_det, x1, y1, x2, y2, signed_area
REAL(r8) :: pts(2,2), dl(2), dn(3), dby(3), eta_p_loc
eq => self%parent_sim%tkmr%gs_equil
dev => self%parent_sim%tkmr%gs_device
F0 = self%parent_sim%F0_node
!---Constrain By to a constant (= By at F0) over the plasma DOFs
ALLOCATE(by_res_plasma(lag_rep%ne))
by_res_plasma = by_weights - by_weights(F0)
CALL fem_dirichlet_vec(lag_rep, by_res_plasma, by_res, self%parent_sim%plasma_flag)
DEALLOCATE(by_res_plasma)
!---Get the total poloidal flux (solved psi + vacuum psi_vac) for the bounds test
NULLIFY(psi_weights, pvac)
CALL a%get_local(psi_weights, 6)
CALL self%parent_sim%mug%psi_vac%get_local(pvac)
F0_res = 0.d0
!---Volume contribution: toroidal flux over the plasma region (region 1)
ALLOCATE(basis_vals_2(lag_rep%nce), psi_weights_loc(lag_rep%nce), cell_dofs_2(lag_rep%nce))
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
      F0_res = F0_res + SQRT(self%parent_sim%tkmr%f_scale*eq%I%f(psi) + by_weights(F0)**2) &
               *jac_det*quad%wts(m)/(coords(1)+gs_epsilon)
    ELSE
      ! in-region but outside the plasma: vacuum F = F0
      F0_res = F0_res + by_weights(F0)*jac_det*quad%wts(m)/(coords(1)+gs_epsilon)
    END IF
  END DO
END DO
DEALLOCATE(basis_vals_2, psi_weights_loc, cell_dofs_2)
!---Boundary contribution: voltage integral around the limiter contour
nlim = dev%nlim_con
signed_area = 0.d0
ALLOCATE(elist(2,nlim), cell_b_dofs(lag_rep%nce), basis_vals(lag_rep%nce), &
         basis_grads(3,lag_rep%nce), by_weights_loc(lag_rep%nce), ff(SIZE(quad%pts,1)))
!---For each contour segment, find the non-plasma cell and its local edge index
DO i=1,nlim
  IF(i < nlim)THEN
    je = ABS(mesh_local_findedge(mesh,[dev%lim_con(i),dev%lim_con(i+1)]))
    x1 = mesh%r(1,dev%lim_con(i));   y1 = mesh%r(2,dev%lim_con(i))
    x2 = mesh%r(1,dev%lim_con(i+1)); y2 = mesh%r(2,dev%lim_con(i+1))
    signed_area = signed_area + (x1*y2 - x2*y1) ! shoelace (for orientation)
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
!---Integrate eta_p*dt*dBy/dn/R along the contour
CALL set_quad_1d(quad_1d, lag_rep%order+2)
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
  dn = [-dl(2), dl(1), 0.d0] ! outward-ish normal*|edge|
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
DEALLOCATE(elist, cell_b_dofs, basis_vals, basis_grads, by_weights_loc, ff)
!---Overwrite the F0-node residual with the physical value
by_res(F0) = F0_res
CALL a%restore_local(psi_weights, 6)
END BLOCK
END IF
CALL b%restore_local(by_res, 7)
CALL a%restore_local(by_weights, 7)
end subroutine add_f_diff

subroutine build_blankettd_jacobian(self, mat, a, update_vac)
class(oft_blanket_td_sim), intent(inout) :: self
class(oft_matrix), pointer, intent(inout) :: mat
class(oft_matrix), pointer:: vac_op, mug_native
class(oft_vector), intent(inout) :: a ! Solution for computing Jacobian
LOGICAL, INTENT(in), optional :: update_vac
CLASS(oft_native_matrix), POINTER :: V => NULL()
CLASS(oft_native_matrix), POINTER :: M => NULL()
INTEGER(4) :: i, n, j, colcount, row_block, col_block, jp, jn
INTEGER(4), ALLOCATABLE :: cols(:)
REAL(r8), ALLOCATABLE :: vals(:)
REAL(r8), pointer :: tmp_vec(:)
INTEGER(4), ALLOCATABLE :: bc_rows(:)
INTEGER(4) :: nbc_rows, k
REAL(r8) :: diag_val(1,1)
class(oft_vector), pointer :: tmp

CALL mat%zero

! CALL fem_dirichlet_diag(lag_rep,mat,self%mug%n_bc,1)
! CALL fem_dirichlet_diag(lag_rep,mat,self%mug%velx_bc,2)
! CALL fem_dirichlet_diag(lag_rep,mat,self%mug%vely_bc,3)
! CALL fem_dirichlet_diag(lag_rep,mat,self%mug%velz_bc,4)
! CALL fem_dirichlet_diag(self%mug%fe_rep%fields(5)%fe,mat,self%mug%T_bc,5)
! CALL fem_dirichlet_diag(lag_rep,mat,self%mug%psi_bc,6)
! CALL fem_dirichlet_diag(lag_rep,mat,self%mug%by_bc,7)

! ! Set block (8,8) as identity matrix for coil currents
! IF (self%tkmr%gs_device%ncoils > 0) THEN
!     diag_val = 1.d0
!     DO i = 1, self%tkmr%gs_device%ncoils
!         CALL mat%add_values([i], [i], diag_val, 1, 1, 8, 8)
!     END DO
! END IF

select type(vac_op => self%tkmr%vac_op)
type is (oft_native_matrix)
    V => vac_op
class default
    call oft_abort("vac_op must be an oft_native_matrix", &
                   "build_blankettd_jacobian", __FILE__)
end select

select type(mug_jac => self%mug%jacobian)
type is (oft_native_matrix)
    M => mug_jac
class default
    call oft_abort("mug jacobian must be oft_native_matrix", &
                   "build_blankettd_jacobian", __FILE__)
end select


write(*,*) "Building individual jacobians"
!Populate MUG and TokaMaker matrices
! build_approx_jacobian scales all implicit terms by mug%jac_dt; blanket never
! calls mug%run_simulation (which normally sets it), so set it here to match the
! timestep used in the residual, else jac_dt stays at its -1.0 default and the
! mug Jacobian is ~1/dt too large and sign-flipped.
self%mug%jac_dt = self%nlfun%dt
CALL build_approx_jacobian(self%mug, a)
IF (PRESENT(update_vac)) THEN
    IF (update_vac) CALL build_vac_op(self%tkmr,self%tkmr%vac_op)
ELSE
    CALL build_vac_op(self%tkmr,self%tkmr%vac_op)
END IF

! Count and collect BC row indices
! nbc_rows = 0
! DO i = 1, SIZE(self%mug%psi_bc)
!     IF (self%mug%psi_bc(i)) nbc_rows = nbc_rows + 1
! END DO

! IF (nbc_rows > 0) THEN
!     ALLOCATE(bc_rows(nbc_rows))
!     bc_rows = 0
!     k = 0
!     DO i = 1, SIZE(self%mug%psi_bc)
!         IF (self%mug%psi_bc(i)) THEN
!             k = k + 1
!             bc_rows(k) = i
!         END IF
!     END DO
!     ! Zero rows in block (6,6) only (psi equation rows)
!     CALL M%zero_rows(nbc_rows, bc_rows, 6)
!     DEALLOCATE(bc_rows)
! END IF


write(*,*) "Combining jacobians"

DO row_block = 1, M%ni
    DO col_block = 1, M%nj
        DO i = 1, M%i_map(row_block)%n
            jp=M%map(row_block,col_block)%ext(1,i)
            jn=M%map(row_block,col_block)%ext(2,i)
            colcount = jn-jp+1
            ALLOCATE(cols(colcount), vals(colcount))
            cols = M%lc(jp:jn)
            vals = M%M(jp:jn)
            ! IF (MAXVAL(ABS(vals)) /=1.d0 .AND. row_block==col_block) THEN
            !     write(*,*) "WEIRD ", row_block, col_block, " row ", i
            !     write(*,*) MAXVAL(ABS(vals))
            ! END IF
            ! IF (row_block ==6 .AND. col_block == 6) THEN
            !     IF (MAXVAL(vals) > 0.d0) write(*,*) i
            ! END IF
            cols = cols - M%j_map(col_block)%offset
            CALL mat%add_values([i], cols, RESHAPE(vals, [1,colcount]), &
                1, colcount, row_block, col_block)
            DEALLOCATE(cols, vals)
        END DO
    END DO 
END DO

!ADD TOKAMAKER JACOBIAN to 6,6 block
! Add only block (1,1) from vac_op (no coil blocks for now)

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
    ! IF(ANY(vals /= 0.d0)) THEN
    !     write(*,*) SIZE(vals)
    ! END IF
    cols = cols - V%j_map(col_block)%offset
    CALL mat%add_values([i], cols, RESHAPE(vals, [1,colcount]), &
        1, colcount, 6, 6)
    DEALLOCATE(cols, vals)
END DO

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

CALL self%aug_vec%new(tmp)

CALL mat%assemble(tmp)
! NULLIFY(tmp_vec)

DEALLOCATE(tmp)
end subroutine build_blankettd_jacobian


subroutine delete_blanket_td(self)
class(oft_blanket_td_sim), intent(inout) :: self !< NL operator object
INTEGER(4) :: i
DEBUG_STACK_PUSH

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

subroutine delete_mfop(self)
class(oft_blanket_td_mfop), intent(inout) :: self !< NL operator object
DEBUG_STACK_PUSH
!
self%dt=-1.d0

!
IF(ASSOCIATED(self%jac_op))THEN
    CALL self%jac_op%delete()
    DEALLOCATE(self%jac_op)
END IF
DEBUG_STACK_POP
end subroutine

SUBROUTINE blanket_mfnk_update(a)
CLASS(oft_vector), TARGET, INTENT(inout) :: a
! CALL active_tMaker_td%mfop%update_lims(a)
CALL current_sim%mfmat%update(a)
! CALL build_jop(active_tMaker_td%mfop,adv_op,a)
!CALL active_tMaker_td%adv_solver%update(.TRUE.)
END SUBROUTINE blanket_mfnk_update

!IMPLEMENTTTTT
! !---------------------------------------------------------------------------
! !> Update matrix-free Jacobian on all levels with new solution
! !---------------------------------------------------------------------------
! subroutine mfnk_update(uin)
! class(oft_vector), target, intent(inout) :: uin !< Current field
! IF(oft_debug_print(1))write(*,*)'Updating 2D MUG MF-Jacobian'
! CALL current_sim%mf_mat%update(uin)
! END SUBROUTINE mfnk_update
! !---------------------------------------------------------------------------
! !> Update Jacobian matrices on all levels with new solution
! !---------------------------------------------------------------------------
! subroutine update_jacobian(uin)
! class(oft_vector), target, intent(inout) :: uin !< Current solution
! IF(oft_debug_print(1))write(*,*)'Updating 2D MUG approximate Jacobian'
! CALL build_approx_jacobian(current_sim,uin)
! END SUBROUTINE update_jacobian


END MODULE oft_blanket_td