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
USE oft_gs, ONLY: gs_epsilon, build_dels, gs_equil, gs_update_bounds, gs_test_bounds, set_bcmat
USE oft_gs_td, ONLY: oft_tmaker_td_mfop, tMaker_td_mfnk_update, build_vac_op, apply_rhs
USE oft_mesh_local_util, ONLY: mesh_local_findedge
USE xmhd_2d, ONLY: oft_xmhd_2d_sim, build_approx_jacobian
USE oft_stitching, ONLY: oft_seam, seam_list
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

contains
    !> Apply the matrix
    procedure :: setup => setup_blanket_td
    !> Needs Docs
    procedure :: delete => delete_blanket_td
    !> Needs Docs
    procedure :: step => step_blanket_td
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
TYPE(oft_graph), TARGET :: dense_graph
type(oft_1d_int), pointer, dimension(:) :: bc_nodes
integer(i4), allocatable :: dense_flag(:)
LOGICAL, ALLOCATABLE :: p_dir_set(:)
INTEGER(i4), POINTER, DIMENSION(:) :: cell_dofs, cell_dofs_p
type(seam_list), pointer, dimension(:) :: stitch_tmp
type(map_list), pointer, dimension(:) :: map_tmp
CLASS(oft_vector), pointer :: tmp_vec

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

CALL mhd_sim%setup(mg_mesh, lag_rep%order, lag_rep)

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

mhd_sim%psi_bc = .TRUE. 
mhd_sim%by_bc = .TRUE.  !-> deal with by later
mhd_sim%velx_bc = .TRUE. 
mhd_sim%vely_bc = .TRUE. 
mhd_sim%velz_bc = .TRUE. 

mhd_sim%n_bc = .TRUE.
mhd_sim%T_bc = .TRUE.

! mhd_sim%psi_bc = .TRUE. 
! mhd_sim%by_bc = .TRUE.  !-> deal with by later
! mhd_sim%velx_bc = .FALSE. 
! mhd_sim%vely_bc = .FALSE. 
! mhd_sim%velz_bc = .FALSE. 

! mhd_sim%n_bc = .TRUE.
! mhd_sim%T_bc = .TRUE.


! ALLOCATE(cell_dofs(lag_rep%nce))
! ALLOCATE(cell_dofs_p(mhd_sim%fe_rep%fields(5)%fe%ne))
! DO i = 1, mesh%nc
!     call mhd_sim%fe_rep%fields(5)%fe%ncdofs(i,cell_dofs_p) ! Get global index of local DOFs
!     call lag_rep%ncdofs(i,cell_dofs) ! Get global index of local DOFs
!     DO j=1, SIZE(cell_dofs)
!         IF (mhd_flag(mesh%reg(i))) THEN
!             mhd_sim%psi_bc(cell_dofs(j)) = .FALSE. !add contributions everywhere in MHD regions, including boundaries
!             mhd_sim%T_bc(cell_dofs_p(j)) = .FALSE. !add contributions everywhere in MHD regions, including boundaries
!             IF (.NOT. mhd_sim%incomp) mhd_sim%n_bc(cell_dofs(j)) = .FALSE. !add contributions everywhere in MHD regions, including boundaries, if incompressible
!         ELSE
!             mhd_sim%velx_bc(cell_dofs(j)) = .TRUE. !pin velocity to zero everywhere in non-MHD regions, including boundaries
!             mhd_sim%vely_bc(cell_dofs(j)) = .TRUE. 
!             mhd_sim%velz_bc(cell_dofs(j)) = .TRUE.
!         END IF 
!     END DO
! END DO

! IF (mhd_sim%incomp) THEN
!     ALLOCATE(p_dir_set(mesh%nreg))
!     p_dir_set = .FALSE.
!     DO i=1, mesh%nc
!         IF (mhd_flag(mesh%reg(i)) .AND. .NOT. p_dir_set(mesh%reg(i))) THEN
!             call mhd_sim%fe_rep%fields(5)%fe%ncdofs(i,cell_dofs_p) ! Get global index of local DOFs
!             mhd_sim%T_bc(cell_dofs_p(1)) = .TRUE.
!             p_dir_set(mesh%reg(i)) = .TRUE.
!             write(*,*) "Pinning node ", cell_dofs_p(1), " in region ", mesh%reg(i)
!         END IF
!     END DO
!     DEALLOCATE(p_dir_set)
! END IF
! DEALLOCATE(cell_dofs, cell_dofs_p)

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
CALL self%u%set(1000.d0, 5)
CALL self%u%set(self%tkmr%gs_equil%I%f_offset, 7)

NULLIFY(tmp_arr, vals_out)
CALL self%tkmr%gs_device%fe_rep%vec_create(tmp_vec)
CALL tmp_vec%add(0.d0,1.d0,self%tkmr%gs_equil%psi)
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
DEALLOCATE(tmp_vec)
DEALLOCATE(fe_graphs)


! Preconditioner should use approximate jacobian
ALLOCATE(self%pre) !CHECKBACK
self%pre%A=>self%nlfun%jac_op 

!------------------------------------------------------------------------------
! Setup matrix free solver
!------------------------------------------------------------------------------
ALLOCATE(self%mfmat) 
self%mfmat%f=>self%nlfun
CALL self%rhs%new(self%mfmat%u0)
CALL self%rhs%new(self%mfmat%f0)
CALL self%rhs%new(self%mfmat%tmp)
CALL self%rhs%new(self%mfmat%utyp)


ALLOCATE(self%mf_solver)
self%mfmat%b0=1.d-5
self%mf_solver%A=>self%mfmat
self%mf_solver%its=1000
self%mf_solver%nrits=20
self%mf_solver%atol=self%lin_tol
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
self%nksolver%atol=self%nl_tol
self%nksolver%rtol=1.d-20 ! Disable relative tolerance
self%nksolver%backtrack=.FALSE.
!self%nksolver%J_update=>gs_mfnk_update
self%nksolver%up_freq=1

!TEST MUG JACOBIAN
CALL mhd_sim%u%set(5.d0)
CALL mhd_sim%u0%set(5.d0)
CALL build_approx_jacobian(mhd_sim, mhd_sim%u)
self%pre%A => mhd_sim%jacobian
! self%pre%A => self%tkmr%vac_op
CALL self%pre%apply(self%u, self%tmp)
CALL mhd_sim%jacobian%apply(mhd_sim%u, mhd_sim%u0)
CALL mhd_sim%u0%get_local(tmp_arr)
write(*,*) MAXVAL(tmp_arr)
write(*,*) MINVAL(tmp_arr)

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
self%pre%A => self%tkmr%vac_op
self%mf_solver%pre => self%pre

CALL self%pre%update(.TRUE.)
CALL self%mf_solver%pre%update(.TRUE.)
! !diagnose preconditioner:
! select type(pre=> self%pre%A)
! type is (oft_native_matrix)
!     P => pre
! class default
!     call oft_abort("mug jacobian must be oft_native_matrix", &
!                    "build_blankettd_jacobian", __FILE__)
! end select
! write(*,*) 'nnz:', P%nnz
! write(*,*) 'nrg:', P%nrg
! write(*,*) 'ncg:', P%ncg
! DO i = 1, P%ni
!     DO j = 1, P%nj
!         write(*,*) "Block (", i, ",", j, ") nnz:", P%map(i,j)%nnz
!     END DO
! END DO

!Build RHS and apply boundary conditions
CALL self%tmp%add(0.d0,1.d0,self%u)
CALL apply_rhs_blanket(self%nlfun,self%u,self%rhs)

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

end subroutine step_blanket_td

subroutine apply_rhs_blanket(self, a, b)
class(oft_blanket_td_mfop), intent(inout) :: self
class(oft_vector), target, intent(inout) :: a !< Source field
class(oft_vector), intent(inout) :: b !< Result of metric function
class(oft_vector), pointer :: tmp_in, tmp_out !< Temporary vectors
REAL(r8), POINTER, DIMENSION(:) :: tmp_arr1, tmp_arr2

self%parent_sim%mug%nlfun%dt = 0.d0
CALL b%set(0.d0)
CALL self%parent_sim%mug%nlfun%apply_real(a,b)
NULLIFY(tmp_arr1)
NULLIFY(tmp_arr2)
CALL b%get_local(tmp_arr1, 6) 
tmp_arr1 = tmp_arr1/self%dt !Divide by dt so form of psi equation matches tokamaker implementation

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
CALL self%parent_sim%tkmr%gs_device%zerob_bc%apply(tmp_out)

! Extract component 1 from tmp_out
CALL tmp_out%get_local(tmp_arr2, 1)
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
NULLIFY(tmp_arr1)
NULLIFY(tmp_arr2)
CALL b%get_local(tmp_arr1, 6)
tmp_arr1 = tmp_arr1/self%dt

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
CALL self%parent_sim%tkmr%gs_device%zerob_bc%apply(tmp_out)

! Extract component 1 from tmp_out
CALL tmp_out%get_local(tmp_arr2, 1)
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

CALL fem_dirichlet_diag(lag_rep,mat,self%mug%n_bc,1)
CALL fem_dirichlet_diag(lag_rep,mat,self%mug%velx_bc,2)
CALL fem_dirichlet_diag(lag_rep,mat,self%mug%vely_bc,3)
CALL fem_dirichlet_diag(lag_rep,mat,self%mug%velz_bc,4)
CALL fem_dirichlet_diag(self%mug%fe_rep%fields(5)%fe,mat,self%mug%T_bc,5)
CALL fem_dirichlet_diag(lag_rep,mat,self%mug%psi_bc,6)
CALL fem_dirichlet_diag(lag_rep,mat,self%mug%by_bc,7)

! Set block (8,8) as identity matrix for coil currents
IF (self%tkmr%gs_device%ncoils > 0) THEN
    diag_val = 1.d0
    DO i = 1, self%tkmr%gs_device%ncoils
        CALL mat%add_values([i], [i], diag_val, 1, 1, 8, 8)
    END DO
END IF

! select type(vac_op => self%tkmr%vac_op)
! type is (oft_native_matrix)
!     V => vac_op
! class default
!     call oft_abort("vac_op must be an oft_native_matrix", &
!                    "build_blankettd_jacobian", __FILE__)
! end select

! select type(mug_jac => self%mug%jacobian)
! type is (oft_native_matrix)
!     M => mug_jac
! class default
!     call oft_abort("mug jacobian must be oft_native_matrix", &
!                    "build_blankettd_jacobian", __FILE__)
! end select


! write(*,*) "Building individual jacobians"
! !Populate MUG and TokaMaker matrices
! CALL build_approx_jacobian(self%mug, a)
! IF (PRESENT(update_vac)) THEN
!     IF (update_vac) CALL build_vac_op(self%tkmr,self%tkmr%vac_op)
! ELSE
!     CALL build_vac_op(self%tkmr,self%tkmr%vac_op)
! END IF

! ! Count and collect BC row indices
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


! write(*,*) "Combining jacobians"

! DO row_block = 1, M%ni
!     DO col_block = 1, M%nj
!         DO i = 1, M%i_map(row_block)%n
!             jp=M%map(row_block,col_block)%ext(1,i)
!             jn=M%map(row_block,col_block)%ext(2,i)
!             colcount = jn-jp+1
!             ALLOCATE(cols(colcount), vals(colcount))
!             cols = M%lc(jp:jn)
!             vals = M%M(jp:jn)
!             ! IF (row_block ==6 .AND. col_block == 6) THEN
!             !     IF (MAXVAL(vals) > 0.d0) write(*,*) i
!             ! END IF
!             cols = cols - M%j_map(col_block)%offset
!             CALL mat%add_values([i], cols, RESHAPE(vals, [1,colcount]), &
!                 1, colcount, row_block, col_block)
!             DEALLOCATE(cols, vals)
!         END DO
!     END DO 
! END DO

! !ADD TOKAMAKER JACOBIAN to 6,6 block
! ! Add only block (1,1) from vac_op (no coil blocks for now)

! row_block = 1
! col_block = 1
! DO i = 1, V%i_map(row_block)%n
!     jp=V%map(row_block,col_block)%ext(1,i)
!     jn=V%map(row_block,col_block)%ext(2,i)
!     colcount = jn-jp+1
!     ALLOCATE(cols(colcount), vals(colcount))
!     cols = V%lc(jp:jn)
!     vals = V%M(jp:jn)
!     cols = cols - V%j_map(col_block)%offset
!     CALL mat%add_values([i], cols, RESHAPE(vals, [1,colcount]), &
!         1, colcount, 6, 6)
!     DEALLOCATE(cols, vals)
! END DO

! row_block = 1
! col_block = 2
! DO i = 1, V%i_map(row_block)%n
!     jp=V%map(row_block,col_block)%ext(1,i)
!     jn=V%map(row_block,col_block)%ext(2,i)
!     colcount = jn-jp+1
!     ALLOCATE(cols(colcount), vals(colcount))
!     cols = V%lc(jp:jn)
!     vals = V%M(jp:jn)
!     cols = cols - V%j_map(col_block)%offset
!     CALL mat%add_values([i], cols, RESHAPE(vals, [1,colcount]), &
!         1, colcount, 6, 8)
!     DEALLOCATE(cols, vals)
! END DO

! row_block = 2
! col_block = 1
! DO i = 1, V%i_map(row_block)%n
!     jp=V%map(row_block,col_block)%ext(1,i)
!     jn=V%map(row_block,col_block)%ext(2,i)
!     colcount = jn-jp+1
!     ALLOCATE(cols(colcount), vals(colcount))
!     cols = V%lc(jp:jn)
!     vals = V%M(jp:jn)
!     cols = cols - V%j_map(col_block)%offset
!     CALL mat%add_values([i], cols, RESHAPE(vals, [1,colcount]), &
!         1, colcount, 8, 6)
!     DEALLOCATE(cols, vals)
! END DO

! row_block = 2
! col_block = 2
! DO i = 1, V%i_map(row_block)%n
!     jp=V%map(row_block,col_block)%ext(1,i)
!     jn=V%map(row_block,col_block)%ext(2,i)
!     colcount = jn-jp+1
!     ALLOCATE(cols(colcount), vals(colcount))
!     cols = V%lc(jp:jn)
!     vals = V%M(jp:jn)
!     cols = cols - V%j_map(col_block)%offset
!     CALL mat%add_values([i], cols, RESHAPE(vals, [1,colcount]), &
!         1, colcount, 8, 8)
!     DEALLOCATE(cols, vals)
! END DO

CALL self%aug_vec%new(tmp)

CALL mat%assemble(tmp)
! NULLIFY(tmp_vec)
! CALL tmp%get_local(tmp_vec)
! DO i = 1, SIZE(tmp_vec)
!     IF(tmp_vec(i) /= 1.d0) THEN
!         write(*,*) i
!     END IF
! END DO

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