PROGRAM hartmann
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
USE fem_utils, ONLY: bfem_map_flag
USE xmhd_2d
IMPLICIT NONE
INTEGER(i4) :: io_unit,ierr
REAL(r8), POINTER :: vec_vals(:)
TYPE(oft_xmhd_2d_sim) :: mhd_sim
TYPE(multigrid_mesh) :: mg_mesh
TYPE(oft_blag_zerob), TARGET :: blag_zerob ! setting boundary vals to zero
!---Mass matrix solver
TYPE(poss_scalar_bfield) :: field_init
CLASS(oft_solver), POINTER :: minv => NULL()
CLASS(oft_matrix), POINTER :: mop => NULL()
CLASS(oft_vector), POINTER :: u,v
!---Runtime options
INTEGER(i4) :: order = 2
INTEGER(i4) :: i
INTEGER(i4) :: nsteps = 100
INTEGER(i4) :: rst_freq = 10
REAL(r8) :: n0 = 1.d0
REAL(r8) :: velx0 = 1.d0
REAL(r8) :: vely0 = 1.d0
REAL(r8) :: velz0 = 1.d0
REAL(r8) :: t0 = 1.d0
REAL(r8) :: psi0 = 1.d0
REAL(r8) :: by0 = 1.d0
REAL(r8) :: chi=1.E-12 !< Needs docs
REAL(r8) :: eta=1.d0 !< Needs docs
REAL(r8) :: nu=1.E-12 !< Needs docs
REAL(r8) :: gamma=1.d0
REAL(r8) :: D_diff=1.E-12
REAL(r8) :: m_i=proton_mass
REAL(r8) :: mu = mu0
REAL(r8) :: den_scale=1.d24
REAL(r8) :: dt = 1.d-3
REAL(r8) :: v_delta = 0.d0
REAL(r8) :: B0
LOGICAL, ALLOCATABLE :: vert_flag(:),edge_flag(:)
LOGICAL :: pm=.TRUE.
LOGICAL :: use_mfnk=.FALSE.
NAMELIST/xmhd_options/order,chi,eta,nu,gamma, D_diff, m_i, mu,&
dt,nsteps,rst_freq,use_mfnk,pm, n0, psi0, velx0,vely0,velz0, t0, by0, den_scale
CALL oft_init
!---Read in options
OPEN(NEWUNIT=io_unit,FILE=oft_env%ifile)
READ(io_unit,xmhd_options,IOSTAT=ierr)
CLOSE(io_unit)
!---------------------------------------------------------------------------
! Setup grid
!---------------------------------------------------------------------------
CALL multigrid_construct_surf(mg_mesh)
CALL mhd_sim%setup(mg_mesh,order)
blag_zerob%ML_lag_rep=>ML_oft_blagrange

!---------------------------------------------------------------------------
! Set intial conditions from analytic functions
!---------------------------------------------------------------------------
!---Generate mass matrix
NULLIFY(u,v,mop,vec_vals) ! Ensure the matrix is unallocated (pointer is NULL)
CALL oft_blag_getmop(ML_oft_blagrange%current_level,mop) ! Construct mass matrix with "none" BC
!---Setup linear solver
CALL create_cg_solver(minv)
minv%A=>mop ! Set matrix to be solved
minv%its=-2 ! Set convergence type (in this case "full" CG convergence)
CALL create_diag_pre(minv%pre) ! Setup Preconditioner
!---Create fields for solver
CALL ML_oft_blagrange%vec_create(u)
CALL ML_oft_blagrange%vec_create(v)

!---Project n initial condition onto scalar Lagrange basis
CALL u%set(0.d0)
CALL u%get_local(vec_vals)
vec_vals = (1.d0 - vec_vals*1.d-4)*n0
! CALL mesh%save_vertex_scalar(vec_vals,mhd_sim%xdmf_plot,'n0')
mhd_sim%den_scale = den_scale
vec_vals = vec_vals / den_scale
CALL mhd_sim%u%restore_local(vec_vals,1)

v_delta = 1.d0
!---Project v_x initial condition onto scalar Lagrange basis
field_init%func=>vx_init
field_init%mesh=>mesh
CALL oft_blag_project(ML_oft_blagrange%current_level,field_init,v)
CALL u%set(0.d0)
CALL minv%apply(u,v)
CALL u%scale(v_delta)
! CALL blag_zerob%apply(u)
CALL u%get_local(vec_vals)
vec_vals = vec_vals - MINVAL(vec_vals) ! Shift to ensure positivity
! vec_vals = vec_vals - MINVAL(vec_vals) ! Shift to ensure positivity
! CALL mesh%save_vertex_scalar(vec_vals,mhd_sim%xdmf_plot,'vx0')
CALL mhd_sim%u%restore_local(vec_vals,2)

!---Project v_y initial condition onto scalar Lagrange basis

CALL u%set(0.d0)
CALL u%get_local(vec_vals)
! CALL mesh%save_vertex_scalar(vec_vals,mhd_sim%xdmf_plot,'vy0')
CALL mhd_sim%u%restore_local(vec_vals,3)

!---Project v_z initial condition onto scalar Lagrange basis
CALL u%set(0.d0)
CALL u%get_local(vec_vals)
! CALL mesh%save_vertex_scalar(vec_vals,mhd_sim%xdmf_plot,'vz0')
CALL mhd_sim%u%restore_local(vec_vals,4)

!---Project T initial condition onto scalar Lagrange basis
CALL u%set(1.d0)
CALL u%scale(t0)
CALL u%get_local(vec_vals)
! CALL mesh%save_vertex_scalar(vec_vals,mhd_sim%xdmf_plot,'T0')
CALL mhd_sim%u%restore_local(vec_vals,5)

!---Project psi initial condition onto scalar Lagrange basis
field_init%func=>psi_alf
CALL oft_blag_project(ML_oft_blagrange%current_level,field_init,v)
CALL u%set(0.d0)
CALL minv%apply(u,v)
CALL u%scale(psi0)
CALL u%get_local(vec_vals)
! CALL mesh%save_vertex_scalar(vec_vals,mhd_sim%xdmf_plot,'psi0')
CALL mhd_sim%u%restore_local(vec_vals,6)

!---Project by initial condition onto scalar Lagrange basis
CALL u%set(1.d0)
CALL u%scale(by0)
CALL u%get_local(vec_vals)
! CALL mesh%save_vertex_scalar(vec_vals,mhd_sim%xdmf_plot,'by0')
CALL mhd_sim%u%restore_local(vec_vals,7)

!---Cleanup objects used for projection
CALL u%delete ! Destroy LHS vector
CALL v%delete ! Destroy RHS vector
CALL mop%delete ! Destroy mass matrix
DEALLOCATE(u,v,mop) ! Deallocate objects
CALL minv%pre%delete ! Destroy preconditioner
DEALLOCATE(minv%pre)
CALL minv%delete ! Destroy solver
DEALLOCATE(minv)

!---------------------------------------------------------------------------
! Run simulation
!---------------------------------------------------------------------------
mhd_sim%chi=chi
mhd_sim%eta=eta
mhd_sim%gamma=gamma
mhd_sim%D_diff=D_diff
mhd_sim%m_i = m_i*proton_mass
mhd_sim%nu=nu/(mhd_sim%m_i*n0)
! mhd_sim%mu = mu
mhd_sim%dt=dt
mhd_sim%nsteps=nsteps
mhd_sim%rst_freq=rst_freq
mhd_sim%mfnk=use_mfnk
oft_env%pm=.FALSE.

!---Set boundary conditions for pipe flow in square mesh
ALLOCATE(vert_flag(mesh%np),edge_flag(mesh%ne))
!---Temperature is fixed on outlet only
NULLIFY(mhd_sim%T_bc) ! Ensure the vector is unallocated (pointer is NULL)
ALLOCATE(mhd_sim%T_bc(oft_blagrange%ne))
vert_flag=.FALSE.; edge_flag=.FALSE.
DO i=1,mesh%nbe
  IF(mesh%bes(i)==2)THEN
    edge_flag(mesh%lbe(i))=.TRUE.
    vert_flag(mesh%le(1,mesh%lbe(i)))=.TRUE.
    vert_flag(mesh%le(2,mesh%lbe(i)))=.TRUE.
  END IF
END DO
CALL bfem_map_flag(oft_blagrange,vert_flag,edge_flag,mhd_sim%T_bc)
!---Density is fixed on inlet/outlet
NULLIFY(mhd_sim%n_bc) ! Ensure the vector is unallocated (pointer is NULL)
ALLOCATE(mhd_sim%n_bc(oft_blagrange%ne))
vert_flag=.FALSE.; edge_flag=.FALSE.
DO i=1,mesh%nbe
  IF(mesh%bes(i)<3)THEN
    edge_flag(mesh%lbe(i))=.TRUE.
    vert_flag(mesh%le(1,mesh%lbe(i)))=.TRUE.
    vert_flag(mesh%le(2,mesh%lbe(i)))=.TRUE.
  END IF
END DO
CALL bfem_map_flag(oft_blagrange,vert_flag,edge_flag,mhd_sim%n_bc)
! !---Velocity outlet needs work
! ALLOCATE(self%velx_bc(oft_blagrange%ne))
! vert_flag=.FALSE.; edge_flag=.FALSE.
! DO i=1,mesh%nbe
!   IF(mesh%bes(i)/=2)THEN
!     edge_flag(mesh%lbe(i))=.TRUE.
!     vert_flag(mesh%le(1,mesh%lbe(i)))=.TRUE.
!     vert_flag(mesh%le(2,mesh%lbe(i)))=.TRUE.
!   END IF
! END DO
! CALL bfem_map_flag(oft_blagrange,vert_flag,edge_flag,self%velx_bc)
! ! self%vely_bc=>self%velx_bc
! ! self%velz_bc=>self%velx_bc
!---Clean up flags
DEALLOCATE(vert_flag,edge_flag)


CALL mhd_sim%run_simulation()
CALL xmhd_2d_plot(mhd_sim)
!---Finalize enviroment
CALL oft_finalize
CONTAINS
!


SUBROUTINE psi_alf(pt, val)
REAL(r8), INTENT(in) :: pt(3)
REAL(r8), INTENT(out) :: val
val=1.d0*pt(1)
END SUBROUTINE psi_alf

SUBROUTINE vx_init(pt, val)
REAL(r8), INTENT(in) :: pt(3)
REAL(r8), INTENT(out) :: val
val=1.d0*EXP(-pt(2)**16)-1.d0*EXP(-1.d0)
END SUBROUTINE vx_init


END PROGRAM hartmann