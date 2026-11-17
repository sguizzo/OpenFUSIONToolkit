!---------------------------------------------------------------------------
! Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
!---------------------------------------------------------------------------
!> @file empump_td.F90
!
!> Time-dependent wrapper around xmhd_2d for electromagnetic induction pumps.
!!
!! Region types:
!!  - Vacuum/air:  a region with a large resistivity. 
!!  - Solid conductors: a conducting region with its real resistivity and the velocity
!!                 pinned to zero.
!!  - Coils:       a prescribed current density per region (`j_coil`), added to
!!                 the psi right-hand side. Coil regions have a large
!!                 resistivity so they carry no induced current on top of it.
!!  - Iron:        nonlinear permeability mu(|B|) from gs_mu_of_B, added as a
!!                 correction dt*(1/mu - 1)*stiffness on top of xmhd_2d's
!!                 mu = mu0 operator. Large resistivity by default.
!!  - Outer boundary: psi Dirichlet (by default psi = 0, a perfect conductor),
!!                 so the boundary must sit far enough out that the field has
!!                 decayed there.
!!
!! Geometry is slab (y ignorable) or cylindrical (phi ignorable), selected at
!! setup. In slab psi = A_y; in cylindrical psi = R*A_phi.
!!
!! Residual for one step:
!!   F(u) = xmhd_2d(u; dt) + iron(u; dt)
!!   rhs  = xmhd_2d(u_old; 0) + dt*mu0*j_coil
!! and the Newton-Krylov solve finds F(u) = rhs.
!!
!! xmhd_2d keeps its FE spaces and mesh in module variables, so only one
!! empump_td object may exist per program.
!---------------------------------------------------------------------------
MODULE empump_td
USE oft_base
USE oft_quadrature, ONLY: oft_quad_type
USE oft_mesh_type, ONLY: oft_bmesh, cell_is_curved
USE multigrid, ONLY: multigrid_mesh
!
USE oft_la_base, ONLY: oft_vector, oft_matrix
USE oft_deriv_matrices, ONLY: oft_noop_matrix, oft_mf_matrix
USE oft_native_solvers, ONLY: oft_nksolver, oft_native_gmres_solver
USE oft_lu, ONLY: oft_lusolver
!
USE oft_lag_basis, ONLY: oft_scalar_bfem, oft_blag_eval, oft_blag_geval
USE mhd_utils, ONLY: mu0
USE oft_gs, ONLY: gs_epsilon, gs_mu_of_B
USE xmhd_2d, ONLY: oft_xmhd_2d_sim, build_approx_jacobian, xmhd_2d_plot, oft_blagrange
IMPLICIT NONE
#include "local.h"
#if !defined(XMHD_RST_LEN)
#define XMHD_RST_LEN 5
#endif
PRIVATE

!---------------------------------------------------------------------------
!> EM device simulation object
!---------------------------------------------------------------------------
TYPE, PUBLIC :: oft_empump_td
  TYPE(oft_xmhd_2d_sim), POINTER :: mug => NULL() !< Underlying xmhd_2d simulation
  CLASS(oft_vector), POINTER :: u => NULL() !< Current solution (aliases mug%u)
  CLASS(oft_vector), POINTER :: u_old => NULL() !< Solution at the start of the step
  CLASS(oft_vector), POINTER :: rhs => NULL() !< Right-hand side of the step
  TYPE(oft_empump_td_nlfun), POINTER :: nlfun => NULL() !< Time-advance operator
  TYPE(oft_mf_matrix), POINTER :: mfmat => NULL() !< Matrix-free Jacobian
  TYPE(oft_native_gmres_solver), POINTER :: mf_solver => NULL() !< Outer linear solver
  TYPE(oft_lusolver), POINTER :: pre => NULL() !< LU of the approximate Jacobian (mug%jacobian)
  TYPE(oft_nksolver) :: nksolver !< Newton-Krylov solver
  LOGICAL :: pm = .FALSE. !< Print solver histories
  LOGICAL :: save_rst = .FALSE. !< Write 'xmhd2d_NNNNN.rst' files (readable by xmhd_2d_plot)
  INTEGER(i4) :: rst_freq = 1 !< Restart output frequency [steps]
  INTEGER(i4) :: rst_count = 0 !< Completed steps, used as the restart index
  LOGICAL :: has_iron = .FALSE. !< Any region flagged in iron_reg
  LOGICAL :: iron_jac_lift = .FALSE. !< Include the (nu(|B|) - 1) stiffness lift in the approximate Jacobian (preconditioner only). Off (default) leaves the iron at mu0 in the preconditioner, as mugtok_td does with mu_relax = 0; on the em_pump test the lift tripled the GMRES iterations (2026-10-23)
  LOGICAL :: iron_jac_deriv = .FALSE. !< Include the d(mu)/d|B| term in the approximate Jacobian (preconditioner only). Off by default: on the em_pump mesh it made Newton fail and halve dt (2026-10-23 comparison), matching mugtok_td's default mu_jac_deriv = F
  LOGICAL, ALLOCATABLE, DIMENSION(:) :: fluid_reg !< True where the velocity and pressure are solved, by region
  LOGICAL, ALLOCATABLE, DIMENSION(:) :: iron_reg !< True where the nonlinear permeability applies, by region
  REAL(r8), ALLOCATABLE, DIMENSION(:) :: j_coil !< Prescribed out-of-plane current density [A/m^2], by region
  LOGICAL :: flow_checked = .FALSE. !< The untested-flow warning has been checked (done once, on the first step)
CONTAINS
  !> Setup simulation
  PROCEDURE :: setup => setup_empump_td
  !> Take one timestep
  PROCEDURE :: step => step_empump_td
  !> Write a restart file for the current state
  PROCEDURE :: rst_save => empump_rst_save
  !> Convert restart files to XDMF (wraps xmhd_2d_plot)
  PROCEDURE :: plot => empump_td_plot
  !> Delete simulation
  PROCEDURE :: delete => delete_empump_td
END TYPE oft_empump_td

!---------------------------------------------------------------------------
!> Time-advance operator: xmhd_2d's residual plus the iron correction
!---------------------------------------------------------------------------
TYPE, EXTENDS(oft_noop_matrix) :: oft_empump_td_nlfun
  REAL(r8) :: dt = -1.d0 !< Timestep [s]
  TYPE(oft_empump_td), POINTER :: parent_sim => NULL() !< Parent simulation
CONTAINS
  !> Apply the nonlinear function
  PROCEDURE :: apply_real => nlfun_apply
END TYPE oft_empump_td_nlfun

TYPE(oft_empump_td), POINTER :: current_sim => NULL() !< Active simulation, for the MF Jacobian update
CLASS(oft_bmesh), POINTER :: mesh => NULL() !< Mesh (xmhd_2d's)
CLASS(oft_scalar_bfem), POINTER :: lag_rep => NULL() !< Order-p Lagrange space (xmhd_2d's oft_blagrange)

CONTAINS

!---------------------------------------------------------------------------
!> Setup the simulation
!!
!! Default boundary conditions, all of which the driver may edit afterwards
!! through self%mug%*_bc:
!!  - velocity: pinned in every non-fluid cell, which makes fluid/solid
!!    interfaces no-slip, and on the outer boundary (v_Z is left free on the
!!    axis). v_y/v_phi is pinned everywhere unless toroidal_flow.
!!  - pressure: free on DOFs of fluid cells, pinned elsewhere, plus one gauge
!!    DOF per fluid region (set pin_gauge=.FALSE. if the driver imposes the
!!    pressure on a face instead)
!!  - psi: Dirichlet on the outer boundary, holding its initial value (a perfect
!!    conductor), or free there (dpsi/dn = 0, an infinitely permeable wall) with
!!    psi_outer_dirichlet=.FALSE.; pinned on the axis either way
!!  - B_y/B_phi: pinned everywhere unless evolve_by (nothing drives it in a
!!    pump)
!!  - axis regularity (cylindrical), via xmhd_2d's apply_axis_bc
!! With freeze_flow=.TRUE. every velocity and pressure DOF is pinned at zero
!! (an EM-only run: the fluid still carries eddy currents but does not move).
!---------------------------------------------------------------------------
subroutine setup_empump_td(self,mg_mesh,order,dt,lin_tol,nl_tol,fluid_reg,dens_reg,visc_reg, &
                           eta_reg,cyl,iron_reg,toroidal_flow,evolve_by,pin_gauge,freeze_flow, &
                           psi_outer_dirichlet)
CLASS(oft_empump_td), TARGET, INTENT(inout) :: self !< Simulation object
CLASS(multigrid_mesh), TARGET, INTENT(in) :: mg_mesh !< Mesh
INTEGER(i4), INTENT(in) :: order !< FE order (pressure uses order-1)
REAL(r8), INTENT(in) :: dt !< Timestep [s]
REAL(r8), INTENT(in) :: lin_tol !< Linear solver tolerance
REAL(r8), INTENT(in) :: nl_tol !< Nonlinear solver tolerance
LOGICAL, INTENT(in) :: fluid_reg(:) !< True for regions where the fluid moves
REAL(r8), INTENT(in) :: dens_reg(:) !< Mass density [kg/m^3], by region (used in fluid regions only)
REAL(r8), INTENT(in) :: visc_reg(:) !< Kinematic viscosity [m^2/s], by region (used in fluid regions only)
REAL(r8), INTENT(in) :: eta_reg(:,:) !< Magnetic diffusivity eta/mu0 [m^2/s], by region, [in plane, out of plane]. Must be positive everywhere; make it large in air, coils and laminated iron
LOGICAL, INTENT(in) :: cyl !< Cylindrical (R,Z) if true, slab (x,z) if false
LOGICAL, OPTIONAL, INTENT(in) :: iron_reg(:) !< True for regions with nonlinear permeability (default none)
LOGICAL, OPTIONAL, INTENT(in) :: toroidal_flow !< Solve the out-of-plane velocity (default false)
LOGICAL, OPTIONAL, INTENT(in) :: evolve_by !< Solve the out-of-plane field (default false)
LOGICAL, OPTIONAL, INTENT(in) :: pin_gauge !< Pin one pressure DOF per fluid region (default true)
LOGICAL, OPTIONAL, INTENT(in) :: freeze_flow !< Pin all velocity and pressure DOFs at zero, for an EM-only run (default false)
LOGICAL, OPTIONAL, INTENT(in) :: psi_outer_dirichlet !< psi fixed on the outer boundary (default true); false leaves it free
TYPE(oft_xmhd_2d_sim), POINTER :: mhd
INTEGER(i4) :: i,j,nreg
INTEGER(i4), ALLOCATABLE :: cell_dofs(:),cell_dofs_p(:)
LOGICAL :: tor_flow,do_by,do_gauge,do_freeze,psi_dir
LOGICAL, ALLOCATABLE :: gauge_set(:)

nreg = mg_mesh%smesh%nreg
IF(SIZE(fluid_reg)/=nreg .OR. SIZE(dens_reg)/=nreg .OR. SIZE(visc_reg)/=nreg .OR. SIZE(eta_reg,1)/=nreg) &
  CALL oft_abort('Region arrays must have one entry per mesh region','setup_empump_td',__FILE__)
IF(.NOT.ANY(fluid_reg)) &
  CALL oft_abort('No fluid region flagged','setup_empump_td',__FILE__)
DO i=1,nreg
  IF(ANY(eta_reg(i,:) <= 0.d0))THEN
    WRITE(*,'(A,I0,A,2ES11.3)')'Region ',i,' has non-positive resistivity: ',eta_reg(i,:)
    CALL oft_abort('Resistivity must be positive in every region','setup_empump_td',__FILE__)
  END IF
  IF(fluid_reg(i) .AND. dens_reg(i) <= 0.d0) &
    CALL oft_abort('Fluid regions need a positive density','setup_empump_td',__FILE__)
END DO
tor_flow = .FALSE.
IF(PRESENT(toroidal_flow)) tor_flow = toroidal_flow
do_by = .FALSE.
IF(PRESENT(evolve_by)) do_by = evolve_by
do_gauge = .TRUE.
IF(PRESENT(pin_gauge)) do_gauge = pin_gauge
do_freeze = .FALSE.
IF(PRESENT(freeze_flow)) do_freeze = freeze_flow
psi_dir = .TRUE.
IF(PRESENT(psi_outer_dirichlet)) psi_dir = psi_outer_dirichlet

ALLOCATE(self%fluid_reg(nreg), self%iron_reg(nreg), self%j_coil(nreg))
self%fluid_reg = fluid_reg
self%iron_reg = .FALSE.
IF(PRESENT(iron_reg))THEN
  IF(SIZE(iron_reg)/=nreg) CALL oft_abort('iron_reg must have one entry per mesh region', &
                                          'setup_empump_td',__FILE__)
  self%iron_reg = iron_reg
END IF
self%has_iron = ANY(self%iron_reg)
self%j_coil = 0.d0

!------------------------------------------------------------------------------
! xmhd_2d object. Material arrays are allocated before setup() so its
! defaults do not overwrite them. Non-fluid regions get placeholder density
! and viscosity: their velocity rows are pinned, but the residual still
! divides by the density there.
!------------------------------------------------------------------------------
ALLOCATE(mhd)
mhd%incomp = .TRUE.
mhd%cyl_flag = cyl
mhd%dt = dt
mhd%jac_dt = dt
mhd%den_scale = 1.d0 ! field 1 holds n = 1, so m_i*n is the mass density
mhd%lin_tol = lin_tol
mhd%nl_tol = nl_tol
ALLOCATE(mhd%m_i(nreg), mhd%nu(nreg), mhd%eta(nreg,2))
mhd%m_i = MERGE(dens_reg, 1.d0, fluid_reg)
mhd%nu = MERGE(visc_reg, 0.d0, fluid_reg)
mhd%eta = eta_reg
CALL mhd%setup(mg_mesh, order)
self%mug => mhd
mesh => mg_mesh%smesh
lag_rep => oft_blagrange

!------------------------------------------------------------------------------
! Boundary conditions. setup() installed xmhd_2d's defaults, some of which
! alias the FE library's own gbe array, so the pointers are dropped (not
! deallocated) and every flag gets its own storage.
!------------------------------------------------------------------------------
NULLIFY(mhd%n_bc, mhd%velx_bc, mhd%vely_bc, mhd%velz_bc, mhd%T_bc, mhd%psi_bc, mhd%by_bc)
ALLOCATE(mhd%n_bc(mhd%fe_rep%fields(1)%fe%ne))
ALLOCATE(mhd%velx_bc(mhd%fe_rep%fields(2)%fe%ne))
ALLOCATE(mhd%vely_bc(mhd%fe_rep%fields(3)%fe%ne))
ALLOCATE(mhd%velz_bc(mhd%fe_rep%fields(4)%fe%ne))
ALLOCATE(mhd%T_bc(mhd%fe_rep%fields(5)%fe%ne))
ALLOCATE(mhd%psi_bc(mhd%fe_rep%fields(6)%fe%ne))
ALLOCATE(mhd%by_bc(mhd%fe_rep%fields(7)%fe%ne))

mhd%n_bc = .TRUE. ! incompressible: density is not evolved
!---Velocity: no-slip on the outer boundary. v_Z is finite on the axis, which is
!   part of the outer boundary, so it is released there.
mhd%velx_bc = lag_rep%global%gbe
mhd%velz_bc = lag_rep%global%gbe
IF(ASSOCIATED(mhd%axis_bc)) mhd%velz_bc = mhd%velz_bc .AND. .NOT.mhd%axis_bc
IF(tor_flow)THEN
  mhd%vely_bc = lag_rep%global%gbe
ELSE
  mhd%vely_bc = .TRUE.
END IF
mhd%T_bc = .TRUE.
mhd%psi_bc = lag_rep%global%gbe
IF(do_by)THEN
  mhd%by_bc = lag_rep%global%gbe
ELSE
  mhd%by_bc = .TRUE.
END IF

ALLOCATE(cell_dofs(lag_rep%nce), cell_dofs_p(mhd%fe_rep%fields(5)%fe%nce))
ALLOCATE(gauge_set(nreg))
gauge_set = .FALSE.
DO i=1,mesh%nc
  CALL lag_rep%ncdofs(i,cell_dofs)
  CALL mhd%fe_rep%fields(5)%fe%ncdofs(i,cell_dofs_p)
  IF(fluid_reg(mesh%reg(i)))THEN
    !---Pressure is solved on every DOF of a fluid cell, interface DOFs included
    DO j=1,SIZE(cell_dofs_p)
      mhd%T_bc(cell_dofs_p(j)) = .FALSE.
    END DO
  ELSE
    !---Solids, coils, iron and air do not move. Pinning every DOF of these
    !   cells includes the interface DOFs, which makes the fluid wall no-slip.
    DO j=1,SIZE(cell_dofs)
      mhd%velx_bc(cell_dofs(j)) = .TRUE.
      mhd%vely_bc(cell_dofs(j)) = .TRUE.
      mhd%velz_bc(cell_dofs(j)) = .TRUE.
    END DO
  END IF
END DO
!---Pressure gauge: the enclosed incompressible problem fixes p only up to a
!   constant, once per fluid region. Done after the loop above so a later cell
!   cannot release the pinned DOF.
IF(do_gauge)THEN
  DO i=1,mesh%nc
    IF(.NOT.fluid_reg(mesh%reg(i)) .OR. gauge_set(mesh%reg(i)))CYCLE
    CALL mhd%fe_rep%fields(5)%fe%ncdofs(i,cell_dofs_p)
    mhd%T_bc(cell_dofs_p(1)) = .TRUE.
    gauge_set(mesh%reg(i)) = .TRUE.
  END DO
END IF
DEALLOCATE(cell_dofs, cell_dofs_p, gauge_set)
!---psi left free on the outer boundary if asked; the axis is re-pinned below
IF(.NOT.psi_dir) mhd%psi_bc = .FALSE.
!---EM-only: the flow is frozen. A pinned row holds its old value, and the
!   initial state below is at rest, so v = p = 0 exactly at every step.
IF(do_freeze)THEN
  mhd%velx_bc = .TRUE.
  mhd%vely_bc = .TRUE.
  mhd%velz_bc = .TRUE.
  mhd%T_bc = .TRUE.
END IF
!---Axis regularity (cylindrical only; a no-op in slab)
CALL mhd%apply_axis_bc()

!------------------------------------------------------------------------------
! Initial state: at rest, no field, n = 1
!------------------------------------------------------------------------------
self%u => mhd%u
CALL self%u%set(0.d0)
CALL self%u%set(1.d0, 1)
CALL mhd%fe_rep%vec_create(self%u_old)
CALL mhd%fe_rep%vec_create(self%rhs)

!------------------------------------------------------------------------------
! Setup linear and nonlinear solvers
!------------------------------------------------------------------------------
ALLOCATE(self%nlfun)
self%nlfun%dt = dt
self%nlfun%parent_sim => self

ALLOCATE(self%pre)
self%pre%A => mhd%jacobian

ALLOCATE(self%mfmat)
CALL self%mfmat%setup(self%u_old, self%nlfun)
self%mfmat%b0 = 1.d-4

ALLOCATE(self%mf_solver)
self%mf_solver%A => self%mfmat
self%mf_solver%its = 1000
self%mf_solver%nrits = 20
self%mf_solver%atol = lin_tol
self%mf_solver%itplot = 1
oft_env%pm = self%pm
self%mf_solver%pm = self%pm
self%mf_solver%pre => self%pre

self%nksolver%A => self%nlfun
self%nksolver%J_inv => self%mf_solver
self%nksolver%its = 20
self%nksolver%atol = nl_tol
self%nksolver%rtol = 1.d-20 ! disable the relative tolerance
self%nksolver%backtrack = .FALSE.
self%nksolver%J_update => empump_mfnk_update
self%nksolver%up_freq = 1

self%mug%rst_freq = self%rst_freq
IF(self%save_rst) CALL self%rst_save(0.d0)
end subroutine setup_empump_td

!---------------------------------------------------------------------------
!> Attempt one timestep of length dt
!---------------------------------------------------------------------------
subroutine step_empump_td(self,time,dt,nl_its,lin_its)
CLASS(oft_empump_td), TARGET, INTENT(inout) :: self !< Simulation object
REAL(r8), INTENT(inout) :: time !< Current time [s]; advanced on success
REAL(r8), INTENT(in) :: dt !< Timestep [s]
INTEGER(i4), INTENT(out) :: nl_its !< Nonlinear iterations (-1 if the solve failed)
INTEGER(i4), INTENT(out) :: lin_its !< Linear iterations
current_sim => self
!---Checked here rather than in setup, because drivers edit the velocity flags after
!   setup (e.g. pinning everything for an EM-only run)
IF(.NOT.self%flow_checked)THEN
  self%flow_checked = .TRUE.
  IF(.NOT.(ALL(self%mug%velx_bc) .AND. ALL(self%mug%vely_bc) .AND. ALL(self%mug%velz_bc))) &
    CALL oft_warn("empump_td: the flow is being solved (some velocity DOFs are free). Only the " // &
                  "frozen-flow (EM-only) configuration has been tested; treat flow results as unvalidated.")
END IF
!---A new timestep has to reach every operator that uses it
IF(dt /= self%nlfun%dt)THEN
  self%nlfun%dt = ABS(dt)
  self%mug%dt = ABS(dt)
  self%mug%jac_dt = ABS(dt)
END IF
CALL self%u_old%add(0.d0,1.d0,self%u)
!---Preconditioner: xmhd_2d's approximate Jacobian, which has the iron at mu0.
!   The iron terms are added only if one of the optional preconditioner
!   terms is on
CALL build_approx_jacobian(self%mug, self%u)
IF(self%has_iron .AND. (self%iron_jac_lift .OR. self%iron_jac_deriv))THEN
  CALL add_iron_jacobian(self, self%mug%jacobian, self%u)
  !---Re-assemble so the cached diagonal includes the iron entries
  BLOCK
  CLASS(oft_vector), POINTER :: tmp
  CALL self%mug%fe_rep%vec_create(tmp)
  CALL self%mug%jacobian%assemble(tmp)
  CALL tmp%delete()
  DEALLOCATE(tmp)
  END BLOCK
END IF
CALL self%pre%update(.TRUE.)
!---Right-hand side: xmhd_2d's dt = 0 operator applied to the old state (the
!   mass terms; pinned rows return their old values), plus the coil source
self%mug%nlfun%dt = 0.d0
CALL self%mug%nlfun%apply_real(self%u_old, self%rhs)
CALL add_source_terms(self, self%rhs, self%nlfun%dt)
CALL self%nksolver%apply(self%u, self%rhs)
lin_its = self%nksolver%lits
IF(self%nksolver%cits < 0)THEN
  !---Failed: restore the state
  CALL self%u%add(0.d0,1.d0,self%u_old)
  nl_its = -1
  RETURN
END IF
time = time + self%nlfun%dt
self%mug%t = time
nl_its = self%nksolver%nlits
self%rst_count = self%rst_count + 1
IF(self%save_rst .AND. MOD(self%rst_count, self%rst_freq)==0) CALL self%rst_save(time)
end subroutine step_empump_td

!---------------------------------------------------------------------------
!> F(a): xmhd_2d's residual at the current dt, plus the iron correction
!---------------------------------------------------------------------------
subroutine nlfun_apply(self,a,b)
CLASS(oft_empump_td_nlfun), INTENT(inout) :: self
CLASS(oft_vector), TARGET, INTENT(inout) :: a !< Current iterate
CLASS(oft_vector), INTENT(inout) :: b !< Result
self%parent_sim%mug%nlfun%dt = self%dt
CALL self%parent_sim%mug%nlfun%apply_real(a,b)
CALL add_iron_terms(self%parent_sim, a, b, self%dt)
end subroutine nlfun_apply

!---------------------------------------------------------------------------
!> Add the prescribed coil current to the psi right-hand side
!!
!! xmhd_2d's psi rows carry dt*grad(phi).grad(psi)/R (cylindrical) or
!! dt*grad(phi).grad(psi) (slab), the weak form of -dt*Delta*psi/R or
!! -dt*grad^2 psi. Both equal dt*mu0*J by Ampere's law, so the source is
!! dt*mu0*J*phi in either geometry.
!---------------------------------------------------------------------------
subroutine add_source_terms(self,b,dt)
CLASS(oft_empump_td), INTENT(inout) :: self
CLASS(oft_vector), INTENT(inout) :: b !< Right-hand side; field 6 is modified
REAL(r8), INTENT(in) :: dt !< Timestep [s]
REAL(r8), POINTER, DIMENSION(:) :: psi_res
TYPE(oft_quad_type), POINTER :: quad
INTEGER(i4) :: i,m,jr
INTEGER(i4), ALLOCATABLE :: cell_dofs(:)
REAL(r8) :: jac_mat(3,4),jac_det,src
REAL(r8), ALLOCATABLE :: basis_vals(:),res_loc(:)
LOGICAL :: curved
IF(ALL(ABS(self%j_coil) < 1.d-30))RETURN
quad => lag_rep%quad
NULLIFY(psi_res)
CALL b%get_local(psi_res, 6)
!$omp parallel private(m,jr,curved,cell_dofs,jac_mat,jac_det,src,basis_vals,res_loc)
ALLOCATE(cell_dofs(lag_rep%nce), basis_vals(lag_rep%nce), res_loc(lag_rep%nce))
!$omp do
DO i=1,mesh%nc
  src = self%j_coil(mesh%reg(i))
  IF(ABS(src) < 1.d-30)CYCLE
  curved = cell_is_curved(mesh,i)
  CALL lag_rep%ncdofs(i,cell_dofs)
  res_loc = 0.d0
  DO m=1,quad%np
    IF(curved.OR.(m==1))CALL mesh%jacobian(i,quad%pts(:,m),jac_mat,jac_det)
    DO jr=1,lag_rep%nce
      CALL oft_blag_eval(lag_rep,i,jr,quad%pts(:,m),basis_vals(jr))
      res_loc(jr) = res_loc(jr) + basis_vals(jr)*jac_det*quad%wts(m)
    END DO
  END DO
  DO jr=1,lag_rep%nce
    IF(self%mug%psi_bc(cell_dofs(jr)))CYCLE ! pinned rows hold their old value
    !$omp atomic
    psi_res(cell_dofs(jr)) = psi_res(cell_dofs(jr)) + dt*mu0*src*res_loc(jr)
  END DO
END DO
DEALLOCATE(cell_dofs, basis_vals, res_loc)
!$omp end parallel
CALL b%restore_local(psi_res, 6)
DEALLOCATE(psi_res)
end subroutine add_source_terms

!---------------------------------------------------------------------------
!> Residual correction for the nonlinear permeability of iron regions
!!
!! xmhd_2d assembles the psi stiffness with mu = mu0 everywhere. In iron this
!! adds dt*(nu(|B|) - 1)*grad(phi).grad(psi)/Rfac, with nu = 1/mu_r and
!! Rfac = R (cylindrical) or 1 (slab), so the total is nu(|B|) times the
!! stiffness. |B| is the in-plane field |grad(psi)|/Rfac; B_y is not included.
!---------------------------------------------------------------------------
subroutine add_iron_terms(self,a,b,dt)
CLASS(oft_empump_td), INTENT(inout) :: self
CLASS(oft_vector), INTENT(inout) :: a !< Current iterate
CLASS(oft_vector), INTENT(inout) :: b !< Residual; field 6 is modified
REAL(r8), INTENT(in) :: dt !< Timestep [s]
REAL(r8), POINTER, DIMENSION(:) :: psi_weights,psi_res
TYPE(oft_quad_type), POINTER :: quad
INTEGER(i4) :: i,m,jr
INTEGER(i4), ALLOCATABLE :: cell_dofs(:)
REAL(r8) :: jac_mat(3,4),jac_det,int_factor,coords(3),dpsi(3),Bmag,nu,Rfac
REAL(r8), ALLOCATABLE :: basis_grads(:,:),psi_loc(:),res_loc(:)
LOGICAL :: curved,cyl
IF(.NOT.self%has_iron)RETURN
IF(dt <= 0.d0)RETURN ! the correction is proportional to dt
cyl = self%mug%cyl_flag
quad => lag_rep%quad
NULLIFY(psi_weights, psi_res)
CALL a%get_local(psi_weights, 6)
CALL b%get_local(psi_res, 6)
!$omp parallel private(m,jr,curved,cell_dofs,jac_mat,jac_det,int_factor,coords, &
!$omp   dpsi,Bmag,nu,Rfac,basis_grads,psi_loc,res_loc)
ALLOCATE(basis_grads(3,lag_rep%nce), cell_dofs(lag_rep%nce), &
         psi_loc(lag_rep%nce), res_loc(lag_rep%nce))
!$omp do schedule(dynamic,16)
DO i=1,mesh%nc
  IF(.NOT.self%iron_reg(mesh%reg(i)))CYCLE
  curved = cell_is_curved(mesh,i)
  CALL lag_rep%ncdofs(i,cell_dofs)
  psi_loc = psi_weights(cell_dofs)
  res_loc = 0.d0
  DO m=1,quad%np
    IF(curved.OR.(m==1))CALL mesh%jacobian(i,quad%pts(:,m),jac_mat,jac_det)
    DO jr=1,lag_rep%nce
      CALL oft_blag_geval(lag_rep,i,jr,quad%pts(:,m),basis_grads(:,jr),jac_mat)
    END DO
    coords = mesh%log2phys(i,quad%pts(:,m))
    Rfac = 1.d0
    IF(cyl) Rfac = coords(1) + gs_epsilon
    int_factor = jac_det*quad%wts(m)
    dpsi = 0.d0
    DO jr=1,lag_rep%nce
      dpsi = dpsi + psi_loc(jr)*basis_grads(:,jr)
    END DO
    Bmag = SQRT(SUM(dpsi(1:2)**2))/Rfac
    nu = 1.d0/gs_mu_of_B(Bmag)
    DO jr=1,lag_rep%nce
      res_loc(jr) = res_loc(jr) + dt*(nu-1.d0) &
        *DOT_PRODUCT(basis_grads(1:2,jr),dpsi(1:2))*int_factor/Rfac
    END DO
  END DO
  DO jr=1,lag_rep%nce
    IF(self%mug%psi_bc(cell_dofs(jr)))CYCLE ! pinned rows are owned by the BC
    !$omp atomic
    psi_res(cell_dofs(jr)) = psi_res(cell_dofs(jr)) + res_loc(jr)
  END DO
END DO
DEALLOCATE(basis_grads, cell_dofs, psi_loc, res_loc)
!$omp end parallel
CALL b%restore_local(psi_res, 6)
DEALLOCATE(psi_res, psi_weights)
end subroutine add_iron_terms

!---------------------------------------------------------------------------
!> Jacobian of add_iron_terms, added to the psi-psi block. Optional, as it may
!! not improve convergence
!! (A) dt*(nu - 1)*grad(phi_r).grad(phi_c)/Rfac, which lifts xmhd_2d's mu0
!!     stiffness to nu(|B|). Included when iron_jac_lift (default off);
!!     without it the preconditioner sees the iron as mu0.
!! (B) the field dependence of mu, rank one in g_r = grad(phi_r).grad(psi):
!!     -dt*(dmu/dB)/mu^2 * g_r*g_c/(|grad psi|*Rfac^2). Included when
!!     iron_jac_deriv (default off)
!! Preconditioner only: the converged solution does not depend on either term.
!---------------------------------------------------------------------------
subroutine add_iron_jacobian(self,mat,a)
CLASS(oft_empump_td), INTENT(inout) :: self
CLASS(oft_matrix), POINTER, INTENT(inout) :: mat !< Matrix to add entries to
CLASS(oft_vector), INTENT(inout) :: a !< State to linearize about
REAL(r8), POINTER, DIMENSION(:) :: psi_weights
TYPE(oft_quad_type), POINTER :: quad
INTEGER(i4) :: i,m,jr,jc
INTEGER(i4), ALLOCATABLE :: cell_dofs(:)
REAL(r8) :: jac_mat(3,4),jac_det,int_factor,coords(3),dpsi(3),Bmag,absP
REAL(r8) :: nu,mu_loc,dmu,coefA,coefB,dt,Rfac
REAL(r8), ALLOCATABLE :: basis_grads(:,:),psi_loc(:),jac_loc(:,:),gdotp(:)
LOGICAL :: curved,cyl
REAL(r8), PARAMETER :: absP_min = 1.d-12 ! below this |grad psi| the rank-1 term is dropped
IF(.NOT.self%has_iron)RETURN
dt = self%nlfun%dt
cyl = self%mug%cyl_flag
quad => lag_rep%quad
NULLIFY(psi_weights)
CALL a%get_local(psi_weights, 6)
!$omp parallel private(m,jr,jc,curved,cell_dofs,jac_mat,jac_det,int_factor,coords, &
!$omp   dpsi,Bmag,absP,nu,mu_loc,dmu,coefA,coefB,Rfac,basis_grads,psi_loc,jac_loc,gdotp)
ALLOCATE(basis_grads(3,lag_rep%nce), cell_dofs(lag_rep%nce), psi_loc(lag_rep%nce), &
         jac_loc(lag_rep%nce,lag_rep%nce), gdotp(lag_rep%nce))
!$omp do schedule(dynamic,16)
DO i=1,mesh%nc
  IF(.NOT.self%iron_reg(mesh%reg(i)))CYCLE
  curved = cell_is_curved(mesh,i)
  CALL lag_rep%ncdofs(i,cell_dofs)
  psi_loc = psi_weights(cell_dofs)
  jac_loc = 0.d0
  DO m=1,quad%np
    IF(curved.OR.(m==1))CALL mesh%jacobian(i,quad%pts(:,m),jac_mat,jac_det)
    DO jr=1,lag_rep%nce
      CALL oft_blag_geval(lag_rep,i,jr,quad%pts(:,m),basis_grads(:,jr),jac_mat)
    END DO
    coords = mesh%log2phys(i,quad%pts(:,m))
    Rfac = 1.d0
    IF(cyl) Rfac = coords(1) + gs_epsilon
    int_factor = jac_det*quad%wts(m)
    dpsi = 0.d0
    DO jr=1,lag_rep%nce
      dpsi = dpsi + psi_loc(jr)*basis_grads(:,jr)
    END DO
    absP = SQRT(SUM(dpsi(1:2)**2))
    Bmag = absP/Rfac
    mu_loc = gs_mu_of_B(Bmag, dmu_dB=dmu)
    nu = 1.d0/mu_loc
    coefA = 0.d0
    IF(self%iron_jac_lift) coefA = dt*(nu-1.d0)*int_factor/Rfac
    coefB = 0.d0
    IF(self%iron_jac_deriv .AND. absP > absP_min) &
      coefB = -dt*dmu/(mu_loc*mu_loc)*int_factor/(absP*Rfac*Rfac)
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
    IF(self%mug%psi_bc(cell_dofs(jr)))CYCLE ! pinned rows keep their identity row
    CALL mat%add_values([cell_dofs(jr)], cell_dofs, RESHAPE(jac_loc(jr,:),[1,lag_rep%nce]), &
                        1, lag_rep%nce, 6, 6)
  END DO
  !$omp end critical
END DO
DEALLOCATE(basis_grads, cell_dofs, psi_loc, jac_loc, gdotp)
!$omp end parallel
DEALLOCATE(psi_weights)
end subroutine add_iron_jacobian

!---------------------------------------------------------------------------
!> Save the current state as 'xmhd2d_NNNNN.rst', indexed by rst_count, so
!! xmhd_2d_plot can read it
!---------------------------------------------------------------------------
subroutine empump_rst_save(self,t)
CLASS(oft_empump_td), INTENT(inout) :: self
REAL(r8), INTENT(in) :: t !< Solution time [s]
CHARACTER(LEN=XMHD_RST_LEN) :: rst_char
INTEGER(i4) :: io_stat,rst_tmp
100 FORMAT (I XMHD_RST_LEN.XMHD_RST_LEN)
WRITE(rst_char,100) self%rst_count
READ(rst_char,100,IOSTAT=io_stat) rst_tmp
IF((io_stat/=0).OR.(rst_tmp/=self%rst_count)) &
  CALL oft_abort('Step count exceeds restart format width','empump_rst_save',__FILE__)
CALL self%mug%rst_save(self%u, t, self%nlfun%dt, 'xmhd2d_'//rst_char//'.rst', 'U')
end subroutine empump_rst_save

!---------------------------------------------------------------------------
!> Convert the restart files to XDMF with xmhd_2d_plot (options group
!! `xmhd_plot_options`). B is written as 'B*R' in cylindrical, 'B' in slab.
!---------------------------------------------------------------------------
subroutine empump_td_plot(self)
CLASS(oft_empump_td), INTENT(inout) :: self
self%mug%rst_freq = self%rst_freq
CALL xmhd_2d_plot(self%mug)
end subroutine empump_td_plot

!---------------------------------------------------------------------------
!> Delete the simulation's solver objects
!---------------------------------------------------------------------------
subroutine delete_empump_td(self)
CLASS(oft_empump_td), INTENT(inout) :: self
IF(ASSOCIATED(self%mfmat))THEN
  CALL self%mfmat%delete()
  DEALLOCATE(self%mfmat)
END IF
IF(ASSOCIATED(self%mf_solver))THEN
  CALL self%mf_solver%delete()
  DEALLOCATE(self%mf_solver)
END IF
IF(ASSOCIATED(self%pre))THEN
  CALL self%pre%delete()
  DEALLOCATE(self%pre)
END IF
CALL self%nksolver%delete()
IF(ASSOCIATED(self%nlfun))DEALLOCATE(self%nlfun)
IF(ASSOCIATED(self%u_old))THEN
  CALL self%u_old%delete()
  CALL self%rhs%delete()
  DEALLOCATE(self%u_old, self%rhs)
END IF
NULLIFY(self%u)
IF(ALLOCATED(self%fluid_reg))DEALLOCATE(self%fluid_reg, self%iron_reg, self%j_coil)
end subroutine delete_empump_td

!---------------------------------------------------------------------------
!> Update the matrix-free Jacobian with the current Newton iterate
!---------------------------------------------------------------------------
subroutine empump_mfnk_update(a)
CLASS(oft_vector), TARGET, INTENT(inout) :: a !< Current iterate
CALL current_sim%mfmat%update(a)
end subroutine empump_mfnk_update

END MODULE empump_td
