!---------------------------------------------------------------------------
! Hartmann flow driven by an out-of-plane body force
!---------------------------------------------------------------------------
!> Steady Hartmann flow in a plane channel, driven by a body force in the
!! out-of-plane (y) direction.
!!
!! The channel occupies -a <= z <= a with no-slip walls, and carries a uniform
!! magnetic field B0 normal to those walls. The flow v_y is out-of-plane, so it
!! feels no pressure gradient (y is the ignorable direction) and is driven
!! purely by the body force `force_y`, balanced by viscosity and by the Lorentz
!! force from the field it induces.
!!
!! In the quasi-static limit the induced out-of-plane field b_y satisfies
!!     0 = B0 dv/dz + eta d^2(b_y)/dz^2
!! and the momentum balance is
!!     0 = F + mu d^2v/dz^2 + (B0/mu0) d(b_y)/dz
!! which combine to the Hartmann equation
!!     mu v'' - sigma B0^2 v = -F
!! with solution
!!     v(z) = (F/(sigma B0^2)) * [1 - cosh(Ha z/a)/cosh(Ha)],
!!     Ha = B0 a sqrt(sigma/mu)
!! Leaving b_y free gives the natural (zero-gradient) wall condition, which is
!! the constant of integration for which the expression above is exact.
!!
!! The field is imposed through the poloidal flux, psi = B0*x, which is frozen
!! for the whole run. That makes psi vary along the streamwise direction, so the
!! domain cannot be periodic there; instead v_y is left free on the streamwise
!! ends, which is the correct traction-free condition for fully developed flow.
!---------------------------------------------------------------------------
MODULE hartmann_helpers
USE oft_local
IMPLICIT NONE
REAL(r8) :: B0_par      !< Applied field strength [T]
REAL(r8) :: a_par       !< Channel half-width [m]
REAL(r8) :: Ha_par      !< Hartmann number
REAL(r8) :: vscale_par  !< F*a^2/mu, the Shercliff velocity scale [m/s]
REAL(r8) :: aspect_par  !< duct half-width across the field, over a
LOGICAL  :: insulating_par !< .TRUE. -> Shercliff (insulating), .FALSE. -> conducting
CONTAINS

!> Poloidal flux giving a uniform field B0 normal to the channel walls
SUBROUTINE psi_init(pt, val)
REAL(r8), INTENT(in) :: pt(3)
REAL(r8), INTENT(out) :: val
val = B0_par*pt(1)
END SUBROUTINE psi_init

!> Shercliff's analytic solution for fully developed flow in a rectangular duct
!! with insulating walls, driven by a uniform body force.
!!
!! Non-dimensionalising by V = v*mu/(F a^2), with zeta along the applied field and
!! xi across it, the steady equations are
!!     grad^2 V + Ha dB/dzeta = -1,  grad^2 B + Ha dV/dzeta = 0
!! Substituting V± = V ± B decouples them into
!!     grad^2 V± ± Ha dV±/dzeta = -1
!! with V± = 0 on every wall (no-slip, and B=0 for insulating walls). Expanding in
!! the cross-field eigenfunctions cos(lam_n xi) reduces each to a two-point ODE
!! whose solution is written below in a form with only non-positive exponents, so
!! it stays finite at large Ha and large n.
FUNCTION shercliff(xi, zeta, Ha, aspect, nmax) RESULT(V)
REAL(r8), INTENT(in) :: xi, zeta, Ha, aspect
INTEGER(i4), INTENT(in) :: nmax
REAL(r8) :: V, lam, Amp, fp, disc, r1, r2, den, fsum, sgn, term
INTEGER(i4) :: n, sgn_i
V = 0.d0
DO n=1,nmax
  lam = (2*n-1)*pi/(2.d0*aspect)
  Amp = 2.d0*((-1)**(n+1))/(aspect*lam)   ! coefficient of 1 in the cos series
  fp = Amp/lam**2                          ! particular solution
  fsum = 0.d0
  DO sgn_i=1,2
    sgn = 1.d0; IF(sgn_i==2)sgn = -1.d0    ! V+ then V-
    disc = SQRT(Ha**2 + 4.d0*lam**2)
    r1 = (-sgn*Ha + disc)/2.d0             ! > 0
    r2 = (-sgn*Ha - disc)/2.d0             ! < 0
    den = 1.d0 - EXP(2.d0*(r2-r1))
    term = fp*(1.d0 &
      + (EXP(2.d0*r2) - 1.d0)*EXP(r1*(zeta-1.d0))/den &
      - (1.d0 - EXP(-2.d0*r1))*EXP(r2*(zeta+1.d0))/den)
    fsum = fsum + term
  END DO
  V = V + 0.5d0*fsum*COS(lam*xi)
END DO
END FUNCTION shercliff

!> Fully-developed duct profile for PERFECTLY CONDUCTING walls.
!!
!! With conducting walls the circuit closes through the walls with no resistance,
!! so E=0 and Ohm's law gives J = sigma*v*B0 pointwise -- there is no <v> term as
!! in the insulating case. The braking is then sigma*B0^2*v everywhere and the
!! momentum equation collapses to a single modified Helmholtz problem,
!!     mu grad^2 v - sigma B0^2 v = -F,   v = 0 on all walls,
!! i.e. non-dimensionally  grad^2 V - Ha^2 V = -1. No V+/V- splitting is needed.
FUNCTION hunt_cond(xi, zeta, Ha, aspect, nmax) RESULT(V)
REAL(r8), INTENT(in) :: xi, zeta, Ha, aspect
INTEGER(i4), INTENT(in) :: nmax
REAL(r8) :: V, lam, Amp, m, ratio
INTEGER(i4) :: n
V = 0.d0
DO n=1,nmax
  lam = (2*n-1)*pi/(2.d0*aspect)
  Amp = 2.d0*((-1)**(n+1))/(aspect*lam)
  m = SQRT(lam**2 + Ha**2)
  !---cosh(m*zeta)/cosh(m) written with only non-positive exponents
  ratio = (EXP(m*(zeta-1.d0)) + EXP(-m*(zeta+1.d0)))/(1.d0 + EXP(-2.d0*m))
  V = V + (Amp/(lam**2 + Ha**2))*(1.d0 - ratio)*COS(lam*xi)
END DO
END FUNCTION hunt_cond

!> Analytic duct profile, in physical units
SUBROUTINE vy_analytic(pt, val)
REAL(r8), INTENT(in) :: pt(3)
REAL(r8), INTENT(out) :: val
!---pt(2) is along the applied field (Hartmann walls), pt(1) across it
IF(insulating_par)THEN
  val = vscale_par*shercliff(pt(1)/a_par, pt(2)/a_par, Ha_par, aspect_par, 200)
ELSE
  val = vscale_par*hunt_cond(pt(1)/a_par, pt(2)/a_par, Ha_par, aspect_par, 200)
END IF
END SUBROUTINE vy_analytic

END MODULE hartmann_helpers

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
USE oft_blag_operators, ONLY: oft_blag_getmop, oft_blag_project, oft_lag_brinterp
USE oft_scalar_inits, ONLY: poss_scalar_bfield
USE mhd_utils, ONLY: proton_mass, mu0
USE diagnostic, ONLY: scal_energy_2d
USE fem_utils, ONLY: diff_interp_2d, bfem_map_flag
USE xmhd_2d
USE hartmann_helpers
IMPLICIT NONE
INTEGER(i4) :: io_unit,ierr,i,j,ed,p1,p2
REAL(r8), POINTER :: vec_vals(:)
TYPE(oft_xmhd_2d_sim) :: mhd_sim
TYPE(multigrid_mesh) :: mg_mesh
TYPE(oft_lag_brinterp), TARGET :: ana_field, num_field
TYPE(diff_interp_2d) :: err_field
TYPE(poss_scalar_bfield) :: field_init
CLASS(oft_solver), POINTER :: minv => NULL()
CLASS(oft_matrix), POINTER :: mop => NULL()
CLASS(oft_vector), POINTER :: u,v,vy_ana,tmp
LOGICAL, ALLOCATABLE :: vert_flag(:),edge_flag(:),wall_flag(:)
REAL(r8) :: err_num,err_den,rho,B0,vmax
!---Runtime options
INTEGER(i4) :: order = 2
INTEGER(i4) :: nsteps = 300
INTEGER(i4) :: rst_freq = 50
INTEGER(i4) :: ittarget = 1000
REAL(r8) :: a_chan = 1.d-2      !< duct half-width ALONG the field (Hartmann walls) [m]
REAL(r8) :: b_chan = 1.d-2      !< duct half-width ACROSS the field (side walls) [m]
REAL(r8) :: Ha = 1.d0           !< Hartmann number
REAL(r8) :: rho_fluid = 1800.d0 !< mass density [kg/m^3]
REAL(r8) :: mu_fluid = 4.d-3    !< dynamic viscosity [Pa-s]
REAL(r8) :: sigma = 200.d0      !< electrical conductivity [S/m]
REAL(r8) :: force_y = 1.d2      !< driving body force density [N/m^3]
REAL(r8) :: m_i = 205.d0        !< ion mass [proton masses]
REAL(r8) :: t0 = 1.d0
REAL(r8) :: dt = 1.d0
REAL(r8) :: den_scale = 1.d27
REAL(r8) :: lin_tol = 1.d-13
REAL(r8) :: nl_tol = 1.d-11
REAL(r8) :: wall_tol = 1.d-8
LOGICAL :: pm = .FALSE.
LOGICAL :: use_mfnk = .TRUE.
LOGICAL :: insulating = .FALSE.  !< .TRUE. pins b_y=0 on the Hartmann walls (insulating); .FALSE. leaves b_y free (perfectly conducting / short-circuited)

NAMELIST/hartmann_options/order,nsteps,rst_freq,ittarget,a_chan,b_chan,Ha,rho_fluid, &
mu_fluid,sigma,force_y,m_i,t0,dt,den_scale,lin_tol,nl_tol,wall_tol,pm,use_mfnk,insulating
CALL oft_init
OPEN(NEWUNIT=io_unit,FILE=oft_env%ifile)
READ(io_unit,hartmann_options,IOSTAT=ierr)
CLOSE(io_unit)
IF(ierr/=0)CALL oft_abort('Error reading "hartmann_options"','hartmann',__FILE__)
!---------------------------------------------------------------------------
! Derived parameters
!---------------------------------------------------------------------------
!---Ha = B0*a*sqrt(sigma/mu) fixes the field for the requested Hartmann number
B0 = Ha/(a_chan*SQRT(sigma/mu_fluid))
rho = m_i*proton_mass*(rho_fluid/(m_i*proton_mass))
B0_par = B0
a_par = a_chan
Ha_par = Ha
aspect_par = b_chan/a_chan
vscale_par = force_y*a_chan**2/mu_fluid          ! Shercliff scaling F*a^2/mu
insulating_par = insulating
IF(insulating)THEN
  vmax = vscale_par*shercliff(0.d0,0.d0,Ha,aspect_par,200)
ELSE
  vmax = vscale_par*hunt_cond(0.d0,0.d0,Ha,aspect_par,200)
END IF
IF(oft_env%head_proc)THEN
  WRITE(*,'(A)')      ' Hartmann flow driven by an out-of-plane body force'
  WRITE(*,'(2X,A,ES14.6)')'Ha              = ',Ha
  WRITE(*,'(2X,A,ES14.6)')'B0        [T]   = ',B0
  WRITE(*,'(2X,A,ES14.6)')'nu      [m^2/s] = ',mu_fluid/rho_fluid
  WRITE(*,'(2X,A,ES14.6)')'eta     [m^2/s] = ',1.d0/(sigma*mu0)
  WRITE(*,'(2X,A,ES14.6)')'damping time[s] = ',rho_fluid/(sigma*B0**2)
  WRITE(*,'(2X,A,ES14.6)')'v_max analytic  = ',vmax
END IF
!---------------------------------------------------------------------------
! Setup grid and simulation
!---------------------------------------------------------------------------
CALL multigrid_construct_surf(mg_mesh)
mhd_sim%incomp = .TRUE.
CALL mhd_sim%setup(mg_mesh, order)

!---------------------------------------------------------------------------
! Initial conditions
!---------------------------------------------------------------------------
NULLIFY(u,v,mop,vec_vals)
CALL oft_blag_getmop(ML_oft_blagrange%current_level,mop)
CALL create_cg_solver(minv)
minv%A=>mop
minv%its=-2
CALL create_diag_pre(minv%pre)
CALL ML_oft_blagrange%vec_create(u)
CALL ML_oft_blagrange%vec_create(v)
CALL ML_oft_blagrange%vec_create(vy_ana)
CALL ML_oft_blagrange%vec_create(tmp)

!---Uniform density; rho = m_i*n is what enters the momentum equation
mhd_sim%den_scale = den_scale
CALL u%set(rho_fluid/(m_i*proton_mass))
CALL u%get_local(vec_vals)
vec_vals = vec_vals/den_scale
CALL mhd_sim%u%restore_local(vec_vals,1)

!---Fluid starts at rest
CALL u%set(0.d0)
CALL u%get_local(vec_vals)
CALL mhd_sim%u%restore_local(vec_vals,2)
CALL mhd_sim%u%restore_local(vec_vals,3)
CALL mhd_sim%u%restore_local(vec_vals,4)
CALL mhd_sim%u%restore_local(vec_vals,7)

!---Uniform temperature (frozen; it does not couple to the out-of-plane flow)
CALL u%set(t0)
CALL u%get_local(vec_vals)
CALL mhd_sim%u%restore_local(vec_vals,5)

!---Applied field, via psi = B0*x. Frozen for the whole run.
field_init%func=>psi_init
field_init%mesh=>mesh
CALL oft_blag_project(ML_oft_blagrange%current_level,field_init,v)
CALL u%set(0.d0)
CALL minv%apply(u,v)
CALL u%get_local(vec_vals)
CALL mhd_sim%u%restore_local(vec_vals,6)

!---Analytic profile, for the error check below
field_init%func=>vy_analytic
CALL oft_blag_project(ML_oft_blagrange%current_level,field_init,v)
CALL u%set(0.d0)
CALL minv%apply(u,v)
CALL u%get_local(vec_vals)
CALL vy_ana%restore_local(vec_vals)

CALL u%delete; CALL v%delete; CALL mop%delete
DEALLOCATE(u,v,mop)
CALL minv%pre%delete; DEALLOCATE(minv%pre)
CALL minv%delete; DEALLOCATE(minv)

!---------------------------------------------------------------------------
! Boundary conditions
!---------------------------------------------------------------------------
!---Everything except v_y and b_y is held fixed: the in-plane velocity stays at
!   zero (so div(v)=0 exactly and the density never moves), the density and
!   temperature are prescribed, and psi holds the applied field.
!   b_y is left free everywhere, which gives the natural zero-gradient wall
!   condition assumed by the analytic profile.
NULLIFY(mhd_sim%n_bc,mhd_sim%velx_bc,mhd_sim%vely_bc,mhd_sim%velz_bc)
NULLIFY(mhd_sim%T_bc,mhd_sim%psi_bc,mhd_sim%by_bc)
ALLOCATE(mhd_sim%n_bc(oft_blagrange%ne))
ALLOCATE(mhd_sim%velx_bc(oft_blagrange%ne))
ALLOCATE(mhd_sim%vely_bc(oft_blagrange%ne))
ALLOCATE(mhd_sim%velz_bc(oft_blagrange%ne))
ALLOCATE(mhd_sim%T_bc(oft_blagrange%ne))
ALLOCATE(mhd_sim%psi_bc(oft_blagrange%ne))
ALLOCATE(mhd_sim%by_bc(oft_blagrange%ne))
mhd_sim%n_bc    = .TRUE.
mhd_sim%velx_bc = .TRUE.
mhd_sim%velz_bc = .TRUE.
mhd_sim%T_bc    = .TRUE.
mhd_sim%psi_bc  = .TRUE.
mhd_sim%by_bc   = .FALSE.

!---v_y: no-slip on the Hartmann walls (|z| = a) only. The streamwise ends are
!   left free, which is the traction-free condition for fully developed flow;
!   pinning them instead would make this a duct rather than a plane channel.
ALLOCATE(vert_flag(mg_mesh%smesh%np),edge_flag(mg_mesh%smesh%ne))
ALLOCATE(wall_flag(oft_blagrange%ne))
vert_flag = .FALSE.
edge_flag = .FALSE.
wall_flag = .FALSE.
DO i=1,mg_mesh%smesh%np
  IF(ABS(ABS(mg_mesh%smesh%r(2,i))-a_chan) < wall_tol)vert_flag(i)=.TRUE.  ! Hartmann walls
  IF(ABS(ABS(mg_mesh%smesh%r(1,i))-b_chan) < wall_tol)vert_flag(i)=.TRUE.  ! side walls
END DO
DO ed=1,mg_mesh%smesh%ne
  p1 = mg_mesh%smesh%le(1,ed)
  p2 = mg_mesh%smesh%le(2,ed)
  IF(vert_flag(p1).AND.vert_flag(p2))edge_flag(ed)=.TRUE.
END DO
CALL bfem_map_flag(oft_blagrange,vert_flag,edge_flag,wall_flag)
mhd_sim%vely_bc = wall_flag
!---Electrical wall condition.
!   insulating (default): b_y = 0 on every wall, so no current enters the walls and
!     the circuit closes entirely within the fluid. This is Shercliff's case, which
!     the analytic profile above assumes.
!   conducting: b_y left free, giving the natural db_y/dn = 0 (zero wall current,
!     hence E=0, short-circuited). That is Hunt's case and the analytic does NOT apply.
IF(insulating)mhd_sim%by_bc = wall_flag
IF(oft_env%head_proc)WRITE(*,'(2X,A,I6,A,I6)')'no-slip v_y DOFs = ', &
  COUNT(wall_flag),'  of ',oft_blagrange%ne
IF(COUNT(wall_flag)==0)CALL oft_abort('No wall DOFs found; check a_chan/wall_tol', &
  'hartmann',__FILE__)

!---------------------------------------------------------------------------
! Run
!---------------------------------------------------------------------------
mhd_sim%nu = mu_fluid/rho_fluid          ! kinematic viscosity [m^2/s]
mhd_sim%eta = 1.d0/(sigma*mu0)           ! magnetic diffusivity [m^2/s]
mhd_sim%force_y = force_y                ! out-of-plane drive [N/m^3]
mhd_sim%chi = 1.d0
mhd_sim%D_diff = 1.d0
mhd_sim%gamma = 5.d0/3.d0
mhd_sim%m_i = m_i*proton_mass
mhd_sim%dt = dt
mhd_sim%nsteps = nsteps
mhd_sim%rst_freq = rst_freq
mhd_sim%mfnk = use_mfnk
mhd_sim%lin_tol = lin_tol
mhd_sim%nl_tol = nl_tol
mhd_sim%ittarget = ittarget
oft_env%pm = pm

IF(oft_env%head_proc)THEN
  WRITE(*,'(2X,A,ES16.8)')'DBG m_i [kg]     = ',mhd_sim%m_i(1)
  WRITE(*,'(2X,A,ES16.8)')'DBG n   [m^-3]   = ',rho_fluid/(m_i*proton_mass)
  WRITE(*,'(2X,A,ES16.8)')'DBG rho = m_i*n  = ',mhd_sim%m_i(1)*(rho_fluid/(m_i*proton_mass))
  WRITE(*,'(2X,A,ES16.8)')'DBG nu           = ',mhd_sim%nu(1)
  WRITE(*,'(2X,A,ES16.8)')'DBG force_y      = ',mhd_sim%force_y(1)
  WRITE(*,'(2X,A,ES16.8)')'DBG a_chan       = ',a_chan
END IF
CALL mhd_sim%run_simulation()

!---------------------------------------------------------------------------
! Compare to the analytic Hartmann profile
!---------------------------------------------------------------------------
BLOCK
REAL(r8), POINTER :: dbg(:)
NULLIFY(dbg); CALL mhd_sim%u%get_local(dbg,6)
WRITE(*,'(2X,A,2ES14.6)')'DBG psi min/max = ',MINVAL(dbg),MAXVAL(dbg)
WRITE(*,'(2X,A,ES14.6)') 'DBG psi expected max = ',B0*b_chan
DEALLOCATE(dbg); NULLIFY(dbg)
CALL mhd_sim%u%get_local(dbg,7)
WRITE(*,'(2X,A,2ES14.6)')'DBG b_y min/max = ',MINVAL(dbg),MAXVAL(dbg)
DEALLOCATE(dbg)
END BLOCK
CALL mhd_sim%u%get_local(vec_vals,3)
CALL tmp%restore_local(vec_vals)
ana_field%u=>vy_ana
num_field%u=>tmp
CALL ana_field%setup(ML_oft_blagrange%current_level)
CALL num_field%setup(ML_oft_blagrange%current_level)
err_field%dim=1
err_field%a=>ana_field
err_field%b=>num_field
err_den = scal_energy_2d(mg_mesh%smesh,ana_field,order*2)
err_num = scal_energy_2d(mg_mesh%smesh,err_field,order*2)
IF(oft_env%head_proc)THEN
  WRITE(*,'(A)')       ' Hartmann profile check'
  WRITE(*,'(2X,A,ES16.8)')'final time [s]   = ',mhd_sim%t
  WRITE(*,'(2X,A,ES16.8)')'v_y max (num)    = ',MAXVAL(vec_vals)
  WRITE(*,'(2X,A,ES16.8)')'v_y max (analytic)= ',vmax
  WRITE(*,'(2X,A,ES16.8)')'relative L2 error= ',SQRT(err_num/err_den)
END IF
!---Profile dump: point DOFs come first in the Lagrange ordering, so entries
!   1..np pair with the mesh vertices.
BLOCK
REAL(r8), POINTER :: byv(:)
NULLIFY(byv); CALL mhd_sim%u%get_local(byv,7)
OPEN(NEWUNIT=io_unit,FILE='profile.txt')
DO i=1,mg_mesh%smesh%np
  WRITE(io_unit,'(4ES20.10)')mg_mesh%smesh%r(1,i),mg_mesh%smesh%r(2,i), &
    vec_vals(i),byv(i)
END DO
CLOSE(io_unit)
DEALLOCATE(byv)
END BLOCK
OPEN(NEWUNIT=io_unit,FILE='hartmann.results')
WRITE(io_unit,*)SQRT(err_num/err_den)
CLOSE(io_unit)

CALL xmhd_2d_plot(mhd_sim)
CALL oft_finalize
END PROGRAM hartmann
