!------------------------------------------------------------------------------
!> Note from Sophia: Still having many issues with flow. Run with flow_on = .FALSE.!
!> 
!> Travelling-wave induction pump on a cylindrical duct.
!!
!! The fluid fills a round pipe on the symmetry axis.  A stack of `n_seg` coil
!! segments outside the duct is driven with a phase step between neighbours, so
!! the applied field is a wave travelling in +Z; the slip between that wave and
!! the fluid drives an axial force.  An iron yoke returns the flux.  Build the
!! matching mesh with make_pump_mesh.py, and update any region indices as needed in oft.in
!!
!! BOUNDARY CONDITIONS (fixed, not options):
!!   inlet  (-Z)  velocity Dirichlet, Hagen-Poiseuille profile of mean v_inlet
!!                pressure free, so the developed head can be read off it
!!   outlet (+Z)  velocity free (Neumann)
!!                pressure Dirichlet p = 0 across the whole face
!!   duct wall    no-slip
!!   axis  (R=0)  installed by apply_axis_bc: v_R, v_phi, psi and
!!                B_phi pinned, v_Z left FREE (it peaks on the axis)
!!
!! The coils are driven as a prescribed current density (mugtok_td's
!! `source_term`) rather than as TokaMaker circuit coils.  That is required by
!! the iron: TokaMaker pre-solves a per-coil vacuum field at mu = mu0 and merely
!! scales it by current, which a field-dependent permeability invalidates.
!------------------------------------------------------------------------------
PROGRAM em_pump_driver
USE oft_base
USE multigrid, ONLY: multigrid_mesh
USE multigrid_build, ONLY: multigrid_construct_surf
USE mhd_utils, ONLY: mu0
USE oft_mesh_type, ONLY: oft_bmesh
USE oft_gs, ONLY: gs_equil, gs_factory, gs_update_bounds, gs_setup_walls, compute_bcmat
USE oft_gs_util, ONLY: gs_profile_alloc
USE oft_lag_basis, ONLY: oft_lag_setup, oft_scalar_bfem, oft_2D_lagrange_cast
USE fem_base, ONLY: oft_ml_fem_type
USE fem_utils, ONLY: bfem_map_flag
USE mugtok_td, ONLY: oft_mugtok_td

IMPLICIT NONE
TYPE(multigrid_mesh) :: mg_mesh
TYPE(oft_ml_fem_type), TARGET :: ML_blagrange
CLASS(oft_scalar_bfem), POINTER :: blagrange
TYPE(gs_factory), TARGET :: machine
TYPE(gs_equil), TARGET :: equil
TYPE(oft_mugtok_td) :: pump_sim
CLASS(oft_bmesh), POINTER :: smesh

!==============================================================================
! Runtime options (override with namelist `em_pump_options` in oft.in)
!==============================================================================
INTEGER(i4) :: order     = 2       !< FE order
INTEGER(i4) :: n_seg     = 24      !< Coil segments (must match the mesh)
INTEGER(i4) :: n_per_lam = 12      !< Segments per wavelength; sets the phase step
INTEGER(i4) :: nsteps    = 2000    !< Timesteps to take
INTEGER(i4) :: rst_freq  = 50      !< Restart output frequency [steps]
REAL(r8) :: dt        = 5.d-4      !< Timestep [s]
REAL(r8) :: freq      = 50.d0      !< Drive frequency [Hz]
REAL(r8) :: I0        = 1.d4       !< Peak current per coil segment [A]
REAL(r8) :: ramp_time = 1.d-1      !< Linear ramp of the drive envelope [s]
REAL(r8) :: rho_fluid = 1.94d3     !< Density [kg/m^3] (FLiBe ~600 C)
REAL(r8) :: eta_fluid = 1.625d-3   !< Fluid resistivity [Ohm-m] (FLiBe)
REAL(r8) :: eta_wall  = 1.2d-6     !< Duct wall resistivity [Ohm-m] (steel)
REAL(r8) :: mf_b0     = 1.d-4      !< Assumed relative error in the residual, sets the matrix-free Jacobian FD step
REAL(r8) :: nu_fluid  = 3.d-5      !< Kinematic viscosity [m^2/s]
REAL(r8) :: v_inlet   = 1.d0       !< Mean axial velocity imposed at the inlet [m/s]
REAL(r8) :: lin_tol   = 1.d-7      !< Linear solver tolerance
REAL(r8) :: nl_tol    = 1.d-6      !< Nonlinear solver tolerance
REAL(r8) :: mu_relax  = 0.5d0      !< Relaxation weight for the cached iron permeability
INTEGER(i4) :: wall_region = 4     !< Mesh region of the duct wall
INTEGER(i4) :: iron_region = 3     !< Mesh region of the iron yoke
INTEGER(i4) :: coil_reg0   = 5     !< Mesh region of coil segment 0; segments are contiguous
INTEGER(i4) :: dt_recover  = 20    !< Clean steps before dt is doubled back toward its namelist value (0 disables)
LOGICAL :: flow_on = .TRUE.        !< Solve the flow. .FALSE. freezes velocity and pressure at zero everywhere, leaving only the electromagnetics: the fluid still carries eddy currents, so this computes the induced currents in every component with no flow physics in the loop
LOGICAL :: mu_jac_deriv = .FALSE.   !< Assemble the d(mu)/d|B| term into the approximate Jacobian. Preconditioner-only, so converged answers must not move; set F to recover the old behaviour
LOGICAL :: mu_knee_damp = .FALSE.   !< Damp the mu relaxation target where the B-H curve is steep, so saturated cells fall back toward mu=1 instead of tracking mu(|B|). Only meaningful with mu_relax > 0
LOGICAL :: solver_pm = .FALSE.      !< Print the nonlinear and linear solver convergence history (diagnostic; very verbose)
NAMELIST/em_pump_options/order,n_seg,n_per_lam,nsteps,rst_freq,dt,freq,I0,ramp_time, &
  rho_fluid,eta_fluid,eta_wall,nu_fluid,v_inlet,lin_tol,nl_tol, &
  mu_relax,wall_region,iron_region,coil_reg0,dt_recover,flow_on,mf_b0,solver_pm, &
  mu_knee_damp,mu_jac_deriv

!==============================================================================
! Working storage
!==============================================================================
INTEGER(i4) :: i,j,k,io_unit,ierr,nreg,nl_its,l_its,nretry,nclean = 0
REAL(r8) :: t = 0.d0, dt_target, t_next
REAL(r8) :: lambda,v_sync,pipe_rout,chan_zmax,pump_len,seg_pitch
REAL(r8), ALLOCATABLE, DIMENSION(:) :: coil_zc
REAL(r8), ALLOCATABLE, DIMENSION(:) :: dens_reg,visc_reg,coil_step,coil_area
REAL(r8), ALLOCATABLE, DIMENSION(:,:) :: eta_reg
LOGICAL, ALLOCATABLE, DIMENSION(:) :: mhd_flag
INTEGER(i4), ALLOCATABLE, DIMENSION(:) :: inlet_pdofs
REAL(r8), PARAMETER :: ztol = 1.d-6, rtol = 1.d-6

!==============================================================================
! Read options
!==============================================================================
CALL oft_init
OPEN(NEWUNIT=io_unit,FILE=oft_env%ifile)
READ(io_unit,em_pump_options,IOSTAT=ierr)
CLOSE(io_unit)

IF(ierr /= 0)THEN
  WRITE(*,'(A,I0)')'Failed to read &em_pump_options from '//TRIM(oft_env%ifile)//', IOSTAT = ',ierr
  CALL oft_abort('Could not read &em_pump_options (stale variable name?)','em_pump_driver',__FILE__)
END IF
IF(mu_relax < 0.d0 .OR. mu_relax > 1.d0) &
  CALL oft_abort('"mu_relax" must lie in [0,1]','em_pump_driver',__FILE__)
dt_target = dt   ! step() may shrink dt; this is the value to climb back to

!==============================================================================
! Mesh and FE space
!==============================================================================
CALL multigrid_construct_surf(mg_mesh)
CALL oft_lag_setup(mg_mesh,order,ML_blag_obj=ML_blagrange,minlev=-1)
IF(.NOT.oft_2D_lagrange_cast(blagrange,ML_blagrange%current_level)) &
  CALL oft_abort("Invalid lagrange FE object","em_pump_driver",__FILE__)
smesh => mg_mesh%smesh
nreg = smesh%nreg
IF(coil_reg0 < 2) &
  CALL oft_abort("coil_reg0 must be at least 2 (region 1 is the fluid)","em_pump_driver",__FILE__)
IF(coil_reg0+n_seg-1 > nreg) &
  CALL oft_abort("Coil regions run past the end of the mesh; check n_seg and coil_reg0", &
                 "em_pump_driver",__FILE__)

!---Pipe radius and half-length, computed from the mesh rather than declared.
!   Both are pure consequences of the geometry, and a mismatch is silent: the
!   end-face search below finds nothing, so the duct comes out sealed and the
!   prescribed inflow never enters.
pipe_rout = 0.d0
chan_zmax = 0.d0
DO i=1,smesh%nc
  IF(smesh%reg(i) /= 1)CYCLE
  DO j=1,smesh%cell_np
    pipe_rout = MAX(pipe_rout,     smesh%r(1,smesh%lc(j,i)))
    chan_zmax = MAX(chan_zmax,ABS(smesh%r(2,smesh%lc(j,i))))
  END DO
END DO

!==============================================================================
! TokaMaker device
!
! MUG is handling all the coils, so every coil array is zero-size and every coil loop runs empty
! The duct wall is a solid conductor
!==============================================================================
CALL machine%setup(ML_blagrange)
machine%save_visit = .FALSE.
machine%ncoil_regs = 0
machine%ncoils = 0
ALLOCATE(machine%coil_regions(0))
machine%ncond_regs = 1
ALLOCATE(machine%cond_regions(1))
machine%cond_regions(1)%id = wall_region
machine%cond_regions(1)%eta = eta_wall/mu0
CALL gs_setup_walls(machine)
ALLOCATE(machine%coil_nturns(nreg,0))

!---Give the yoke nonlinear-permeability. mag_suscep defaults to
!   "not magnetic" (0) everywhere. The value set here is only a flag, since
!   mu itself comes from the B-H table inside gs_update_mu.
machine%mag_suscep(iron_region) = 1.d0

!---Suppress the O-/X-point search, which is not necessary without a plasma
!---and is computationally expensive
ALLOCATE(machine%saddle_rmask(nreg))
machine%saddle_rmask = .TRUE.
machine%free = .TRUE.
CALL machine%init()
CALL compute_bcmat(machine)
ALLOCATE(machine%coil_vcont(0))

!---Coil region areas, needed to compute current density
ALLOCATE(coil_area(n_seg),coil_zc(n_seg))
coil_area = 0.d0
coil_zc = 0.d0
BLOCK
REAL(r8) :: f3(3),gop3(3,4),vol3,pt3(3)
INTEGER(i4) :: jreg
f3 = 1.d0/3.d0
DO j=1,smesh%nc
  jreg = smesh%reg(j)
  IF(jreg >= coil_reg0 .AND. jreg <= coil_reg0+n_seg-1)THEN
    CALL smesh%jacobian(j,f3,gop3,vol3)
    pt3 = smesh%log2phys(j,f3)
    coil_area(jreg-coil_reg0+1) = coil_area(jreg-coil_reg0+1) + vol3
    coil_zc(jreg-coil_reg0+1)   = coil_zc(jreg-coil_reg0+1)   + vol3*pt3(2)
  END IF
END DO
END BLOCK
IF(ANY(coil_area <= 0.d0)) &
  CALL oft_abort("A coil region has zero area; check n_seg and coil_reg0","em_pump_driver",__FILE__)

!---Coil stack length, computed from the mesh 
coil_zc = coil_zc/coil_area
seg_pitch = (coil_zc(n_seg)-coil_zc(1))/REAL(n_seg-1,8)
IF(seg_pitch <= 0.d0) &
  CALL oft_abort("Coil regions are not ordered by increasing Z; check coil_reg0 and n_seg", &
                 "em_pump_driver",__FILE__)
DO j=1,n_seg-1
  IF(ABS((coil_zc(j+1)-coil_zc(j))/seg_pitch - 1.d0) > 0.05d0) &
    CALL oft_abort("Coil segment spacing is not uniform; the mesh does not match n_seg", &
                   "em_pump_driver",__FILE__)
END DO
pump_len = REAL(n_seg,8)*seg_pitch

!==============================================================================
! Plasma-free equilibrium
!
! Everything starts at rest with zero flux and zero current, matching the ramped
! drive. The FF' and P' profiles are zeroed because there is no plasma
!==============================================================================
CALL equil%new(machine)
CALL equil%psi%set(0.d0)
equil%has_plasma = .FALSE.
equil%Ip_target = 0.d0
equil%p_scale = 0.d0
equil%ffp_scale = 0.d0
equil%mode = 0
equil%vcontrol_val = 0.d0
CALL gs_profile_alloc('zero',equil%I)
CALL gs_profile_alloc('zero',equil%P)
CALL gs_update_bounds(equil)

!==============================================================================
! Set MHD region fluid properties and wall resistivity, and set up coupled simulation
!==============================================================================
ALLOCATE(mhd_flag(nreg),dens_reg(nreg),visc_reg(nreg),eta_reg(nreg,2))
mhd_flag = .FALSE.  ; mhd_flag(1) = .TRUE.
dens_reg = -1.d0    ; dens_reg(1) = rho_fluid
visc_reg = -1.d0    ; visc_reg(1) = nu_fluid
eta_reg  = -1.d0    ; eta_reg(1,:) = eta_fluid/mu0
eta_reg(wall_region,:) = eta_wall/mu0   ! not MHD, but add_f_terms wants a real eta

lambda = REAL(n_per_lam,8)*seg_pitch
v_sync = lambda*freq
WRITE(*,'(A)')         '=== EM pump drive ==='
WRITE(*,'(A,ES11.3,A)')'  pipe radius        = ',pipe_rout,' m'
WRITE(*,'(A,ES11.3,A)')'  duct length        = ',2.d0*chan_zmax,' m'
WRITE(*,'(A,ES11.3,A)')'  pump length        = ',pump_len,' m (from mesh)'
WRITE(*,'(A,ES11.3,A)')'  segment pitch      = ',seg_pitch,' m (from mesh)'
WRITE(*,'(A,ES11.3,A)')'  wavelength         = ',lambda,' m'
WRITE(*,'(A,ES11.3,A)')'  drive frequency    = ',freq,' Hz'
WRITE(*,'(A,ES11.3,A)')'  synchronous speed  = ',v_sync,' m/s'

pump_sim%mu_relax = mu_relax   ! step() refreshes the cached iron permeability
pump_sim%mu_jac_deriv = mu_jac_deriv
pump_sim%mu_knee_damp = mu_knee_damp
pump_sim%pm = solver_pm        ! setup() copies this into oft_env%pm, which gates the solver histories
pump_sim%save_rst = .TRUE.
pump_sim%rst_freq = rst_freq
CALL pump_sim%setup(equil,dt,lin_tol,nl_tol,mhd_flag,dens_reg,visc_reg, &
                    eta_reg=eta_reg,incomp=.TRUE.,toroidal_flow=.FALSE.,plasma_reg=0)
pump_sim%mfmat%b0 = mf_b0   ! setup() hard-codes 1e-4; override after it runs

!==============================================================================
! Boundary conditions
!
! setup() pins velocity at every DOF touching a non-MHD cell, which seals both
! end faces, and pins one interior pressure DOF to remove the null space of the
! enclosed problem. Both are edited below. The end faces are interior
! region-1/region-2 interfaces rather than domain boundaries, so mesh%bes does
! not see them and they have to be found geometrically.
!
!==============================================================================
BLOCK
LOGICAL, ALLOCATABLE, DIMENSION(:) :: vert_flag,edge_flag,in_flag,outlet_flag
INTEGER(i4), POINTER, DIMENSION(:) :: pdofs
REAL(r8), POINTER, DIMENSION(:) :: vals
REAL(r8), ALLOCATABLE, DIMENSION(:) :: dof_r
INTEGER(i4) :: ic,je,ed,p1,p2,nin,nout,nface
IF(blagrange%ne /= SIZE(pump_sim%mug%velz_bc)) &
  CALL oft_abort("Velocity BC array does not match the FE representation","em_pump_driver",__FILE__)
ALLOCATE(vert_flag(smesh%np),edge_flag(smesh%ne))
ALLOCATE(outlet_flag(blagrange%ne),in_flag(blagrange%ne))

!---Radius of every velocity DOF. Order-2 Lagrange lays out the np vertices
!   first, then the edge midpoints, so an edge DOF sits at the mean of its ends.
ALLOCATE(dof_r(blagrange%ne))
IF(blagrange%ne /= smesh%np + smesh%ne) &
  CALL oft_abort("Unexpected velocity DOF layout","em_pump_driver",__FILE__)
DO k=1,smesh%np
  dof_r(k) = smesh%r(1,k)
END DO
DO k=1,smesh%ne
  dof_r(smesh%np+k) = 0.5d0*(smesh%r(1,smesh%le(1,k)) + smesh%r(1,smesh%le(2,k)))
END DO

!------------------------------------------------------------------------------
! Outlet (+Z): release the velocity, leaving the Neumann outflow that is
! the natural condition of the weak form.
!------------------------------------------------------------------------------
vert_flag = .FALSE.; edge_flag = .FALSE.
DO ic=1,smesh%nc
  IF(smesh%reg(ic) /= 1)CYCLE
  DO je=1,smesh%cell_ne
    ed = ABS(smesh%lce(je,ic)); p1 = smesh%le(1,ed); p2 = smesh%le(2,ed)
    IF(ABS(smesh%r(2,p1)-chan_zmax) < ztol .AND. ABS(smesh%r(2,p2)-chan_zmax) < ztol)THEN
      edge_flag(ed) = .TRUE.
      IF(ABS(smesh%r(1,p1)-pipe_rout) > rtol) vert_flag(p1) = .TRUE.
      IF(ABS(smesh%r(1,p2)-pipe_rout) > rtol) vert_flag(p2) = .TRUE.
    END IF
  END DO
END DO
CALL bfem_map_flag(blagrange,vert_flag,edge_flag,outlet_flag)
pump_sim%mug%velx_bc = pump_sim%mug%velx_bc .AND. .NOT.outlet_flag
pump_sim%mug%velz_bc = pump_sim%mug%velz_bc .AND. .NOT.outlet_flag
nout = COUNT(edge_flag)

!------------------------------------------------------------------------------
! Inlet (-Z): keep the Dirichlet flags setup() installed and write the inlet profile
! value into u. Hagen-Poiseuille, v_Z = 2*v_inlet*(1-(R/a)^2): peak on the axis,
! zero at the wall, area-weighted mean exactly v_inlet
!------------------------------------------------------------------------------
vert_flag = .FALSE.; edge_flag = .FALSE.
DO ic=1,smesh%nc
  IF(smesh%reg(ic) /= 1)CYCLE
  DO je=1,smesh%cell_ne
    ed = ABS(smesh%lce(je,ic)); p1 = smesh%le(1,ed); p2 = smesh%le(2,ed)
    IF(ABS(smesh%r(2,p1)+chan_zmax) < ztol .AND. ABS(smesh%r(2,p2)+chan_zmax) < ztol)THEN
      edge_flag(ed) = .TRUE.
      vert_flag(p1) = .TRUE.
      vert_flag(p2) = .TRUE.
    END IF
  END DO
END DO
CALL bfem_map_flag(blagrange,vert_flag,edge_flag,in_flag)
nin = COUNT(edge_flag)

!---Set up the whole pipe with the correct IC, so the solve starts from a
!   divergence-free state
NULLIFY(vals)
CALL pump_sim%u%get_local(vals,4)   ! field 4 = axial velocity
WHERE(.NOT.pump_sim%mug%velz_bc) vals = 2.d0*v_inlet*(1.d0-(dof_r/pipe_rout)**2)
WHERE(in_flag)                   vals = 2.d0*v_inlet*(1.d0-(dof_r/pipe_rout)**2)
CALL pump_sim%u%restore_local(vals,4)
DEALLOCATE(vals)
NULLIFY(vals)
CALL pump_sim%u%get_local(vals,2)   ! field 2 = radial velocity
WHERE(in_flag) vals = 0.d0
CALL pump_sim%u%restore_local(vals,2)
DEALLOCATE(vals)

!------------------------------------------------------------------------------
! Pressure: free at the inlet, p = 0 across the outlet face.
!
! Prescribing p over a whole face already removes the null space, so setup()'s
! interior pin must be cleared rather than kept -- holding both over-determines
! the pressure
!------------------------------------------------------------------------------
ALLOCATE(pdofs(pump_sim%mug%fe_rep%fields(5)%fe%nce))
DO ic=1,smesh%nc
  IF(smesh%reg(ic) /= 1)CYCLE
  CALL pump_sim%mug%fe_rep%fields(5)%fe%ncdofs(ic,pdofs)
  DO k=1,SIZE(pdofs)
    pump_sim%mug%T_bc(pdofs(k)) = .FALSE.
  END DO
END DO
DEALLOCATE(pdofs)
!---Pressure is order-1, so its DOFs sit on mesh vertices and the DOF index maps
!   straight to the vertex index.
NULLIFY(vals)
CALL pump_sim%u%get_local(vals,5)
nface = 0
DO k=1,MIN(SIZE(pump_sim%mug%T_bc),smesh%np)
  IF(ABS(smesh%r(2,k)-chan_zmax) > ztol)CYCLE
  IF(smesh%r(1,k) > pipe_rout+rtol)CYCLE
  pump_sim%mug%T_bc(k) = .TRUE.
  vals(k) = 0.d0
  nface = nface + 1
END DO
CALL pump_sim%u%restore_local(vals,5)
DEALLOCATE(vals)
IF(nface == 0) CALL oft_abort("No outlet-face pressure DOF found; check chan_zmax", &
                              "em_pump_driver",__FILE__)

!---Inlet pressure DOFs for the head diagnostic
ALLOCATE(inlet_pdofs(COUNT( &
  ABS(smesh%r(2,1:MIN(SIZE(pump_sim%mug%T_bc),smesh%np))+chan_zmax) <= ztol .AND. &
  smesh%r(1,1:MIN(SIZE(pump_sim%mug%T_bc),smesh%np)) <= pipe_rout+rtol)))
k = 0
DO ic=1,MIN(SIZE(pump_sim%mug%T_bc),smesh%np)
  IF(ABS(smesh%r(2,ic)+chan_zmax) > ztol)CYCLE
  IF(smesh%r(1,ic) > pipe_rout+rtol)CYCLE
  k = k + 1
  inlet_pdofs(k) = ic
END DO

IF(nin == 0 .OR. nout == 0) CALL oft_warn("End-face edges not found; check chan_zmax")
DEALLOCATE(vert_flag,edge_flag,in_flag,outlet_flag,dof_r)
END BLOCK

!==============================================================================
! EM-only mode: freeze the flow
!
! A pinned DOF's residual row reads u_new = u_old, so it holds whatever value it
! carries here for the whole run. Flagging every velocity component and zeroing
! the field therefore gives v == 0 exactly, at every step, without touching the
! momentum equation: advection, the JxB acceleration and viscous transport all
! still assemble, but their rows are overwritten.
!
! The pressure has to be pinned too, because it's equation is meaningless without velocity
!
! This must come after the boundary-condition block, which it deliberately
! overrides.
!==============================================================================
IF(.NOT.flow_on)THEN
BLOCK
REAL(r8), POINTER, DIMENSION(:) :: fld_vals
INTEGER(i4) :: ifld
pump_sim%mug%velx_bc = .TRUE.
pump_sim%mug%vely_bc = .TRUE.
pump_sim%mug%velz_bc = .TRUE.
pump_sim%mug%T_bc    = .TRUE.
DO ifld=2,5   ! 2 = v_R, 3 = v_phi, 4 = v_Z, 5 = pressure
  NULLIFY(fld_vals)
  CALL pump_sim%u%get_local(fld_vals,ifld)
  fld_vals = 0.d0
  CALL pump_sim%u%restore_local(fld_vals,ifld)
  DEALLOCATE(fld_vals)
END DO
WRITE(*,'(A)')'=== EM-only mode (flow_on = F): velocity and pressure pinned at zero ==='
END BLOCK
END IF

!==============================================================================
! Time advance
!==============================================================================
ALLOCATE(coil_step(n_seg))
DO i=1,nsteps
  !---step() prescribes the coil currents at the END of the step
  t_next = t + dt
  CALL pump_waveform(t_next,coil_step)
  pump_sim%source_term = 0.d0
  DO j=1,n_seg
    pump_sim%source_term(coil_reg0+j-1) = coil_step(j)/coil_area(j)
  END DO
  CALL pump_sim%step(coil_step(1:0),t,dt,nl_its,l_its,nretry)


  CALL report(i,t,dt,nl_its,l_its,nretry)

  IF(nretry == 0)THEN
    nclean = nclean + 1
  ELSE
    nclean = 0
  END IF
  !if dt_recover, reset dt to the target value
  IF(dt_recover > 0 .AND. nclean >= dt_recover .AND. dt < dt_target)THEN
    dt = MIN(dt*2.d0,dt_target)
    nclean = 0
    WRITE(*,'(A,ES10.3)')'   dt recovered to ',dt
  END IF

END DO

CALL pump_sim%plot()
CALL oft_finalize

CONTAINS
!------------------------------------------------------------------------------
!> Travelling-wave coil currents at a given time.
!!
!! Segment k is driven with a phase lag k*2*pi/n_per_lam, so the applied field is
!! a wave travelling in +Z. Currents come back in TokaMaker's normalized units
!! (Amps*mu0). The envelope ramps from zero so t=0 -- no field, no flow -- is a
!! consistent initial condition rather than a step transient.
!------------------------------------------------------------------------------
SUBROUTINE pump_waveform(time,currs)
REAL(r8), INTENT(in) :: time      !< Time [s]
REAL(r8), INTENT(out) :: currs(:) !< Coil currents [A*mu0]
INTEGER(i4) :: jc
REAL(r8) :: omega,dphi,env
omega = 2.d0*pi*freq
dphi = 2.d0*pi/REAL(n_per_lam,8)
env = 1.d0
IF(ramp_time > 0.d0) env = MIN(1.d0,time/ramp_time)
DO jc=1,SIZE(currs)
  currs(jc) = env*I0*COS(omega*time - REAL(jc-1,8)*dphi)*mu0
END DO
END SUBROUTINE pump_waveform
!------------------------------------------------------------------------------
!> One line per step: solver health, peak velocity, and the developed head
!>  (meaningless with no flow, ignore)
!------------------------------------------------------------------------------
SUBROUTINE report(istep,time,dt_now,nl,nlin,nret)
INTEGER(i4), INTENT(in) :: istep,nl,nlin,nret
REAL(r8), INTENT(in) :: time,dt_now
REAL(r8), POINTER, DIMENSION(:) :: vz,pv
REAL(r8) :: vmax,p_in
INTEGER(i4) :: kk
!---With the flow frozen both of these are identically zero, so they are noise
IF(.NOT.flow_on)THEN
  WRITE(*,'(A,I6,A,ES12.5,A,ES10.3,A,I3,A,I5,A,I2)') &
    'Step ',istep,'  t = ',time,'  dt = ',dt_now,'  NL = ',nl,'  lin = ',nlin, &
    '  retry = ',nret
  RETURN
END IF
NULLIFY(vz,pv)
CALL pump_sim%u%get_local(vz,4)
vmax = 0.d0
DO kk=1,SIZE(vz)
  IF(.NOT.pump_sim%mug%velz_bc(kk)) vmax = MAX(vmax,ABS(vz(kk)))
END DO
DEALLOCATE(vz)
CALL pump_sim%u%get_local(pv,5)
p_in = 0.d0
DO kk=1,SIZE(inlet_pdofs)
  p_in = p_in + pv(inlet_pdofs(kk))
END DO
IF(SIZE(inlet_pdofs) > 0) p_in = p_in/REAL(SIZE(inlet_pdofs),8)
DEALLOCATE(pv)
p_in = -p_in   ! developed head: positive means the pump is raising pressure
WRITE(*,'(A,I6,A,ES12.5,A,ES10.3,A,I3,A,I5,A,I2,A,ES11.4,A,ES13.5)') &
  'Step ',istep,'  t = ',time,'  dt = ',dt_now,'  NL = ',nl,'  lin = ',nlin, &
  '  retry = ',nret,'  max|v_z| = ',vmax,'  head = ',p_in
END SUBROUTINE report
END PROGRAM em_pump_driver
