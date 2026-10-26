!------------------------------------------------------------------------------
!> Travelling-wave induction pump driver via empump_td.F90, in a cylindrical or slab duct
!! Still having many issues with flow, so run as pure EM simulation for now (flow_on = F)
!!
!! Mesh and inputs are read from oft.in
!! Run as ./em_pump oft.in
!------------------------------------------------------------------------------
PROGRAM em_pump
USE oft_base
USE multigrid, ONLY: multigrid_mesh
USE multigrid_build, ONLY: multigrid_construct_surf
USE mhd_utils, ONLY: mu0
USE oft_mesh_type, ONLY: oft_bmesh
USE fem_utils, ONLY: bfem_map_flag
USE xmhd_2d, ONLY: oft_blagrange
USE empump_td, ONLY: oft_empump_td
IMPLICIT NONE
TYPE(multigrid_mesh) :: mg_mesh
TYPE(oft_empump_td) :: pump_sim
CLASS(oft_bmesh), POINTER :: smesh

!==============================================================================
! Runtime options (defaults overwritten by namelist `em_pump_options` in oft.in)
!==============================================================================
INTEGER(i4) :: order     = 2       !< FE order
INTEGER(i4) :: n_seg     = 24      !< Coil segments (must match the mesh)
INTEGER(i4) :: n_per_lam = 12      !< Segments per wavelength; sets the phase step
INTEGER(i4) :: coil_reg0 = 5       !< Mesh region of coil segment 0; segments are contiguous
REAL(r8) :: freq      = 50.d0      !< Drive frequency [Hz]
REAL(r8) :: I0        = 1.d4       !< Peak current per coil segment [A]
REAL(r8) :: ramp_time = 1.d-1      !< Linear ramp of the drive envelope [s]
INTEGER(i4) :: nsteps    = 2000    !< Timesteps to take
INTEGER(i4) :: rst_freq  = 50      !< Restart output frequency [steps]
REAL(r8) :: dt        = 5.d-4      !< Timestep [s]
REAL(r8) :: rho_fluid = 1.94d3     !< Density [kg/m^3] (FLiBe ~600 C)
REAL(r8) :: eta_fluid = 1.625d-3   !< Fluid resistivity [Ohm-m] (FLiBe)
REAL(r8) :: eta_wall  = 1.2d-6     !< Duct wall resistivity [Ohm-m] (steel)
REAL(r8) :: eta_vac   = 1.d7       !< Magnetic diffusivity of air, coils and iron [m^2/s]; vacuum-like when >> L^2*omega
REAL(r8) :: mf_b0     = 1.d-4      !< Matrix-free Jacobian FD step parameter
REAL(r8) :: nu_fluid  = 3.d-5      !< Kinematic viscosity [m^2/s]
REAL(r8) :: v_inlet   = 1.d0       !< Mean axial velocity imposed at the inlet [m/s]
REAL(r8) :: lin_tol   = 1.d-7      !< Linear solver tolerance
REAL(r8) :: nl_tol    = 1.d-6      !< Nonlinear solver tolerance
INTEGER(i4) :: wall_region = 4     !< Mesh region of the duct wall
INTEGER(i4) :: iron_region = 3     !< Mesh region of the iron yoke (0 = no iron)
LOGICAL :: cyl = .TRUE.            !< Cylindrical (R,Z) geometry; .FALSE. for a slab (x,z)
LOGICAL :: iron_jac_lift = .FALSE. !< Include the (nu - 1) iron lift in the preconditioner (empump_td%iron_jac_lift); off by default, it tripled the GMRES iterations on this mesh
LOGICAL :: iron_jac_deriv = .FALSE. !< Include the d(mu)/d|B| term in the preconditioner (empump_td%iron_jac_deriv); off by default, it caused failed Newton solves on this mesh
LOGICAL :: flow_on = .FALSE.        !< Solve the flow. .FALSE. pins velocity and pressure at zero everywhere
LOGICAL :: solver_pm = .FALSE.     !< Print solver convergence histories
LOGICAL :: psi_outer_dirichlet = .TRUE. !< psi = 0 on the outer mesh boundary (perfect conductor); .FALSE. makesdpsi/dn = 0 (perfect insulator)
NAMELIST/em_pump_options/order,n_seg,n_per_lam,nsteps,rst_freq,dt,freq,I0,ramp_time, &
  rho_fluid,eta_fluid,eta_wall,eta_vac,nu_fluid,v_inlet,lin_tol,nl_tol,wall_region,iron_region, &
  coil_reg0,cyl,iron_jac_lift,iron_jac_deriv,flow_on,mf_b0,solver_pm,psi_outer_dirichlet

!==============================================================================
! Other variables we need
!==============================================================================
INTEGER(i4) :: i,j,k,io_unit,ierr,nreg,nl_its,l_its,nretry
REAL(r8) :: t = 0.d0
REAL(r8) :: lambda,v_sync,chan_zmax,pump_len,seg_pitch,chan_xmin,chan_xmax
LOGICAL :: solid_pipe
REAL(r8), ALLOCATABLE, DIMENSION(:) :: coil_area,coil_zc,dens_reg,visc_reg
REAL(r8), ALLOCATABLE, DIMENSION(:,:) :: eta_reg
LOGICAL, ALLOCATABLE, DIMENSION(:) :: fluid_reg,iron_reg
INTEGER(i4), ALLOCATABLE, DIMENSION(:) :: inlet_pdofs
REAL(r8), PARAMETER :: ztol = 1.d-6, rtol = 1.d-6

!==============================================================================
! Read options from oft.in
!==============================================================================
CALL oft_init
OPEN(NEWUNIT=io_unit,FILE=oft_env%ifile)
READ(io_unit,em_pump_options,IOSTAT=ierr)
CLOSE(io_unit)
IF(ierr /= 0)THEN
  WRITE(*,'(A,I0)')'Failed to read &em_pump_options from '//TRIM(oft_env%ifile)//', IOSTAT = ',ierr
  CALL oft_abort('Could not read &em_pump_options (stale variable name?)','em_pump',__FILE__)
END IF

!==============================================================================
! Set up mesh and determine duct geometry 
!==============================================================================
CALL multigrid_construct_surf(mg_mesh)
smesh => mg_mesh%smesh
nreg = smesh%nreg
IF(coil_reg0 < 2) &
  CALL oft_abort("coil_reg0 must be at least 2 (region 1 is the fluid)","em_pump",__FILE__)
IF(coil_reg0+n_seg-1 > nreg) &
  CALL oft_abort("Coil regions run past the end of the mesh region list; check n_seg and coil_reg0", &
                 "em_pump",__FILE__)

!---Compute channel geometry from min/max x values and max z values. 
chan_xmin = HUGE(1.d0)
chan_xmax = -HUGE(1.d0)
chan_zmax = 0.d0
DO i=1,smesh%nc
  IF(smesh%reg(i) /= 1)CYCLE
  DO j=1,smesh%cell_np
    chan_xmin = MIN(chan_xmin,smesh%r(1,smesh%lc(j,i)))
    chan_xmax = MAX(chan_xmax,smesh%r(1,smesh%lc(j,i)))
    chan_zmax = MAX(chan_zmax,ABS(smesh%r(2,smesh%lc(j,i))))
  END DO
END DO
solid_pipe = cyl .AND. chan_xmin < 1.d-8

!---Compute coil areas and centroids to apply appropriate current density
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
  CALL oft_abort("A coil region has zero area; check n_seg and coil_reg0","em_pump",__FILE__)
coil_zc = coil_zc/coil_area
seg_pitch = (coil_zc(n_seg)-coil_zc(1))/REAL(n_seg-1,8)
IF(seg_pitch <= 0.d0) &
  CALL oft_abort("Coil regions are not ordered by increasing Z; check coil_reg0 and n_seg", &
                 "em_pump",__FILE__)
DO j=1,n_seg-1
  IF(ABS((coil_zc(j+1)-coil_zc(j))/seg_pitch - 1.d0) > 0.05d0) &
    CALL oft_abort("Coil segment spacing is not uniform; the mesh does not match n_seg", &
                   "em_pump",__FILE__)
END DO
pump_len = REAL(n_seg,8)*seg_pitch

!==============================================================================
! Set up region properties
!==============================================================================
ALLOCATE(fluid_reg(nreg),iron_reg(nreg),dens_reg(nreg),visc_reg(nreg),eta_reg(nreg,2))
fluid_reg = .FALSE. ; fluid_reg(1) = .TRUE.
dens_reg = 1.d0     ; dens_reg(1) = rho_fluid
visc_reg = 0.d0     ; visc_reg(1) = nu_fluid
eta_reg = eta_vac
eta_reg(1,:) = eta_fluid/mu0
eta_reg(wall_region,:) = eta_wall/mu0
iron_reg = .FALSE.
IF(iron_region > 0) iron_reg(iron_region) = .TRUE.   ! nonlinear B-H permeability

lambda = REAL(n_per_lam,8)*seg_pitch
v_sync = lambda*freq
WRITE(*,'(A)')         '=== EM pump drive (xmhd_2d / empump_td) ==='
WRITE(*,'(A,L2)')      '  cylindrical        = ',cyl
IF(solid_pipe)THEN
  WRITE(*,'(A,ES11.3,A)')'  pipe radius        = ',chan_xmax,' m'
ELSE IF(cyl)THEN
  WRITE(*,'(A,2ES11.3,A)')'  annulus R in/out   = ',chan_xmin,chan_xmax,' m'
ELSE
  WRITE(*,'(A,2ES11.3,A)')'  channel x from/to  = ',chan_xmin,chan_xmax,' m'
END IF
WRITE(*,'(A,ES11.3,A)')'  duct length        = ',2.d0*chan_zmax,' m'
WRITE(*,'(A,ES11.3,A)')'  pump length        = ',pump_len,' m (from mesh)'
WRITE(*,'(A,ES11.3,A)')'  segment pitch      = ',seg_pitch,' m (from mesh)'
WRITE(*,'(A,ES11.3,A)')'  wavelength         = ',lambda,' m'
WRITE(*,'(A,ES11.3,A)')'  drive frequency    = ',freq,' Hz'
WRITE(*,'(A,ES11.3,A)')'  synchronous speed  = ',v_sync,' m/s'
WRITE(*,'(A,ES11.3,A)')'  eta_vac            = ',eta_vac,' m^2/s'

pump_sim%pm = solver_pm
pump_sim%iron_jac_deriv = iron_jac_deriv
pump_sim%iron_jac_lift = iron_jac_lift
pump_sim%save_rst = .TRUE.
pump_sim%rst_freq = rst_freq
CALL pump_sim%setup(mg_mesh,order,dt,lin_tol,nl_tol,fluid_reg,dens_reg,visc_reg,eta_reg, &
                    cyl,iron_reg=iron_reg,pin_gauge=.FALSE.,freeze_flow=.NOT.flow_on, &
                    psi_outer_dirichlet=psi_outer_dirichlet)
pump_sim%mfmat%b0 = mf_b0
IF(.NOT.flow_on) WRITE(*,'(A)')'=== EM-only mode (flow_on = F): velocity and pressure pinned at zero ==='

!==============================================================================
! Flow boundary conditions (flow_on only). empump_td has already pinned the
! velocity on every DOF of a non-fluid cell, which seals both end faces; the
! outlet is released here and the inlet profile written. The end faces are
! interior region interfaces, so they are found geometrically.
!==============================================================================
IF(flow_on)THEN
BLOCK
LOGICAL, ALLOCATABLE, DIMENSION(:) :: vert_flag,edge_flag,in_flag,outlet_flag
REAL(r8), POINTER, DIMENSION(:) :: vals
REAL(r8), ALLOCATABLE, DIMENSION(:) :: dof_r,v_prof
INTEGER(i4) :: ic,je,ed,p1,p2,nin,nout,nface
ALLOCATE(vert_flag(smesh%np),edge_flag(smesh%ne))
ALLOCATE(outlet_flag(oft_blagrange%ne),in_flag(oft_blagrange%ne))

!---Transverse coordinate (R or x) of every velocity DOF: order-2 Lagrange lays
!   out the vertices first, then the edge midpoints.
IF(oft_blagrange%ne /= smesh%np + smesh%ne) &
  CALL oft_abort("Unexpected velocity DOF layout (order 2 expected)","em_pump",__FILE__)
ALLOCATE(dof_r(oft_blagrange%ne),v_prof(oft_blagrange%ne))
DO k=1,smesh%np
  dof_r(k) = smesh%r(1,k)
END DO
DO k=1,smesh%ne
  dof_r(smesh%np+k) = 0.5d0*(smesh%r(1,smesh%le(1,k)) + smesh%r(1,smesh%le(2,k)))
END DO
!---Developed laminar profile with area-mean v_inlet, zero on every wall
IF(solid_pipe)THEN
  !---Hagen-Poiseuille
  v_prof = 2.d0*v_inlet*(1.d0-(dof_r/chan_xmax)**2)
ELSE IF(cyl)THEN
  !---Annular Poiseuille, u ~ (ro^2 - r^2) + (ro^2 - ri^2) ln(r/ro)/ln(ro/ri), divided
  !   by its area mean (ro^2 + ri^2)/2 - (ro^2 - ri^2)/(2 ln(ro/ri)). Tends to
  !   Hagen-Poiseuille as ri -> 0. MAX() only keeps the log finite at DOFs
  !   inside the inner wall, which are pinned and never use the value.
  BLOCK
  REAL(r8) :: ri,ro,lnr,umean
  ri = chan_xmin; ro = chan_xmax; lnr = LOG(ro/ri)
  umean = 0.5d0*(ro**2 + ri**2) - 0.5d0*(ro**2 - ri**2)/lnr
  v_prof = v_inlet*((ro**2 - dof_r**2) + (ro**2 - ri**2)*LOG(MAX(dof_r,ri)/ro)/lnr)/umean
  END BLOCK
ELSE
  !---Plane Poiseuille across the gap
  v_prof = 1.5d0*v_inlet*(1.d0-((dof_r - 0.5d0*(chan_xmin + chan_xmax)) &
                                /(0.5d0*(chan_xmax - chan_xmin)))**2)
END IF

!---Outlet (+Z): release the velocity, leaving the natural (Neumann) outflow.
!   The wall corners stay no-slip.
vert_flag = .FALSE.; edge_flag = .FALSE.
DO ic=1,smesh%nc
  IF(smesh%reg(ic) /= 1)CYCLE
  DO je=1,smesh%cell_ne
    ed = ABS(smesh%lce(je,ic)); p1 = smesh%le(1,ed); p2 = smesh%le(2,ed)
    IF(ABS(smesh%r(2,p1)-chan_zmax) < ztol .AND. ABS(smesh%r(2,p2)-chan_zmax) < ztol)THEN
      edge_flag(ed) = .TRUE.
      IF(.NOT.on_wall(smesh%r(1,p1))) vert_flag(p1) = .TRUE.
      IF(.NOT.on_wall(smesh%r(1,p2))) vert_flag(p2) = .TRUE.
    END IF
  END DO
END DO
CALL bfem_map_flag(oft_blagrange,vert_flag,edge_flag,outlet_flag)
pump_sim%mug%velx_bc = pump_sim%mug%velx_bc .AND. .NOT.outlet_flag
pump_sim%mug%velz_bc = pump_sim%mug%velz_bc .AND. .NOT.outlet_flag
IF(ASSOCIATED(pump_sim%mug%axis_bc)) &   ! v_R stays pinned on the axis
  pump_sim%mug%velx_bc = pump_sim%mug%velx_bc .OR. pump_sim%mug%axis_bc
nout = COUNT(edge_flag)

!---Inlet (-Z): keep the velocity pinned and write the profile into it
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
CALL bfem_map_flag(oft_blagrange,vert_flag,edge_flag,in_flag)
nin = COUNT(edge_flag)

!---Whole duct starts on the developed profile, so the first solve begins
!   divergence-free
NULLIFY(vals)
CALL pump_sim%u%get_local(vals,4)
WHERE(.NOT.pump_sim%mug%velz_bc) vals = v_prof
WHERE(in_flag)                   vals = v_prof
CALL pump_sim%u%restore_local(vals,4)
DEALLOCATE(vals)

!---Pressure p = 0 across the outlet face (order-1 DOFs = mesh vertices). This
!   removes the null space, so setup() was told not to pin a gauge DOF.
NULLIFY(vals)
CALL pump_sim%u%get_local(vals,5)
nface = 0
DO k=1,MIN(SIZE(pump_sim%mug%T_bc),smesh%np)
  IF(ABS(smesh%r(2,k)-chan_zmax) > ztol)CYCLE
  IF(.NOT.in_channel(smesh%r(1,k)))CYCLE
  pump_sim%mug%T_bc(k) = .TRUE.
  vals(k) = 0.d0
  nface = nface + 1
END DO
CALL pump_sim%u%restore_local(vals,5)
DEALLOCATE(vals)
IF(nface == 0) CALL oft_abort("No outlet-face pressure DOF found","em_pump",__FILE__)

!---Inlet pressure DOFs for the head diagnostic
k = 0
DO ic=1,MIN(SIZE(pump_sim%mug%T_bc),smesh%np)
  IF(ABS(smesh%r(2,ic)+chan_zmax) <= ztol .AND. in_channel(smesh%r(1,ic))) k = k + 1
END DO
ALLOCATE(inlet_pdofs(k))
k = 0
DO ic=1,MIN(SIZE(pump_sim%mug%T_bc),smesh%np)
  IF(ABS(smesh%r(2,ic)+chan_zmax) > ztol)CYCLE
  IF(.NOT.in_channel(smesh%r(1,ic)))CYCLE
  k = k + 1
  inlet_pdofs(k) = ic
END DO

IF(nin == 0 .OR. nout == 0) CALL oft_warn("End-face edges not found; check the mesh")
DEALLOCATE(vert_flag,edge_flag,in_flag,outlet_flag,dof_r,v_prof)
END BLOCK
!---Restart 0 was written by setup(), before the inlet profile existed
CALL pump_sim%rst_save(0.d0)
END IF

!==============================================================================
! Time advance, cuts timestep in half if solve fails
!==============================================================================
DO i=1,nsteps
  nretry = 0
  DO
    CALL set_coils(t+dt)
    CALL pump_sim%step(t,dt,nl_its,l_its)
    IF(nl_its >= 0)EXIT
    IF(nretry >= 4)THEN
      CALL pump_sim%rst_save(t)
      CALL oft_abort('Nonlinear solve failed after all retries','em_pump',__FILE__)
    END IF
    nretry = nretry + 1
    dt = dt/2.d0
  END DO
  CALL report(i,t,dt,nl_its,l_its,nretry)
END DO

CALL pump_sim%plot()
CALL oft_finalize

CONTAINS
!------------------------------------------------------------------------------
!> Travelling-wave coil current density [A/m^2] at a given time, into j_coil.
!!
!! Segment k carries I0*cos(omega*t - k*2*pi/n_per_lam) spread uniformly over
!! its area, so the applied field is a wave travelling in +Z. The envelope
!! ramps from zero so that t=0 (no field, no flow) is a consistent start.
!------------------------------------------------------------------------------
SUBROUTINE set_coils(time)
REAL(r8), INTENT(in) :: time
INTEGER(i4) :: jc
REAL(r8) :: omega,dphi,env
omega = 2.d0*pi*freq
dphi = 2.d0*pi/REAL(n_per_lam,8)
env = 1.d0
IF(ramp_time > 0.d0) env = MIN(1.d0,time/ramp_time)
pump_sim%j_coil = 0.d0
DO jc=1,n_seg
  pump_sim%j_coil(coil_reg0+jc-1) = env*I0*COS(omega*time - REAL(jc-1,8)*dphi)/coil_area(jc)
END DO
END SUBROUTINE set_coils
!------------------------------------------------------------------------------
!> True if transverse coordinate r lies on a duct wall
!------------------------------------------------------------------------------
LOGICAL FUNCTION on_wall(r)
REAL(r8), INTENT(in) :: r
on_wall = ABS(r - chan_xmax) < rtol .OR. (.NOT.solid_pipe .AND. ABS(r - chan_xmin) < rtol)
END FUNCTION on_wall
!------------------------------------------------------------------------------
!> True if transverse coordinate r lies within the fluid, walls included
!------------------------------------------------------------------------------
LOGICAL FUNCTION in_channel(r)
REAL(r8), INTENT(in) :: r
in_channel = r >= chan_xmin - rtol .AND. r <= chan_xmax + rtol
END FUNCTION in_channel
!------------------------------------------------------------------------------
!> One line per step: solver health and, with flow, peak velocity and head
!------------------------------------------------------------------------------
SUBROUTINE report(istep,time,dt_now,nl,nlin,nret)
INTEGER(i4), INTENT(in) :: istep,nl,nlin,nret
REAL(r8), INTENT(in) :: time,dt_now
REAL(r8), POINTER, DIMENSION(:) :: vz,pv
REAL(r8) :: vmax,p_in
INTEGER(i4) :: kk
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
WRITE(*,'(A,I6,A,ES12.5,A,ES10.3,A,I3,A,I5,A,I2,A,ES11.4,A,ES13.5)') &
  'Step ',istep,'  t = ',time,'  dt = ',dt_now,'  NL = ',nl,'  lin = ',nlin, &
  '  retry = ',nret,'  max|v_z| = ',vmax,'  head = ',-p_in
END SUBROUTINE report
END PROGRAM em_pump
