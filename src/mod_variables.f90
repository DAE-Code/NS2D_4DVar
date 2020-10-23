!-----------------------------------------------------------------------
module mod_variables
  implicit none
!
  integer             :: modef              ! mode (1:TLM&ADJ Check, 2:4DV)
  integer             :: n_prb              ! problem (1:Karman vortex, 2:vortex advection)
  integer             :: n_mes              ! size of measurement vector
  integer             :: ncycl              ! number of 4DV period
!
  type variances
    real(8)           :: ref                ! Error variance of pseudo measurement
    real(8)           :: mes                ! Measurement error variance
    real(8)           :: bkg                ! Coefficit of background term
  end type variances
  type(variances) vars
!
  integer             :: iskip              ! Every iskip-th mesh point as measuremnt 
  integer             :: jskip              ! Every jskip-th mesh point as measuremnt 
  integer             :: tskip              ! Every tskip-th time step as measuremnt 
  real(8)             :: mbot               ! Bottom of measurement area
  real(8)             :: mtop               ! Top of measurement area
  real(8)             :: mlft               ! Left corner of measurement area
  real(8)             :: mrht               ! Right corner of measurement area
!
  real(8)             :: u_inf = 1.d0
  real(8)             :: v_inf = 0.d0
  real(8)             :: p_inf = 0.d0
! 
  integer             :: icycl              ! index  of 4DV period
  integer             :: itrp               ! # of iteration
  integer             :: ITRMAX
  real(8)             :: divmax
!
  real(8)             :: dt
  real(8)             :: dx
  real(8)             :: re
!
  real(8)             :: cfl
  integer             :: nstep              ! total # of step (DA window)
  integer             :: nsall              ! total # of step (for evaluation)
  real(8)             :: tstep              ! current time
  real(8)             :: tstep_ini          ! initial time
  integer             :: istep              ! time step
  integer             :: ostep              ! current step
  integer             :: ostep_ini          ! initial step
  integer             :: irestart           ! restart flag
  integer             :: ioutput            ! output flag
  integer             :: iskip_plot
  integer             :: iskip_clcd
  integer             :: iskip_hist
  real(8)             :: sumenergy
! 
  integer             :: iobj               ! object size
  integer             :: jobj               !
  integer             :: imax               ! # of mesh points
  integer             :: jmax               !
  integer             :: imx1               !
  integer             :: jmx1               !
  integer             :: corner(2,4)
  integer             :: iprb(4)
  integer             :: jprb(4)
!
  real(8),allocatable :: xcen(:)
  real(8),allocatable :: ycen(:)
!
  real(8),allocatable :: ustg(:,:)
  real(8),allocatable :: vstg(:,:)
  real(8),allocatable :: udlt(:,:)
  real(8),allocatable :: vdlt(:,:)
  real(8),allocatable :: pcnt(:,:)
  real(8),allocatable :: divg(:,:)
!
  real(8),allocatable :: ustg_b(:,:)
  real(8),allocatable :: vstg_b(:,:)
  real(8),allocatable :: udlt_b(:,:)
  real(8),allocatable :: vdlt_b(:,:)
  real(8),allocatable :: pcnt_b(:,:)
  real(8),allocatable :: divg_b(:,:)
!
  real(8),allocatable :: ustg_t(:,:)
  real(8),allocatable :: vstg_t(:,:)
  real(8),allocatable :: udlt_t(:,:)
  real(8),allocatable :: vdlt_t(:,:)
  real(8),allocatable :: pcnt_t(:,:)
!
  real(8),allocatable :: ustg_a(:,:)
  real(8),allocatable :: vstg_a(:,:)
  real(8),allocatable :: udlt_a(:,:)
  real(8),allocatable :: vdlt_a(:,:)
  real(8),allocatable :: pcnt_a(:,:)
!
  integer,allocatable :: iblk(:,:)
  real(8),allocatable :: uinf(:)
  real(8),allocatable :: prob(:,:)
!
! Measurement
  integer,allocatable :: ih (:,:)           ! Index for measurement points (1D to 1D)
  integer,allocatable :: ihx(:,:,:)         ! Index for measurement points (2D mesh)
  integer,allocatable :: ihy(:,:,:)         ! Index for measurement points (2D mesh)
!
end module mod_variables
!-----------------------------------------------------------------------
