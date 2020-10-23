!-----------------------------------------------------------------------
!
!   4D Variational Method (4D-Var) coupled with 2D Navier-Stokes code 
!   prepared for the book "Data Assimilation Fluid Science"
!
!   Further information is avairable on: https://github.com/DAE-Code
!
!-----------------------------------------------------------------------
program NS2D_4DVar
  use mod_variables
  implicit none
!
  real(4) :: time1,tnow
!-----------------------------------------------------------------------
!
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Reading initial parameters
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  open(28,file='4DVar.inp',form='formatted')
  read(28,*) modef               ! Mode (1:TLM&ADJ Check, 2:4DV, 3:4DV restart)
  read(28,*) n_prb               ! Problem 1:Karman vortex, 2:Vortex advection
  read(28,*) ncycl,ITRMAX        ! Number of 4DV cycles, iteration on each cycle
  read(28,*) vars%ref            ! Error variance of pseudo measurement
  read(28,*) vars%mes            ! Measurement error variance
  read(28,*) vars%bkg            ! Coefficit of background term
  read(28,*) iskip,jskip,tskip   ! Every iskip & jskip & tskip for measuremnt 
  read(28,*) mlft,mrht           ! Left and right of measurement area
  read(28,*) mbot,mtop           ! Bottom and top of measurement area
  read(28,*)
  read(28,*) nstep,nsall         ! DA window steps, total # of time steps
  read(28,*) Re                  ! Reynolds number 
  read(28,*) jmax                ! Number of mesh: jmax (imax x jmax, imax = 3*jmax)
  read(28,*) dt                  ! Time step
  read(28,*) irestart            ! 0:initial, 1:restart
  read(28,*) ioutput             ! 0:No output to screen,1:otherwise
  read(28,*) iskip_plot          ! Output plot3d interval
  read(28,*) iskip_hist          ! Output history interval
  close(28)
  nsall = max0(nstep,nsall)
!
! Specify the number of points per one side of square cylinder
  jobj = int(0.100*jmax)  ! 8 points per one side when jmax = 80
! jobj = int(0.064*jmax)  ! jmax = 320: the original setting 2018/12/26
  iobj = jobj
  imax = 3*jmax
!
  imx1 = imax+1
  jmx1 = jmax+1
!
  if(ioutput==1) call system('mkdir output')
  call cpu_time(time1)
!
!
!-----------------------------------------------------------------------
!-Initial setup
!-----------------------------------------------------------------------
  call sub_initial
!
  if(irestart==1) then
    if(ioutput==1) write(*,*) 'Restart'
    call sub_restart
    tstep_ini = tstep
    ostep_ini = ostep
  else
    if(ioutput==1) write(*,*) 'Impulse start'
    tstep = 0.0
    ostep = 0
  endif
  if(ioutput==1) call sub_p3dwrite(ostep,tstep,'res')
!-----------------------------------------------------------------------
!-Initial setup end
!-----------------------------------------------------------------------
!
!
  write(*,*) 
  write(*,*) '----------------------------------------------------------'
  write(*,*) '1:ADJ Chk, 2:4DV, 3:4DV-R = ',modef
  write(*,*) 'Number of 4DV cycle       = ',ncycl
  write(*,*) 'Number of time-step       = ',nstep
  write(*,*) 'Mesh points               = ',imax,jmax
  write(*,*) 'Time-step size            = ',dt
  write(*,*) 'Grid spacing              = ',dx
  write(*,*) 'Reynolds number           = ',Re
  write(*,*) '----------------------------------------------------------'
  write(*,*)
!
!
!=======================================================================
!-TLM & ADJ check
!=======================================================================
  if(modef==1) then
    call sub_check_foa
!
!=======================================================================
!-4D-Var analysis
!=======================================================================
  elseif(modef==2 .or. modef==3) then
    call sub_4dvar
!
  endif
!
  call cpu_time(tnow)
  if(ioutput==1) write(*,'("Program end",f15.3," [sec]")') tnow-time1
!
end program NS2D_4DVar
!-----------------------------------------------------------------------
