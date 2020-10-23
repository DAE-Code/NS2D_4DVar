!-----------------------------------------------------------------------
subroutine sub_initial
!-----------------------------------------------------------------------
  use mod_variables
  implicit none
!
  integer :: ic,jc,kc,icen,jcen,ib,nblock,ierr,jmax_tmp,ios
  real(8) :: rblank
  real(4) :: xout(imax),yout(jmax)
!-----------------------------------------------------------------------
!
  ostep     = 0
  tstep     = 0.d0
!
!-Variables
  allocate(ustg(imx1,jmax),stat=ierr)
  allocate(vstg(imax,jmx1),stat=ierr)
  allocate(udlt(imx1,jmax),stat=ierr)
  allocate(vdlt(imax,jmx1),stat=ierr)
!
  allocate(pcnt(imax,jmax),stat=ierr)
  allocate(divg(imax,jmax),stat=ierr)
!
  allocate(iblk(imax,jmax),stat=ierr)
!
  allocate(uinf(jmax)     ,stat=ierr)
!
  allocate(xcen(imx1)     ,stat=ierr)
  allocate(ycen(jmx1)     ,stat=ierr)
!
  allocate(ustg_b(imx1,jmax),stat=ierr)
  allocate(vstg_b(imax,jmx1),stat=ierr)
  allocate(udlt_b(imx1,jmax),stat=ierr)
  allocate(vdlt_b(imax,jmx1),stat=ierr)
  allocate(pcnt_b(imax,jmax),stat=ierr)
!
  allocate(ustg_t(imx1,jmax),stat=ierr)
  allocate(vstg_t(imax,jmx1),stat=ierr)
  allocate(udlt_t(imx1,jmax),stat=ierr)
  allocate(vdlt_t(imax,jmx1),stat=ierr)
  allocate(pcnt_t(imax,jmax),stat=ierr)
!
  allocate(ustg_a(imx1,jmax),stat=ierr)
  allocate(vstg_a(imax,jmx1),stat=ierr)
  allocate(udlt_a(imx1,jmax),stat=ierr)
  allocate(vdlt_a(imax,jmx1),stat=ierr)
  allocate(pcnt_a(imax,jmax),stat=ierr)
!

  dx = 1.d0/dble(jobj)
!
  do ic=1,imx1
    xcen(ic) = dx* dble(ic-1-0.25*imax)
  enddo
  do ic=1,imax
    xout(ic) = dx*(dble(ic-0.25*imax)-0.5d0)
  enddo
  do jc=1,jmx1
    ycen(jc) = dx* dble(jc-1-0.50*jmax)
  enddo
  do jc=1,jmax
    yout(jc) = dx*(dble(jc-0.50*jmax)-0.5d0)
  enddo
!
!-Definition of object (a rectangular cylinder)
  icen = int(0.25d0*dble(imax))
  jcen = int(0.50d0*dble(jmax))
!
  do jc=1,jmax
  do ic=1,imax
    if(abs(xout(ic))<=0.5d0 .and. abs(yout(jc))<=0.5d0) then
      iblk(ic,jc) = 0
    else
      iblk(ic,jc) = 1
    endif
  enddo
  enddo
!
  do jc=2,jmax-1
  do ic=2,imax-1
    if(iblk(ic,jc)==0) then
      if(iblk(ic-1,jc)==1 .and. iblk(ic,jc-1)==1) then
        corner(1,1) = ic
        corner(2,1) = jc
      endif
      if(iblk(ic+1,jc)==1 .and. iblk(ic,jc-1)==1) then
        corner(1,2) = ic
        corner(2,2) = jc
      endif
      if(iblk(ic-1,jc)==1 .and. iblk(ic,jc+1)==1) then
        corner(1,3) = ic
        corner(2,3) = jc
      endif
      if(iblk(ic+1,jc)==1 .and. iblk(ic,jc+1)==1) then
        corner(1,4) = ic
        corner(2,4) = jc
      endif
    endif
  enddo
  enddo
!
  do jc=1,jmax
  do ic=1,imax
    rblank = dble(iblk(ic,jc))
    ustg(ic  ,jc) = u_inf*rblank
    vstg(ic,jc  ) = v_inf*rblank
    ustg(ic+1,jc) = u_inf*rblank
    vstg(ic,jc+1) = v_inf*rblank
    udlt(ic  ,jc) = 0.0
    vdlt(ic,jc  ) = 0.0
    udlt(ic+1,jc) = 0.0
    vdlt(ic,jc+1) = 0.0
  enddo
  enddo
!
! Inflow velocity
  uinf(:) = u_inf
!
  do jc=1,jmax
  do ic=1,imax
    pcnt(ic,jc) = 0.0
    divg(ic,jc) = 0.0
  enddo
  enddo
!
  cfl = 0.02
!!call sub_timestep(dt)
  call sub_cfl(cfl)
  call sub_bc_wall(ustg,vstg,pcnt)
  call sub_bc_outer(ustg,vstg,pcnt)
!
!
  if(ioutput==1) then
! Output BCM mesh >>>>>>
  open(30,file='./output/mesh.g',form='unformatted',status='replace',action='write')
  nblock = 1
  write(30) nblock
  write(30) (imax,jmax,ib=1,nblock)
  do ib=1,nblock
    write(30) ((xout(ic)   ,ic=1,imax),jc=1,jmax), &
              ((yout(jc)   ,ic=1,imax),jc=1,jmax), &
              ((iblk(ic,jc),ic=1,imax),jc=1,jmax)
  enddo
  close(30)
! Output BCM mesh <<<<<<
  endif
!
  return
end subroutine sub_initial
!-----------------------------------------------------------------------
!
!
!
! ----------------------------------------------------------------------
subroutine sub_restart
!-----------------------------------------------------------------------
  use mod_variables
  implicit none

  integer :: ic,jc,kc,m,imx,jmx,ib,nblock
  real(4) :: re4,step4,time4
  real(4) :: q4(imax,jmax,4)
!-----------------------------------------------------------------------
!
!
!-Read an instantaneous field
  open(15,file='./restart.q',form='unformatted',action='read',status='old')
  nblock = 1
  read(15) nblock
  read(15) (imx,jmx,ib=1,nblock)

  ! consistency check >>>>>>
  if(imx/=imax .or. jmx/=jmax) stop 'stop imx/=imax or jmx/=jmax'

  do ib=1,nblock
    read(15) re4,step4,re4,time4
    read(15) (((q4(ic,jc,m),ic=1,imax),jc=1,jmax),m=1,4)
  enddo
  close(15)
!
  ostep = int(step4)
  tstep = dble(time4)
!
  do jc=1,jmax
    do ic=2,imax
      ustg(ic,jc) = dble(q4(ic-1,jc,2)+q4(ic,jc,2))*0.5d0
    enddo
    ustg(1   ,jc) = dble(q4(1   ,jc,2))
    ustg(imx1,jc) = dble(q4(imax,jc,2))
  enddo
  do ic=1,imax
    do jc=2,jmax
      vstg(ic,jc) = dble(q4(ic,jc-1,3)+q4(ic,jc,3))*0.5d0
    enddo
    vstg(ic,1   ) = dble(q4(ic,1   ,3))
    vstg(ic,jmx1) = dble(q4(ic,jmax,3))
  enddo
  do jc=1,jmax
  do ic=1,imax
    pcnt(ic,jc) = dble(q4(ic,jc,4))
  enddo
  enddo
!
! Inlet boundary condition
  do jc=1,jmax
    uinf(jc) = ustg(1,jc)
  enddo
!
  return
end subroutine sub_restart
!-----------------------------------------------------------------------
!
!
!
! ----------------------------------------------------------------------
subroutine sub_restart_4DV(ios)
!-----------------------------------------------------------------------
  use mod_variables
  implicit none

  integer :: ic,jc,kc,m,imx,jmx,ib,nblock,ios
  real(4) :: re4,step4,time4
  real(4) :: q4(imax,jmax,4)
!-----------------------------------------------------------------------
!
!
!-Read an instantaneous field
  open(15,file='./initial.q',form='unformatted',action='read',status='old',iostat=ios)
  if(ios==0) then
  nblock = 1
  read(15) nblock
  read(15) (imx,jmx,ib=1,nblock)

  ! consistency check >>>>>>
  if(imx/=imax .or. jmx/=jmax) stop 'stop imx/=imax or jmx/=jmax'

  do ib=1,nblock
    read(15) re4,step4,re4,time4
    read(15) (((q4(ic,jc,m),ic=1,imax),jc=1,jmax),m=1,4)
  enddo
  close(15)
  endif
!
  ostep = int(step4)
  tstep = dble(time4)
!
  do jc=1,jmax
    do ic=2,imax
      ustg(ic,jc) = dble(q4(ic-1,jc,2)+q4(ic,jc,2))*0.5d0
    enddo
    ustg(1   ,jc) = dble(q4(1   ,jc,2))
    ustg(imx1,jc) = dble(q4(imax,jc,2))
  enddo
  do ic=1,imax
    do jc=2,jmax
      vstg(ic,jc) = dble(q4(ic,jc-1,3)+q4(ic,jc,3))*0.5d0
    enddo
    vstg(ic,1   ) = dble(q4(ic,1   ,3))
    vstg(ic,jmx1) = dble(q4(ic,jmax,3))
  enddo
  do jc=1,jmax
  do ic=1,imax
    pcnt(ic,jc) = dble(q4(ic,jc,4))
  enddo
  enddo
!
! Inlet boundary condition
  do jc=1,jmax
    uinf(jc) = ustg(1,jc)
  enddo
!
  return
end subroutine sub_restart_4DV
!-----------------------------------------------------------------------
!
!
!
!-----------------------------------------------------------------------
subroutine sub_addvtx(xc,yc)
!-----------------------------------------------------------------------
  use mod_variables
  implicit none
!
  integer :: ic,jc,icen,jcen,ib,nblock,ierr,jmax_tmp,ios
  real(8),parameter :: rc = 0.5d0
  real(8),parameter :: ga = 1.0d0 
  real(8) :: rblank,gp,xc,yc,rr,pi,yl
!-----------------------------------------------------------------------
!
!
  pi = 4.d0*datan(1.d0)
  gp = ga/(2.d0*pi)
  yl = ycen(jmx1)-ycen(1)
!
  do jc=1,jmax
  do ic=1,imax
    rblank = real(iblk(ic,jc))
!   Main vortex
    rr     = (xcen(ic)-xc)**2+(ycen(jc)-yc)**2+rc**2
    ustg(ic,jc) = ustg(ic,jc)-(ga/rr)*(ycen(jc)-yc)*rblank
    vstg(ic,jc) = vstg(ic,jc)+(ga/rr)*(xcen(ic)-xc)*rblank
!   Lower mirror vortex
    rr     = (xcen(ic)-xc)**2+(ycen(jc)-yc-yl)**2+rc**2
    ustg(ic,jc) = ustg(ic,jc)+(ga/rr)*(ycen(jc)-yc)*rblank
    vstg(ic,jc) = vstg(ic,jc)-(ga/rr)*(xcen(ic)-xc)*rblank
!   Upper mirror vortex
    rr     = (xcen(ic)-xc)**2+(ycen(jc)-yc+yl)**2+rc**2
    ustg(ic,jc) = ustg(ic,jc)+(ga/rr)*(ycen(jc)-yc)*rblank
    vstg(ic,jc) = vstg(ic,jc)-(ga/rr)*(xcen(ic)-xc)*rblank
  enddo
  enddo
!
  return
end subroutine sub_addvtx
!-----------------------------------------------------------------------
