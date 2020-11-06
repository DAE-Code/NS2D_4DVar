!-----------------------------------------------------------------------
! Calculate summation of energy in entire flow field
!-----------------------------------------------------------------------
subroutine sub_energy(sum)
!-----------------------------------------------------------------------
  use mod_variables,only : imax,jmax,dx,ustg,vstg
  implicit none
!
  integer :: ic,jc
  real(8) :: vol
  real(8) :: sum 
!-----------------------------------------------------------------------
!
  sum = 0.d0
  vol = dx**2
!
  do jc = 1,jmax
  do ic = 1,imax
    sum = sum+vol*0.5d0*(ustg(ic,jc)**2 &
                        +vstg(ic,jc)**2)
  enddo
  enddo
!
  return
end subroutine sub_energy
!-----------------------------------------------------------------------
!
!
!
!-----------------------------------------------------------------------
! Define time step size based on CFL number
!-----------------------------------------------------------------------
subroutine sub_timestep(dt)
!-----------------------------------------------------------------------
  use mod_variables,only : imax,jmax,ustg,vstg,dx,cfl
  implicit none
!
  integer :: ic,jc
  real(8) :: dt
!-----------------------------------------------------------------------
!
  dt  = 1.d0
!
  do jc = 1,jmax
  do ic = 1,imax
    dt = min(dt,cfl*dx/(abs(ustg(ic,jc)) &
                       +abs(vstg(ic,jc))))
  enddo
  enddo
!
  return
end subroutine sub_timestep
!-----------------------------------------------------------------------
!
!
!
!-----------------------------------------------------------------------
! Calculate maximum value of local CFL number
!-----------------------------------------------------------------------
subroutine sub_cfl(cfl)
!-----------------------------------------------------------------------
  use mod_variables,only : imax,jmax,ustg,vstg,dx,dt
  implicit none
!
  integer :: ic,jc
  real(8) :: cfl
!-----------------------------------------------------------------------
!
  cfl  = 0.d0
!
  do jc = 1,jmax
  do ic = 1,imax
    cfl = dmax1(cfl,(dabs(ustg(ic,jc)) &
                    +dabs(vstg(ic,jc)))*dt/dx)
  enddo
  enddo
!
  return
end subroutine sub_cfl
!-----------------------------------------------------------------------
!
!
!
!-----------------------------------------------------------------------
! Probing u-velocity
!-----------------------------------------------------------------------
subroutine sub_probe(istep)
!-----------------------------------------------------------------------
  use mod_variables,only : imax,jmax,xcen,ycen,ustg,vstg,nstep,iprb,jprb,prob
  implicit none
!
  real(8),parameter :: u_inf = 1.d0

  integer :: ic,jc,istep,ierr
!-----------------------------------------------------------------------
!
!
  if(istep==1) then
!
   allocate(prob(3,nstep),stat=ierr)
   prob(:,:) = 0.d0
!
   do ic=1,imax-1
      if(xcen(ic)<3.d0 .and. 3.d0<xcen(ic+1)) then
        iprb(1) = ic
        iprb(2) = ic
        iprb(3) = ic
        exit
      endif
    enddo
    do jc=1,jmax-1
      if(ycen(jc)<0.d0 .and. 0.d0<ycen(jc+1)) then
        jprb(1) = jc
      endif
      if(ycen(jc)<1.d0 .and. 1.d0<ycen(jc+1)) then
        jprb(2) = jc
      endif
      if(ycen(jc)<2.d0 .and. 2.d0<ycen(jc+1)) then
        jprb(3) = jc
      endif
    enddo
! 
  endif
!
  prob(1,istep) = dsqrt(ustg(iprb(1),jprb(1))**2+vstg(iprb(1),jprb(1))**2)
  prob(2,istep) = dsqrt(ustg(iprb(2),jprb(2))**2+vstg(iprb(2),jprb(2))**2)
  prob(3,istep) = dsqrt(ustg(iprb(3),jprb(3))**2+vstg(iprb(3),jprb(3))**2)
!
  if(istep==nstep) then
    open(60,file='probe.dat',form='formatted',status='replace',action='write')
    ic = 0
    write(60,'(i6,6f12.7)') ic,xcen(iprb(1)),xcen(iprb(2)),xcen(iprb(3)), &
                               ycen(jprb(1)),ycen(jprb(2)),ycen(jprb(3))
    do ic=1,nstep,8
      write(60,'(i6,3f12.7)') ic,prob(1,ic),prob(2,ic),prob(3,ic)
    enddo
    close(60)
  endif
!
  return
end subroutine sub_probe
!-----------------------------------------------------------------------
!
!
!
!-----------------------------------------------------------------------
! Probing u-velocity
!-----------------------------------------------------------------------
subroutine sub_probe2(istep,iter,icycl,char)
!-----------------------------------------------------------------------
  use mod_variables,only : jmax,imx1,xcen,ycen,ustg, &
                           iprb,jprb,prob,nstep,ncycl,nsall,ihx
  implicit none
!
  integer :: ic,jc,istep,ierr,iter,icycl,ip,icc
  real(8) :: xx,yy,xt,yt,dist,dmin
  character(len=3) :: char
!-----------------------------------------------------------------------
!
!
  if(istep==1 .and. iter==1 .and. icycl==1) then
!
    if(allocated(prob)) deallocate(prob)
    allocate(prob(4,nsall*ncycl),stat=ierr)
    prob(:,:) = 0.d0
!
!
!   Position 1
    xt   =-2.5d0 
    yt   =-0.5d0
    dmin = 9.d0
    do jc=1,jmax
    do ic=1,imx1
!     if(ihx(ic,jc,1)==1) then      
        xx   = xcen(ic)
        yy   = 0.5d0*(ycen(jc)+ycen(jc+1))
        dist = (xx-xt)**2+(yy-yt)**2
        if(dist<=dmin) then
          dmin    = dist
          iprb(1) = ic
          jprb(1) = jc
        endif
!     endif
    enddo
    enddo
!   Position 2
    xt   = 0.0d0 
    yt   =-2.5d0
    dmin = 9.d0
    do jc=1,jmax
    do ic=1,imx1
!     if(ihx(ic,jc,1)==1) then      
        xx   = xcen(ic)
        yy   = 0.5d0*(ycen(jc)+ycen(jc+1))
        dist = (xx-xt)**2+(yy-yt)**2
        if(dist<=dmin) then
          dmin    = dist
          iprb(2) = ic
          jprb(2) = jc
        endif
!     endif
    enddo
    enddo
!   Position 3
    xt   = 2.5d0 
    yt   = 0.0d0
    dmin = 9.d0
    do jc=1,jmax
    do ic=1,imx1
!     if(ihx(ic,jc,1)==1) then      
        xx   = xcen(ic)
        yy   = 0.5d0*(ycen(jc)+ycen(jc+1))
        dist = (xx-xt)**2+(yy-yt)**2
        if(dist<=dmin) then
          dmin    = dist
          iprb(3) = ic
          jprb(3) = jc
        endif
!     endif
    enddo
    enddo
!   Position 4
    xt   = 7.5d0 
    yt   = 0.0d0
    dmin = 9.d0
    do jc=1,jmax
    do ic=1,imx1
!     if(ihx(ic,jc,1)==1) then      
        xx   = xcen(ic)
        yy   = 0.5d0*(ycen(jc)+ycen(jc+1))
        dist = (xx-xt)**2+(yy-yt)**2
        if(dist<dmin) then
          dmin    = dist
          iprb(4) = ic
          jprb(4) = jc
        endif
!     endif
    enddo
    enddo
    jprb(4) = jprb(3)
!
    if(char=="ref") then
      open(60,file='probe_ref.dat',form='formatted',status='replace',action='write')
      write(60,'(" pos(x,y)",8f12.7)') (xcen(iprb(ip)),ycen(jprb(ip)),ip=1,4)
      write(60,'(" iter   step      ref1        ref2        ref3        ref4")')
      close(60)
      open(60,file='probe_mes.dat',form='formatted',status='replace',action='write')
      write(60,'(" pos(x,y)",8f12.7)') (xcen(iprb(ip)),ycen(jprb(ip)),ip=1,4)
      write(60,'(" iter   step      mes1        mes2        mes3        mes4")')
      close(60)
    elseif(char=="est") then
      open(60,file='probe_est.dat',form='formatted',status='replace',action='write')
      write(60,'(" pos(x,y)",8f12.7)') (xcen(iprb(ip)),ycen(jprb(ip)),ip=1,4)
      write(60,'(" iter   step      est1        est2        est3        est4")')
      close(60)
    elseif(char=="bst") then
      open(60,file='probe_bst.dat',form='formatted',status='replace',action='write')
      write(60,'(" pos(x,y)",8f12.7)') (xcen(iprb(ip)),ycen(jprb(ip)),ip=1,4)
      write(60,'(" iter   step      bst1        bst2        bst3        bst4")')
      close(60)
    endif
!
  endif
!
  prob(1,istep) = ustg(iprb(1),jprb(1))
  prob(2,istep) = ustg(iprb(2),jprb(2))
  prob(3,istep) = ustg(iprb(3),jprb(3))
  prob(4,istep) = ustg(iprb(4),jprb(4))
!
  if(char=="ref" .and. istep==icycl*nsall) then
    open(60,file='probe_ref.dat',form='formatted',position='append')
    do ic=1+(icycl-1)*nsall,(icycl)*nsall
      write(60,'(2i6,4f12.7)') iter,ic,prob(1,ic),prob(2,ic),prob(3,ic),prob(4,ic)
    enddo
    close(60)
    open(60,file='probe_mes.dat',form='formatted',position='append')
    do ic=1,nstep
      icc = ic+(icycl-1)*nstep
      if(maxval(ihx(:,:,ic))>0) then
        write(60,'(2i6,4f12.7)') iter,icc,prob(1,icc),prob(2,icc),prob(3,icc),prob(4,icc)
      endif
    enddo
    close(60)
  elseif(char=="est" .and. istep==icycl*nstep) then
    open(60,file='probe_est.dat',form='formatted',position='append')
    do ic=1+(icycl-1)*nstep,(icycl)*nstep
      write(60,'(2i6,4f12.7)') iter,ic,prob(1,ic),prob(2,ic),prob(3,ic),prob(4,ic)
    enddo
    write(60,'(2i6,4f12.7)')
    close(60)
  elseif(char=="bst" .and. istep==icycl*nsall) then
    open(60,file='probe_bst.dat',form='formatted',position='append')
    do ic=1+(icycl-1)*nsall,(icycl)*nsall
      write(60,'(2i6,4f12.7)') iter,ic,prob(1,ic),prob(2,ic),prob(3,ic),prob(4,ic)
    enddo
    close(60)
  endif
!
  return
end subroutine sub_probe2
!-----------------------------------------------------------------------
