!-----------------------------------------------------------------------
subroutine sub_filter(ustg,vstg,Jdf)
!-----------------------------------------------------------------------
  use mod_variables,only : imax,jmax,imx1,jmx1,vars,xcen,ycen,ihx,ihy
  implicit none
!
  integer :: ic,jc
  real(8) :: coef
  real(8) :: hflag,Jdf
  real(8) :: ustg(imx1,jmax),vstg(imax,jmx1)
  real(8) :: utp1(imx1,jmax),vtp1(imax,jmx1)
  real(8) :: utp2(imx1,jmax),vtp2(imax,jmx1)
!-----------------------------------------------------------------------
!
!
  coef = vars%bkg
!
!---> Filtering
  do jc=1,jmax
  do ic=2,imax-1
    utp1(ic,jc) = 0.25d0*(ustg(ic-1,jc)+2.d0*ustg(ic,jc)+ustg(ic+1,jc))
    vtp1(ic,jc) = 0.25d0*(vstg(ic-1,jc)+2.d0*vstg(ic,jc)+vstg(ic+1,jc))
  enddo
  enddo
  do jc=2,jmax-1
  do ic=1,imax
    utp2(ic,jc) = 0.25d0*(utp1(ic,jc-1)+2.d0*utp1(ic,jc)+utp1(ic,jc+1))
    vtp2(ic,jc) = 0.25d0*(vtp1(ic,jc-1)+2.d0*vtp1(ic,jc)+vtp1(ic,jc+1))
  enddo
  enddo
!
!---> Jdf
  Jdf = 0.d0
  do jc=3,jmax-2
  do ic=3,imx1-2
    if(-7.0d0<xcen(ic) .and. xcen(ic)<-2.d0 .and. -3.0d0<ycen(jc) .and. ycen(jc)<3.0d0) then   ! Transient vortex
      hflag = dble(ihx(ic,jc,1))
      Jdf = Jdf+0.5d0*coef*(ustg(ic,jc)-utp2(ic,jc))**2 !*hflag
    endif
  enddo
  enddo
  do jc=3,jmx1-2
  do ic=3,imax-2
    if(-7.0d0<xcen(ic) .and. xcen(ic)<-2.d0 .and. -3.0d0<ycen(jc) .and. ycen(jc)<3.0d0) then   ! Transient vortex
      hflag = dble(ihy(ic,jc,1))
      Jdf = Jdf+0.5d0*coef*(vstg(ic,jc)-vtp2(ic,jc))**2 !*hflag
    endif
  enddo
  enddo
!
  return
end subroutine sub_filter
!-----------------------------------------------------------------------
!
!
!
!-----------------------------------------------------------------------
subroutine sub_filter_adj(ustg,vstg,ustg_b,vstg_b)
!-----------------------------------------------------------------------
  use mod_variables,only : imax,jmax,imx1,jmx1,vars,xcen,ycen,ihx,ihy
  implicit none
!
  integer :: ic,jc
  real(8) :: coef
  real(8) :: hflag
  real(8) :: ustg_b(imx1,jmax),vstg_b(imax,jmx1)
  real(8) :: ustg(imx1,jmax),vstg(imax,jmx1)
  real(8) :: utp1(imx1,jmax),vtp1(imax,jmx1)
  real(8) :: utp2(imx1,jmax),vtp2(imax,jmx1)
!-----------------------------------------------------------------------
!
!
  coef = vars%bkg
!
!---> Filtering
  do jc=1,jmax
  do ic=2,imax-1
    utp1(ic,jc) = 0.25d0*(ustg_b(ic-1,jc)+2.d0*ustg_b(ic,jc)+ustg_b(ic+1,jc))
    vtp1(ic,jc) = 0.25d0*(vstg_b(ic-1,jc)+2.d0*vstg_b(ic,jc)+vstg_b(ic+1,jc))
  enddo
  enddo
  do jc=2,jmax-1
  do ic=1,imax
    utp2(ic,jc) = 0.25d0*(utp1(ic,jc-1)+2.d0*utp1(ic,jc)+utp1(ic,jc+1))
    vtp2(ic,jc) = 0.25d0*(vtp1(ic,jc-1)+2.d0*vtp1(ic,jc)+vtp1(ic,jc+1))
  enddo
  enddo
!
!---> Gradient of Jdf
  do jc=3,jmax-2
  do ic=3,imx1-2
    if(-7.0d0<xcen(ic) .and. xcen(ic)<-2.d0 .and. -3.0d0<ycen(jc) .and. ycen(jc)<3.0d0) then   ! Transient vortex
      hflag = dble(ihx(ic,jc,1))
      ustg(ic,jc) = ustg(ic,jc)+coef*(ustg_b(ic,jc)-utp2(ic,jc)) !*hflag
    endif
  enddo
  enddo
  do jc=3,jmx1-2
  do ic=3,imax-2
    if(-7.0d0<xcen(ic) .and. xcen(ic)<-2.d0 .and. -3.0d0<ycen(jc) .and. ycen(jc)<3.0d0) then   ! Transient vortex
      hflag = dble(ihy(ic,jc,1))
      vstg(ic,jc) = vstg(ic,jc)+coef*(vstg_b(ic,jc)-vtp2(ic,jc)) !*hflag
    endif
  enddo
  enddo
!
  return
end subroutine sub_filter_adj
!**********************************************************************
