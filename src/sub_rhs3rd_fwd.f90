!-----------------------------------------------------------------------
subroutine sub_rhs_fwd
!-----------------------------------------------------------------------
  use mod_variables
  implicit none
!
  integer :: ic,jc
  real(8) :: dci,dci2,r12,dc12,dc12al
  real(8) :: xibl,yibl
!
  real(8) :: uxc,uxm,uxp,uym,uyp
  real(8) :: vyc,vxm,vxp,vym,vyp
  real(8) :: uxmm,uxpp,uymm,uyc,uypp
  real(8) :: vxmm,vxc,vxpp,vymm,vypp
!
  real(8) :: dudx1,dudy1
  real(8) :: dvdx1,dvdy1
  real(8) :: dudx2,dudy2
  real(8) :: dvdx2,dvdy2
!
  real(8) :: abs_uxc,abs_uyc
  real(8) :: abs_vxc,abs_vyc
!
  real(8) :: uold(imx1,jmax)
  real(8) :: vold(imax,jmx1)
!-----------------------------------------------------------------------
!
!
!---> Coefficients
  r12    = 1.d0/12.d0
  dci    = 1.d0/dx
  dci2   = dci**2/Re
  dc12   = dci*r12
  dc12al = dc12*3.d0
!
!
  udlt(1:imx1,1:jmax) = 0.d0
  vdlt(1:imax,1:jmx1) = 0.d0
  uold(1:imx1,1:jmax) = 0.d0
  vold(1:imax,1:jmx1) = 0.d0
!
!-----------------------------------------------------------------------
!---- x-direction
!-----------------------------------------------------------------------
    do jc=3,jmax-2
    do ic=3,imx1-2
!
      uxm  = ustg(ic-1,jc)
      uxp  = ustg(ic+1,jc)
      uym  = ustg(ic,jc-1)
      uyp  = ustg(ic,jc+1)
!
      uxmm = ustg(ic-2,jc)
      uxpp = ustg(ic+2,jc)
      uymm = ustg(ic,jc-2)
      uypp = ustg(ic,jc+2)
!
      uxc  = ustg(ic,jc)
      vxc  = 0.25d0*(vstg(ic  ,jc+1)+vstg(ic  ,jc) &
                    +vstg(ic-1,jc+1)+vstg(ic-1,jc))
!
      dudx1 = dc12  *(uxmm-8.d0*uxm         +8.d0*uxp-uxpp)
      dudy1 = dc12  *(uymm-8.d0*uym         +8.d0*uyp-uypp)
      dudx2 = dc12al*(uxmm-4.d0*uxm+6.d0*uxc-4.d0*uxp+uxpp)
      dudy2 = dc12al*(uymm-4.d0*uym+6.d0*uxc-4.d0*uyp+uypp)
!
      abs_uxc = dabs(uxc)
      abs_vxc = dabs(vxc)
!
!---- Inviscid term
      udlt(ic,jc) = (uxc*dudx1+abs_uxc*dudx2) &
                   +(vxc*dudy1+abs_vxc*dudy2)
!---- Viscous term
      udlt(ic,jc) = udlt(ic,jc)-(uxp-2.d0*uxc+uxm)*dci2 &
                               -(uyp-2.d0*uxc+uym)*dci2
!---- Pressure gradient term
      udlt(ic,jc) = udlt(ic,jc)+(pcnt(ic,jc)-pcnt(ic-1,jc))*dci
!
    enddo
    enddo
!#endif
!
!
!-----------------------------------------------------------------------
!---- y-direction
!-----------------------------------------------------------------------
    do jc=3,jmx1-2
    do ic=3,imax-2
!
      vxm  = vstg(ic-1,jc)
      vxp  = vstg(ic+1,jc)
      vym  = vstg(ic,jc-1)
      vyp  = vstg(ic,jc+1)
!
      vxmm = vstg(ic-2,jc)
      vxpp = vstg(ic+2,jc)
      vymm = vstg(ic,jc-2)
      vypp = vstg(ic,jc+2)
!
      uyc  = 0.25d0*(ustg(ic+1,jc  )+ustg(ic,jc  ) &
                    +ustg(ic+1,jc-1)+ustg(ic,jc-1))
      vyc  = vstg(ic,jc)
!
      dvdx1 = dc12  *(vxmm-8.d0*vxm         +8.d0*vxp-vxpp)
      dvdy1 = dc12  *(vymm-8.d0*vym         +8.d0*vyp-vypp)
      dvdx2 = dc12al*(vxmm-4.d0*vxm+6.d0*vyc-4.d0*vxp+vxpp)
      dvdy2 = dc12al*(vymm-4.d0*vym+6.d0*vyc-4.d0*vyp+vypp)
!
      abs_uyc = dabs(uyc)
      abs_vyc = dabs(vyc)
!
!---- Inviscid term
      vdlt(ic,jc) = (uyc*dvdx1+abs_uyc*dvdx2) &
                   +(vyc*dvdy1+abs_vyc*dvdy2)
!---- Viscous term
      vdlt(ic,jc) = vdlt(ic,jc)-(vxp-2.d0*vyc+vxm)*dci2 &
                               -(vyp-2.d0*vyc+vym)*dci2 
!---- Pressure gradient term
      vdlt(ic,jc) = vdlt(ic,jc)+(pcnt(ic,jc)-pcnt(ic,jc-1))*dci
!
  enddo
  enddo
!
!
  uold(1:imx1,1:jmax) = ustg(1:imx1,1:jmax)
  vold(1:imax,1:jmx1) = vstg(1:imax,1:jmx1)
  ustg(1:imx1,1:jmax) = 0.d0
  vstg(1:imax,1:jmx1) = 0.d0
!
  do jc=3,jmx1-2
  do ic=3,imx1-2
    if(jc==jmx1-2) goto 10
    xibl = dble(min0(iblk(ic-1,jc),iblk(ic,jc)))
    ustg(ic,jc) = uold(ic,jc)-dt*udlt(ic,jc)*xibl
 10 continue
    if(ic==imx1-2) goto 20
    yibl = dble(min0(iblk(ic,jc-1),iblk(ic,jc)))
    vstg(ic,jc) = vold(ic,jc)-dt*vdlt(ic,jc)*yibl
 20 continue
  enddo
  enddo
!
  return
end subroutine sub_rhs_fwd
!-----------------------------------------------------------------------

