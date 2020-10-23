!-----------------------------------------------------------------------
subroutine sub_rhs_tlm
!-----------------------------------------------------------------------
  use mod_variables
  implicit none
!
  integer :: ic,jc
  real(8) :: dci,dci2,r12,dc12,dc12al
  real(8) :: xibl,yibl
!
!=====Base flow
  real(8) :: uxc_b,uxm_b,uxp_b,uym_b,uyp_b
  real(8) :: vyc_b,vxm_b,vxp_b,vym_b,vyp_b
  real(8) :: uxmm_b,uxpp_b,uymm_b,uyc_b,uypp_b
  real(8) :: vxmm_b,vxc_b,vxpp_b,vymm_b,vypp_b
!
  real(8) :: dudx1_b,dudy1_b
  real(8) :: dvdx1_b,dvdy1_b
  real(8) :: dudx2_b,dudy2_b
  real(8) :: dvdx2_b,dvdy2_b
!
  real(8) :: abs_uxc_b,abs_uyc_b
  real(8) :: abs_vxc_b,abs_vyc_b
!
!=====Disturbance
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
!-----------------------------------------------------------------------
!---- x-direction
!-----------------------------------------------------------------------
    do jc=3,jmax-2
    do ic=3,imx1-2
!
!=====Base flow
      uxm_b  = ustg_b(ic-1,jc)
      uxp_b  = ustg_b(ic+1,jc)
      uym_b  = ustg_b(ic,jc-1)
      uyp_b  = ustg_b(ic,jc+1)
!
      uxmm_b = ustg_b(ic-2,jc)
      uxpp_b = ustg_b(ic+2,jc)
      uymm_b = ustg_b(ic,jc-2)
      uypp_b = ustg_b(ic,jc+2)
!
      uxc_b  = ustg_b(ic,jc)
      vxc_b  = 0.25d0*(vstg_b(ic  ,jc+1)+vstg_b(ic  ,jc) &
                      +vstg_b(ic-1,jc+1)+vstg_b(ic-1,jc))
!
      dudx1_b = dc12  *(uxmm_b-8.d0*uxm_b           +8.d0*uxp_b-uxpp_b)
      dudy1_b = dc12  *(uymm_b-8.d0*uym_b           +8.d0*uyp_b-uypp_b)
      dudx2_b = dc12al*(uxmm_b-4.d0*uxm_b+6.d0*uxc_b-4.d0*uxp_b+uxpp_b)
      dudy2_b = dc12al*(uymm_b-4.d0*uym_b+6.d0*uxc_b-4.d0*uyp_b+uypp_b)
!
      abs_uxc_b = dabs(uxc_b)
      abs_vxc_b = dabs(vxc_b)
!
!---- Inviscid term
      udlt_b(ic,jc) = (uxc_b*dudx1_b+abs_uxc_b*dudx2_b) &
                     +(vxc_b*dudy1_b+abs_vxc_b*dudy2_b)
!---- Viscous term
      udlt_b(ic,jc) = udlt_b(ic,jc)-(uxp_b-2.d0*uxc_b+uxm_b)*dci2 &
                                   -(uyp_b-2.d0*uxc_b+uym_b)*dci2
!---- Pressure gradient term
      udlt_b(ic,jc) = udlt_b(ic,jc)+(pcnt_b(ic,jc)-pcnt_b(ic-1,jc))*dci
!
!
!=====Disturbance
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
      abs_uxc   = dsign(1.d0,uxc_b)*uxc
      abs_vxc   = dsign(1.d0,vxc_b)*vxc
!
!---- Inviscid term
      udlt(ic,jc) = (uxc*dudx1_b+uxc_b*dudx1+abs_uxc*dudx2_b+abs_uxc_b*dudx2) &
                   +(vxc*dudy1_b+vxc_b*dudy1+abs_vxc*dudy2_b+abs_vxc_b*dudy2)
!---- Viscous term
      udlt(ic,jc) = udlt(ic,jc)-(uxp-2.d0*uxc+uxm)*dci2 &
                               -(uyp-2.d0*uxc+uym)*dci2
!---- Pressure gradient term
      udlt(ic,jc) = udlt(ic,jc)+(pcnt(ic,jc)-pcnt(ic-1,jc))*dci
!
    enddo
    enddo
!
!
!-----------------------------------------------------------------------
!---- y-direction
!-----------------------------------------------------------------------
    do jc=3,jmx1-2
    do ic=3,imax-2
!
!=====Base flow
      vxm_b  = vstg_b(ic-1,jc)
      vxp_b  = vstg_b(ic+1,jc)
      vym_b  = vstg_b(ic,jc-1)
      vyp_b  = vstg_b(ic,jc+1)
!
      vxmm_b = vstg_b(ic-2,jc)
      vxpp_b = vstg_b(ic+2,jc)
      vymm_b = vstg_b(ic,jc-2)
      vypp_b = vstg_b(ic,jc+2)
!
      uyc_b  = 0.25d0*(ustg_b(ic+1,jc  )+ustg_b(ic,jc  ) &
                      +ustg_b(ic+1,jc-1)+ustg_b(ic,jc-1))
      vyc_b  = vstg_b(ic,jc)
!
      dvdx1_b = dc12  *(vxmm_b-8.d0*vxm_b           +8.d0*vxp_b-vxpp_b)
      dvdy1_b = dc12  *(vymm_b-8.d0*vym_b           +8.d0*vyp_b-vypp_b)
      dvdx2_b = dc12al*(vxmm_b-4.d0*vxm_b+6.d0*vyc_b-4.d0*vxp_b+vxpp_b)
      dvdy2_b = dc12al*(vymm_b-4.d0*vym_b+6.d0*vyc_b-4.d0*vyp_b+vypp_b)
!
      abs_uyc_b = dabs(uyc_b)
      abs_vyc_b = dabs(vyc_b)
!
!---- Inviscid term
      vdlt_b(ic,jc) = (uyc_b*dvdx1_b+abs_uyc_b*dvdx2_b) &
                     +(vyc_b*dvdy1_b+abs_vyc_b*dvdy2_b)
!---- Viscous term
      vdlt_b(ic,jc) = vdlt_b(ic,jc)-(vxp_b-2.d0*vyc_b+vxm_b)*dci2 &
                                   -(vyp_b-2.d0*vyc_b+vym_b)*dci2 
!---- Pressure gradient term
      vdlt_b(ic,jc) = vdlt_b(ic,jc)+(pcnt_b(ic,jc)-pcnt_b(ic,jc-1))*dci
!
!
!=====Disturbance
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
      abs_uyc = dsign(1.d0,uyc_b)*uyc
      abs_vyc = dsign(1.d0,vyc_b)*vyc
!
!---- Inviscid term
      vdlt(ic,jc) = (uyc*dvdx1_b+uyc_b*dvdx1+abs_uyc*dvdx2_b+abs_uyc_b*dvdx2) &
                   +(vyc*dvdy1_b+vyc_b*dvdy1+abs_vyc*dvdy2_b+abs_vyc_b*dvdy2)
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
!
  do jc=3,jmx1-2
  do ic=3,imx1-2
    if(jc==jmx1-2) goto 10
    xibl = dble(min0(iblk(ic-1,jc),iblk(ic,jc)))
    ustg_b(ic,jc) = ustg_b(ic,jc)-dt*udlt_b(ic,jc)*xibl
    ustg  (ic,jc) = uold  (ic,jc)-dt*udlt  (ic,jc)*xibl
 10 continue
    if(ic==imx1-2) goto 20
    yibl = dble(min0(iblk(ic,jc-1),iblk(ic,jc)))
    vstg_b(ic,jc) = vstg_b(ic,jc)-dt*vdlt_b(ic,jc)*yibl
    vstg  (ic,jc) = vold  (ic,jc)-dt*vdlt  (ic,jc)*yibl
 20 continue
  enddo
  enddo
!
  return
end subroutine sub_rhs_tlm
!-----------------------------------------------------------------------
