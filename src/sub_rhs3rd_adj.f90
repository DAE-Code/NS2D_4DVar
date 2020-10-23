!-----------------------------------------------------------------------
subroutine sub_rhs_adj
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
  uold(1:imx1,1:jmax) = 0.d0  ! remove edges
  vold(1:imax,1:jmx1) = 0.d0
!
  do jc=3,jmx1-2
  do ic=3,imx1-2
    if(jc==jmx1-2) goto 10
    xibl = dble(min0(iblk(ic-1,jc),iblk(ic,jc)))
    udlt(ic,jc) =-dt*ustg(ic,jc)*xibl
    uold(ic,jc) = ustg(ic,jc)
 10 continue
    if(ic==imx1-2) goto 20
    yibl = dble(min0(iblk(ic,jc-1),iblk(ic,jc)))
    vdlt(ic,jc) =-dt*vstg(ic,jc)*yibl
    vold(ic,jc) = vstg(ic,jc)
 20 continue
  enddo
  enddo
!
! ustg(1:imx1,1:jmax) = 0.d0
! vstg(1:imax,1:jmx1) = 0.d0
  ustg(1:imx1,1:jmax) = uold(1:imx1,1:jmax)
  vstg(1:imax,1:jmx1) = vold(1:imax,1:jmx1)
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
      dvdx1_b = dc12  *(vxmm_b-8.d0*vxm_b          +8.0*vxp_b-vxpp_b)
      dvdy1_b = dc12  *(vymm_b-8.d0*vym_b          +8.0*vyp_b-vypp_b)
      dvdx2_b = dc12al*(vxmm_b-4.d0*vxm_b+6.0*vyc_b-4.0*vxp_b+vxpp_b)
      dvdy2_b = dc12al*(vymm_b-4.d0*vym_b+6.0*vyc_b-4.0*vyp_b+vypp_b)
!
      abs_uyc_b = dabs(uyc_b)
      abs_vyc_b = dabs(vyc_b)
!
!---- Inviscid term
!     vdlt_b(ic,jc) = (uyc_b*dvdx1_b+abs(uyc_b)*dvdx2_b) &
!                    +(vyc_b*dvdy1_b+abs(vyc_b)*dvdy2_b)
!---- Viscous term
!     vdlt_b(ic,jc) = vdlt_b(ic,jc)-(vxp_b-2.0*vyc_b+vxm_b)*dci2 &
!                                  -(vyp_b-2.0*vyc_b+vym_b)*dci2 
!---- Pressure gradient term
!     vdlt_b(ic,jc) = vdlt_b(ic,jc)+(pcnt_b(ic,jc)-pcnt_b(ic,jc-1))*dci
!
!
!=====Adjoint
!---- Pressure gradient term
      pcnt(ic,jc  ) = (pcnt(ic,jc  ))+vdlt(ic,jc)*dci
      pcnt(ic,jc-1) = (pcnt(ic,jc-1))-vdlt(ic,jc)*dci
!---- Viscous term
      uyc = 0.d0
      vyc = 0.d0
      vxp =      -vdlt(ic,jc)*dci2
      vyc = (vyc)+vdlt(ic,jc)*dci2*2.d0
      vxm =      -vdlt(ic,jc)*dci2
      vyp =      -vdlt(ic,jc)*dci2
      vyc = (vyc)+vdlt(ic,jc)*dci2*2.d0
      vym =      -vdlt(ic,jc)*dci2
!
!---- Inviscid term
      uyc     =       dvdx1_b  *vdlt(ic,jc)
      dvdx1   =       uyc_b    *vdlt(ic,jc)
      abs_uyc =       dvdx2_b  *vdlt(ic,jc)
      dvdx2   =       abs_uyc_b*vdlt(ic,jc)
      vyc     = (vyc)+dvdy1_b  *vdlt(ic,jc)
      dvdy1   =       vyc_b    *vdlt(ic,jc)
      abs_vyc =       dvdy2_b  *vdlt(ic,jc)
      dvdy2   =       abs_vyc_b*vdlt(ic,jc)
!
      uyc   = (uyc)+dsign(1.d0,uyc_b)*abs_uyc
      vyc   = (vyc)+dsign(1.d0,vyc_b)*abs_vyc
!
      vxmm =            dc12*dvdx1     +dc12al*dvdx2
      vymm =            dc12*dvdy1     +dc12al*dvdy2
      vxm  = (vxm)-8.d0*dc12*dvdx1-4.d0*dc12al*dvdx2
      vym  = (vym)-8.d0*dc12*dvdy1-4.d0*dc12al*dvdy2
      vyc  = (vyc)                +6.d0*dc12al*(dvdx2+dvdy2)
      vxp  = (vxp)+8.d0*dc12*dvdx1-4.d0*dc12al*dvdx2
      vyp  = (vyp)+8.d0*dc12*dvdy1-4.d0*dc12al*dvdy2
      vxpp =           -dc12*dvdx1     +dc12al*dvdx2
      vypp =           -dc12*dvdy1     +dc12al*dvdy2
!
      ustg(ic+1,jc  ) = (ustg(ic+1,jc  ))+0.25d0*uyc
      ustg(ic  ,jc  ) = (ustg(ic  ,jc  ))+0.25d0*uyc
      ustg(ic+1,jc-1) = (ustg(ic+1,jc-1))+0.25d0*uyc
      ustg(ic  ,jc-1) = (ustg(ic  ,jc-1))+0.25d0*uyc
      vstg(ic  ,jc  ) = (vstg(ic  ,jc  ))       +vyc
!
      vstg(ic-2,jc) = (vstg(ic-2,jc))+vxmm 
      vstg(ic+2,jc) = (vstg(ic+2,jc))+vxpp 
      vstg(ic,jc-2) = (vstg(ic,jc-2))+vymm  
      vstg(ic,jc+2) = (vstg(ic,jc+2))+vypp 
!
      vstg(ic-1,jc) = (vstg(ic-1,jc))+vxm 
      vstg(ic+1,jc) = (vstg(ic+1,jc))+vxp 
      vstg(ic,jc-1) = (vstg(ic,jc-1))+vym 
      vstg(ic,jc+1) = (vstg(ic,jc+1))+vyp 
!
  enddo
  enddo
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
!     udlt_b(ic,jc) = (uxc_b*dudx1_b+abs(uxc_b)*dudx2_b) &
!                    +(vxc_b*dudy1_b+abs(vxc_b)*dudy2_b)
!---- Viscous term
!     udlt_b(ic,jc) = udlt_b(ic,jc)-(uxp_b-2.0*uxc_b+uxm_b)*dci2 &
!                                  -(uyp_b-2.0*uxc_b+uym_b)*dci2
!---- Pressure gradient term
!     udlt_b(ic,jc) = udlt_b(ic,jc)+(pcnt_b(ic,jc)-pcnt_b(ic-1,jc))*dci
!
!
!=====Adjoint
!---- Pressure gradient term
      pcnt(ic  ,jc) = (pcnt(ic  ,jc))+udlt(ic,jc)*dci
      pcnt(ic-1,jc) = (pcnt(ic-1,jc))-udlt(ic,jc)*dci
!---- Viscous term
      uxc = 0.d0
      vxc = 0.d0
      uxp =      -udlt(ic,jc)*dci2
      uxc = (uxc)+udlt(ic,jc)*dci2*2.d0
      uxm =      -udlt(ic,jc)*dci2
      uyp =      -udlt(ic,jc)*dci2
      uxc = (uxc)+udlt(ic,jc)*dci2*2.d0
      uym =      -udlt(ic,jc)*dci2
!---- Inviscid term
      uxc     = (uxc)+dudx1_b  *udlt(ic,jc)
      dudx1   =       uxc_b    *udlt(ic,jc)
      abs_uxc =       dudx2_b  *udlt(ic,jc)
      dudx2   =       abs_uxc_b*udlt(ic,jc)
      vxc     =       dudy1_b  *udlt(ic,jc)
      dudy1   =       vxc_b    *udlt(ic,jc)
      abs_vxc =       dudy2_b  *udlt(ic,jc)
      dudy2   =       abs_vxc_b*udlt(ic,jc)
!
      uxc   = (uxc)+dsign(1.d0,uxc_b)*abs_uxc
      vxc   = (vxc)+dsign(1.d0,vxc_b)*abs_vxc
!
      uxmm =            dc12*dudx1     +dc12al*dudx2
      uymm =            dc12*dudy1     +dc12al*dudy2
      uxm  = (uxm)-8.d0*dc12*dudx1-4.d0*dc12al*dudx2
      uym  = (uym)-8.d0*dc12*dudy1-4.d0*dc12al*dudy2
      uxc  = (uxc)                +6.d0*dc12al*(dudx2+dudy2)
      uxp  = (uxp)+8.d0*dc12*dudx1-4.d0*dc12al*dudx2
      uyp  = (uyp)+8.d0*dc12*dudy1-4.d0*dc12al*dudy2
      uxpp =           -dc12*dudx1     +dc12al*dudx2
      uypp =           -dc12*dudy1     +dc12al*dudy2
!
      ustg(ic  ,jc  ) = (ustg(ic  ,jc  ))       +uxc
      vstg(ic  ,jc+1) = (vstg(ic  ,jc+1))+0.25d0*vxc
      vstg(ic  ,jc  ) = (vstg(ic  ,jc  ))+0.25d0*vxc
      vstg(ic-1,jc+1) = (vstg(ic-1,jc+1))+0.25d0*vxc
      vstg(ic-1,jc  ) = (vstg(ic-1,jc  ))+0.25d0*vxc
!
      ustg(ic-2,jc)=(ustg(ic-2,jc))+uxmm
      ustg(ic+2,jc)=(ustg(ic+2,jc))+uxpp
      ustg(ic,jc-2)=(ustg(ic,jc-2))+uymm
      ustg(ic,jc+2)=(ustg(ic,jc+2))+uypp
!
      ustg(ic-1,jc)=(ustg(ic-1,jc))+uxm
      ustg(ic+1,jc)=(ustg(ic+1,jc))+uxp
      ustg(ic,jc-1)=(ustg(ic,jc-1))+uym
      ustg(ic,jc+1)=(ustg(ic,jc+1))+uyp
!
    enddo
    enddo
!
  return
end subroutine sub_rhs_adj
!-----------------------------------------------------------------------

