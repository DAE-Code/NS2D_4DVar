!-----------------------------------------------------------------------
subroutine sub_HSMAC_tlm(ustg,vstg,pcnt)
!-----------------------------------------------------------------------
  use mod_variables,only : imax,jmax,dx,divmax,itrp,divg,iblk,dt
  implicit none
!
!-----------------------------------------------------------------------
  integer :: itrmax = 10     ! total # of iteration
  real(8) :: epsmax = 1.d-12 ! max tolerance
  real(8) :: omg    = 1.5d0  ! accel. coef.
!-----------------------------------------------------------------------
!
  integer :: ic,jc,im,jm,ip,jp
  real(8) :: ibl,dp,dive,dtdx
  real(8) :: ibl_xm,ibl_xp,ibl_ym,ibl_yp
  real(8) :: UD1X_TMP,VD1Y_TMP
  real(8) :: C0,dxi
  real(8) :: ustg(imax+1,jmax  )
  real(8) :: vstg(imax  ,jmax+1)
  real(8) :: pcnt(imax  ,jmax  )
!-----------------------------------------------------------------------
!
  C0  =-omg/(2.d0*dt*(2.d0/(dx**2)))
  dxi = 1.d0/dx
!
  do itrp=1,itrmax
!
    call sub_bc_wall     (ustg,vstg,pcnt)
    call sub_bc_outer_tlm(ustg,vstg,pcnt)
!
!   divmax = 0.0
!
    do jc=2,jmax-1
    do ic=2,imax-1
!
      im             = ic-1
      jm             = jc-1
      ip             = ic+1
      jp             = jc+1
!
      ibl            = dble(iblk(ic,jc))
      ibl_xm         = dble(iblk(im,jc))
      ibl_xp         = dble(iblk(ip,jc))
      ibl_ym         = dble(iblk(ic,jm))
      ibl_yp         = dble(iblk(ic,jp))
!
!---- U
      UD1X_TMP       = (ustg(ip,jc)-ustg(ic,jc))
!---- V
      VD1Y_TMP       = (vstg(ic,jp)-vstg(ic,jc))
!
!---> Divergence
      dive           = (UD1X_TMP+VD1Y_TMP)*dxi*ibl
      dp             = C0*dive
      pcnt(ic,jc)    = pcnt(ic,jc)+dp
      divg(ic,jc)    = dive
!
      dp             = dp*dt*dxi
      ustg(ic,jc)    = ustg(ic,jc)-dp*ibl_xm
      ustg(ip,jc)    = ustg(ip,jc)+dp*ibl_xp
      vstg(ic,jc)    = vstg(ic,jc)-dp*ibl_ym
      vstg(ic,jp)    = vstg(ic,jp)+dp*ibl_yp
!
      divmax         = max(divmax,abs(divg(ic,jc))) 
!
    enddo
    enddo
!
!   if(divmax<epsmax) exit
!
  enddo
!
  return
end subroutine sub_HSMAC_tlm
!-----------------------------------------------------------------------
