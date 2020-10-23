!-----------------------------------------------------------------------
subroutine sub_HSMAC_adj(ustg,vstg,pcnt)
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
!   divmax = 0.d0
!
    do jc=jmax-1,2,-1
    do ic=imax-1,2,-1
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
!-----velocity correction
      dp          =     -ustg(ic,jc)*ibl_xm
      dp          = (dp)+ustg(ip,jc)*ibl_xp
      dp          = (dp)-vstg(ic,jc)*ibl_ym
      dp          = (dp)+vstg(ic,jp)*ibl_yp
!
!-----pressure correction
      dp          = dp*dt*dxi
!
!-----pressure update
      dp          = (dp)+pcnt(ic,jc)
      dive        = C0*dp
!
!---> Divergence
      UD1X_TMP    = dive*dxi*ibl
      VD1Y_TMP    = dive*dxi*ibl
!
!---- U
      ustg(ic,jc) = (ustg(ic,jc))-UD1X_TMP
      ustg(ip,jc) = (ustg(ip,jc))+UD1X_TMP
!---- V
      vstg(ic,jc) = (vstg(ic,jc))-VD1Y_TMP
      vstg(ic,jp) = (vstg(ic,jp))+VD1Y_TMP
!
!     divmax         = max(divmax,abs(divg(ic,jc))) 
!
    enddo
    enddo
!
    call sub_bc_outer_adj(ustg,vstg,pcnt)
    call sub_bc_wall_adj (ustg,vstg,pcnt)
!
!   if(divmax<epsmax) exit
!
  enddo
!
  return
end subroutine sub_HSMAC_adj
!-----------------------------------------------------------------------
