!-----------------------------------------------------------------------
subroutine sub_bc_outer(ustg,vstg,pcnt)
!-----------------------------------------------------------------------
  use mod_variables,only : imax,jmax,imx1,jmx1,     &
                           u_inf,v_inf,p_inf,uinf,  &
                           dt,dx   
  implicit none
!
  integer :: ic,jc,kc,is,js,ks
  real(8) :: uave,dtdxi
  real(8) :: ustg(imax+1,jmax  )
  real(8) :: vstg(imax  ,jmax+1)
  real(8) :: pcnt(imax  ,jmax  )
!-----------------------------------------------------------------------
!
!
! Boundaries in y-direction: Slip boundary
  do ic=3,imax-1
    ustg(ic,1:2        ) = ustg(ic,3     )
    ustg(ic,jmax-1:jmax) = ustg(ic,jmax-2)
  enddo
  do ic=3,imax-2
    vstg(ic,1          ) =-vstg(ic,4     )
    vstg(ic,2          ) =-vstg(ic,3     )
    vstg(ic,jmx1       ) =-vstg(ic,jmx1-3)
    vstg(ic,jmx1-1     ) =-vstg(ic,jmx1-2)
  enddo
  do ic=2,imax-1
    pcnt(ic,1          ) = pcnt(ic,2     )
    pcnt(ic,jmax       ) = pcnt(ic,jmax-1)
  enddo
!
!
! Boundaries in x-direction
  do jc=1,jmax
    ustg(1:2        ,jc) = uinf(jc)
    ustg(imx1-1:imx1,jc) = ustg(imx1-2,jc)
    pcnt(1          ,jc) = pcnt(2     ,jc)
    pcnt(imax       ,jc) = pcnt(imax-1,jc)
  enddo
  do jc=1,jmx1
    vstg(1:2        ,jc) = 0.d0
    vstg(imax-1:imax,jc) = vstg(imax-2,jc)
  enddo
!
  return
end subroutine sub_bc_outer
!-----------------------------------------------------------------------
!
!
!
!-----------------------------------------------------------------------
subroutine sub_bc_outer_tlm(ustg,vstg,pcnt)
!-----------------------------------------------------------------------
  use mod_variables,only : imax,jmax,imx1,jmx1, &
                           u_inf,v_inf,p_inf,uinf,  &
                           dt,dx   
  implicit none
!
  integer :: ic,jc,kc,is,js,ks
  real(8) :: uave,dtdxi
  real(8) :: ustg(imax+1,jmax  )
  real(8) :: vstg(imax  ,jmax+1)
  real(8) :: pcnt(imax  ,jmax  )
!-----------------------------------------------------------------------
!
!
! Boundaries in y-direction: Slip boundary
  do ic=3,imax-1
    ustg(ic,1:2        ) = ustg(ic,3     )
    ustg(ic,jmax-1:jmax) = ustg(ic,jmax-2)
  enddo
  do ic=3,imax-2
    vstg(ic,1          ) =-vstg(ic,4     )
    vstg(ic,2          ) =-vstg(ic,3     )
    vstg(ic,jmx1       ) =-vstg(ic,jmx1-3)
    vstg(ic,jmx1-1     ) =-vstg(ic,jmx1-2)
  enddo
  do ic=2,imax-1
    pcnt(ic,1          ) = pcnt(ic,2     )
    pcnt(ic,jmax       ) = pcnt(ic,jmax-1)
  enddo
!
!
! Boundaries in x-direction
  do jc=1,jmax
    ustg(1:2        ,jc) = 0.d0
    ustg(imx1-1:imx1,jc) = ustg(imx1-2,jc)
    pcnt(1          ,jc) = pcnt(2     ,jc)
    pcnt(imax       ,jc) = pcnt(imax-1,jc)
  enddo
  do jc=1,jmx1
    vstg(1:2        ,jc) = 0.d0
    vstg(imax-1:imax,jc) = vstg(imax-2,jc)
  enddo
!!
  return
end subroutine sub_bc_outer_tlm
!-----------------------------------------------------------------------
!
!
!
!-----------------------------------------------------------------------
subroutine sub_bc_outer_adj(ustg,vstg,pcnt)
!-----------------------------------------------------------------------
  use mod_variables,only : imax,jmax,imx1,jmx1, &
                           u_inf,v_inf,p_inf,uinf,  &
                           dt,dx   
  implicit none
!
  integer :: ic,jc,kc,is,js,ks
  real(8) :: uave,dtdxi
  real(8) :: ustg(imax+1,jmax  )
  real(8) :: vstg(imax  ,jmax+1)
  real(8) :: pcnt(imax  ,jmax  )
!-----------------------------------------------------------------------
!
!
! Boundaries in x-direction
  do jc=1,jmax
    ustg(1:2   ,jc) = 0.d0
    ustg(imx1-2,jc) = (ustg(imx1-2,jc))+ustg(imx1-1,jc) 
    ustg(imx1-2,jc) = (ustg(imx1-2,jc))+ustg(imx1  ,jc) 
    ustg(imx1  ,jc) = 0.d0
    ustg(imx1-1,jc) = 0.d0
    pcnt(2     ,jc) = (pcnt(2     ,jc))+pcnt(1     ,jc)
    pcnt(1     ,jc) = 0.d0
    pcnt(imax-1,jc) = (pcnt(imax-1,jc))+pcnt(imax  ,jc)
    pcnt(imax  ,jc) = 0.d0
  enddo
  do jc=1,jmx1
    vstg(1:2   ,jc) = 0.d0
    vstg(imax-2,jc) = (vstg(imax-2,jc))+vstg(imax-1,jc) 
    vstg(imax-2,jc) = (vstg(imax-2,jc))+vstg(imax  ,jc) 
    vstg(imax-1,jc) = 0.d0
    vstg(imax  ,jc) = 0.d0
  enddo
!
!
! Boundaries in y-direction: Slip boundary
  do ic=3,imax-1
    ustg(ic,3          ) = (ustg(ic,3     ))+ustg(ic,1     )
    ustg(ic,3          ) = (ustg(ic,3     ))+ustg(ic,2     )
    ustg(ic,1:2        ) = 0.d0
    ustg(ic,jmax-2     ) = (ustg(ic,jmax-2))+ustg(ic,jmax-1) 
    ustg(ic,jmax-2     ) = (ustg(ic,jmax-2))+ustg(ic,jmax  ) 
    ustg(ic,jmax-1:jmax) = 0.d0
  enddo
  do ic=3,imax-2
    vstg(ic,3          ) = (vstg(ic,3     ))-vstg(ic,2     )
    vstg(ic,4          ) = (vstg(ic,4     ))-vstg(ic,1     )
    vstg(ic,1:2        ) = 0.d0
    vstg(ic,jmx1-3     ) = (vstg(ic,jmx1-3))-vstg(ic,jmx1  )
    vstg(ic,jmx1-2     ) = (vstg(ic,jmx1-2))-vstg(ic,jmx1-1)
    vstg(ic,jmx1-1:jmx1) = 0.d0
  enddo
  do ic=2,imax-1
    pcnt(ic,2          ) = (pcnt(ic,2     ))+pcnt(ic,1     ) 
    pcnt(ic,1          ) = 0.d0
    pcnt(ic,jmax-1     ) = (pcnt(ic,jmax-1))+pcnt(ic,jmax  ) 
    pcnt(ic,jmax       ) = 0.d0
  enddo
!
  return
end subroutine sub_bc_outer_adj
!-----------------------------------------------------------------------