!-----------------------------------------------------------------------
subroutine sub_bc_wall(ustg,vstg,pcnt)
!-----------------------------------------------------------------------
  use mod_variables,only : imax,jmax,imx1,jmx1,corner,iblk
  implicit none
!
  integer :: i,ic,jc,kc,ico,jco
  real(8) :: ustg(imax+1,jmax  )
  real(8) :: vstg(imax  ,jmax+1)
  real(8) :: pcnt(imax  ,jmax  )
!-----------------------------------------------------------------------
!
!
! Wall
  do jc=3,jmax-2
  do ic=3,imax-2
    if(iblk(ic,jc)==0) then
      if(iblk(ic-1,jc)==1 .or. iblk(ic-1,jc)==0) then
        ustg(ic  ,jc) = 0.d0
      endif
      if(iblk(ic+1,jc)==1 .or. iblk(ic+1,jc)==0) then
        ustg(ic+1,jc) = 0.d0
      endif
      if(iblk(ic,jc-1)==1 .or. iblk(ic,jc-1)==0) then
        vstg(ic,jc  ) = 0.d0
      endif
      if(iblk(ic,jc+1)==1 .or. iblk(ic,jc+1)==0) then
        vstg(ic,jc+1) = 0.d0
      endif
    endif
  enddo
  enddo
!
  do jc=3,jmax-2
  do ic=3,imax-2
    if(iblk(ic,jc)==0) then
      if(iblk(ic-1,jc)==0 .and. iblk(ic+1,jc)==0 .and. &
         iblk(ic,jc-1)==0 .and. iblk(ic,jc+1)==0) then
        pcnt(ic,jc) = 0.d0
      endif
    endif
  enddo
  enddo
!
! Wall boundaries in y-direction
  do jc=3,jmax-2
  do ic=3,imax-2
    if(iblk(ic,jc)==0) then
      if(iblk(ic,jc-1)==1) then
        ustg(ic,jc  ) =-ustg(ic,jc-1)
        vstg(ic,jc+1) =-vstg(ic,jc-1)
        ustg(ic,jc+1) =-ustg(ic,jc-2)
        vstg(ic,jc+2) =-vstg(ic,jc-2)
        pcnt(ic,jc  ) = pcnt(ic,jc-1)
      endif
      if(iblk(ic,jc+1)==1) then
        ustg(ic,jc  ) =-ustg(ic,jc+1)
        vstg(ic,jc  ) =-vstg(ic,jc+2)
        ustg(ic,jc-1) =-ustg(ic,jc+2)
        vstg(ic,jc-1) =-vstg(ic,jc+3)
        pcnt(ic,jc  ) = pcnt(ic,jc+1)
      endif
    endif
  enddo
  enddo
!
! Wall boundaries in x-direction
  do jc=3,jmax-2
  do ic=3,imax-2
    if(iblk(ic,jc)==0 .and. iblk(ic,jc-1)/=1 .and. iblk(ic,jc+1)/=1) then
      if(iblk(ic-1,jc)==1) then
        ustg(ic+1,jc) =-ustg(ic-1,jc)
        vstg(ic  ,jc) =-vstg(ic-1,jc)
        ustg(ic+2,jc) =-ustg(ic-2,jc)
        vstg(ic+1,jc) =-vstg(ic-2,jc)
        pcnt(ic  ,jc) = pcnt(ic-1,jc)
      endif
      if(iblk(ic+1,jc)==1) then
        ustg(ic  ,jc) =-ustg(ic+2,jc)
        vstg(ic  ,jc) =-vstg(ic+1,jc)
        ustg(ic-1,jc) =-ustg(ic+3,jc)
        vstg(ic-1,jc) =-vstg(ic+2,jc)
        pcnt(ic  ,jc) = pcnt(ic+1,jc)
      endif
    endif
  enddo
  enddo
!
!-Corners
  ico = corner(1,1)
  jco = corner(2,1)
  pcnt(ico,jco) = pcnt(ico-1,jco-1)
  ico = corner(1,2)
  jco = corner(2,2)
  pcnt(ico,jco) = pcnt(ico+1,jco-1)
  ico = corner(1,3)
  jco = corner(2,3)
  pcnt(ico,jco) = pcnt(ico-1,jco+1)
  ico = corner(1,4)
  jco = corner(2,4)
  pcnt(ico,jco) = pcnt(ico+1,jco+1)
!
  return
end subroutine sub_bc_wall
!-----------------------------------------------------------------------
!
!
!
!-----------------------------------------------------------------------
subroutine sub_bc_wall_adj(ustg,vstg,pcnt)
!-----------------------------------------------------------------------
  use mod_variables,only : imax,jmax,imx1,jmx1,corner,iblk
  implicit none
!
  integer :: i,ic,jc,kc,ico,jco
  real(8) :: ustg(imax+1,jmax  )
  real(8) :: vstg(imax  ,jmax+1)
  real(8) :: pcnt(imax  ,jmax  )
!-----------------------------------------------------------------------
!
!
!-Corners
  ico = corner(1,4)
  jco = corner(2,4)
  pcnt(ico+1,jco+1) = (pcnt(ico+1,jco+1))+pcnt(ico,jco)
  pcnt(ico  ,jco  ) = 0.d0
  ico = corner(1,3)
  jco = corner(2,3)
  pcnt(ico-1,jco+1) = (pcnt(ico-1,jco+1))+pcnt(ico,jco)
  pcnt(ico  ,jco  ) = 0.d0
  ico = corner(1,2)
  jco = corner(2,2)
  pcnt(ico+1,jco-1) = (pcnt(ico+1,jco-1))+pcnt(ico,jco)
  pcnt(ico  ,jco  ) = 0.d0
  ico = corner(1,1)
  jco = corner(2,1)
  pcnt(ico-1,jco-1) = (pcnt(ico-1,jco-1))+pcnt(ico,jco)
  pcnt(ico  ,jco  ) = 0.d0
!
!
! Wall boundaries in x-direction
  do jc=3,jmax-2
  do ic=3,imax-2
    if(iblk(ic,jc)==0 .and. iblk(ic,jc-1)/=1 .and. iblk(ic,jc+1)/=1) then
      if(iblk(ic-1,jc)==1) then
        ustg(ic-1,jc) = (ustg(ic-1,jc))-ustg(ic+1,jc)
        vstg(ic-1,jc) = (vstg(ic-1,jc))-vstg(ic  ,jc)
        ustg(ic-2,jc) = (ustg(ic-2,jc))-ustg(ic+2,jc)
        vstg(ic-2,jc) = (vstg(ic-2,jc))-vstg(ic+1,jc)
        pcnt(ic-1,jc) = (pcnt(ic-1,jc))+pcnt(ic  ,jc)
        ustg(ic+1,jc) = 0.d0
        vstg(ic  ,jc) = 0.d0
        ustg(ic+2,jc) = 0.d0
        vstg(ic+1,jc) = 0.d0
        pcnt(ic  ,jc) = 0.d0
      endif
      if(iblk(ic+1,jc)==1) then
        ustg(ic+2,jc) = (ustg(ic+2,jc))-ustg(ic  ,jc)
        vstg(ic+1,jc) = (vstg(ic+1,jc))-vstg(ic  ,jc)
        ustg(ic+3,jc) = (ustg(ic+3,jc))-ustg(ic-1,jc)
        vstg(ic+2,jc) = (vstg(ic+2,jc))-vstg(ic-1,jc)
        pcnt(ic+1,jc) = (pcnt(ic+1,jc))+pcnt(ic  ,jc)
        ustg(ic  ,jc) = 0.d0
        vstg(ic  ,jc) = 0.d0
        ustg(ic-1,jc) = 0.d0
        vstg(ic-1,jc) = 0.d0
        pcnt(ic  ,jc) = 0.d0
      endif
    endif
  enddo
  enddo
!
!
!-Wall boundaries in y-direction
  do jc=3,jmax-2
  do ic=3,imax-2
    if(iblk(ic,jc)==0) then
      if(iblk(ic,jc-1)==1) then
        ustg(ic,jc-1) = (ustg(ic,jc-1))-ustg(ic,jc  )
        vstg(ic,jc-1) = (vstg(ic,jc-1))-vstg(ic,jc+1)
        ustg(ic,jc-2) = (ustg(ic,jc-2))-ustg(ic,jc+1)
        vstg(ic,jc-2) = (vstg(ic,jc-2))-vstg(ic,jc+2)
        pcnt(ic,jc-1) = (pcnt(ic,jc-1))+pcnt(ic,jc  )
        ustg(ic,jc  ) = 0.d0
        vstg(ic,jc+1) = 0.d0
        ustg(ic,jc+1) = 0.d0
        vstg(ic,jc+2) = 0.d0
        pcnt(ic,jc  ) = 0.d0
      endif
      if(iblk(ic,jc+1)==1) then
        ustg(ic,jc+1) = (ustg(ic,jc+1))-ustg(ic,jc  )
        vstg(ic,jc+2) = (vstg(ic,jc+2))-vstg(ic,jc  )
        ustg(ic,jc+2) = (ustg(ic,jc+2))-ustg(ic,jc-1)
        vstg(ic,jc+3) = (vstg(ic,jc+3))-vstg(ic,jc-1)
        pcnt(ic,jc+1) = (pcnt(ic,jc+1))+pcnt(ic,jc  )
        ustg(ic,jc  ) = 0.d0
        vstg(ic,jc  ) = 0.d0
        ustg(ic,jc-1) = 0.d0
        vstg(ic,jc-1) = 0.d0
        pcnt(ic,jc  ) = 0.d0
      endif
    endif
  enddo
  enddo
!
!
! Wall
  do jc=3,jmax-2
  do ic=3,imax-2
    if(iblk(ic,jc)==0) then
      if(iblk(ic-1,jc)==1 .or. iblk(ic-1,jc)==0) then
        ustg(ic  ,jc) = 0.d0
      endif
      if(iblk(ic+1,jc)==1 .or. iblk(ic+1,jc)==0) then
        ustg(ic+1,jc) = 0.d0
      endif
      if(iblk(ic,jc-1)==1 .or. iblk(ic,jc-1)==0) then
        vstg(ic,jc  ) = 0.d0
      endif
      if(iblk(ic,jc+1)==1 .or. iblk(ic,jc+1)==0) then
        vstg(ic,jc+1) = 0.d0
      endif
    endif
  enddo
  enddo
!
  do jc=3,jmax-2
  do ic=3,imax-2
    if(iblk(ic,jc)==0) then
      if(iblk(ic-1,jc)==0 .and. iblk(ic+1,jc)==0 .and. &
         iblk(ic,jc-1)==0 .and. iblk(ic,jc+1)==0) then
        pcnt(ic,jc) = 0.d0
      endif
    endif
  enddo
  enddo
!
  return
end subroutine sub_bc_wall_adj
!-----------------------------------------------------------------------
