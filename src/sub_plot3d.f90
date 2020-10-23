!-----------------------------------------------------------------------
! output result file (plot3d format)
!-----------------------------------------------------------------------
subroutine sub_p3dwrite(ostep,tstep,header)
!-----------------------------------------------------------------------
  use mod_variables,only : imax,jmax,ustg,vstg,pcnt,divg, &
                           Re,istep,nstep,iskip_plot
  implicit none
!
  integer :: ic,jc,kc,m,ostep,ib,nblock
  real(8) :: tstep
  real(4) :: re4,step4,time4
  real(4) :: q4(imax,jmax,4)
  character(len=6) :: cicount = '000000'
  character(len=3) :: header
!-----------------------------------------------------------------------
!
!
  do jc=1,jmax
  do ic=1,imax
    q4(ic,jc,1) = real(divg(ic,jc))
    q4(ic,jc,2) = real(ustg(ic,jc)+ustg(ic+1,jc))*0.5d0
    q4(ic,jc,3) = real(vstg(ic,jc)+vstg(ic,jc+1))*0.5d0
    q4(ic,jc,4) = real(pcnt(ic,jc))
  enddo
  enddo
!
  write(cicount,'(i6.6)') ostep
  re4   = re
  step4 = ostep
  time4 = tstep 
!
  if(header(1:3)=='rst') then
    open(30,file='./output/restart.q',&
         form='unformatted',status='replace',action='write')
  else
    open(30,file='./output/'//trim(header)//'_'//cicount//'.q',&
         form='unformatted',status='replace',action='write')
  endif
  nblock = 1
  write(30) nblock
  write(30) (imax,jmax,ib=1,nblock)
  do ib = 1,nblock
    write(30) re4,step4,re4,time4
    write(30) (((q4(ic,jc,m),ic=1,imax),jc=1,jmax),m=1,4)
  enddo
  close(30)
!
  return
end subroutine sub_p3dwrite
!-----------------------------------------------------------------------
