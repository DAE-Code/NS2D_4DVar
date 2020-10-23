!-----------------------------------------------------------------------
subroutine sub_measure
!-----------------------------------------------------------------------
  use mod_variables,only : imax,jmax,imx1,jmx1,nstep,n_prb,n_mes,variances,  &
                           xcen,ycen,ih,ihx,ihy,iskip,jskip,tskip,mbot,mtop, &
                           mlft,mrht
  implicit none
!
  integer :: ic,jc,istep
  integer :: iht(imax,jmax)
!-----------------------------------------------------------------------
!
!
  iht(:,:) = 0
!
  do jc=1,jmax,iskip
  do ic=1,imax,jskip
    if(mlft<xcen(ic) .and. xcen(ic)<mrht .and. mbot<ycen(jc) .and. ycen(jc)<mtop) then
      iht(ic,jc) = 1
    endif
  enddo
  enddo
!
  allocate(ihx(imx1,jmax,nstep),ihy(imax,jmx1,nstep))
  ihx(:,:,:) = 0
  ihy(:,:,:) = 0
!
  do istep=1,nstep
    if(mod(istep,tskip)==0 .or. istep==1) then
      do jc=1,jmax
      do ic=1,imax
        if(iht(ic,jc)==1) then
          ihx(ic  ,jc,istep) = 1
          ihx(ic+1,jc,istep) = 1
          ihy(ic,jc  ,istep) = 1
          ihy(ic,jc+1,istep) = 1
        endif
      enddo
      enddo
      n_mes = 0
      do jc=1,jmax
      do ic=1,imx1
        if(ihx(ic,jc,istep)==1) then
          n_mes = n_mes+1
        endif
      enddo
      enddo
      do jc=1,jmx1
      do ic=1,imax
        if(ihy(ic,jc,istep)==1) then
          n_mes = n_mes+1
        endif
      enddo
      enddo
    endif
  enddo
  write(*,*)
  write(*,*) "Number of measurement points:",n_mes
!
  return
end subroutine sub_measure
!-----------------------------------------------------------------------
