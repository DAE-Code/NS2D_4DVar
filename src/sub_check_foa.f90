!-----------------------------------------------------------------------
  module mod_temp
!-----------------------------------------------------------------------
!
    integer,parameter :: iaverage = 10
    integer,parameter :: ialpha   = 14
    integer,parameter :: nva      =  3
    integer           :: ialp,iave,nc,icc,jcc
    real(8)           :: alp,dif,den,ang,random,fLHS,fRHS,fJ1,fJ2,foa,soa
    real(8)           :: tlv_max,tlv_min,tlv_sum
    real(8)           :: tla_max,tla_min,tla_sum
    real(8)           :: adj_max,adj_min,adj_sum
    real(8)           :: xmin,xmax,ymin,ymax
    real(8)           :: ibl,xibl,yibl,zibl
    real(8),dimension(:,:,:),allocatable :: urec,vrec,prec
    real(8),dimension(:,:,:),allocatable :: uref,vref,pref
    real(8),dimension(:,:,:),allocatable :: ufrc,vfrc,pfrc
    real(8),dimension(:,:,:),allocatable :: utlm,vtlm,ptlm
    real(8),dimension(:,:,:),allocatable :: uadj,vadj,padj
    real(8),dimension(:,:,:),allocatable :: uwns,vwns
!
!-- Adjoint check
    real(8),dimension(:,:,:),allocatable :: PQ1,PQ2,DLT,DSP,QIN
    real(8),dimension(:,:,:),allocatable :: rslt
!
  end module mod_temp
!-----------------------------------------------------------------------
!
!
!
!-----------------------------------------------------------------------
subroutine sub_check_foa
!-----------------------------------------------------------------------
  use mod_variables
  use mod_temp
  implicit none
!
  integer :: ic,jc,kc
!-----------------------------------------------------------------------
!
!
!-Domain bounds
  xmin = minval(xcen(:))
  xmax = maxval(xcen(:))
  ymin = minval(ycen(:))
  ymax = maxval(ycen(:))
!
!-Sine wave forms for a base flow field
  do jc=1,jmax
  do ic=1,imx1
    icc = min0(ic,imax)
    ustg(ic,jc) = ustg(ic,jc)                                  &
                 +(dsin(2.d0*3.141592d0*xcen(icc)/(xmax-xmin)) &
                  +dsin(2.d0*3.141592d0*ycen(jc )/(ymax-ymin)))*0.01d0
  enddo
  enddo
  do jc=1,jmx1
  do ic=1,imax
    jcc = min0(jc,jmax)
    vstg(ic,jc) = vstg(ic,jc)                                  &
                 +(dsin(2.d0*3.141592d0*xcen(jc )/(xmax-xmin)) &
                  +dsin(2.d0*3.141592d0*ycen(jcc)/(ymax-ymin)))*0.01d0
  enddo
  enddo
  do jc=1,jmax
  do ic=1,imax
    pcnt(ic,jc) = 0.d0
  enddo
  enddo
!
!-Variables for adjoint development
  allocate(PQ1(imx1,jmx1,nva),PQ2(imx1,jmx1,nva), &
           DLT(imx1,jmx1,nva),DSP(imx1,jmx1,nva), &
           QIN(imx1,jmx1,nva))
!
  do nc=1,nva
  do jc=1,jmx1
  do ic=1,imx1
    PQ1(ic,jc,nc) = 0.d0
    PQ2(ic,jc,nc) = 0.d0
    DLT(ic,jc,nc) = 0.d0
    DSP(ic,jc,nc) = 0.d0
    QIN(ic,jc,nc) = 0.d0
  enddo
  enddo
  enddo
  allocate(rslt(iaverage,ialpha,10))
  rslt(1:iaverage,1:ialpha,1:10) = 0.d0
!
!-Save intial field
  QIN(1:imx1,1:jmax,1) = ustg(1:imx1,1:jmax)
  QIN(1:imax,1:jmx1,2) = vstg(1:imax,1:jmx1)
  QIN(1:imax,1:jmax,3) = pcnt(1:imax,1:jmax)
!
! goto 1000
! goto 2000
!
!======================================================================
!===== Tangent linear check ===========================================
!======================================================================
  do iave=1,iaverage
!
!---Generate randum displacement
    do nc=1,nva
    do jc=1,jmx1
    do ic=1,imx1
      call random_number(random)
!     DSP(ic,jc,nc) = (2.d0*random-1.d0)*1.d0
      DSP(ic,jc,nc) = random*100.d0
    enddo
    enddo
    enddo
!   do jc=3,jmx1-2
!   do ic=3,imx1-2
!     call random_number(random)
!     if(jc==jmx1-2) goto 15
!     xibl = dble(min0(iblk(ic-1,jc),iblk(ic,jc)))
!     DSP(ic,jc,1) = (2.d0*random-1.d0)*0.1d0*xibl
!     QIN(ic,jc,1) = QIN(ic,jc,1)*xibl
!  15 continue
!     if(ic==imx1-2) goto 25
!     yibl = dble(min0(iblk(ic,jc-1),iblk(ic,jc)))
!     DSP(ic,jc,2) = (2.d0*random-1.d0)*0.1d0*yibl
!     QIN(ic,jc,2) = QIN(ic,jc,2)*yibl
!  25 continue
!   enddo
!   enddo
!
!---Sine wave forms
!   DSP(:,:,:) = 0.d0
!   do nc=1,nva
!   do jc=1,jmx1
!   do ic=1,imx1
!     DSP(ic,jc,nc) = (dsin(12.d0*3.141592d0*(xcen(ic))/(xmax-xmin)) &
!                     +dsin(12.d0*3.141592d0*(ycen(jc))/(ymax-ymin)))/2.d0
!   enddo
!   enddo
!   enddo
!
! 
  alp = 1.d+2
  do ialp=1,ialpha
    alp = alp*0.1d0
!
!
!---Vector 1
    ustg(1:imx1,1:jmax) = QIN(1:imx1,1:jmax,1)
    vstg(1:imax,1:jmx1) = QIN(1:imax,1:jmx1,2)
    pcnt(1:imax,1:jmax) = 0.d0
!
!rhs
!   call sub_rhs_fwd
!prs
!   call sub_HSMAC_fwd(ustg,vstg,pcnt)
!full
    call sub_rhs_fwd
    call sub_bc_outer(ustg,vstg,pcnt)
    call sub_bc_wall (ustg,vstg,pcnt)
    call sub_HSMAC_fwd(ustg,vstg,pcnt)
!
    PQ1(1:imx1,1:jmax,1) = ustg(1:imx1,1:jmax)
    PQ1(1:imax,1:jmx1,2) = vstg(1:imax,1:jmx1)
    PQ1(1:imax,1:jmax,3) = pcnt(1:imax,1:jmax)
!
!
!---Vector 2
    ustg(1:imx1,1:jmax) = QIN(1:imx1,1:jmax,1)+DSP(1:imx1,1:jmax,1)*alp
    vstg(1:imax,1:jmx1) = QIN(1:imax,1:jmx1,2)+DSP(1:imax,1:jmx1,2)*alp
    pcnt(1:imax,1:jmax) =                      DSP(1:imax,1:jmax,3)*alp
!
!rhs
!   call sub_rhs_fwd
!prs
!   call sub_HSMAC_fwd(ustg,vstg,pcnt)
!full
    call sub_rhs_fwd
    call sub_bc_outer(ustg,vstg,pcnt)
    call sub_bc_wall (ustg,vstg,pcnt)
    call sub_HSMAC_fwd(ustg,vstg,pcnt)
!
    PQ2(1:imx1,1:jmax,1) = ustg(1:imx1,1:jmax)
    PQ2(1:imax,1:jmx1,2) = vstg(1:imax,1:jmx1)
    PQ2(1:imax,1:jmax,3) = pcnt(1:imax,1:jmax)
!
!
!---dy  
    dif = 0.d0
    do jc=3,jmax-2
    do ic=3,imx1-2
      xibl = dble(min0(iblk(ic-1,jc),iblk(ic,jc)))
      DLT(ic,jc,1) = PQ2(ic,jc,1)-PQ1(ic,jc,1)
      dif = dif+DLT(ic,jc,1)**2 !*xibl
    enddo
    enddo
    do jc=3,jmx1-2
    do ic=3,imax-2
      yibl = dble(min0(iblk(ic,jc-1),iblk(ic,jc)))
      DLT(ic,jc,2) = PQ2(ic,jc,2)-PQ1(ic,jc,2)
      dif = dif+DLT(ic,jc,2)**2 !*yibl
    enddo
    enddo
    do jc=3,jmax-2
    do ic=3,imax-2
      ibl             = dble(iblk(ic,jc))
      DLT(ic,jc,3) = PQ2(ic,jc,3)-PQ1(ic,jc,3)
      dif = dif+DLT(ic,jc,3)**2 !*ibl
    enddo
    enddo
    dif = dsqrt(dif)
!
!
!---dq
    ustg_b(1:imx1,1:jmax) = QIN(1:imx1,1:jmax,1)
    vstg_b(1:imax,1:jmx1) = QIN(1:imax,1:jmx1,2)
    pcnt_b(1:imax,1:jmax) = 0.d0
    ustg  (1:imx1,1:jmax) = DSP(1:imx1,1:jmax,1)
    vstg  (1:imax,1:jmx1) = DSP(1:imax,1:jmx1,2)
    pcnt  (1:imax,1:jmax) = DSP(1:imax,1:jmax,3)
!
!rhs
!   call sub_rhs_tlm
!prs
!   call sub_HSMAC_fwd(ustg,vstg,pcnt)
!full -- the order may be important for IBM 
!     -- (the base field should not be altered before the usage??)
    call sub_rhs_tlm
    call sub_bc_outer_tlm(ustg  ,vstg  ,pcnt  )
!!  call sub_bc_outer    (ustg_b,vstg_b,pcnt_b)  ! unchange the base 
    call sub_bc_wall(ustg  ,vstg  ,pcnt  )
!!  call sub_bc_wall(ustg_b,vstg_b,pcnt_b)       ! unchange the base 
    call sub_HSMAC_tlm(ustg,vstg,pcnt)
!
    den = 0.d0
    do jc=3,jmax-2
    do ic=3,imx1-2
      xibl = dble(min0(iblk(ic-1,jc),iblk(ic,jc)))
      den  = den+ustg(ic,jc)**2 !*xibl
    enddo
    enddo
    do jc=3,jmx1-2
    do ic=3,imax-2
      yibl = dble(min0(iblk(ic,jc-1),iblk(ic,jc)))
      den  = den+vstg(ic,jc)**2 !*yibl
    enddo
    enddo
    do jc=3,jmax-2
    do ic=3,imax-2
      ibl  = dble(iblk(ic,jc))
      den  = den+pcnt(ic,jc)**2 !*ibl
    enddo
    enddo
    den = dsqrt(den)
!
    ang = 0.d0
    do jc=3,jmax-2
    do ic=3,imx1-2
      xibl = dble(min0(iblk(ic-1,jc),iblk(ic,jc)))
      ang = ang+ustg(ic,jc)*DLT(ic,jc,1) !*xibl
    enddo
    enddo
    do jc=3,jmx1-2
    do ic=3,imax-2
      yibl = dble(min0(iblk(ic,jc-1),iblk(ic,jc)))
      ang = ang+vstg(ic,jc)*DLT(ic,jc,2) !*yibl
    enddo
    enddo
    do jc=3,jmax-2
    do ic=3,imax-2
      ibl  = dble(iblk(ic,jc))
      ang = ang+pcnt(ic,jc)*DLT(ic,jc,3) !*ibl
    enddo
    enddo
    ang = dabs(ang)
!
    rslt(iave,ialp,1) = dif/(alp*den)-1.d0
    rslt(iave,ialp,2) = ang/(dif*den)-1.d0
    rslt(iave,ialp,4) = alp
!
    write(*,'(2i7,6e20.12)') ialp,iave,dif,alp*den,ang,dif*den, &
                             rslt(iave,ialp,1),rslt(iave,ialp,2)
!
    enddo
  enddo
!
  write(*,*) ' -- Norm --'
  do ialp=1,ialpha
    tlv_max = -1.d3
    tlv_min =  1.d3
    tlv_sum =  0.d0
    tla_max = -1.d3
    tla_min =  1.d3
    tla_sum =  0.d0
    do iave=1,iaverage
      tlv_sum = tlv_sum+rslt(iave,ialp,1)
      tlv_max = dmax1(rslt(iave,ialp,1),tlv_max)
      tlv_min = dmin1(rslt(iave,ialp,1),tlv_min)
    enddo
    write(*  ,'(3e15.7,e32.25)') rslt(1,ialp,4),tlv_max,tlv_min, &
                                 tlv_sum/dble(iaverage)
    write(100,'(3e15.7,e32.25)') rslt(1,ialp,4),tlv_max,tlv_min, &
                                 tlv_sum/dble(iaverage)
  enddo
  write(*,*) ' -- Angle --'
  do ialp=1,ialpha
    tlv_max = -1.d3
    tlv_min =  1.d3
    tlv_sum =  0.d0
    tla_max = -1.d3
    tla_min =  1.d3
    tla_sum =  0.d0
    do iave=1,iaverage
      tla_sum = tla_sum+rslt(iave,ialp,2)
      tla_max = dmax1(rslt(iave,ialp,2),tla_max)
      tla_min = dmin1(rslt(iave,ialp,2),tla_min)
    enddo
    write(*  ,'(3e15.7,e32.25)') rslt(1,ialp,4),tla_max,tla_min, &
                                 tla_sum/dble(iaverage)
    write(200,'(3e15.7,e32.25)') rslt(1,ialp,4),tla_max,tla_min, &
                                 tla_sum/dble(iaverage)
  enddo
!
!
1000 continue
!======================================================================
!===== adjoint check ==================================================
!======================================================================
  write(*  ,*  )
  write(*  ,*  ) "adjoint inner product"
  allocate(urec(1:imx1,1:jmx1,nstep), &
           vrec(1:imx1,1:jmx1,nstep), &
           prec(1:imx1,1:jmx1,nstep))
  urec(1:imx1,1:jmx1,:) = 0.d0
  vrec(1:imx1,1:jmx1,:) = 0.d0
  prec(1:imx1,1:jmx1,:) = 0.d0
  do iave=1,iaverage
!
!---Generate randum displacement
    do nc=1,nva
    do jc=1,jmx1
    do ic=1,imx1
      call random_number(random)
      DSP(ic,jc,nc) = (2.d0*random-1.d0)*0.1d0
    enddo
    enddo
    enddo
!    do jc=3,jmx1-2
!    do ic=3,imx1-2
!      call random_number(random)
!      if(jc==jmx1-2) goto 10
!      xibl = real(min0(iblk(ic-1,jc),iblk(ic,jc)))
!      DSP(ic,jc,1) = (2.d0*random-1.d0)*0.1d0*xibl
!      QIN(ic,jc,1) = QIN(ic,jc,1)*xibl
!   10 continue
!      if(ic==imx1-2) goto 20
!      yibl = real(min0(iblk(ic,jc-1),iblk(ic,jc)))
!      DSP(ic,jc,2) = (2.d0*random-1.d0)*0.1d0*yibl
!      QIN(ic,jc,2) = QIN(ic,jc,2)*yibl
!   20 continue
!      if(ic==imx1 .or. jc==jmx1) goto 30
!      ibl = real(iblk(ic,jc))
!      DSP(ic,jc,4) = (2.d0*random-1.d0)*0.1d0*ibl
!      QIN(ic,jc,4) = 0.d0
!   30 continue
!    enddo
!    enddo
!
!
!---Tangent linear
    ustg_b(1:imx1,1:jmax) = QIN(1:imx1,1:jmax,1)
    vstg_b(1:imax,1:jmx1) = QIN(1:imax,1:jmx1,2)
    pcnt_b(1:imax,1:jmax) = 0.d0
!   udlt  (1:imx1,1:jmax) = 0.d0
!   vdlt  (1:imax,1:jmx1) = 0.d0
!   uold  (1:imx1,1:jmax) = 0.d0
!   vold  (1:imax,1:jmx1) = 0.d0
    ustg  (1:imx1,1:jmax) = DSP(1:imx1,1:jmax,1)
    vstg  (1:imax,1:jmx1) = DSP(1:imax,1:jmx1,2)
    pcnt  (1:imax,1:jmax) = DSP(1:imax,1:jmax,3)
!
    do istep=1,nstep
      urec(1:imx1,1:jmax,istep) = ustg_b(1:imx1,1:jmax)
      vrec(1:imax,1:jmx1,istep) = vstg_b(1:imax,1:jmx1)
      prec(1:imax,1:jmax,istep) = pcnt_b(1:imax,1:jmax)
      call sub_bc_outer    (ustg_b,vstg_b,pcnt_b)
      call sub_bc_outer_tlm(ustg  ,vstg  ,pcnt  )
      call sub_bc_wall (ustg_b,vstg_b,pcnt_b)
      call sub_bc_wall (ustg  ,vstg  ,pcnt  )
      call sub_rhs_tlm
      call sub_HSMAC_fwd(ustg_b,vstg_b,pcnt_b)
      call sub_HSMAC_tlm(ustg  ,vstg  ,pcnt  )
    enddo
!   call sub_HSMAC_fwd(ustg,vstg,pcnt)
!
    fLHS = 0.d0
    do jc=1,jmax
    do ic=1,imx1
      xibl = dble(min0(iblk(max0(ic-1,1),jc),iblk(min0(ic,imax),jc)))
      fLHS = fLHS+ustg(ic,jc)**2 !*xibl
!     ustg(ic,jc) = ustg(ic,jc)**2
    enddo
    enddo
    do jc=1,jmx1
    do ic=1,imax
      yibl = dble(min0(iblk(ic,max0(jc-1,1)),iblk(ic,min0(jc,jmax))))
      fLHS = fLHS+vstg(ic,jc)**2 !*yibl
!     vstg(ic,jc) = vstg(ic,jc)**2
    enddo
    enddo
    do jc=1,jmax
    do ic=1,imax
      ibl  = dble(iblk(ic,jc))
      fLHS = fLHS+pcnt(ic,jc)**2 !*ibl
!     pcnt(ic,jc) = pcnt(ic,jc)**2
    enddo
    enddo
!
!   call sub_p3dwrite(1,1.d0,'chk')
!
!---Adjoint
    ustg_b(1:imx1,1:jmax) = QIN(1:imx1,1:jmax,1)
    vstg_b(1:imax,1:jmx1) = QIN(1:imax,1:jmx1,2)
    pcnt_b(1:imax,1:jmax) = 0.d0
!   udlt  (1:imx1,1:jmax) = 0.d0
!   vdlt  (1:imax,1:jmx1) = 0.d0
!   uold  (1:imx1,1:jmax) = 0.d0
!   vold  (1:imax,1:jmx1) = 0.d0
!   ustg  (1:imx1,1:jmax) = 0.d0
!   vstg  (1:imax,1:jmx1) = 0.d0
!   pcnt  (1:imax,1:jmax) = 0.d0
!
    do istep=nstep,1,-1
      ustg_b(1:imx1,1:jmax) = urec(1:imx1,1:jmax,istep) 
      vstg_b(1:imax,1:jmx1) = vrec(1:imax,1:jmx1,istep) 
      pcnt_b(1:imax,1:jmax) = prec(1:imax,1:jmax,istep) 
      call sub_bc_wall     (ustg_b,vstg_b,pcnt_b)
      call sub_bc_outer    (ustg_b,vstg_b,pcnt_b)
      call sub_HSMAC_adj   (ustg  ,vstg  ,pcnt  )
      call sub_rhs_adj
      call sub_bc_wall_adj (ustg  ,vstg  ,pcnt  )
      call sub_bc_outer_adj(ustg  ,vstg  ,pcnt  )
    enddo
!   call sub_HSMAC_adj(ustg,vstg,pcnt)
!
    fRHS = 0.d0
    do jc=1,jmax
    do ic=1,imx1
      xibl = dble(min0(iblk(max0(ic-1,1),jc),iblk(min0(ic,imax),jc)))
      fRHS = fRHS+ustg(ic,jc)*DSP(ic,jc,1) !*xibl
      ustg(ic,jc) = ustg(ic,jc)*DSP(ic,jc,1)
    enddo
    enddo
    do jc=1,jmx1
    do ic=1,imax
      yibl = dble(min0(iblk(ic,max0(jc-1,1)),iblk(ic,min0(jc,jmax))))
      fRHS = fRHS+vstg(ic,jc)*DSP(ic,jc,2) !*yibl
      vstg(ic,jc) = vstg(ic,jc)*DSP(ic,jc,2)
    enddo
    enddo
    do jc=1,jmax
    do ic=1,imax
      ibl  = dble(iblk(ic,jc))
      fRHS = fRHS+pcnt(ic,jc)*DSP(ic,jc,3) !*ibl
      pcnt(ic,jc) = pcnt(ic,jc)*DSP(ic,jc,3)
    enddo
    enddo
!
!   call sub_p3dwrite(2,2.d0,'chk')
!
    write(*  ,'(3e32.25)') (fRHS-fLHS)/fLHS,fRHS,fLHS
    write(300,'(3e32.25)') (fRHS-fLHS)/fLHS,fRHS,fLHS
!
  enddo
!
!
 2000 continue
!======================================================================
!===== gradient check =================================================
!======================================================================
  write(*  ,*  )
  write(*  ,*  ) "FOA gradient"
  if(.not.allocated(urec)) then
    allocate(urec(1:imx1,1:jmx1,nstep), &
             vrec(1:imx1,1:jmx1,nstep), &
             prec(1:imx1,1:jmx1,nstep))
  endif
  urec(1:imx1,1:jmx1,:) = 0.d0
  vrec(1:imx1,1:jmx1,:) = 0.d0
  prec(1:imx1,1:jmx1,:) = 0.d0
  allocate(uref(1:imx1,1:jmx1,nstep), &
           vref(1:imx1,1:jmx1,nstep), &
           pref(1:imx1,1:jmx1,nstep))
  uref(1:imx1,1:jmx1,:) = 0.d0
  vref(1:imx1,1:jmx1,:) = 0.d0
  pref(1:imx1,1:jmx1,:) = 0.d0
  allocate(ufrc(1:imx1,1:jmx1,nstep), &
           vfrc(1:imx1,1:jmx1,nstep), &
           pfrc(1:imx1,1:jmx1,nstep))
  ufrc(1:imx1,1:jmx1,:) = 0.d0
  vfrc(1:imx1,1:jmx1,:) = 0.d0
  pfrc(1:imx1,1:jmx1,:) = 0.d0
!
!-
! do jc=3,jmx1-2
! do ic=3,imx1-2
!   if(jc==jmx1-2) goto 11
!   xibl = real(min0(iblk(ic-1,jc),iblk(ic,jc)))
!   QIN(ic,jc,1) = QIN(ic,jc,1)*xibl
!11 continue
!   if(ic==imx1-2) goto 21
!   yibl = real(min0(iblk(ic,jc-1),iblk(ic,jc)))
!   QIN(ic,jc,2) = QIN(ic,jc,2)*yibl
!21 continue
! enddo
! enddo
!
!-Reference flow field (as measurement data)
  ustg(1:imx1,1:jmax) = QIN(1:imx1,1:jmax,1)-1.d-3
  vstg(1:imax,1:jmx1) = QIN(1:imax,1:jmx1,2)
  pcnt(1:imax,1:jmax) = 0.d0
  do istep=1,nstep
    uref(1:imx1,1:jmax,istep) = ustg(1:imx1,1:jmax)
    vref(1:imax,1:jmx1,istep) = vstg(1:imax,1:jmx1)
    pref(1:imax,1:jmax,istep) = pcnt(1:imax,1:jmax)
    call sub_bc_outer (ustg,vstg,pcnt)
    call sub_bc_wall  (ustg,vstg,pcnt)
    call sub_rhs_fwd
    call sub_HSMAC_fwd(ustg,vstg,pcnt)
  enddo
!
!-
  do iave=1,iaverage
!
!---Generate random displacement
!   do nc=1,nva
!   do jc=1,jmx1
!   do ic=1,imx1
!     call random_number(random)
!!    DSP(ic,jc,nc) = (2.d0*random-1.d0)*0.1d0
!     DSP(ic,jc,nc) = random
!   enddo
!   enddo
!   enddo
    do jc=1,jmax
    do ic=1,imax
!     call random_number(random)
      DSP(ic,jc,1) =((dsin(8.d0*3.141592d0*xcen(ic)/(xmax-xmin)) &
                     +dsin(8.d0*3.141592d0*ycen(jc)/(ymax-ymin)))+1.d0)*0.0001d0 !0.1d0
      DSP(ic,jc,2) =((dsin(8.d0*3.141592d0*xcen(ic)/(xmax-xmin)) &
                     +dsin(8.d0*3.141592d0*ycen(jc)/(ymax-ymin)))+1.d0)*0.0001d0 !0.1d0
      DSP(ic,jc,3) = 0.d0
    enddo
    enddo
    do jc=3,jmx1-2
    do ic=3,imx1-2
!     call random_number(random)
      if(jc==jmx1-2) goto 12
      xibl = real(min0(iblk(ic-1,jc),iblk(ic,jc)))
      DSP(ic,jc,1) = DSP(ic,jc,1)*xibl  !!(2.d0*random-1.d0)*0.1d0*xibl
   12 continue
      if(ic==imx1-2) goto 22
      yibl = real(min0(iblk(ic,jc-1),iblk(ic,jc)))
      DSP(ic,jc,2) = DSP(ic,jc,2)*yibl  !!(2.d0*random-1.d0)*0.1d0*yibl
   22 continue
    enddo
    enddo
!
    alp = 1.d+5
    do ialp=1,ialpha
      alp = alp*0.1d0
!
!
!-----J(qin)
      ustg(1:imx1,1:jmax) = QIN(1:imx1,1:jmax,1)
      vstg(1:imax,1:jmx1) = QIN(1:imax,1:jmx1,2)
      pcnt(1:imax,1:jmax) = 0.d0
      fJ1 = 0.d0
      do istep=1,nstep
        urec(1:imx1,1:jmax,istep) = ustg(1:imx1,1:jmax)
        vrec(1:imax,1:jmx1,istep) = vstg(1:imax,1:jmx1)
        prec(1:imax,1:jmax,istep) = pcnt(1:imax,1:jmax)
        do jc=1,jmax
        do ic=1,imx1
          ufrc(ic,jc,istep) = ustg(ic,jc)-uref(ic,jc,istep) 
          fJ1 = fJ1+0.5d0*ufrc(ic,jc,istep)**2
        enddo
        enddo
        do jc=1,jmx1
        do ic=1,imax
          vfrc(ic,jc,istep) = vstg(ic,jc)-vref(ic,jc,istep) 
          fJ1 = fJ1+0.5d0*vfrc(ic,jc,istep)**2
        enddo
        enddo
        call sub_bc_outer (ustg,vstg,pcnt)
        call sub_bc_wall  (ustg,vstg,pcnt)
        call sub_rhs_fwd
        call sub_HSMAC_fwd(ustg,vstg,pcnt)
      enddo
!
!
!-----J(qin+alpha*dsp)
      ustg(1:imx1,1:jmax) = QIN(1:imx1,1:jmax,1)+DSP(1:imx1,1:jmax,1)*alp
      vstg(1:imax,1:jmx1) = QIN(1:imax,1:jmx1,2)+DSP(1:imax,1:jmx1,2)*alp
      pcnt(1:imax,1:jmax) = 0.d0               !+DSP(1:imax,1:jmax,3)*alp
      fJ2 = 0.d0
      do istep=1,nstep
        do jc=1,jmax
        do ic=1,imx1
          fJ2 = fJ2+0.5d0*(ustg(ic,jc)-uref(ic,jc,istep))**2
        enddo
        enddo
        do jc=1,jmx1
        do ic=1,imax
          fJ2 = fJ2+0.5d0*(vstg(ic,jc)-vref(ic,jc,istep))**2
        enddo
        enddo
        call sub_bc_outer (ustg,vstg,pcnt)
        call sub_bc_wall  (ustg,vstg,pcnt)
        call sub_rhs_fwd
        call sub_HSMAC_fwd(ustg,vstg,pcnt)
      enddo
!
!
!-----Adjoint
      udlt(1:imx1,1:jmax) = 0.d0
      vdlt(1:imax,1:jmx1) = 0.d0
      ustg(1:imx1,1:jmax) = 0.d0
      vstg(1:imax,1:jmx1) = 0.d0
      pcnt(1:imax,1:jmax) = 0.d0
      do istep=nstep,1,-1
        ustg_b(1:imx1,1:jmax) = urec(1:imx1,1:jmax,istep) 
        vstg_b(1:imax,1:jmx1) = vrec(1:imax,1:jmx1,istep) 
        pcnt_b(1:imax,1:jmax) = prec(1:imax,1:jmax,istep) 
        call sub_bc_outer    (ustg_b,vstg_b,pcnt_b)
        call sub_bc_wall     (ustg_b,vstg_b,pcnt_b)
        call sub_HSMAC_adj   (ustg  ,vstg  ,pcnt  )
        call sub_rhs_adj
        call sub_bc_wall_adj (ustg  ,vstg  ,pcnt  )
        call sub_bc_outer_adj(ustg  ,vstg  ,pcnt  )
        ustg(1:imx1,1:jmax) = (ustg(1:imx1,1:jmax))+ufrc(1:imx1,1:jmax,istep) 
        vstg(1:imax,1:jmx1) = (vstg(1:imax,1:jmx1))+vfrc(1:imax,1:jmx1,istep) 
        pcnt(1:imax,1:jmax) = (pcnt(1:imax,1:jmax)) !+pfrc(1:imax,1:jmax,istep) 
      enddo
!
      den = 0.d0
      do jc=1,jmax
      do ic=1,imx1
        den = den+ustg(ic,jc)*DSP(ic,jc,1)
      enddo
      enddo
      do jc=1,jmx1
      do ic=1,imax
        den = den+vstg(ic,jc)*DSP(ic,jc,2)
      enddo
      enddo
!
!
      rslt(iave,ialp,1) = (fJ2-fJ1)/(alp*den)-1.d0
      rslt(iave,ialp,4) = alp
      write(*  ,'(4e25.17,2x,e32.25)') alp,(fJ2-fJ1),(alp*den),rslt(iave,ialp,1),rslt(iave,ialp,1) !+1.d0
      write(300,'(4e25.17,2x,e32.25)') alp,(fJ2-fJ1),(alp*den),rslt(iave,ialp,1),rslt(iave,ialp,1) !+1.d0
!
    enddo
  enddo
!
!
  write(*,*)
  do ialp=1,ialpha
    tlv_max = -1.d9
    tlv_min =  1.d9
    tlv_sum =  0.d0
    do iave =1,iaverage 
      tlv_sum = tlv_sum+dabs(rslt(iave,ialp,1))
      tlv_max = dmax1(dabs(rslt(iave,ialp,1)),tlv_max) 
      tlv_min = dmin1(dabs(rslt(iave,ialp,1)),tlv_min) 
    enddo 
    iave = iaverage
    write(*  ,'(3e25.17,e32.25)') rslt(1,ialp,4),tlv_max,tlv_min,tlv_sum/dble(iave)
    write(400,'(3e25.17,e32.25)') rslt(1,ialp,4),tlv_max,tlv_min,tlv_sum/dble(iave)
  enddo
!
  return
end subroutine sub_check_foa
!-----------------------------------------------------------------------
