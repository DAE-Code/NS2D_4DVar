!**********************************************************************
subroutine sub_4dvar
!**********************************************************************
  use mod_variables
  use mod_temp
  use m_random3
  implicit none
!
  integer :: ic,jc,kc,i,j,n
!
!----------------------------------------------------------------------
  character(len=6)  :: ioid = '000000'
  character(len=99) :: fname1,fname2
  integer iter,nmax,mmax,indx,iistep,rec1,rec2,iolength,ios
  real(8) eps1,eps2,eps3,eta 
  real(8) funcJ,oldfJ,newfJ,Jdf,hflag
  real(8) gnorm,snorm,ynorm,sy,sigma_b
  data    eps1,eps2,eps3,eta  / 1.d-4,1.d-7,1.d-8,1.d-4 /
!
  real(8),allocatable :: gJ(:),oldgJ(:),ns(:),ny(:),v(:)
  real(8),allocatable :: ss(:,:),yy(:,:),oldss(:,:),oldyy(:,:)
  real(8),allocatable :: ssyy(:),oldssyy(:),vv(:),sigma(:)
  real(8),allocatable :: nx(:),newx(:),oldx(:),direc(:)
  real(8),allocatable :: Jy2D(:,:),Jy2Dg(:,:),Jy2Dl(:,:)
!----------------------------------------------------------------------
!
!
  write(*,*) 
  write(*,'(//," ===== 4D-Var iteration =====")')
!
  mmax = 20
  nmax = imx1*jmax+imax*jmx1
! ITRMAX = 50
  allocate(gJ(nmax),oldgJ(nmax),ns(nmax),ny(nmax),v(nmax))
  allocate(ss(nmax,mmax),yy(nmax,mmax),oldss(nmax,mmax),oldyy(nmax,mmax))
  allocate(ssyy(mmax),oldssyy(mmax),vv(nmax),sigma(mmax))
  allocate(nx(nmax),newx(nmax),oldx(nmax),direc(nmax))
  gJ(:)   = 0.d0 ; oldgJ(:)   = 0.d0 ; ns(:)      = 0.d0 ; ny(:)      = 0.d0 ; v(:) = 0.d0
  ss(:,:) = 0.d0 ; yy(:,:)    = 0.d0 ; oldss(:,:) = 0.d0 ; oldyy(:,:) = 0.d0
  ssyy(:) = 0.d0 ; oldssyy(:) = 0.d0 ; vv(:)      = 0.d0 ; sigma(:)   = 0.d0
  nx(:)   = 0.d0 ; newx(:)    = 0.d0 ; oldx(:)    = 0.d0 ; direc(:)   = 0.d0
!
!
!-Measurement
  call sub_measure
!
!
!-Temporary variables
  allocate(PQ1(imx1,jmx1,nva),PQ2(imx1,jmx1,nva), &
           DLT(imx1,jmx1,nva),DSP(imx1,jmx1,nva), &
           QIN(imx1,jmx1,nva))
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
!
!
!-Save intial field
  QIN(1:imx1,1:jmax,1) = ustg(1:imx1,1:jmax)
  QIN(1:imax,1:jmx1,2) = vstg(1:imax,1:jmx1)
  QIN(1:imax,1:jmax,3) = pcnt(1:imax,1:jmax)
!
!
  allocate(urec(1:imx1,1:jmx1,nstep), &
           vrec(1:imx1,1:jmx1,nstep), &
           prec(1:imx1,1:jmx1,nstep))
  urec(1:imx1,1:jmx1,:) = 0.d0
  vrec(1:imx1,1:jmx1,:) = 0.d0
  prec(1:imx1,1:jmx1,:) = 0.d0
  allocate(uref(1:imx1,1:jmx1,nsall), &
           vref(1:imx1,1:jmx1,nsall), &
           pref(1:imx1,1:jmx1,nsall))
  uref(1:imx1,1:jmx1,:) = 0.d0
  vref(1:imx1,1:jmx1,:) = 0.d0
  pref(1:imx1,1:jmx1,:) = 0.d0
  allocate(ufrc(1:imx1,1:jmx1,nstep), &
           vfrc(1:imx1,1:jmx1,nstep), &
           pfrc(1:imx1,1:jmx1,nstep))
  ufrc(1:imx1,1:jmx1,:) = 0.d0
  vfrc(1:imx1,1:jmx1,:) = 0.d0
  pfrc(1:imx1,1:jmx1,:) = 0.d0
  allocate(uwns(1:imx1,1:jmx1,nstep), &
           vwns(1:imx1,1:jmx1,nstep))
  do istep=1,nstep
    call random3(uwns(:,:,istep),imx1,jmx1)
    call random3(vwns(:,:,istep),imx1,jmx1)
    uwns(:,:,istep) = dsqrt(vars%ref)*uwns(:,:,istep)
    vwns(:,:,istep) = dsqrt(vars%ref)*vwns(:,:,istep)
  enddo
!
!
!-> Forward (reference) -------------------------------------------------------
  ustg(1:imx1,1:jmax) = QIN(1:imx1,1:jmax,1)
  vstg(1:imax,1:jmx1) = QIN(1:imax,1:jmx1,2)
  pcnt(1:imax,1:jmax) = QIN(1:imax,1:jmax,3)
  if(n_prb==2) then
!   Set vortex
    call sub_addvtx(-4.d0,0.d0)
  else
!   Forward (shift in time)
    do istep=1,int(7.2d0/dt)/2   ! Shift half period of vortex shedding
      call sub_bc_outer (ustg,vstg,pcnt)
      call sub_bc_wall  (ustg,vstg,pcnt)
      call sub_rhs_fwd
      call sub_HSMAC_fwd(ustg,vstg,pcnt)
    enddo
!------------------------------------------------------------------------------
  endif
!
  open(10,file='RMSE.dat')
  write(10,'("  Step       RMSE")') 
  close(10)
!
!==============================================================================
!====Sequential loop comes back here===========================================
!==============================================================================
  icycl = 1
1234 continue
  call sub_p3dwrite(ostep,tstep,'ref')
  do istep=1,nsall
    tstep = tstep+dt ! current time
    ostep = ostep+1  ! current time-step
    uref(1:imx1,1:jmax,istep) = ustg(1:imx1,1:jmax)
    vref(1:imax,1:jmx1,istep) = vstg(1:imax,1:jmx1)
    pref(1:imax,1:jmax,istep) = pcnt(1:imax,1:jmax)
    call sub_bc_outer (ustg,vstg,pcnt)
    call sub_bc_wall  (ustg,vstg,pcnt)
    call sub_rhs_fwd
    call sub_HSMAC_fwd(ustg,vstg,pcnt)
    call sub_probe2(istep+(icycl-1)*nstep,1,icycl,"ref")
    if(mod(istep,iskip_plot)==0) call sub_p3dwrite(ostep,tstep,'ref')
  enddo
  PQ1(1:imx1,1:jmax,1) = ustg(1:imx1,1:jmax) ! Save as initial considion of ref case in the next period
  PQ1(1:imax,1:jmx1,2) = vstg(1:imax,1:jmx1)
  PQ1(1:imax,1:jmax,3) = pcnt(1:imax,1:jmax)
!------------------------------------------------------------------------------
!
!
!-> Base flow state (initial value of cost function)
  if(modef==3) then
    call sub_restart_4DV(ios)
    if(ios==0) then
      write(*,*) "Reading 4DV restart file"
      QIN(1:imx1,1:jmax,1) = ustg(1:imx1,1:jmax)
      QIN(1:imax,1:jmx1,2) = vstg(1:imax,1:jmx1)
      QIN(1:imax,1:jmax,3) = pcnt(1:imax,1:jmax)
    endif
  endif
  ustg(1:imx1,1:jmax) = QIN(1:imx1,1:jmax,1)
  vstg(1:imax,1:jmx1) = QIN(1:imax,1:jmx1,2)
  pcnt(1:imax,1:jmax) = QIN(1:imax,1:jmax,3)
  call sub_p3dwrite(0,dble(0),'ini')
!
!-- array conversion
  indx = 0
  do jc=1,jmax
  do ic=1,imx1
    indx = indx+1
    nx(indx) = ustg(ic,jc)
  enddo
  enddo
  do jc=1,jmx1
  do ic=1,imax
    indx = indx+1
    nx(indx) = vstg(ic,jc)
  enddo
  enddo
!
!-> Forward -------------------------------------------------------------------
  funcJ = 0.d0
  if(vars%bkg>0.d0) then
    call sub_filter(ustg,vstg,Jdf)  ! Add Jdf to funcJ
    write(*,*) "Jdf:",Jdf
    funcJ = funcJ+Jdf
  endif
  do istep=1,nstep
    urec(1:imx1,1:jmax,istep) = ustg(1:imx1,1:jmax)
    vrec(1:imax,1:jmx1,istep) = vstg(1:imax,1:jmx1)
    prec(1:imax,1:jmax,istep) = pcnt(1:imax,1:jmax)
    do jc=1,jmax
    do ic=1,imx1
      hflag = dble(ihx(ic,jc,istep))
      ufrc(ic,jc,istep) = (ustg(ic,jc)-uref(ic,jc,istep)-uwns(ic,jc,istep))*hflag
      funcJ = funcJ+0.5d0*ufrc(ic,jc,istep)**2
    enddo
    enddo
    do jc=1,jmx1
    do ic=1,imax
      hflag = dble(ihy(ic,jc,istep))
      vfrc(ic,jc,istep) = (vstg(ic,jc)-vref(ic,jc,istep)-vwns(ic,jc,istep))*hflag
      funcJ = funcJ+0.5d0*vfrc(ic,jc,istep)**2
    enddo
    enddo
    call sub_bc_outer (ustg,vstg,pcnt)
    call sub_bc_wall  (ustg,vstg,pcnt)
    call sub_rhs_fwd
    call sub_HSMAC_fwd(ustg,vstg,pcnt)
    call sub_probe2(istep+(icycl-1)*nstep,1,icycl,"est")
  enddo

!------------------------------------------------------------------------------
!
  iter = 0
  write(*  ,'(i3,2e15.7)') iter,funcJ
  write(300,'(i3,2e15.7)') iter,funcJ
!
!-> Gradient ------------------------------------------------------------------
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
    do jc=1,jmax
    do ic=1,imx1
      hflag = dble(ihx(ic,jc,istep))
      ustg(ic,jc) = (ustg(ic,jc))+ufrc(ic,jc,istep)*hflag
    enddo
    enddo
    do jc=1,jmx1
    do ic=1,imax
      hflag = dble(ihy(ic,jc,istep))
      vstg(ic,jc) = (vstg(ic,jc))+vfrc(ic,jc,istep)*hflag
    enddo
    enddo
  enddo
  if(vars%bkg>0.d0) then
    call sub_filter_adj(ustg,vstg,ustg_b,vstg_b)  ! Add grad of Jdf
  endif
!------------------------------------------------------------------------------
!
!-- array conversion
  indx = 0
  do jc=1,jmax
  do ic=1,imx1
    indx = indx+1
    gJ(indx) = ustg(ic,jc)
  enddo
  enddo
  do jc=1,jmx1
  do ic=1,imax
    indx = indx+1
    gJ(indx) = vstg(ic,jc)
  enddo
  enddo
!
!-> Initial search direction
  do i=1,nmax
    direc(i) = -gJ(i)
  enddo
!
!******************************************************
!-> Newton iteration
  iter = 1
  do while(.true.)
!
!---- array conversion
    indx = 0
    do jc=1,jmax
    do ic=1,imx1
      indx = indx+1
      ustg(ic,jc) = direc(indx)
    enddo
    enddo
    do jc=1,jmx1
    do ic=1,imax
      indx = indx+1
      vstg(ic,jc) = direc(indx)
    enddo
    enddo
    call sub_p3dwrite(iter,dble(iter),'dir')
!
!---> linear search
    write(*  ,'(//," **** Linear Search ****")') 
    write(100,'(//," **** Linear Search ****")') 
    write(400,'(//," **** Linear Search ****")')
    call linear_search(nmax,funcJ,newfJ,gJ,nx,newx,direc)
!
!---> evaluate gradient at the new iterate
    oldfJ = funcJ
    funcJ = newfJ
    do i=1,nmax
      oldx(i)  = nx(i)
      nx(i)    = newx(i)
      ns(i)    = nx(i)-oldx(i)
      oldgJ(i) = gJ(i) 
    enddo
    write(*  ,'(i3,2e15.7)') iter,funcJ
    write(300,'(i3,2e15.7)') iter,funcJ
!
!
!---- array conversion
    indx = 0
    do jc=1,jmax
    do ic=1,imx1
      indx = indx+1
      ustg(ic,jc) = nx(indx)
    enddo
    enddo
    do jc=1,jmx1
    do ic=1,imax
      indx = indx+1
      vstg(ic,jc) = nx(indx)
    enddo
    enddo
    pcnt(:,:) = QIN(1:imax,1:jmax,3)
!
    call sub_p3dwrite(iter,dble(iter),'ini')
!
!---> Forward -----------------------------------------------------------------
    do istep=1,nstep
      urec(1:imx1,1:jmax,istep) = ustg(1:imx1,1:jmax)
      vrec(1:imax,1:jmx1,istep) = vstg(1:imax,1:jmx1)
      prec(1:imax,1:jmax,istep) = pcnt(1:imax,1:jmax)
      do jc=1,jmax
      do ic=1,imx1
        hflag = dble(ihx(ic,jc,istep))
        ufrc(ic,jc,istep) = (ustg(ic,jc)-uref(ic,jc,istep)-uwns(ic,jc,istep))*hflag
      enddo
      enddo
      do jc=1,jmx1
      do ic=1,imax
        hflag = dble(ihy(ic,jc,istep))
        vfrc(ic,jc,istep) = (vstg(ic,jc)-vref(ic,jc,istep)-vwns(ic,jc,istep))*hflag
      enddo
      enddo
      call sub_bc_outer (ustg,vstg,pcnt)
      call sub_bc_wall  (ustg,vstg,pcnt)
      call sub_rhs_fwd
      call sub_HSMAC_fwd(ustg,vstg,pcnt)
!     call sub_measurement(1,Jy2Dg,Jy2Dl)
      call sub_probe2(istep+(icycl-1)*nstep,iter+1,icycl,"est")
    enddo
!------------------------------------------------------------------------------
!
!---> Gradient ----------------------------------------------------------------
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
      do jc=1,jmax
      do ic=1,imx1
        hflag = dble(ihx(ic,jc,istep))
        ustg(ic,jc) = (ustg(ic,jc))+ufrc(ic,jc,istep)*hflag
      enddo
      enddo
      do jc=1,jmx1
      do ic=1,imax
        hflag = dble(ihy(ic,jc,istep))
        vstg(ic,jc) = (vstg(ic,jc))+vfrc(ic,jc,istep)*hflag
      enddo
      enddo
    enddo
    if(vars%bkg>0.d0) then
      call sub_filter_adj(ustg,vstg,ustg_b,vstg_b)  ! Add grad of Jdf
    endif
!------------------------------------------------------------------------------
!
    call sub_p3dwrite(iter,dble(iter),'grd')
!
!---- array conversion
    indx = 0
    do jc=1,jmax
    do ic=1,imx1
      indx = indx+1
      gJ(indx) = ustg(ic,jc)
    enddo
    enddo
    do jc=1,jmx1
    do ic=1,imax
      indx = indx+1
      gJ(indx) = vstg(ic,jc)
    enddo
    enddo
!
!
    do n=1,nmax
      ny(n) = gJ(n)-oldgJ(n)
    enddo
!
!---> Check termination criteria
    do i=1,NMAX
      if(nx(i)>1.d+20) exit
    enddo
!
    gnorm = 0.d0
    do i=1,NMAX
      gnorm = gnorm+gJ(i)**2
    enddo
    gnorm = dsqrt(gnorm)
    if((oldfJ-funcJ)<(eps2*(1.d0+dabs(oldfJ))) .and. gnorm<eps1) exit
!
    write(*  ,'(i3,"  gJ norm:",e15.7,"  fJ diff:",e15.7)') iter,gnorm,oldfJ-funcJ
    write(100,'(i3,"  gJ norm:",e15.7,"  fJ diff:",e15.7)') iter,gnorm,oldfJ-funcJ
    write(200,'(i3,"  gJ norm:",e15.7,"  fJ diff:",e15.7)') iter,gnorm,oldfJ-funcJ
!
!---> Update inverse Hessian approximation
!
    snorm = 0.d0
    ynorm = 0.d0
    sy    = 0.d0
    do i=1,NMAX
      snorm = snorm+ns(i)**2
      ynorm = ynorm+ny(i)**2
      sy    = sy+ns(i)*ny(i)
    enddo
    snorm = dsqrt(snorm)
    ynorm = dsqrt(ynorm)
!
    if(sy < eta*snorm*ynorm) then
      write(*,*) "Hessian update is skipped at iteration: ",iter
      write(*,*) "To preserve positive definitenss"
      goto 300
    endif
!
!
!---> Store s(i), y(i), rho(i)=1/(s^t y)
    if(iter > mmax) then
      do j=1,mmax 
         do n=1,nmax 
            oldss(n,j) = ss(n,j) 
            oldyy(n,j) = yy(n,j) 
         enddo 
         oldssyy(j) = ssyy(j) 
      enddo 
      do j=1,mmax-1 
         do n=1,nmax 
            ss(n,j) = oldss(n,j+1) 
            yy(n,j) = oldyy(n,j+1) 
         enddo 
         ssyy(j) = oldssyy(j+1) 
      enddo 
    endif
!
    do n=1,nmax 
       ss(n,min0(iter,mmax)) = ns(n) 
       yy(n,min0(iter,mmax)) = ny(n) 
    end do 
    ssyy(min0(iter,mmax)) = 1.d0/sy 
!
!---> Hessian update based on limited memory BFGS
! step 1
    do n=1,nmax 
       v(n) = gJ(n) 
    enddo 
! step 2
    do i=min0(iter,mmax),1,-1
      sigma(i) = 0.d0
      do j=1,nmax 
         sigma(i) = sigma(i)+ssyy(i)*ss(j,i)*v(j) 
      enddo 
      do n=1,nmax 
         v(n) = v(n)-sigma(i)*yy(n,i) 
      enddo 
    enddo
! step 3
!   v(:) = v(:)
! step 4
    do i=1,min0(iter,MMAX)
      sigma_b = 0.d0
      do j=1,nmax 
         sigma_b = sigma_b+ssyy(i)*yy(j,i)*v(j) 
      enddo 
      do n=1,nmax 
         v(n) = v(n)+(sigma(i)-sigma_b)*ss(n,i) 
      enddo 
    enddo
!---> new search direction
    do n=1,nmax
      direc(n) = -v(n)
    enddo
!
300 continue
    if(iter>=ITRMAX) exit
    iter = iter+1
!
  enddo
!******************************************************
!
!+Evaluate RMSE++++++++++++++++++++++++++++++++++++++++
!---- array conversion
  indx = 0
  do jc=1,jmax
  do ic=1,imx1
    indx = indx+1
    ustg(ic,jc) = nx(indx)
  enddo
  enddo
  do jc=1,jmx1
  do ic=1,imax
    indx = indx+1
    vstg(ic,jc) = nx(indx)
  enddo
  enddo
  pcnt(:,:) = QIN(1:imax,1:jmax,3)
!
  open(10,file='RMSE.dat',position='append')
!
  tstep = tstep_ini
  ostep = ostep_ini
  newfJ = 0.d0
  do istep=1,nsall
    funcJ = 0.d0
    tstep = tstep+dt ! current time
    ostep = ostep+1  ! current time-step
    do jc=1,jmax
    do ic=1,imx1
      funcJ = funcJ+(ustg(ic,jc)-uref(ic,jc,istep))**2  ! RMSE without measurement error 
    enddo
    enddo
    do jc=1,jmx1
    do ic=1,imax
      funcJ = funcJ+(vstg(ic,jc)-vref(ic,jc,istep))**2
    enddo
    enddo
    call sub_bc_outer (ustg,vstg,pcnt)
    call sub_bc_wall  (ustg,vstg,pcnt)
    call sub_rhs_fwd
    call sub_HSMAC_fwd(ustg,vstg,pcnt)
    call sub_probe2(istep+(icycl-1)*nstep,1,icycl,"bst")
    if(mod(istep+(icycl-1)*nsall,iskip_plot)==0) call sub_p3dwrite(ostep,tstep,'bst')
!   if(mod(istep+(icycl-1)*nsall,tskip)==0) then
      write(10,'(i6,e15.7)') istep+(icycl-1)*nsall,dsqrt(funcJ/dble(imx1*jmax+imax*jmx1))
!   endif
    newfJ = newfJ+funcJ
  enddo
  tstep_ini = tstep
  ostep_ini = ostep
  newfJ = dsqrt(newfJ/dble(nsall*(imx1*jmax+imax*jmx1)))
  close(10)
  write(*,*)
  open(10,file='RMSE_all.dat')
  write(10,'("RMSE:",e15.7)') newfJ
  write(* ,'("RMSE:",e15.7)') newfJ
  close(10)
!++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!==============================================================================
! Data assimilation cycle
! --use best trajectry as initial conditions of nest period--
!==============================================================================
  QIN(1:imx1,1:jmax,1) = ustg(1:imx1,1:jmax) ! Analysis
  QIN(1:imax,1:jmx1,2) = vstg(1:imax,1:jmx1)
  QIN(1:imax,1:jmax,3) = pcnt(1:imax,1:jmax)
  ustg(1:imx1,1:jmax) = PQ1(1:imx1,1:jmax,1) ! Restart of reference case
  vstg(1:imax,1:jmx1) = PQ1(1:imax,1:jmx1,2)
  pcnt(1:imax,1:jmax) = PQ1(1:imax,1:jmax,3)
  if(icycl>=ncycl) return
  icycl = icycl+1
  goto 1234
!
  return
end subroutine sub_4dvar
!**********************************************************************
!
!
!
!**********************************************************************
subroutine linear_search(nmax,funcJ,newfJ,gJ,nx,newx,direc)
!**********************************************************************
  use mod_variables
  use mod_temp
  implicit none
!
  logical :: lskip
  integer i,ic,jc,k,m,nmax,ierr,indx,iistep,nmes
  real(8) beta,lambda,mu,eps,gd,newgd,hflag
  real(8) c,st1,st2,st,dirder,dnorm,xnorm,minst
  real(8) funcJ,newfJ,Jdf,gJ(NMAX),newgJ(NMAX)
  real(8) nx(NMAX),newx(NMAX),direc(NMAX)
  data    beta,lambda,mu,eps  / 1.d-3,0.1d0,0.5d0,1.d-8 /
!
!
!---> Directional derivatives
  dirder = 0.d0
  dnorm  = 0.d0
  xnorm  = 0.d0
  do i=1,NMAX
    dirder = dirder+gJ(i)   *direc(i)
    dnorm  = dnorm +direc(i)*direc(i)
    xnorm  = xnorm +nx(i)   *nx(i)
  enddo
  dnorm = dsqrt(dnorm)
  xnorm = dsqrt(xnorm)
  minst = dmin1(1.d0,eps/dnorm)
!
!-> set initial step size
  st = 1.d0 !!LSini
!
  do while(.true.)
    write(*  ,'(" ===============================")') 
    write(*  ,'(" st =",e15.7," minst =",e15.7)')  st,minst
    write(100,'(" ===============================")') 
    write(100,'(" st =",e15.7," minst =",e15.7)')  st,minst
    write(400,'(" ===============================")') 
    write(400,'(" st =",e15.7," minst =",e15.7)')  st,minst
    if(st<minst) return
!
    do i=1,nmax
      newx(i) = nx(i)+st*direc(i)
    enddo
!
!-- array conversion
    indx = 0
    do jc=1,jmax
    do ic=1,imx1
      indx = indx+1
      ustg(ic,jc) = newx(indx)
    enddo
    enddo
    do jc=1,jmax
    do ic=1,imax
      indx = indx+1
      vstg(ic,jc) = newx(indx)
    enddo
    enddo
    pcnt(:,:) = QIN(1:imax,1:jmax,3)
!
!-- Check CFL
    call sub_cfl(cfl)
    if(cfl>0.5d0) then
      st = 0.5d0*st
      cycle
    endif
!
!-> Forward -------------------------------------------------------------------
    newfJ = 0.d0
    if(vars%bkg>0.d0) then
      call sub_filter(ustg,vstg,Jdf)   ! Add Jdf to newfJ
      write(*,*) "Jdf:",Jdf
      newfJ = newfJ+Jdf
    endif
    do istep=1,nstep
      do jc=1,jmax
      do ic=1,imx1
        hflag = dble(ihx(ic,jc,istep))
        ufrc(ic,jc,istep) = (ustg(ic,jc)-uref(ic,jc,istep)-uwns(ic,jc,istep))*hflag
        newfJ = newfJ+0.5d0*ufrc(ic,jc,istep)**2
      enddo
      enddo
      do jc=1,jmx1
      do ic=1,imax
        hflag = dble(ihy(ic,jc,istep))
        vfrc(ic,jc,istep) = (vstg(ic,jc)-vref(ic,jc,istep)-vwns(ic,jc,istep))*hflag
        newfJ = newfJ+0.5d0*vfrc(ic,jc,istep)**2
      enddo
      enddo
      call sub_bc_outer (ustg,vstg,pcnt)
      call sub_bc_wall  (ustg,vstg,pcnt)
      call sub_rhs_fwd
      call sub_HSMAC_fwd(ustg,vstg,pcnt)
    enddo
!------------------------------------------------------------------------------
!
!-> Wolef conditions
    if(newfJ<funcJ+beta*dirder*st) return
!
!-> calculate new step size
    c   = (newfJ-funcJ-dirder*st)/st**2
    st1 =-0.5d0*dirder/c
    st2 = dmin1(st1,mu*st)   
    st  = dmax1(lambda*st,st2)
!
  enddo
!
  return
end subroutine linear_search
!**********************************************************************

