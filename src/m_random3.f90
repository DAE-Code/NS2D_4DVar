module m_random3
contains
subroutine random3(work1,n_var,n_ens)
! Returns a vector of normal (Gaussian) random numbers N(mean=0,variance=1)
!
  implicit none
  integer             :: i,j
  integer,intent(in)  :: n_var,n_ens
  real(8),intent(out) :: work1(n_var,n_ens)
  real(8)             :: work2(n_var,n_ens)
  real(8),parameter   :: pi = 3.141592653589d0
!
  do j=1,n_ens
    call random_number(work1(1:n_var,j))  ! Uniform random number between 0 and 1
    call random_number(work2(1:n_var,j))
  enddo
!
! Box-Muller's method
  do j=1,n_ens
  do i=1,n_var
    work1(i,j)= dsqrt(-2.d0*dlog(work1(i,j)))*dcos(2.d0*pi*work2(i,j))
  enddo
  enddo
!
end subroutine random3
end module m_random3
