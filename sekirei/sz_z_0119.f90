subroutine sz(n,dim,DoF,spin,sztot)
  use variables
  implicit none
  integer :: n,dim,Dof
  integer*8 :: i,j,k,m
  double precision :: spin,sztot
  double precision :: litot


  k=0
  szdim=0

  if (mod(Dof-1,2).eq.1) then
     spinmax=2*spin
     write(*,*) "half-integer spin", spinmax, "/2"
  else
     spinmax=spin
     write(*,*) "integer spin", spinmax
  end if

  do i=1,dim
     k=i-1
     litot=0.0d0

     do j=n,1,-1
        litot = litot + (mod(k,Dof)-spin)
        k=k/DoF
     end do
     if(litot.eq.sztot) then
        szdim=szdim+1
     end if
  end do


  allocate(outpt(szdim),outpt_inv(0:dim))
  allocate(list(szdim,n+1))

  m=0
  outpt_inv=0
  outpt_inv(0)=0
  do i=1,dim
     k=i-1
     litot=0.0d0
     do j=n,1,-1
        litot = litot + (mod(k,Dof)-spin)
        k=k/DoF
     end do
     if(litot.eq.sztot) then
        m=m+1
        outpt(m)=i
        outpt_inv(i)=m
     end if
  end do


  if (mod(Dof-1,2) .eq. 0) then
     call sz_int(n,DoF,spin)
  else
     call sz_half(n,DoF,spin)
  end if

end subroutine sz

subroutine sz_int(n,Dof,spin)
  use variables
  implicit none
  integer :: n,Dof!,szdim
  integer :: i,j,k
  double precision :: spin
!  integer :: list(szdim,n+1),outpt(szdim)

  k=0
  !$omp parallel
  !$omp do private(k,j)
  do i=1,szdim
     k=outpt(i)-1
     list(i,n+1)=outpt(i)
     do j=n,1,-1
        list(i,j) = mod(k,DoF)-spin
        k=k/DoF
     end do
  end do
  !$omp end do
  !$omp end parallel

end subroutine sz_int

subroutine sz_half(n,Dof,spin)
  use variables
  implicit none
  integer :: n,Dof!,szdim
  integer :: i,j,k
  double precision :: spin
!  integer :: list(szdim,n+1),outpt(szdim)

  k=0
  !$omp parallel
  !$omp do private(k,j)  
  do i=1,szdim
     k=outpt(i)-1
     list(i,n+1)=outpt(i)
     do j=n,1,-1
        list(i,j) = (mod(k,DoF)-spin)*2
        k=k/DoF
     end do
  end do
  !$omp end do
  !$omp end parallel


end subroutine sz_half


