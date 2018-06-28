!lanczos method
!if GS is degenerate resuts is not right
subroutine lanczos_main(subspace,v_first,energy,ene_num,eigen_vec,eigen_vec_2,info)
  use variables
  implicit none
  integer ::  i,m,subspace,info,lpnum,ene_num!,incx=1,incy=1
  double precision :: abstol
!  complex(kind(0d0)) zdotc,zres
  double precision energy(ene_num),exp_val(2)
  complex(kind(0d0)) eigen_vec(szdim),v_first(szdim),eigen_vec_2(szdim)

  integer,allocatable :: iwork(:),ifail(:)
  double precision,allocatable :: work(:)
  double precision ,allocatable :: alpha(:),beta(:),tri_vec(:,:)
  complex(kind(0d0)), allocatable :: v(:,:),eigen_vectors(:,:)

  allocate(eigen_vectors(szdim,2))
  allocate(alpha(subspace),beta(subspace))
  allocate(tri_vec(subspace,ene_num))
  allocate(iwork(5*subspace),work(5*subspace),ifail(subspace))
  allocate(v(szdim,subspace+1))


  eigen_vec=(0.0d0,0.0d0)
  eigen_vec_2=(0.0d0,0.0d0)
  eigen_vectors=(0.0d0,0.0d0)
  tri_vec=(0.0d0,0.0d0)
  exp_val=0.0d0
  energy=0.0d0


!  do i=1,2
  call lanczos_tri(subspace,ene_num,tri_vec,energy,v_first,v,alpha,beta,lpnum,info)
!  write(*,*)lpnum
  
  call eigen_vector(lpnum,tri_vec(:,1:2),v(:,1:lpnum),eigen_vectors,exp_val)

  
  if(info.ne.0) write(*,*) "info error lanczos 39",info
  if(abs(energy(1)-exp_val(1))>1.0E-8) write(*,300) "eigen value error 1",energy(1),exp_val(1),energy(1)-exp_val(1)
  if(abs(energy(2)-exp_val(2))>1.0E-8) write(*,300) "eigen value error 2",energy(2),exp_val(2),energy(2)-exp_val(2)


  
  ! tri_vec(:,1)=tri_vec(:,2)
  ! call eigen_vector(lpnum,tri_vec,v(:,1:lpnum),eigen_vectors(:,2),exp_val)
  ! if (info .ne. 0 .or. abs(energy(2)-exp_val)>1.0E-8) then
  !    if(info.ne.0) write(*,*) "info error",info
  !    if(abs(energy(2)-exp_val)>1.0E-8) write(*,*) "eigen value error 2",abs(energy(2)-exp_val(2))
  !    !       stop
  ! end if  
  ! do j=1,szdim
  !    eigen_vectors(j,i)=eigen_vec(j)
  ! end do

!  write(*,200) lpnum,energy(1),energy(2),energy(2)-energy(1)
200 format(I3,f20.13,f20.13,f20.13)
300 format(A,F16.8,F16.8,F16.8)
  


!  end do

  eigen_vec(:)=eigen_vectors(:,1)
  eigen_vec_2(:)=eigen_vectors(:,2)

end subroutine lanczos_main


subroutine lanczos_tri(subspace,ene_num,tri_vec,energy,v_first,v,alpha,beta,lpnum,info)
  use variables
  implicit none
  double precision,parameter :: eps=1.0E-8 !if eps are very small, the excited energy maybe foult,but too large..
  integer ::  j,k,break_num,subspace,incx,incy,m,lpnum,info
  integer :: check_ene_num,ene_num
  double precision :: a_beta,abstol,ebefor
  complex(kind(0d0)) :: ei,zdotc
  character(len=1) :: uplo
  
!  complex(kind(0d0)) ::  eigen_vec(szdim,2)
  double precision :: alpha(subspace),beta(subspace),tri_vec(subspace,ene_num),energy(ene_num)
  complex(kind(0d0)) :: v(szdim,subspace+1),v_first(szdim)

  !   double precision,allocatable :: v_tmp_s(:)
  complex(kind(0d0)),allocatable :: v_tmp(:),y_tmp(:),w_i(:)

  integer,allocatable :: iwork(:),ifail(:)
  double precision,allocatable :: work(:)
  
  allocate(iwork(5*subspace),work(5*subspace),ifail(subspace))
  

  !   integer :: seedsize
  !   integer, allocatable :: seed(:)

  check_ene_num=ene_num
  m=check_ene_num
  allocate(v_tmp(szdim))
  !   allocate(v_tmp_s(2*szdim))
  allocate(y_tmp(szdim),w_i(szdim))
!  eigen_vec=0.0d0


  uplo="U"

  k=0
  break_num=1
  ei=(0.0d0,1.0d0)
  incx=1
  incy=1

  v_tmp=v_first


     alpha=0.0d0
     beta=0.0d0
     v=(0.0d0,0.0d0)

     v(:,1)=v_tmp(:)

     y_tmp=(0.0d0,0.0d0)
     w_i=(0.0d0,0.0d0)

     call zcsrsymv(v(:,1),y_tmp)

     alpha(1)=zdotc(szdim,v(:,1),incx,y_tmp,incy)

     w_i(:)=y_tmp(:)-alpha(1)*v(:,1)

     beta(1)=sqrt(dble(zdotc(szdim,w_i,incx,w_i,incy)))

     a_beta=1/beta(1)
     v_tmp=w_i*a_beta

     v(:,2)=v_tmp(:)

     do j=2,subspace
        y_tmp=(0.0d0,0.0d0)
        w_i=(0.0d0,0.0d0)

        call zcsrsymv(v(:,j),y_tmp)

        alpha(j)=zdotc(szdim,v(:,j),incx,y_tmp,incy)

        w_i(:)=y_tmp(:)-alpha(j)*v(:,j)-beta(j-1)*v(:,j-1)

        beta(j)=sqrt(dble(zdotc(szdim,w_i,incx,w_i,incy)))
        if (j>20 .or. j>szdim/2)then 
           if((mod(j,10)==0))then
              tri_vec=0.0d0
              !alpha(1:j),beta(1:j)
              call dstevx('N', 'I', j, alpha, beta, 1.0, 2.0, 1, check_ene_num, abstol, m, energy, tri_vec, j, work, iwork, ifail, info)
              if(abs((energy(check_ene_num)-ebefor)/energy(check_ene_num))<eps) then
                 !           write(*,*) 'comverged loop',j
                 lpnum=j
!                 write(*,*) lpnum
                 call dstevx('V', 'I', j, alpha, beta, 1.0, 2.0, 1, check_ene_num, abstol, m, energy, tri_vec, j, work, iwork, ifail, info)

                 return
              end if
              ebefor=energy(check_ene_num)

           end if
        end if
        if(j==20.and.j<=szdim/2)then
           tri_vec=0.0d0
           !alpha(1:j),beta(1:j)
           call dstevx('N', 'I', j, alpha, beta, 1.0, 2.0, 1, check_ene_num, abstol,m, energy, tri_vec, j, work, iwork, ifail, info)
           ebefor=energy(check_ene_num)
        end if
        if(j<=20.and.j==szdim/2)then
           tri_vec=0.0d0
           !alpha(1:j),beta(1:j)
           call dstevx('N', 'I', j, alpha, beta, 1.0, 2.0, 1, check_ene_num, abstol,m, energy, tri_vec, j, work, iwork, ifail, info)
           ebefor=energy(check_ene_num)
        end if
        
        a_beta=1.0d0/beta(j)

        v_tmp=w_i*a_beta
        v(:,j+1)=v_tmp(:)
     end do

     v_tmp(:)=v(:,subspace+1)

!  end do

end subroutine lanczos_tri


subroutine eigen_vector(subspace,tri_vec,large_v,eigen_vectors,exp_val)
  use variables
  implicit none

  integer :: i,j,subspace,incx,incy
  double precision :: s_tmp
  complex(kind(0d0)) :: zdotc
!  character(len=1) :: uplo='U'

  complex(kind(0d0)) :: large_v(szdim,subspace),eigen_vectors(szdim,2)
  complex(kind(0d0)),allocatable :: y_tmp(:)
  double precision :: tri_vec(subspace,4),exp_val(2)

  allocate(y_tmp(szdim))
  incx=1
  incy=1
  y_tmp=(0.0d0,0.0d0)
  eigen_vectors=(0.0d0,0.0d0)

  do i=1,szdim
     do j=1,subspace
        eigen_vectors(i,1)=eigen_vectors(i,1)+large_v(i,j)*tri_vec(j,1)
        eigen_vectors(i,2)=eigen_vectors(i,2)+large_v(i,j)*tri_vec(j,2)        
     end do
  end do

  y_tmp=(0.0d0,0.0d0)
  s_tmp=0.0d0
  s_tmp=sqrt(dble(zdotc(szdim,eigen_vectors(:,1),incx,eigen_vectors(:,1),incy)))
  eigen_vectors(:,1)=eigen_vectors(:,1)/s_tmp

  call zcsrsymv(eigen_vectors(:,1),y_tmp)

  exp_val(1)=dble(zdotc(szdim,eigen_vectors(:,1),incx,y_tmp,incy))

  
  y_tmp=(0.0d0,0.0d0)
  s_tmp=0.0d0
  s_tmp=sqrt(dble(zdotc(szdim,eigen_vectors(:,2),incx,eigen_vectors(:,2),incy)))
  eigen_vectors(:,2)=eigen_vectors(:,2)/s_tmp

  call zcsrsymv(eigen_vectors(:,2),y_tmp)

  exp_val(2)=dble(zdotc(szdim,eigen_vectors(:,2),incx,y_tmp,incy))


end subroutine eigen_vector

subroutine zcsrsymv(x,y)
  use variables
  implicit none
  integer i,j,tmp
  complex(kind(0d0)) x(szdim),y(szdim)

  do i=1,szdim
    y(i)=y(i)+sp_ele(ia(i))*x(ja(ia(i)))
    do j=1,ia(i+1)-ia(i)-1
      tmp=ia(i)+j
      y(i)=y(i)+sp_ele(tmp)*x(ja(tmp))
      y(ja(tmp))=y(ja(tmp))+conjg(sp_ele(tmp))*x(i)
    end do
 end do

end subroutine

! subroutine zcsrsymv(uplo,x,y)
!   use variables
!   implicit none
!   integer i,j,m,tmp
!   complex(kind(0d0)) x(m),y(m)
!   complex(kind(0d0)) a(*)
!   integer ia(m/2+1),ja(*)
!   character(len=1) uplo

!   tmp=0
!   do i=1,m/2
!      tmp=tmp+1
!      y(i)=y(i)+a(tmp)*x(ja(tmp))
!      y(m-ja(tmp)+1)=y(m-ja(tmp)+1)+a(tmp)*x(m-i+1) !szdim/2 ver
!      do j=1,ia(i+1)-ia(i)-1
!         tmp=tmp+1
!         y(i)=y(i)+a(tmp)*x(ja(tmp))
!         y(ja(tmp))=y(ja(tmp))+conjg(a(tmp))*x(i)

!         y(m-ja(tmp)+1)=y(m-ja(tmp)+1)+a(tmp)*x(m-i+1) !szdim/2 ver
!         y(m-i+1)=y(m-i+1)+conjg(a(tmp))*x(m-ja(tmp)+1) !szdim/2 ver
!      end do
!   end do

! end subroutine zcsrsymv

subroutine first_vector(v_first,szdim)
  implicit none
  integer :: szdim,i,incx=1,incy=1
  double precision norm_tmp
  complex(kind(0d0)) :: zdotc!,ei=(0.0d0,1.0d0)
  complex(kind(0d0)) v_first(szdim)
  double precision,allocatable :: v_rand(:)
  integer :: seedsize
  integer, allocatable :: seed(:)


  allocate(v_rand(2*szdim))

  call random_seed(size=seedsize) !初期値のサイズを取得
  allocate(seed(seedsize)) !配列の割り当て
  do i = 1, seedsize
     call system_clock(count=seed(i)) !時間を取得
  end do
  call random_seed(put=seed(:)) !初期値を与える

  call random_number(v_rand)

  do i=1,szdim
!     v_first(i)=v_rand(i)+ei*v_rand(i+szdim)
     v_first(i)=v_rand(i)
  end do

  norm_tmp=sqrt(dble(zdotc(szdim,v_first,incx,v_first,incy)))

  v_first=v_first/norm_tmp

  ! write(*,*) v_first
  ! write(*,*) sqrt(dble(zdotc(szdim,v_first,incx,v_first,incy)))

  !   v_first(szdim/2)=(1.0d0,0.0d0)


  end subroutine
