subroutine phase(n,spin,del,subspace,ene_num,phi,twist_type)
  use variables
  implicit none
  integer i,del,subspace,info,incx,incy,n,ene_num
!  integer *8 :: dim
  double precision :: del_theta,theta,spin,phi,total_sz
  complex(kind(0d0)) zdotc!,zdotu

  complex(kind(0d0)) prod(del)
  complex(kind(0d0)),allocatable :: eigen_vec(:),eigen_vec_2(:)
  complex(kind(0d0)),allocatable :: tmp_eigen_vec_0(:),tmp_eigen_vec(:)
  complex(kind(0d0)),allocatable :: v_first(:)
  double precision energy(ene_num)
!  double precision s2_xy(n),s2_z(n)
  character :: fn_bp*20,fn_ene*20,fn_spin*20,fn_n*20,fn_add*20='test'
  character(10) :: twist_type


!  write(*,*)"[start lanczos method]"
  
  allocate(eigen_vec(szdim),eigen_vec_2(szdim),tmp_eigen_vec_0(szdim),tmp_eigen_vec(szdim))
  allocate(v_first(szdim))


  if(del.ne.0) then
     del_theta=2.0d0*pi/dble(del)
  end if

  tmp_eigen_vec=(0.0d0,0.0d0)
  eigen_vec=(0.0d0,0.0d0)
  eigen_vec_2=(0.0d0,0.0d0)
  v_first=(0.0d0,0.0d0)
  incx=1
  incy=1

  call first_vector(v_first,szdim)
!  v_first(szdim/3)=(1.0d0,0.0d0)

!  write(*,*)
  call lanczos_main(subspace,v_first,energy,ene_num,eigen_vec,eigen_vec_2,info)
  !****************************************************************************************************
  ! total Sz
  ! do i=1,n
  !    write(*,*)i, total_sz(eigen_vec,i)
  ! end do
  ! stop
  !****************************************************************************************************
  
  ! call s2(n,spin,eigen_vec,s2_xy,s2_z)

  ! write(*,*) 'SS'
  ! do i=1,n/2
  !    write(*,fmt='(f10.5)',advance='no') s2_xy(i) + s2_z(i)
  ! end do
  
  ! call s2(n,spin,eigen_vec_2,s2_xy,s2_z)
  ! write(*,*)
  ! write(*,*) 'SS excited'
  ! do i=1,n/2
  !    write(*,fmt='(f10.5)',advance='no') s2_xy(i) + s2_z(i)
  ! end do
  ! write(*,*)  


!******************************************************************************************
  ! call phase_change()    
  ! do i=1,5
  !    v_first=(0.0d0,0.0d0)
  !    tmp_eigen_vec=eigen_vec
  !    eigen_vec_2=(0.0d0,0.0d0)
  !    call first_vector(v_first,szdim)

  !    call lanczos_main(subspace,v_first,energy,ene_num,eigen_vec,eigen_vec_2,info)
     
  !    write(*,fmt='(F10.5,F10.5)',advance='no')energy

  !    write(*,*)zdotc(szdim,eigen_vec,incx,tmp_eigen_vec,incy)
  ! end do
  ! stop
!******************************************************************************************

  
!  call correlation(eigen_vec,n,spin)
  tmp_eigen_vec_0=eigen_vec
  energy_theta(:,0)=energy(:)
  theta=0.0d0

      ! open(11,file='energy-int-theta10',position='append')
      ! write(11,100) theta,energy_theta(1,0),energy_theta(2,0)
      ! close(11)
  arg=0.0d0

  ! del_theta=2*pi/10.0d0
  ! call phase_change(dim,del_theta)
  ! write(*,*) sp_ele
  ! stop
  
  do i=1,del
!     write(*,*) 'del=',i
     theta=dble(i)*del_theta
     tmp_eigen_vec=eigen_vec

     call phase_change()

     call lanczos_main(subspace,v_first,energy,ene_num,eigen_vec,eigen_vec_2,info)

     prod(i)=zdotc(szdim,tmp_eigen_vec,incx,eigen_vec,incy)
     arg= arg+atan2(aimag(prod(i)),dble(prod(i)))
     
     ! open(61,file='berry_phase_del_spin1-7',position='append')
     ! if(mod(arg/pi,2.0d0)>0.0d0)then
     !    write(61,'(F8.5,F16.10)') theta,(arg/pi)
     ! else
     !    write(61,'(F8.5,F16.10)') theta,(arg/pi)
     ! end if
     ! close(61)

!     write(*,*) atan2(aimag(prod(i)),dble(prod(i)))

     energy_theta(:,i)=energy(:)

     ! open(11,file='energy-int-theta10',position='append')
     ! write(11,100) theta,energy_theta(1,i),energy_theta(2,i)
     ! close(11)

  end do

!  write(*,*) "arg",arg,mod(arg/pi,2.0d0)
  !  stop

  !*******************************************************************************************************************
  !************************************************** out put ************************************************************

  ! if (fn_add .ne. 'none') then
  !    write(fn_spin,*) int(spin)
  !    write(fn_n,*) int(n)
  !    if(fn_add .ne. '0') then
  !       fn_bp='spin-'//trim(adjustl(fn_spin))//'-n-'//trim(adjustl(fn_n))//'-bp-'//trim(adjustl(twist_type))//'-'//trim(adjustl(fn_add))//''
  !       fn_ene='spin-'//trim(adjustl(fn_spin))//'-n-'//trim(adjustl(fn_n))//'-ene-'//trim(adjustl(twist_type))//'-'//trim(adjustl(fn_add))//''
  !    else
  !       fn_bp='spin-'//trim(adjustl(fn_spin))//'-n-'//trim(adjustl(fn_n))//'-bp-'//trim(adjustl(twist_type))//''
  !       fn_ene='spin-'//trim(adjustl(fn_spin))//'-n-'//trim(adjustl(fn_n))//'-ene-'//trim(adjustl(twist_type))//''
  !    end if
     !  write(*,*) fn_bp,fn_ene
     !  stop

!     open(21,file=fn_bp,position='append')
!     open(32,file=fn_ene,position='append')
     ! if(abs(mod(arg/pi,2.0d0))<1.8d0)then
     !    write(21,210) phi,arg,abs(mod(arg/pi,2.0d0))
     ! else
     !    write(21,210) phi,arg,2.0d0-abs(mod(arg/pi,2.0d0))
     ! end if

     ! write(32,130) '#',spin,n,phi
     ! if(ene_num==3)then
     !    do i=0,del
     !       write(32,140) phi,dble(i)*2.0d0*pi/dble(del),energy_theta(1,i),energy_theta(2,i),energy_theta(3,i)
     !    end do
     ! else if(ene_num==2)then
     !    do i=0,del
     !       write(32,140) phi,dble(i)*2.0d0*pi/dble(del),energy_theta(1,i),energy_theta(2,i),-energy_theta(1,i)+energy_theta(2,i)
     !    end do
     ! end if


     ! write(32,*)
     ! write(32,*)

!     close(21)
!     close(32)
!end if

  !*******************************************************************************************************************
  !************************************************** out put ************************************************************

  
  if(info.ne.0) write(*,*) "info",info

!  write(*,*)"[end lanczos method]"
  
  
100 format(f16.10,f20.12,f20.12)
110 format(f6.2,i5,f20.12,f20.12)
130 format(a3,f5.2,i3,f7.3)
210 format(f8.4,f20.12,f20.12)
140 format(f8.4,f13.8,f25.15,f25.15,f25.15)
150 format(f8.4,f13.8,f25.15,f25.15)

  
end subroutine phase
