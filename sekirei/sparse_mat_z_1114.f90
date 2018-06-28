!nearest neighbor
subroutine sparse_matrix(n,spin,ene_num,del,twist_type,subspace,b_para_loop_max,b_para_min,b_para_max,b_para_del)
  use para
  use variables
  implicit none
  integer :: n,Dof,ind,b_para_loop_max,bond_para,subspace,del,ene_num
  integer*8 :: dupl
  integer :: i,j,ii,tmp_pair1,tmp_pair2,para_ii
  double precision :: spin,half_c,spin_flip,b_para_min,b_para_max,b_para_del,tmp_arg,ele_tmp
  double precision :: phi,del_phi
  double precision,parameter :: eps=1.0E-4

  complex(kind(0d0)) :: ei
!  complex(kind(0d0)) :: const(0:10)

  integer,allocatable :: listud(:,:,:)
  integer,allocatable :: listud_ii(:,:)
  integer,allocatable :: listud_sort(:,:),tmp_sort(:)
  integer,allocatable :: tmp_col(:)
  complex(kind(0d0)),allocatable :: matrix_tmp(:)
  character(10) twist_type

  !
  integer nbun,ista,iend,ierr

  write(*,*) "[START matrix term]"

  allocate(matrix_tmp(2*n_pair),tmp_col(szdim))

  Dof=2*spin+1

  if (mod(Dof-1,2).eq.1) then
     half_c=0.5d0
  else
     half_c=1.0d0
  end if

  allocate(listud(szdim,n_pair,2))
  allocate(listud_sort(szdim,n_pair*2),tmp_sort(n_pair*2))
  allocate(listud_ii(szdim,n_pair))
  allocate(phase_i(szdim*n_pair),phase_j(szdim*n_pair),phase_twist(szdim*n_pair))
!  allocate(phase_i(szdim*n_pair*2/n),phase_j(szdim*n_pair*2/n),phase_twist(szdim*n_pair*2/n))


  listud=0
  tmp_pair1=0
  tmp_pair2=0
  tmp_col=0
  listud_sort=szdim+1
  phase_num=0
  ind=0
  phase_twist=0

  do i=1,szdim
     do ii=1,n_pair
        tmp_pair1=i_pair(2*ii-1)
        tmp_pair2=i_pair(2*ii)

        if(list(i,tmp_pair1).ne.spinmax .and. list(i,tmp_pair2).ne. (-1)*spinmax) then
           listud(i,ii,1)=list(i,n+1) +Dof**(n-tmp_pair1) -Dof**(n-tmp_pair2)
           if(outpt_inv(listud(i,ii,1))>i) then

              tmp_col(i)=tmp_col(i)+1
              listud_sort(i,tmp_col(i))=outpt_inv(listud(i,ii,1))
              listud_ii(i,tmp_col(i))=ii
              ind=ind+1
           end if
        end if
        if(list(i,tmp_pair2).ne.spinmax .and. list(i,tmp_pair1).ne. (-1)*spinmax) then
           listud(i,ii,2)= list(i,n+1) +Dof**(n-tmp_pair2) -Dof**(n-tmp_pair1)
!           if(twi_pair(ii) .ne. 0) then
           if(outpt_inv(listud(i,ii,2))>i) then

              tmp_col(i)=tmp_col(i)+1
              listud_sort(i,tmp_col(i))=outpt_inv(listud(i,ii,2))
              listud_ii(i,tmp_col(i))=ii
              phase_num=phase_num+1
              phase_i(phase_num)=i
              phase_j(phase_num)=outpt_inv(listud(i,ii,2))
              phase_twist(phase_num)=twi_pair(ii)
              ind=ind+1
           end if
        end if
     end do
  end do

!  write(*,*)'phase change number', phase_num

  ind=0
  !$omp parallel
  !$omp do private(j,tmp_sort,dupl)
  do i=1,szdim
     dupl=0

     do j=1,n_pair*2
        tmp_sort(j)=listud_sort(i,j)
     end do

     call quicksort(tmp_sort,1,n_pair*2)

     listud_sort(i,:)=szdim+1
     j=1
     do while(tmp_sort(j).ne.szdim+1)
        listud_sort(i,j)=tmp_sort(j)

        if(tmp_sort(j)==tmp_sort(j+1)) dupl=dupl+1 ! minus duplication
        j=j+1
     end do
!     if(dupl.ne.0) write(*,*) dupl
     tmp_col(i)=tmp_col(i)-dupl
     ind=ind+tmp_col(i)
  end do
  !$omp end do
  !$omp end parallel


  allocate(sp_ele(ind+szdim),ja(ind+szdim),ia(szdim+1))
  sp_size=ind+szdim
!  write(*,*)"num of matrix elements",sp_size


  ind=0
  do i=1,szdim
     ia(i)=ind+i
     ind=ind+tmp_col(i)
  end do
  ia(szdim+1)=ind+szdim+1

!  write(*,*)b_para_loop_max
  del_phi=0.0d0
  if(b_para_loop_max>0)then
     del_phi=(atan2(b_para_max,1.0d0)-atan2(b_para_min,1.0d0))/dble(b_para_loop_max)
  end if
  if(b_para_loop_max<0) then
     b_para_loop_max=-1*b_para_loop_max
     del_phi=pi/dble(b_para_loop_max)/2.0d0
     b_para_min=0.0d0
  end if
!  write(*,*)b_para_loop_max
  !******************************************
  !************** mpi ****************************
  allocate(arg_para(0:b_para_loop_max),energy_para(0:b_para_loop_max,ene_num,0:del))
  call mpi_initialize
  nbun=(b_para_loop_max+1)/nprocs
  ista=myrank*nbun
  iend=ista+nbun-1
  para_ii=0
  allocate(arg_atempt(nbun),energy_atempt(nbun,ene_num,0:del))
  !******************************************
  !******************************************
  loop:do bond_para=ista,iend
     para_ii=para_ii+1
!  do bond_para=0,b_para_loop_max
     phi=del_phi*real(bond_para)+atan2(b_para_min,1.0d0)
     if(myrank==0) write(*,*) para_ii,phi
!     phi=atan2(b_para_min+dble(bond_para)*b_para_del,1.0d0)
     if(n_pair == 2*n) then
        bondj(1:n)     = cos(phi)
        bondj(n+1:2*n) = sin(phi)
!        write(*,'(A,f8.4,A,f8.4,A,f8.4)') ' J1= ',bondj(1),'  j2= ',bondj(n+1),'  ratio(j2/j1) = ',bondj(n+1)/bondj(1)
     end if

     ind=0
     sp_ele=(0.0d0,0.0d0)
     ja=0
     !$omp parallel
     !$omp do private(ii,ele_tmp,matrix_tmp,j,ind,tmp_pair1,tmp_pair2,dupl) schedule(static,1)

     do i=1,szdim
        ! diagonal terms-------------------------------------------------------
        ja(ia(i))=i
        do ii=1,n_pair
           ele_tmp=0.0d0
           tmp_pair1=i_pair(2*ii-1)
           tmp_pair2=i_pair(2*ii)

           ele_tmp=list(i,tmp_pair1)*list(i,tmp_pair2) *half_c*half_c &
                *delta(ii) *bondj(ii) !+ei*(0.0d0)

           sp_ele(ia(i))=sp_ele(ia(i))+ele_tmp
        end do
        !non diagonal terms-----------------------------------------------------
        matrix_tmp=0.0d0
        j=1
        do while(listud_sort(i,j).ne.szdim+1)! j=1,tmp_col(i)
           !        ii=listud_ii(i,j)
           do ii=1,n_pair
              tmp_pair1=i_pair(2*ii-1)
              tmp_pair2=i_pair(2*ii)
              select case(listud_sort(i,j)-outpt_inv(listud(i,ii,1)))
              case(0)
                 matrix_tmp(j) = matrix_tmp(j) + (bondj(ii)) * 0.5d0 * spin_flip(i,spin,tmp_pair1,tmp_pair2,half_c)
              end select

              select case(listud_sort(i,j)-outpt_inv(listud(i,ii,2)))
              case(0)
                 matrix_tmp(j) = matrix_tmp(j) + (bondj(ii)) * 0.5d0 * spin_flip(i,spin,tmp_pair2,tmp_pair1,half_c)
              end select
           end do
           j=j+1
        end do

        j=1
        dupl=0
        do while(listud_sort(i,j).ne.szdim+1)!j=1,tmp_col(i)
           ind=ia(i)+j-dupl
           sp_ele(ind)=matrix_tmp(j)
           ja(ind)=listud_sort(i,j)

           if(listud_sort(i,j).eq.listud_sort(i,j+1)) dupl=dupl+1
           j=j+1
        end do
     end do
     !$omp end do
     !$omp end parallel

     ! write(*,*) "[END matrix term]"
     ! write(*,*)
!      call real_matrix()

     call phase(n,spin,del,subspace,ene_num,phi,twist_type)

!*******************nibunnho**********************************************************************************

     tmp_arg=0.0d0
     if(abs(mod(arg/pi,2.0d0))<1.8d0)then
        tmp_arg = abs(mod(arg/pi,2.0d0))
     else
        tmp_arg = 2.0d0-abs(mod(arg/pi,2.0d0))
     end if

     write(*,'(1X,I5,f16.12,A,f16.12)') bond_para,phi,' Berry phase',tmp_arg
     ! write(*,*) energy_theta(1,0),energy_theta(2,0),-energy_theta(1,0)+energy_theta(2,0)
      write(*,*)

     arg_atempt(para_ii)=arg
!     arg_para(para_ii-1+ista)=arg
     energy_atempt(para_ii,:,:)=energy_theta(:,:)
!     energy_para(para_ii+ista-1,1:ene_num,0:del) = energy_theta(1:ene_num,0:del)


     ! if(tmp_arg > 0.9) then
     !    b_para_min=bondj_2
     ! else
     !    b_para_max=bondj_2
     ! end if

     ! if(abs(bondj_2-(b_para_min+b_para_max)/2.0d0)/bondj_2 < eps)then
     !    return
     ! else
     !    bondj_2 = (b_para_min+b_para_max)/2.0d0
     ! end if

!************************************************************************************************************

  end do loop
  call mpi_allgather(arg_atempt,nbun,MPI_DOUBLE_precision,arg_para,nbun,MPI_DOUBLE_precision,MPI_COMM_WORLD,ierr)
!  call mpi_allgather(energy_atempt,ene_num*(nbun)*(del+1),MPI_DOUBLE_precision,energy_para,ene_num*(nbun)*(del+1),MPI_DOUBLE_precision,MPI_COMM_WORLD,ierr)

  do i = 1,ene_num
     do j = 0,del
        call mpi_allgather(energy_atempt(:,i,j),nbun,MPI_DOUBLE_precision,energy_para(:,i,j),nbun,MPI_DOUBLE_precision,MPI_COMM_WORLD,ierr)
     end do
  end do

  call mpi_fin
  ! if(myrank==0) then
  !    write(*,*)
  !    do i=0,b_para_loop_max
  !       do j = 0,del
  !          write(*,*) energy_para(i,1,j),energy_para(i,2,j)
  !       end do
  !       write(*,*)
  !    end do
  ! end if

  ! stop


  if(myrank==0) then
!     write(*,'(F16.8)') arg_para
     call out_put(b_para_loop_max,b_para_min,twist_type,n,spin,del,ene_num,del_phi)
  end if



    ! write(*,*) ia
    ! write(*,*)
    ! write(*,*) ja
    ! write(*,*)
    ! write(*,*) dble(sp_ele)

end subroutine sparse_matrix

subroutine out_put(b_para_loop_max,b_para_min,twist_type,n,spin,del,ene_num,del_phi)
  use variables
  use para
  implicit none
  integer b_para_loop_max,n,del,ene_num,i,j
  integer,parameter :: ene = 0
  double precision spin,phi,del_phi,b_para_min

  character :: fn_bp*40,fn_ene*40,fn_spin*20,fn_n*20,fn_add*20='2'
  character(10) :: twist_type



  if (fn_add .ne. 'none') then
     write(fn_spin,*) int(spin)
     write(fn_n,*) int(n)
     if(fn_add .ne. '0') then
        fn_bp='spin-'//trim(adjustl(fn_spin))//'-n-'//trim(adjustl(fn_n))//'-bp-'//trim(adjustl(twist_type))//'-'//trim(adjustl(fn_add))//''
        fn_ene='spin-'//trim(adjustl(fn_spin))//'-n-'//trim(adjustl(fn_n))//'-ene-'//trim(adjustl(twist_type))//'-'//trim(adjustl(fn_add))//''
     else
        fn_bp='spin-'//trim(adjustl(fn_spin))//'-n-'//trim(adjustl(fn_n))//'-bp-'//trim(adjustl(twist_type))//''
        fn_ene='spin-'//trim(adjustl(fn_spin))//'-n-'//trim(adjustl(fn_n))//'-ene-'//trim(adjustl(twist_type))//''
     end if
     !  write(*,*) fn_bp,fn_ene
     !  stop

     open(21,file=fn_bp,position='append')
     if(ene==1) open(32,file=fn_ene,position='append')
     do i=0,b_para_loop_max
        phi=del_phi*real(i)+atan2(b_para_min,1.0d0)
        arg=arg_para(i)
        energy_theta(:,:)=energy_para(i,:,:)

        if(abs(mod(arg/pi,2.0d0))<1.8d0)then
           write(21,210) i,phi,arg,abs(mod(arg/pi,2.0d0))
        else
           write(21,210) i,phi,arg,2.0d0-abs(mod(arg/pi,2.0d0))
        end if

        if(ene==1) then
           write(32,130) '#',spin,n,phi
           if(ene_num==3)then
              do j=0,del
                 write(32,140) phi,dble(j)*2.0d0*pi/dble(del),energy_theta(1,j),energy_theta(2,j),energy_theta(3,j)
              end do
           else if(ene_num==2)then
              do j=0,del
                 write(32,140) phi,dble(j)*2.0d0*pi/dble(del),energy_theta(1,j),energy_theta(2,j),-energy_theta(1,j)+energy_theta(2,j)
              end do
           end if


           write(32,*)
           write(32,*)
        end if

     end do

     close(21)
     if(ene==1) close(32)

  end if


100 format(f16.10,f20.12,f20.12)
110 format(f6.2,i5,f20.12,f20.12)
130 format(a3,f5.2,i3,f7.3)
210 format(I5,f8.4,f20.12,f20.12)
140 format(f8.4,f13.8,f25.15,f25.15,f25.15)
150 format(f8.4,f13.8,f25.15,f25.15)
end subroutine out_put



subroutine phase_change()
  use variables
  implicit none
  integer i,j
  integer*8 :: ia_un,ia_up!,dim
  complex(kind(0d0)) ei,del_phase


  ei=(0.0d0,1.0d0)

  do i=1,phase_num
     ia_un=ia(phase_i(i))
     ia_up=ia(phase_i(i)+1)-1
     do j=ia_un,ia_up
        select case(ja(j)-phase_j(i))
        case(0)
           del_phase=exp(ei*twi_unit(phase_twist(i)))
           sp_ele(j)=sp_ele(j)*del_phase
           exit
        end select
     end do
  end do

end subroutine phase_change

function spin_flip(i,spin,site1,site2,half_c)
  use variables
  implicit none
  double precision :: spin_flip,tmp1,tmp2
  integer :: site1,site2,i
  double precision :: spin,half_c

  select case(int(spin*2.0d0))
  case(1) !spin 1/2
     spin_flip=1.0d0

  case(2) !spin 1
     spin_flip=2.0d0

  case(3) !spin 3/2
     select case(list(i,site1))
     case(-1)
        tmp1=2.0d0
     case(1)
        tmp1=sqrt(3.0d0)
     case(-3)
        tmp1=sqrt(3.0d0)
     end select
     select case(list(i,site2))
     case(3)
        tmp2=sqrt(3.0d0)
     case(1)
        tmp2=2.0d0
     case(-1)
        tmp2=sqrt(3.0d0)
     end select
     spin_flip=tmp1*tmp2
  case(4:)
     spin_flip=sqrt((spin-half_c*list(i,site1))    &
                   *(spin+half_c*list(i,site1)+1)  &
                   *(spin+half_c*list(i,site2))    &
                   *(spin-half_c*list(i,site2)+1))
  end select

end function spin_flip

subroutine real_matrix()
  use variables
  implicit none

  integer i,j,tmp
  complex(kind(0d0)),allocatable :: real_mat(:,:)

  allocate(real_mat(szdim,szdim))
  real_mat=(0.0d0,0.0d0)
  tmp=0
  do i=1,szdim
!        write(*,*) tmp_col(i)+1
     do j=ia(i),ia(i+1)-1! i,tmp_col(i)+1
        tmp=tmp+1
        real_mat(i,ja(tmp))=sp_ele(tmp)
        real_mat(ja(tmp),i)=conjg(sp_ele(tmp))
     end do
  end do


  open(33,file='check_mat_j1j2',status='replace')
!   do i=1,szdim
! !     do j=1,szdim
!         write(33,310,advance='no') dble(real_mat(i,i))
! !     end do
! !     write(33,*)
!      end do

  do i = 1, szdim
     write(33, '(100f8.2)')  dble(real_mat(i, 1:szdim))
  enddo
  close(33)
  stop
310 format(f5.1)

end subroutine real_matrix
