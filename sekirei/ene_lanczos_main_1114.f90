
! XXZ Model Hamiltonian (delta,bondj)
module variables
  implicit none
  integer :: spinmax,phase_num
  integer*8 :: sp_size,szdim
  double precision :: arg
  integer,allocatable :: list(:,:),outpt(:),outpt_inv(:)
  integer,allocatable :: ia(:),ja(:),phase_i(:),phase_j(:),phase_twist(:)

  double precision,allocatable :: energy_theta(:,:)
  double precision,allocatable :: twi_unit(:)

  integer :: n_pair
  integer,allocatable :: i_pair(:),twi_pair(:)
  double precision,allocatable :: bondj(:),delta(:)

  complex(kind(0d0)),allocatable :: sp_ele(:)

  double precision,parameter :: pi=atan(1.0d0)*4.0d0
end module variables

module para
  implicit none
  include 'mpif.h'
  integer :: myrank,nprocs
  double precision, allocatable :: arg_para(:),arg_atempt(:)
  double precision, allocatable :: energy_para(:,:,:),energy_atempt(:,:,:)
contains

  subroutine mpi_initialize
    implicit none
    integer :: ierr
    call mpi_init(ierr)
    call mpi_comm_rank(mpi_comm_world,myrank,ierr)
    call mpi_comm_size(mpi_comm_world,nprocs,ierr)

    if(myrank==0) then
       write(*,*)'Num of Parallel calculations: ',nprocs
    end if

  end subroutine mpi_initialize

  subroutine mpi_fin
    implicit none
    integer::ierr

    call mpi_finalize(ierr)
  end subroutine mpi_fin
end module para

program main
  use variables
  implicit none
  integer :: b_para_loop_max
!  integer :: i,j,k
  double precision :: b_para_del,b_para_min,b_para_max
  complex(kind(0d0)),parameter :: ei=(0.0d0,1.0d0)
  !--------------------parameter--------------------------

  !**** if n is large, result would be over-flow ***
  integer,parameter :: del=6
  double precision,parameter :: spin=2.0d0,sztot=0.0d0
  !   double precision,parameter :: pi=atan(1.0d0)*4.0d0

  integer,parameter :: n=11,DoF=2*spin+1,ene_num=2
!  integer,parameter :: twi_num=2
  integer*8,parameter :: dim=DoF**n
  integer,parameter :: subspace=1000
!  character(10) :: twist_type = 'mila1'
!  character(10) :: twist_type = 'j1'
!  character(10) :: twist_type = 'moebius3-3'
  character(10) :: twist_type = 'j2'
!  character(10) :: twist_type = 'moebius3-4'
  character(20) :: model = 'j1j2'

  !************************************************************

  !----------------------------array--------------------------------------

!  integer,allocatable :: i_pair(:),twi_pair(:)
  !When changing the argument of pair, the value of subroutine(matrix_term)
  ! must also be changed
!  double precision,allocatable :: bondj(:),delta(:)
!  double precision :: corr_z(n-1),corr_xy(n-1),energy(ene_num)

  !*****************************MKL*********************************

  integer info

  !********************system clock**************************

  integer t1, t2, t_rate, t_max, diff
  call system_clock(t1)

  allocate(energy_theta(ene_num,0:del))
  energy_theta=0.0d0

  call pair_format(model,twist_type,n,del)

  !**************************** sz *****************************************************
  call sz(n,dim,DoF,spin,sztot)

!  write(*,*) outpt

  write(*,'(A,I4,A,f8.4)') " site= ",n,"  Sz_total= ",sztot
  write(*,'(A,I16)')" Sz_dim= ",szdim

  !************************ sparse matrix elements ***********************************

  b_para_loop_max=9
  b_para_min=tan(0.65d0)
  b_para_max=tan(0.68d0)
  b_para_del=0.0d0
  if (b_para_loop_max .ne. 0)  b_para_del=(b_para_max-b_para_min)/b_para_loop_max
  ! write(*,*)
  ! write(*,*) 'berry phase - split number',del
  ! write(*,310) 'bond parameter - split number ',b_para_loop_max,b_para_min,' ~ ',b_para_max
  ! write(*,*)
310  format(1X,A,I6,2X,f6.3,A,f6.3)

!  write(*,*) "[START matrix term]"


  call sparse_matrix(n,spin,ene_num,del,twist_type,subspace,b_para_loop_max,b_para_min,b_para_max,b_para_del)

!  write(*,*) "[END matrix term]"

!  write(*,*)"num of matrix elements",sp_size

  !***************************** out put ***********************************************
!  if(info.ne.0) write(*,*) "info",info

  !***************************** run time **********************************************
  call system_clock(t2, t_rate, t_max)   ! 終了時を記録
  if ( t2 < t1 ) then
     diff = (t_max - t1) + t2 + 1
  else
     diff = t2 - t1
  endif
  write(*,*)"[run time]",diff/dble(t_rate)

120 format(f5.2,i3,f25.15,f25.15,f15.5)
130 format(a3,f5.2,i3,f7.3)
140 format(f7.3,f13.8,f25.15,f25.15,f25.15)
150 format(f7.3,f13.8,f25.15,f25.15)

end program main
