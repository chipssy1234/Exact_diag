subroutine pair_format(model,twist_type,n,del)
  use variables
  implicit none
  integer i,j,k
  integer n,twi_num,del
  double precision :: bondj_tmp
  
  character model*20,twist_type*10
  
  !************************* chain *****************************************
  if(model == 'chain')then
     
     n_pair=n
     allocate(i_pair(n_pair*2),twi_pair(n_pair), bondj(n_pair),delta(n_pair))

     delta=1.0d0
     bondj=0.0d0
     bondj_tmp=1.0d0
     twi_num=1

     do i=1,n
        i_pair(2*i-1)=i
        if (i.ne.n) then
           i_pair(2*i)=i+1
           bondj(i)=bondj_tmp
        else
           i_pair(2*i)=1
           bondj(i)=bondj_tmp
        end if
     end do

     allocate(twi_unit(0:twi_num))
     twi_unit(0)=0.0d0
     if(del.ne.0) then
        twi_unit(1)=2.0d0*pi/dble(del)
     end if

     twi_pair = 0

     if(twist_type == 'j1') then
        twi_pair(n) = 1
     end if


     do i=1,n_pair
        if(twi_pair(i) .ne. 0) then
           write(*,*) i,'site is chaged for',twi_pair(i)
           j=i_pair(2*i-1)
           k=i_pair(2*i  )
           i_pair(2*i-1)=max(j,k)
           i_pair(2*i  )=min(j,k)
        end if
     end do
  end if

  !************************* chain *****************************************
  if(model == 'open_ch')then
     
     n_pair=n-1
     allocate(i_pair(n_pair*2),twi_pair(n_pair), bondj(n_pair),delta(n_pair))

     delta=1.0d0
     bondj=0.0d0
     bondj_tmp=1.0d0
     twi_num=1

     do i=1,n_pair
        i_pair(2*i-1)=i
        if (i.ne.n) then
           i_pair(2*i)=i+1
           bondj(i)=bondj_tmp
        else
           i_pair(2*i)=1
           bondj(i)=bondj_tmp
        end if
     end do

  end if
  
  !************************* j1_j2 *****************************************
  if(model == 'j1j2')then

     n_pair=2*n
     allocate(i_pair(n_pair*2),twi_pair(n_pair), bondj(n_pair),delta(n_pair))

     delta=1.0d0
     bondj_tmp=1.0d0
     twi_num=2
     
     do i=1,n
        if (i.ne.n) then
           i_pair(2*i-1)=i
           i_pair(2*i)=i+1
           bondj(i)=bondj_tmp
        else
           i_pair(2*i-1)=1
           i_pair(2*i)=i
           bondj(i)=bondj_tmp
        end if
     end do

     do i=1,n
        if (i.le.n-2) then
           i_pair(2*n+2*i-1)=i
           i_pair(2*n+2*i)=i+2
        else
           i_pair(2*n+2*i-1)=i+2-n
           i_pair(2*n+2*i)=i
        end if
     end do
   
     allocate(twi_unit(0:twi_num))
     twi_unit(0)=0.0d0
     if(del.ne.0) then
        twi_unit(1)=2.0d0*pi/dble(del)
        twi_unit(2)=-2.0d0*pi/dble(del)
     end if

     twi_pair = 0
     if(twist_type == 'j1') then
        twi_pair(n) = 1
     end if
     if(twist_type == 'j2') then
        twi_pair(2*n) = 1
!        twi_pair(2*n-1) =1
     end if
     if(twist_type == 'cycle') then
        twi_pair(2*n-1) = 2
        twi_pair(n) = 1
        twi_pair(n-1) = 1
     end if
     if(twist_type == 'cycle2') then
        twi_pair(2*n-1) = 1
        twi_pair(n) = 1
        twi_pair(n-1) = 1
        twi_pair(2*n) = 1
        twi_pair(2*n-2) = 1
     end if
     if(twist_type == '3hon') then
        twi_pair(n-1) = 1
        twi_pair(n) = 1
        twi_pair(2*n) = 1
     end if
     if(twist_type == 'mila1') then
        twi_pair(2*n-2) = 1
        twi_pair(n-1) = 2
        twi_pair(2*n-1) = 1
     end if
     if(twist_type == 'mila2') then
        twi_pair(n) = 1
        twi_pair(2*n) = 1
     end if
     if(twist_type == 'moebius1') then
        twi_pair(1) = 1
        twi_pair(2) = 2
        twi_pair(3) = 1
        twi_pair(n) = 2
        twi_pair(2*n-2) = 1
     end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!! all j2 twist !!!!!!!!!!!!     
     if(twist_type == 'moebius2') then
        twi_pair(n+1:2*n)=1
     end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     
     if(twist_type == 'moebius3') then

        twi_pair(2*n-2) = 1
        do i=1,n-1
           if(mod(i,2)==1)then
              twi_pair(i) =1
           else
              twi_pair(i)=2
           end if
        end do
     end if
     if(twist_type == 'moebius3-2') then
        twi_pair(2*n) = 1
        do i=1,n-1
           if(mod(i,2)==1)then
              twi_pair(i) =1
           else
              twi_pair(i)=2
           end if
        end do
     end if

     if(twist_type == 'moebius3-3') then
        twi_pair(2*n-2) = 1
        twi_pair(n) = 2
        do i=1,n-2
           if(mod(i,2)==1)then
              twi_pair(i) =1
           else
              twi_pair(i)=2
           end if
        end do
     end if

     if(twist_type == 'moebius3-4') then
        twi_pair(1) = 1
        twi_pair(n) = 2
        twi_pair(n+1) =1
     end if     
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!! all j1 twist!!!!!!!!
     if(twist_type == 'moebius4') then
         twi_pair(1:n) = 1
     end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     
     do i=1,n_pair
        if(twi_pair(i) .ne. 0) then
           write(*,*) i,'site is chaged for',twi_pair(i)
           j=i_pair(2*i-1)
           k=i_pair(2*i  )
           i_pair(2*i-1)=max(j,k)
           i_pair(2*i  )=min(j,k)
        end if
     end do
  end if

  !************************* j1-j2 model open B.C.  *****************************************
  if(model == 'j1j2_open')then

     n_pair=2*n-3
     allocate(i_pair(n_pair*2),twi_pair(n_pair), bondj(n_pair),delta(n_pair))

     delta=1.0d0
     bondj_tmp=1.0d0
!NN
     do i=1,n-1
           i_pair(2*i-1)=i
           i_pair(2*i)=i+1
           bondj(i)=bondj_tmp
           bondj(i)=bondj_tmp
     end do
     !NNN
     do i=1,n-2
           i_pair(2*n-2+2*i-1)=i
           i_pair(2*n-2+2*i)=i+2
     end do
  end if
  
  !************************* ladder *****************************************
  !ladder

  ! do i=1,n/2
  !    if (i.ne.n/2) then        
  !       i_pair(2*i-1)=i
  !       i_pair(2*i)=i+1
  !       i_pair(n+2*i-1)=n/2+i
  !       i_pair(n+2*i)=n/2+i+1
  !       ! bondj(n/2+i   , n/2+i+1)=0.0d0
  !       ! bondj(n/2+i+1 , n/2+i  )=0.0d0
  !    else
  !       i_pair(2*i-1)=i
  !       i_pair(2*i)=1
  !       i_pair(n+2*i-1)=n/2+i
  !       i_pair(n+2*i)=n/2+1
  !       ! bondj(n/2+1 , n/2+i)=0.0d0
  !       ! bondj(n/2+i , n/2+1)=0.0d0
  !    end if
  !    i_pair(2*n+2*i-1)=i
  !    i_pair(2*n+2*i  )=n/2+i
  !    bondj(i     , n/2+i)=-1.0d0
  !    bondj(n/2+i , i    )=-1.0d0
  ! end do

  ! i_pair(n-1)=1
  ! i_pair(n)=n/2
  ! i_pair(2*n-1)=n/2+1
  ! i_pair(2*n)=n

!lacked ladder
  ! do i=1,n/2
  !    if (i.ne.n/2) then        
  !       i_pair(2*i-1)=i
  !       i_pair(2*i)=i+1
  !    else
  !       i_pair(2*i-1)=i
  !       i_pair(2*i)=1
  !    end if
  !    i_pair(n+2*i-1)=i
  !    i_pair(n+2*i  )=n/2+i
  !    bondj(i     , n/2+i)=-1.5E-0
  !    bondj(n/2+i , i    )=-1.5E-0
! end do

  !*************************  spin1 j1-j2  made by half spin  *****************************************
  if(model == 's1-j1j2-half')then ! i cant now

     n_pair=2*n
     allocate(i_pair(n_pair*2),twi_pair(n_pair), bondj(n_pair),delta(n_pair))

     delta=1.0d0
     bondj_tmp=1.0d0
     twi_num=2
     
     do i=1,n
        if (i.ne.n) then
           i_pair(2*i-1)=i
           i_pair(2*i)=i+1
           bondj(i)=bondj_tmp
        else
           i_pair(2*i-1)=1
           i_pair(2*i)=i
           bondj(i)=bondj_tmp
        end if
     end do

     do i=1,n
        if (i.le.n-2) then
           i_pair(2*n+2*i-1)=i
           i_pair(2*n+2*i)=i+2
        else
           i_pair(2*n+2*i-1)=i+2-n
           i_pair(2*n+2*i)=i
        end if
     end do
   
     allocate(twi_unit(0:twi_num))
     twi_unit(0)=0.0d0
     if(del.ne.0) then
        twi_unit(1)=2.0d0*pi/dble(del)
        twi_unit(2)=-2.0d0*pi/dble(del)
     end if

     twi_pair = 0
     if(twist_type == 'j1') then
        twi_pair(n) = 1
     end if
     if(twist_type == 'j2') then
        twi_pair(2*n) = 1
     end if
     if(twist_type == 'cycle') then
        twi_pair(1) = 1
        twi_pair(n) = 2
        twi_pair(2*n) = 1
     end if


     do i=1,n_pair
        if(twi_pair(i) .ne. 0) then
           write(*,*) i,'site is chaged for',twi_pair(i)
           j=i_pair(2*i-1)
           k=i_pair(2*i  )
           i_pair(2*i-1)=max(j,k)
           i_pair(2*i  )=min(j,k)
        end if
     end do
end if
!************************************************************************************************






write(*,*)i_pair
end subroutine pair_format
