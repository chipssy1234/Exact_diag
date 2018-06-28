function total_sz(eigen_vec,site)
  use variables
  implicit none
  integer *8 i
  integer site
  complex(kind(0d0)) eigen_vec(szdim)
  double precision total_sz,tmp_i, tmp_n


  total_sz=0.0d0

  do i=1,szdim
     tmp_i=dble(eigen_vec(i))*dble(eigen_vec(i))!+aimag(eigen_vec(i))*aimag(eigen_vec(i))     
     tmp_n=list(i,site)
     total_sz=total_sz+tmp_i*tmp_n
  end do

end function total_sz
