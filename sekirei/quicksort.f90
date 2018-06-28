
recursive subroutine quicksort(a, first, last)
  implicit none
  integer a(*)
  integer first, last
  integer i, j
  integer x, t

  x = a( (first+last) / 2 )
  i = first
  j = last
  do
     do while (a(i) < x)
        i=i+1
     end do
     do while (x < a(j))
        j=j-1
     end do
     if (i >= j) exit
     t = a(i);  a(i) = a(j);  a(j) = t
     i=i+1
     j=j-1
  end do
  if (first < i - 1) call quicksort(a, first, i - 1)
  if (j + 1 < last)  call quicksort(a, j + 1, last)
  end subroutine quicksort
