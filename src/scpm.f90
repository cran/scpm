subroutine h2d(n, s, a, r, h)

!*****************************************************************************80
!
!! h2d subroutine
! compute the distance matrix
! n : number of locations
! s : n*2 array of location (each row is a location)
! a : angle in radians (0-360) for anisotropy
! r : rate of shrinking (0-1) for anisotropy
! h : n*n output matrix with distances
  implicit none

  integer :: n, i, j
  double precision :: s(n,2), h(n,n)
  double precision :: a, r, Hm(2,2), Dm(2,2), Am(2,2)
  double precision :: ds(2), dsp(2)

  Hm(1,1) = cos(a)
  Hm(1,2) = sin(a)
  Hm(2,1) = -sin(a)
  Hm(2,2) = cos(a)
  Dm(1,1) = 1
  Dm(1,2) = 0
  Dm(2,1) = 0
  Dm(2,2) = 1/r
  Am = matmul(Dm,Hm)

  do i = 1, n
   do j = 1, n
	ds(1:2) = s(i,1:2) - s(j,1:2)
	dsp(1:2) = matmul(Am,ds(1:2))
	h(i,j) = sqrt(dsp(1)**2 + dsp(2)**2)
   end do
  end do

end
