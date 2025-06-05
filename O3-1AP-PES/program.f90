program O3surfTest
use callpot1
implicit none
real*8, dimension(3) :: r
real*8, dimension(3) :: dvdr
real*8 :: anener


!      3
!      O
!     / \
!  r3/   \r2
!   /     \
!  /       \
! O---------O
! 1   r1    2


r(1)=2.24d0
r(2)=2.24d0
r(3)=2.24d0
call o31appes(r, anener, dvdr)
print*, r(1), r(2), r(3), anener, dvdr

end program
