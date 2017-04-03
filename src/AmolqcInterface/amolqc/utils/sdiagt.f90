program sdiagt
use utilsmodule
implicit none

integer i,j,k,n,istat,ierr
real*8, allocatable:: A(:,:),ev(:,:),evt(:,:),lambda(:)
real*8 tol,maxerr

n = 3
allocate(A(n,n),ev(n,n),evt(3,3),lambda(n),STAT=istat)
call assert(istat==0," sdiagt: allocation failed ")

A = 0
A(1,1) = -1.d0
A(2,2) = -1.d0
A(3,3) = -1.d0
A(1,3) = -0.1d0
A(3,1) = -0.1d0

call sdiag(A,n,n,lambda,n,ev,n,n,n,ierr)

tol = 1.d-4
call assertEqualAbsolute(lambda(1),-1.1d0,tol,' 1st eigenvalue failed')
call assertEqualAbsolute(lambda(2),-1.0d0,tol,' 2nd eigenvalue failed')
call assertEqualAbsolute(lambda(3),-0.9d0,tol,' 3rd eigenvalue failed')

if (ev(1,1)<0.d0) ev(:,1) = -ev(:,1)
if (ev(2,2)<0.d0) ev(:,2) = -ev(:,2)
if (ev(1,3)<0.d0) ev(:,3) = -ev(:,3)
!print*, ev(:,1)
!print*, ev(:,2)
!print*, ev(:,3)
evt = 0
evt(1,1) = 1.d0/sqrt(2.d0)
evt(3,1) = 1.d0/sqrt(2.d0)
evt(2,2) = 1.d0
evt(1,3) = 1.d0/sqrt(2.d0)
evt(3,3) = -1.d0/sqrt(2.d0)

evt = evt - ev
maxerr = maxval(abs(evt))
!print*, maxerr
call assert(maxerr < tol,' sdiagt: error in eigenvector matrix')
print*, 'all tests passed'

end
   

