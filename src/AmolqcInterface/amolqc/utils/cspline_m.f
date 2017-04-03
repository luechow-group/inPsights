!
! F90 module cubicspline. contains array of 1D cubic splines over [0:infinity[ 
!
! $Id: cspline_m.f,v 1.2 2008/02/07 17:00:53 luechow Exp $
!
! contains:
!    csplinit(...)              :  allocates arrays, set max # of spline points 
!    csplfinal()                :  deallocates arrays
!    csplxarray(np)             :  constructs array of spline knots
!                               :  (corresponds to spline evaluation on
!                               :   csplint)
!    cspline(ispl,...)          :  construct 'ispl'-th cubic spline, 
!                                  using 'n' (x,y) pairs with x points 
!                                  stored in array 'csplx(n)'
!    csplint(ispl,x)            :  evaluate 'ispl'-th spline at x
!
!
!     a cubic spline object is defined by number of points (knots) 
!     and the arrays for y and y".
!
!     a static array of cubic spline objects is defined in the header file,
!     each objects is accessed by an index number (ispl).
!
!     subroutines adapted from Numerical Recipes. This version assumes
!     equispaced values for x in [0,1[ mapped to [0,inf[ 
!     with x -> alpha*x/(1-x). 
!     More precisely:
!      h = 1.d0/(n-1)                         
!      do i=1,n-1
!         x(i) = alpha*(i-1)*h / (1.d0 - (i-1)*h))
!      enddo
!      x(n) = 1.d300                     ! better system parameter 'HUGE'
!
!     The alpha parameter (csalpha) scales the mapping.
!
!     The smaller alpha the closer are the spline knots for small x
!     (and the farther apart for large x) which is important for the 
!     exponential decay.
!
! this version saves the actual coefficients a_k,b_k,c_k,d_k 
! for the cubic polynomial to get faster spline interpolation.
!
!     09.09.1999 SM
!     Cusp-Korrektur der 1S-Funktionen ueber STO's wurde eingefuegt. 

c=========================================================

      MODULE cubicspline
      
      implicit none
!      PRIVATE
!      PUBLIC:: csnpmax,csnsplmax,csalpha,csplx,csplnpnt,csplint,
!     &         csplxarray,cspline,csplinit,csplfinal

      integer :: csnpmax=0                      ! max. # of spline points (knots)
      integer :: csnsplmax=0                    ! max. # of spline objects
      real*8  :: csalpha=1                      ! alpha param of mapping function
      integer :: csplnpnt=0                     ! # of spline points

      real*8, allocatable :: csplx(:)      ! array of x-values for all splines
      real*8, allocatable :: cspla(:,:),     ! array of A coefficients
     .                       csplb(:,:),     ! array of B coefficients
     .                       csplc(:,:),     ! array of C coefficients
     .                       cspld(:,:)      ! array of D coefficients


      CONTAINS

!===========================================================

!     ----------------------------------
      SUBROUTINE csplinit(nspl,np,alpha)
!     ----------------------------------

! csplinit allocates the spline arrays.

      implicit none
      
      integer nspl   ! (max) # spline functions
      integer np     ! (max) # spline points
      real*8 alpha   ! alpha parameter for mapping [0:1] -> [0:infinity[
      integer ierr

      csnsplmax = nspl
      csnpmax = np
      csalpha = alpha

      allocate(cspla(csnsplmax,csnpmax),csplb(csnsplmax,csnpmax),
     .         csplc(csnsplmax,csnpmax),cspld(csnsplmax,csnpmax),
     .         csplx(csnpmax),stat=ierr)
      if (ierr /= 0) then
         if (allocated(cspla)) then
            call abortp("(csplinit): cspl arrays already allocated")
         else
            call abortp("(csplinit):allocation failed")
         endif
      endif

      END SUBROUTINE csplinit

!===========================================================

!     ----------------------
      SUBROUTINE csplfinal()
!     ----------------------

      implicit none
      
      integer ierr

      deallocate(cspla,csplb,csplc,cspld,csplx,stat=ierr)
      if (ierr /= 0) then
         call abortp("csplfinal:deallocation failed")
      endif

      END SUBROUTINE csplfinal

c===========================================================

c     -------------------------
      SUBROUTINE csplxarray(np)
c     -------------------------

c set the array 'csplx(j)' of np spline points 
c for all functions (mapping of [0,1[ to [0,infty[ )

      implicit none

c input parameter:
      INTEGER np              ! # of spline knots
c variables:
      INTEGER i
      REAL*8 h

      h = 1.d0/(np-1)
      do i=1,np-1
         csplx(i) = csalpha*(i-1)*h/(1.d0-(i-1)*h)
      enddo
      csplx(np) = 1.d300          ! "infinity"
      csplnpnt = np

      END SUBROUTINE csplxarray

c=========================================================

c     -------------------------------------
      SUBROUTINE cspline(ispl,y,n,yp1,y2p1)
c     -------------------------------------

c requires: # of spline points 'csplnpnt' and the array csplx(j), 
c i.e. the points x_j, j=1,..,csplnpnt to be set.

c modified 'spline' routine from Numerical Recipes. 

      implicit none

c input arguments:
c   ispl: index of spline object
c      n: # of knots
c    yp1: y'(x_1)
c   y2p1: y''(x_1)
c    ypn: y'(x_n)
c    y: y(x_i)          function values at knots

c input parameter:
      integer ispl      ! index of spline object
      integer n         ! size of array y
      real*8 y(n)       ! function values at spline points
      real*8 yp1,y2p1   ! y'(x_1) and y''(x'1)
      real*8 ypn        ! y'(x_n)
c variables:
      integer i,k,ierr
      real*8 p,qn,sig,un
      real*8 deltax,deltay,deltay2
      real*8 y2(n),u(n)
      
      if (ispl < 1 .or. ispl > csnsplmax) 
     .     call abortp("(cspline): illegal spline index")
      
      if (n /= csplnpnt) then
         call abortp('(cplsine): inconsistent use of cspline')
      endif

      ypn = 0d0
      if (abs(yp1).gt.1d300) then
        y2(1)=0.d0
        u(1) =0.d0
      else
        y2(1)= y2p1
        u(1) =(3.d0/(csplx(2)-csplx(1)))*((y(2)-y(1))
     .      /(csplx(2)-csplx(1))-yp1)
      endif
      do 11 i=2,n-1
        sig=(csplx(i)-csplx(i-1))/(csplx(i+1)-csplx(i-1))
        p=sig*y2(i-1)+2.d0
        y2(i)=(sig-1.d0)/p
        u(i)=(6.d0*((y(i+1)-y(i))/(csplx(i+1)-csplx(i))-(y(i)-y(i-1))
     .      /(csplx(i)-csplx(i-1)))/(csplx(i+1)-csplx(i-1))
     .      -sig*u(i-1))/p
11    continue
      if (abs(ypn).gt.1d300) then
        qn=0.d0
        un=0.d0
      else
        qn=0.5d0
        un=(3.d0/(csplx(n)-csplx(n-1)))*(ypn-(y(n)-y(n-1))
     .    /(csplx(n)-csplx(n-1)))
      endif
      y2(n)=(un-qn*u(n-1))/(qn*y2(n-1)+1.d0)
      do 12 k=n-1,1,-1
        y2(k)=y2(k)*y2(k+1)+u(k)
12    continue

c construct coefficients of cubic spline polynomial and save as 
c spline object 'ispl' 

      do 13 k=1,n-1
         cspla(ispl,k)  = y(k)

         deltax = csplx(k+1) - csplx(k)
         deltay = y(k+1) - y(k)
         deltay2 = y2(k+1) - y2(k)
         csplb(ispl,k)  = deltay/deltax 
     .                  - deltax*(deltay2/6.d0 + y2(k)/2.d0)
         csplc(ispl,k)  = y2(k)/2.d0
         cspld(ispl,k)  = deltay2 / (6.d0*deltax)
13    continue

      END SUBROUTINE cspline

c==========================================

c     -------------------------------
      REAL*8 FUNCTION csplint(ispl,x)
c     -------------------------------

c requires: previous creation of spline points csplx() 
c and # of spline points csplnpt
c requires: previous call to cspline to set 
c coeff-arrays cspla,csplb,csplc,cspld
          
c cspline interpolation: assumes the knots at 
c     x_i = alpha*(i-1)*h / (1 - (i-1)*h), i=1,..,n-1, x_n = HUGE (1.d300)
c (calculated with 'csplxarray(np)'
c inversion gives x in [x_k,x_k+1] with k = int((n-1)*x/(x+alpha))
c

      implicit none

      INTEGER ispl
      REAL*8 x
      INTEGER j
      REAL*8 dx
      
      j = (csplnpnt-1)*x/(csalpha+x)  + 1

      dx = x - csplx(j)
      csplint = cspla(ispl,j)+dx*(csplb(ispl,j)+dx*(csplc(ispl,j)
     .          +dx*cspld(ispl,j)))
      END FUNCTION csplint

      END MODULE cubicspline
