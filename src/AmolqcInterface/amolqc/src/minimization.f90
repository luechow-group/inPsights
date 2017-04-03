MODULE MINIMIZATION
IMPLICIT NONE

private
public :: dfpmin,minimization_setGrad, minimization_getGrad, &
  minimization_setIterMax,minimization_getIterMax,minimization_setStepMax, &
  minimization_getStepMax
logical :: numgrad=.false.
logical :: constraint =.false.
logical :: mPrint = .false.
real*8  :: STPMX = 1d0
integer :: ITMAX=200

CONTAINS

   subroutine minimization_setGrad(g)
      logical, intent(in) :: g
      numgrad = g
   end subroutine minimization_setGrad

   logical function minimization_getGrad()
      minimization_getGrad = numgrad
   end function minimization_getGrad

   subroutine minimization_setIterMax(it)
      integer, intent(in) :: it
      if (it>0) ITMAX = it
   end subroutine minimization_setIterMax

   integer function minimization_getIterMax()
      minimization_getIterMax = ITMAX
   end function minimization_getIterMax

   subroutine minimization_setStepMax(st)
      real*8, intent(in) :: st
      if (st>0) STPMX = st
   end subroutine minimization_setStepMax

   real*8 function minimization_getStepMax()
      minimization_getStepMax = STPMX
   end function minimization_getStepMax



SUBROUTINE dfpmin(theta,n,gtol,func,verbose,ierr)
REAL*8, INTENT(IN)    :: gtol
INTEGER, INTENT(IN)   :: n	        ! number of elements in the parameter vector
REAL*8, INTENT(INOUT) :: theta(n)	! parameter vector over which you which to maximize
integer, intent(in)   :: verbose   ! verbose level
integer, intent(out)  :: ierr      ! 0=grad converged, 1=x converged, 2=# of iter. exceeded, 3=f is NaN
INTERFACE	                        ! This lets the subroutine know that func is not a variable, it is a function
	REAL*8 FUNCTION func(theta,n)   ! This is not the likelihood itself, it just let the program know
    IMPLICIT NONE		  
    INTEGER, INTENT(IN) :: n
    REAL*8, INTENT(IN) :: theta(n)
    END FUNCTION func
END INTERFACE
!INTEGER, PARAMETER :: ITMAX=1000	! Maximum number of iterations (in case it does not converge)
REAL*8, PARAMETER :: EPS=EPSILON(theta),TOLD=EPS ! Tolerance criterions
INTEGER :: i,j,its ! indices
LOGICAL :: check 
REAL*8 :: den,fac,fad,fae,ftheta,sumdg,sumd ! working matrices
REAL*8 :: dg(n),g(n),hdg(n),thetanew(n),d(n)
REAL*8 :: B(n,n),fret,stpmax

ierr = 2
ftheta=func(theta,n)	        ! initial value of the function
g=gradient(theta,n,ftheta,func)	! initial value of the gradient, if analytical gradient known here

if (verbose>=2) then
   write(998,*) ' --- dfpmin start ---'
   write(998,*) '   f=',ftheta,' n=',n
   do i=1,n/3
      write(998,'(i5,3f12.5)') i,theta(i),theta(n/3+i),theta(2*n/3+i)
   end do
   do i=1,n/3
      write(998,'(i5,3f12.5)') i,g(i),g(n/3+i),g(2*n/3+i)
   end do
end if

if (isnan(ftheta)) then
   ierr = 3
   return
end if

B=0.0d0				! initialize the matrix to a negative definite matrix
DO j = 1, n
	B(j,j)=1.0d0
END DO
d=-MATMUL(B,g)	! initial line direction
stpmax=STPMX*MAX(vabs(theta),REAL(SIZE(theta),8))
DO its=1,ITMAX	! Main loop over iterations
	CALL lnsrch(theta,ftheta,g,d,thetanew,fret,stpmax,check,func)
	ftheta=fret	! update the value of the fuction
	if (isnan(ftheta)) then
	   ierr = 3
	   exit
	end if
    d=thetanew-theta	! update the direction of movement
    theta=thetanew	! update the point
    if (verbose>=2) then
       write(998,*) ' --- dfpmin iter ',its
       write(998,*) '   f=',ftheta
       do i=1,n/3
          write(998,'(i5,3f12.5)') i,theta(i),theta(n/3+i),theta(2*n/3+i)
       end do
       write(998,*) '   max(abs(d))/max(abs(f),1)=',maxval(abs(d)/max(abs(theta),1.d0))
    end if
    if (MAXVAL(ABS(d)/MAX(ABS(theta),1.0D0)) < TOLD) then ! Check convergence in Delta-d
       ierr = 1
       exit
    end if
    dg=g				! save the old gradient
    g=gradient(theta,n,ftheta,func)	! and get the new gradient
    den=MAX(ftheta,1.0d0)    !!! MUSS DAS NICHT ABS(FTHETA) HEISSEN???
    if (verbose>=2) then
       write(998,*) '   grad:'
       do i=1,n/3
          write(998,'(i5,3f12.5)') i,g(i),g(n/3+i),g(2*n/3+i)
       end do
       write(998,*) ' grad. criterion = ',maxval(abs(g)*max(abs(theta),1.d0)/den)
    end if
    if (MAXVAL(ABS(g)*MAX(ABS(theta),1.0d0)/den) < gtol) then ! Check for convergence on zero gradient
       ierr = 0
       exit
    end if
    dg=g-dg			! compute difference in gradients
    hdg=MATMUL(B,dg)	        ! and difference times current matrix
    fac=DOT_PRODUCT(dg,d)	! Dot products for denominators
    fae=DOT_PRODUCT(dg,hdg)
    sumdg=DOT_PRODUCT(dg,dg)
    sumd=DOT_PRODUCT(d,d)
    IF (fac > DSQRT(EPS*sumdg*sumd)) THEN	! Skip update if fac not positive enough	
		fac=1.0d0/fac
        fad=1.0d0/fae
        dg=fac*d-fad*hdg			! Vector that makes BFGS different from DFP
        ! BFGS updating formula
        B=B+fac*outerprod(d,d)-fad*outerprod(hdg,hdg)+fae*outerprod(dg,dg)
    END IF
    d=-MATMUL(B,g)	! Next direction to go
END DO		! Go back for another iteration
END SUBROUTINE dfpmin

FUNCTION outerprod(a,b) ! Gets the outerproduct of two vectors
IMPLICIT NONE	
REAL*8, DIMENSION(:), INTENT(IN) :: a(:),b(:)
REAL*8 :: outerprod(SIZE(a),SIZE(b))
outerprod = SPREAD(a,dim=2,ncopies=SIZE(b)) * &
SPREAD(b,dim=1,ncopies=size(a))
END FUNCTION outerprod

FUNCTION GRADIENT(x0,k,f0,func)

! Routine to get the numerical gradient of a function using forward differences
IMPLICIT NONE
INTEGER, INTENT(IN) :: k
REAL*8, INTENT(IN) :: x0(k),f0
INTERFACE	                       ! This lets the subroutine know that func is not a variable, it is a function
	REAL*8 FUNCTION func(theta,k)  ! This is not the likelihood itself, it just let the program know
    IMPLICIT NONE		       ! the arguments and type of arguments it accepts
    INTEGER, INTENT(IN) :: k
    REAL*8, INTENT(IN) :: theta(k)
    END FUNCTION func
END INTERFACE

INTERFACE	                       
SUBROUTINE grad_max(vec,grad,k)          
    IMPLICIT NONE		    
    INTEGER, INTENT(IN)  :: k
    REAL*8, INTENT(IN)  :: vec(k)
    REAL*8, INTENT(OUT) :: grad(k)
    END SUBROUTINE grad_max
END INTERFACE

REAL*8 :: GRADIENT(k),grdd(k),dh(k),ax0(k),xdh(k)
REAL*8 :: arg(k,k),dax0(k),one
INTEGER :: i
if(numgrad)then
!write(*,*) "numerical grad in bfgs used"
one  = 1.0D0
grdd = 0.0D0
! Computation of stepsize (dh) for gradient
ax0  = ABS(x0)
dax0 = 1.0D0
WHERE (x0.ne.0.0D0) dax0 = x0/ax0
dh = (1.0D-6)*(1.0D-2)*one*dax0
WHERE (ax0.gt.(1.0D-2)*one) dh = (1.0D-8)*ax0*dax0
xdh = x0+dh
arg = SPREAD(x0,DIM=2,NCOPIES=k)
DO i = 1, k
	arg(i,i)=xdh(i)
    grdd(i)=func(arg(:,i),k)
END DO
GRADIENT = (grdd-f0)/dh
else
 call grad_max(x0,GRADIENT,k)
endif
END FUNCTION GRADIENT
  
SUBROUTINE lnsrch(xold,fold,g,p,x,f,stpmax,check,func)
! Given an N-dimensional point xold, the value of the function and gradient there, fold
! and g, and a direction p, finds a new point x along the direction p from xold where the
! function func has decreased sufficiently. xold, g, p, and x are all arrays of length N.
! The new function value is returned in f. stpmax is an input quantity that limits the length
! of the steps so that you do not try to evaluate the function in regions where it is undefined
! or subject to overflow. p is usually the Newton direction. The output quantity check is
! false on a normal exit. It is true when x is too close to xold. In a minimization algorithm,
! this usually signals convergence and can be ignored. However, in a zero-finding algorithm
! the calling program should check whether the convergence is spurious.
! Parameters: ALF ensures su.cient decrease in function value; TOLX is the convergence
! criterion on .x.
IMPLICIT NONE
REAL*8, DIMENSION(:), INTENT(IN) :: xold,g
REAL*8, DIMENSION(:), INTENT(INOUT) :: p(SIZE(xold))
REAL*8, INTENT(IN) :: fold,stpmax
REAL*8, INTENT(OUT) :: x(SIZE(xold))
REAL*8, INTENT(OUT) :: f
LOGICAL, INTENT(OUT) :: check
INTERFACE	! This lets the subroutine know that func is not a variable, it is a function
	REAL*8 FUNCTION func(theta,k)  ! This is not the likelihood itself, it just let the program know
    IMPLICIT NONE		
    INTEGER, INTENT(IN) :: k
    REAL*8, INTENT(IN) :: theta(k)
    END FUNCTION func
END INTERFACE
REAL*8, PARAMETER :: ALF=1.0e-4,TOLX=epsilon(x)
INTEGER :: ndum
REAL*8 :: a,alam,alam2,alamin,b,disc,f2,pabs,rhs1,rhs2,slope,tmplam,zero
zero=0.0D0
ndum=SIZE(xold)
IF (SIZE(g).NE.size(xold)) THEN
    call abortp('***FATAL ERROR: Dimensions do not agree (LNSRCH)***')
END IF
check=.false.
pabs=vabs(p(:))
IF (pabs > stpmax) p(:)=p(:)*stpmax/pabs ! Scale if attempted step is too big.
slope=DOT_PRODUCT(g,p)
IF (slope >= zero) THEN
    call abortp('***FATAL ERROR: roundoff problem in lnsrch***')
END IF
alamin=TOLX/MAXVAL(ABS(p(:))/MAX(ABS(xold(:)),1.0D0)) !Compute .min.
alam=1D0 ! Always try full Newton step first.
DO ! Start of iteration loop.
	x(:)=xold(:)+alam*p(:)
    F=FUNC(x,size(x))
    IF (alam < alamin) THEN                 ! Convergence on .x. For zero finding,
		!the calling program should verify the convergence.
		x(:)=xold(:)
		check=.true.
		RETURN
    ELSE IF (f <= fold+ALF*alam*slope) THEN ! Sufficient function decrease.
		RETURN
	ELSE !Backtrack.
		IF (alam == 1D0) THEN !First time.
			tmplam=-slope/(2.0D0*(f-fold-slope))
        ELSE !Subsequent backtracks.
			rhs1=f-fold-alam*slope
            rhs2=f2-fold-alam2*slope
            a=(rhs1/alam**2-rhs2/alam2**2)/(alam-alam2)
            b=(-alam2*rhs1/alam**2+alam*rhs2/alam2**2)/&
              (alam-alam2)
            IF (a == 0.0D0) THEN
				tmplam=-slope/(2.0D0*b)
            ELSE
				disc=b*b-3.0*a*slope
                IF (disc < 0.0D0) THEN
					tmplam=0.5D0*alam
                ELSE IF (b <= 0.0D0) THEN
					tmplam=(-b+SQRT(disc))/(3.0D0*a)
                ELSE
					tmplam=-slope/(b+SQRT(disc))
                END IF
             END IF
             IF (tmplam > 0.5D0*alam) tmplam=0.5D0*alam
		END IF
	END IF
    alam2=alam
    f2=f
    alam=MAX(tmplam,0.1D0*alam) 
END DO !Try again.
END SUBROUTINE lnsrch

FUNCTION vabs(v)
! Return the length (ordinary L2 norm) of a vector.
REAL*8, DIMENSION(:), INTENT(IN) :: v
REAL*8 :: vabs
vabs=sqrt(dot_product(v,v))
END FUNCTION vabs

END MODULE MINIMIZATION


