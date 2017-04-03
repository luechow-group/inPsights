      subroutine dnmtr(n,x,f,g,H,ldH,frtol,fatol,fmin,task,
     +                 delta,diag,scale,isave,dsave,xc,gc,s,gs,wa)
      character*60 task
      logical scale
      integer n, ldH
      integer isave(4)
      real*8 f, frtol, fatol, fmin, delta
      real*8 x(n), g(n), H(ldH,n), diag(n)
      real*8 dsave(7)
      real*8 xc(n), gc(n), s(n), gs(n), wa(n)
c     **********
c
c     Subroutine dnmtr
c
c     This subroutine implements a trust region Newton method for the
c     solution of unconstrained optimization problems
c
c           min { f(x) }
c
c     The user must evaluate the function, gradient, and the Hessian matrix.
c
c     This subroutine uses reverse communication.
c     The user must choose an initial approximation x to the minimizer,
c     and make an initial call with task set to 'START'.
c     On exit task indicates the required action.
c
c     A typical invocation has the following outline:
c
c     Compute a starting vector x.
c
c     task = 'START'
c     search = .true.
c
c     do while (search) 
c
c        if (task .eq. 'F' .or. task .eq. 'START') then
c           Evaluate the function at x and store in f.
c        end if
c        if (task .eq. 'GH' .or. task .eq. 'START') then
c           Evaluate the gradient at x and store in g.
c           Evaluate the Hessian at x and store in Hs.
c        end if
c
c        call dnmtr(n,x,f,g,H,ldH,frtol,fatol,fmin,task,
c    +              delta,diag,scale,isave,dsave,
c    +              s,xc,gs,gc,wa)
c
c        if (task(1:4) .eq. 'CONV') search = .false.
c
c      end do
c
c     NOTE: The user must not alter work arrays between calls.
c
c     The subroutine statement is
c
c       subroutine dnmtr(n,x,f,g,H,ldH,frtol,fatol,fmin,task,
c                        delta,diag,scale,isave,dsave,
c                        xc,gc,s,gs,wa)
c     where
c
c       n is an integer variable.
c         On entry n is the number of variables.
c         On exit n is unchanged.
c
c       x is a double precision array of dimension n.
c          On entry x is an approximation to the solution.
c          On exit x is the current approximation.
c
c       f is a double precision variable.
c         On entry f is the value of the function at x.
c         On final exit f is the value of the function at x.
c
c       g is a double precision array of dimension n.
c         On entry g is the value of the gradient at x.
c         On final exit g is the value of the gradient at x.
c
c       H is a double precision array of dimension (ldH,n).
c         On entry H is the the Hessian matrix at x.
c         On exit H is the Hessian matrix at the previous value
c            of x.
c
c       ldH is an integer variable.
c         On entry ldH is the leading dimension of array H.
c         On exit ldH is unchanged.
c
c       frtol is a double precision variable.
c         On entry frtol specifies the relative error desired in the
c            function. Convergence occurs if the estimate of the
c            relative error between f(x) and f(xsol), where xsol
c            is a local minimizer, is less than frtol.
c         On exit frtol is unchanged.
c
c       fatol is a double precision variable.
c         On entry fatol specifies the absolute error desired in the
c            function. Convergence occurs if the estimate of the
c            absolute error between f(x) and f(xsol), where xsol
c            is a local minimizer, is less than fatol.
c         On exit fatol is unchanged.
c
c       fmin is a double precision variable.
c         On entry fmin specifies a lower bound for the function.
c            The subroutine exits with a warning if f < fmin.
c         On exit fmin is unchanged.
c
c       task is a character variable of length at least 60.
c         On initial entry task must be set to 'START'.
c         On exit task indicates the required action:
c
c            If task(1:1) = 'F' then evaluate the function at x
c            and call dnmtr again.
c
c            If task(1:2) = 'GH' then evaluate the gradient and the
c            Hessian matrix at x and call dnmtr again.
c
c            If task(1:4) = 'NEWX' then a new iterate has been
c            computed. The approximation x, function f, gradient g,
c            and Hessian matrix are available for examination.
c
c            If task(1:4) = 'CONV' then the search is successful.
c
c            If task(1:4) = 'WARN' then the subroutine is not able
c            to satisfy the convergence conditions. The exit value
c            of x contains the best approximation found.
c
c            If task(1:5) = 'ERROR' then there is an error in the
c            input arguments.
c
c         On exit with convergence, a warning or an error, the
c            variable task contains additional information.
c
c       delta is a double precision variable.
c         On entry delta specifies the initial trust region radius.
c            If delta <= 0, delta is initialized internally.
c         On exit delta is the current trust region radius.
c
c       diag is a double precision array of dimension n.
c         On entry diag need not be specified.
c         On exit diag specifies the scaling of the variables.
c            See scale below.
c
c       scale is a logical variable.
c         On entry scale specifies the scaling of the variables.
c            If scale is true, the variables are scaled internally.
c            If scale is false, diag is set to the identity and
c            the variables are not scaled.
c         On exit scale is unchanged.
c
c       isave is an integer work array of dimension 4.
c
c       dsave is a double precision work array of dimension 7.
c
c       xc is a double precision work array of dimension n.
c
c       gc is a double precision work array of dimension n.
c
c       s is a double precision work array of dimension n.
c
c       gs is a double precision work array of dimension n.
c
c       wa is a double precision work array of dimension n.
c
c     Subprograms called
c
c       MINPACK-2  ......  dgqt
c
c       Level 1 BLAS  ...  dcopy, ddot, dnrm2
c
c     MINPACK-2 Project. October 1993.
c     Argonne National Laboratory and University of Minnesota.
c     Brett M. Averick and Jorge J. More'.
c
c     March 2000. Coding improved and simplified.
c     Jorge J. More'
c
c     **********
      double precision zero, p5, one
      parameter (zero=0.0d0,p5=0.5d0,one=1.0d0)

c     Parameters for updating the iterates.

      double precision eta0, eta1, eta2
      parameter(eta0=1d-4,eta1=0.25d0,eta2=0.75d0)

c     Parameters for updating the trust region size delta.

      double precision sigma1, sigma2, sigma3
      parameter(sigma1=0.25d0,sigma2=0.5d0,sigma3=2.0d0)

c     Parameters that govern the trust region subproblem.

      double precision deltaf, stol
      parameter (deltaf=0.1d0,stol=0.1d0)

      logical search
      character*30 work
      integer i, info, iter, j
      double precision actred, fc, par, prered, snorm
      double precision alpha, g0, gnorm

      double precision ddot, dnrm2
      external dcopy, ddot, dgqt, dnrm2

c     Initialialization section.

      if (task(1:5) .eq. 'START') then

c        Check the input arguments for errors.

         if (n .le. 0) task = 'ERROR: N .LE. 0'
         if (ldH .lt. n) task = 'ERROR: LDHES .LT. N'
         if (frtol .le. zero) task = 'ERROR: FRTOL .LE. 0'
         if (fatol .le. zero) task = 'ERROR: FATOL .LE. 0'

c        Exit if there are errors, or initialize local variables.

         if (task(1:5) .eq. 'ERROR') then

            work = 'EXIT'

         else

            work = 'START SEARCH'
            iter = 1
            par = zero
            
c           Compute the initial scaling matrix and
c           scale the initial gradient.
            
            call dcopy(n,g,1,gs,1)
            if (scale) then
               do j = 1, n
                  diag(j) = sqrt(dnrm2(n,H(1,j),1))
                  gs(j) = gs(j)/diag(j)
               end do
            end if
            gnorm = dnrm2(n,gs,1)
            
c           Initialize the step bound delta.
            
            if (delta .le. zero) delta = deltaf*gnorm
            
         end if

      else

c        Restore local variables.

         if (isave(1) .eq. 1) then
            work = 'START SEARCH'
         else if (isave(1) .eq. 2) then
            work = 'COMPUTE'
         else if (isave(1) .eq. 3) then
            work = 'EVALUATE'
         else if (isave(1) .eq. 4) then
            work = 'SCALE'
         end if
         iter = isave(2)

         par = dsave(1)
         actred = dsave(2)
         fc = dsave(3)
         prered = dsave(5)
         snorm = dsave(6)
         prered = dsave(7)

      end if

      if (work .eq. 'START SEARCH') then

c        Scale the gradient and Hessian matrix.

         call dcopy(n,g,1,gs,1)
         if (scale) then
            do j = 1, n
               gs(j) = gs(j)/diag(j)
               do i = 1, n
                  H(i,j) = (H(i,j)/diag(i))/diag(j)
               end do
            end do
         end if

c        Set work to compute the step.

         work = 'COMPUTE'

      end if

c     Search for a lower function value.

      search = .true.
      do while (search)

c        Compute a step and evaluate the function at the trial point.

         if (work .eq. 'COMPUTE') then
         
c           Compute the trust region step.
         
            call dgqt(n,H,ldH,gs,delta,stol,zero,20,par,prered,s,
     +                info,xc,gc,wa)

c           Compute the predicted reduction and check for errors.

            prered = -prered
            if (info .eq. 3) 
     +         task = 'WARNING: ROUNDING ERRORS'

c           Save the best function value, iterate, and gradient.
         
            fc = f
            call dcopy(n,x,1,xc,1)
            call dcopy(n,g,1,gc,1)

c           Compute the new iterate.
         
            snorm = dnrm2(n,s,1)
            if (scale) then
               do i = 1, n
                  s(i) = s(i)/diag(i)
               end do
            end if
            do i = 1, n
               x(i) = x(i) + s(i)
            end do
         
c           Set task to compute the function.
         
            task = 'F'
         
         end if

c        Evaluate the step and determine if the step is successful.

         if (work .eq. 'EVALUATE') then
         
c           Compute the actual reduction.
         
            actred = fc - f
         
c           On the first iteration, adjust the initial step bound.
      
            if (iter .eq. 1)  delta = min(delta,snorm)
         
c           Update the trust region bound.
         
            g0 = ddot(n,gc,1,s,1)
            if (f-fc-g0 .le. zero) then
               alpha = sigma3
            else
               alpha = -p5*(g0/(f-fc-g0))
            end if

c           The trust region bound depends on the ratio of the 
c           actual to predicted reduction.
         
            if (actred .lt. eta0*prered) then
                delta =  min(max(alpha,sigma1)*snorm,sigma2*delta)
             else if (actred .lt. eta1*prered) then
                delta = max(sigma1*delta,min(alpha*snorm,sigma2*delta))
            else if (actred .lt. eta2*prered) then
                delta = max(sigma2*delta,min(alpha*snorm,delta))
            else
                delta = max(delta,min(alpha*snorm,sigma3*delta))
            end if

c           Update the iterate.
         
            if (actred .ge. eta0*prered) then
               task = 'GH'
            else
               task = 'F'
               f = fc
               call dcopy(n,xc,1,x,1)
               call dcopy(n,gc,1,g,1)
            end if
         
c           Test for convergence.
         
            if (f .lt. fmin) task = 'WARNING: F .LT. FMIN'
            if (abs(actred) .le. fatol .and. prered .le. fatol) task =
     +          'CONVERGENCE: FATOL TEST SATISFIED'
            if (abs(actred) .le. frtol*abs(f) .and.
     +          prered .le. frtol*abs(f)) task =
     +          'CONVERGENCE: FRTOL TEST SATISFIED'
         
         end if

c        Test for continuation of search

        if (task .eq. 'F' .and. work .eq. 'EVALUATE') then
           search = .true.
           work = 'COMPUTE'
        else
           search = .false.
        end if

      end do

c     Update the scaling matrix.

      if (work .eq. 'SCALE') then
         iter = iter + 1
         task = 'NEWX'
         if (scale) then
            do j = 1, n
               diag(j) = max(diag(j),sqrt(dnrm2(n,H(1,j),1)))
            end do
         end if
      end if

c     Decide on what work to perform on the enxt iteration.

      if (task .eq. 'F' .and. work .eq. 'COMPUTE') then
         work = 'EVALUATE'
      else if (task .eq. 'F' .and. work .eq. 'EVALUATE') then
         work = 'COMPUTE'
      else if (task .eq. 'GH') then
         work = 'SCALE'
      else if (task .eq. 'NEWX') then
         work = 'START SEARCH'
      end if

c     Save local variables.

      if (work .eq. 'START SEARCH') then
         isave(1) = 1
      else if (work .eq. 'COMPUTE') then
         isave(1) = 2
      else if (work .eq. 'EVALUATE') then
         isave(1) = 3
      else if (work .eq. 'SCALE') then
         isave(1) = 4
      end if
      isave(2) = iter

      dsave(1) = par
      dsave(2) = actred
      dsave(3) = fc
      dsave(5) = prered
      dsave(6) = snorm
      dsave(7) = prered

      end
