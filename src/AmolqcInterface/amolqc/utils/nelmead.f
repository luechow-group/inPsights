c
c standard implementation of the adaptive simplex algorithm
c according to Nelder and Mead.
c

c $Id: nelmead.f,v 1.1.1.1 2007/04/25 13:42:20 luechow Exp $


c     ------------------------------------------------------------------
      subroutine nelmead(x0,n0,maxiter,delta,eps,fmin,func1,iul,logmode)
c     ------------------------------------------------------------------

c dependencies: calls abortp

      implicit none

c parameter
      integer n0                    ! tatsaechliche Dimension
      real*8 x0(n0)                 ! zu optimierender Vektor
      integer maxiter               ! max Anzahl an Iterationen
                                    ! rueckgabe: benoetigte Iter.Zahl
      real*8 delta                  ! anfaengliche Auslenkung
      real*8 eps                    ! Zielgenauigkeit
      real*8 fmin                   ! minimierter Funktionswert
      interface                     ! externe Funktion
        real*8 function func1(x,n)
          integer n
          real*8 x(n)
        end function func1
      end interface
      integer iul                   ! open file unit for log messages
      integer logmode               ! verboseness level

c local constants
      integer nmax
      parameter(nmax=300)
      real*8 alpha,beta,gamma
      parameter(alpha=0.99d0,beta=0.5d0,gamma=2d0)

c local variables
      integer i,j,balt,s,bb,z,iter
      real*8 x(nmax+1,nmax),xrefl(nmax),xquer(nmax),
     &       xneu(nmax)
      real*8 f(nmax+1),falt,frefl,fneu

      if (n0 > nmax) then
         call abortp ('(nelmead): nmax too small in nelmead')
      endif

      do i=1,n0+1
         do j=1,n0
            x(i,j) = x0(j)
         enddo
      enddo
      do j=1,n0
         x(j+1,j) = x(j+1,j) + delta
      enddo

      do i=1,n0+1
        do j=1,n0
          xneu(j) = x(i,j)
        enddo
        f(i) = func1(xneu,n0)
      enddo

      falt = 0d0
      balt = 0

      do iter=1,maxiter
         ! besten b und schlechtesten punkt s feststellen
         s = 1
         bb = 1
         do i=1,n0+1
            if (f(i) > f(s)) s=i
            if (f(i) < f(bb)) bb=i
         enddo
         if (s == 1) then
            z = 2
         else
            z = 1
         endif
         do j=1,n0
            xquer(j) = 0.0
         enddo
         do i=1,n0+1
            if (i /= s) then
               if (f(i) > f(z)) z=i
               do j=1,n0
                  xquer(j) = xquer(j) + x(i,j)/n0
               enddo
            endif
         enddo

         if (bb /= balt) then
            if ( abs(f(bb)-falt) < eps ) goto 105
            falt = f(bb)
            if (logmode.gt.2) then
              write(iul,'(g18.10,i5)') falt,iter
              write(iul,'(5x,3g12.6)') (x(bb,j),j=1,n0)
            endif
            balt = bb
         endif
         do j=1,n0
            xrefl(j) = (1.0+alpha)*xquer(j) - alpha*x(s,j)
         enddo

         frefl = func1(xrefl,n0)

         if (frefl < f(bb)) then
            ! expandiere
            do j=1,n0
               xneu(j) = (1.0-gamma)*xquer(j) + gamma*xrefl(j)
            enddo
            fneu = func1(xneu,n0)
            if (fneu < f(bb)) then
               do j=1,n0
                  x(s,j) = xneu(j)
               enddo
               f(s) = fneu
            else
               do j=1,n0
                  x(s,j) = xrefl(j)
               enddo
               f(s) = frefl
            endif
         else if (frefl < f(z)) then
            ! reflektiere
            do j=1,n0
               x(s,j) = xrefl(j)
            enddo
            f(s) = frefl
         else if (frefl > f(s)) then
            ! (kontrahiere)
            do j=1,n0
               xneu(j) = (1.0-beta)*xquer(j) + beta*x(s,j)
            enddo
            fneu = func1(xneu,n0)
            if (fneu < f(s)) then
               do j=1,n0
                  x(s,j) = xneu(j)
               enddo
               f(s) = fneu
            else
               ! alle neu
               do i=1,n0+1
                  do j=1,n0
                     x(i,j) = ( x(i,j) + x(bb,j) )/2.0
                     xneu(j) = x(i,j)
                  enddo
                  f(i) = func1(xneu,n0)
               enddo
            endif
         else
            ! (f(z) < frefl < f(s): auch kontraktion)
            do j=1,n0
               xneu(j) = (1.0-beta)*xquer(j) + beta*xrefl(j)
            enddo
            fneu = func1(xneu,n0)
            if (fneu < f(s)) then
               do j=1,n0
                  x(s,j) = xneu(j)
               enddo
               f(s) = fneu
            else
               ! alles neu!!
               do i=1,n0+1
                  do j=1,n0
                     x(i,j) = ( x(i,j) + x(bb,j) )/2.0
                     xneu(j) = x(i,j)
                  enddo
                  f(i) = func1(xneu,n0)
               enddo
            endif
         endif

      enddo

 105  continue
      do i=1,n0
         x0(i) = x(bb,i)
      enddo
      fmin = f(bb)

      return
      end

