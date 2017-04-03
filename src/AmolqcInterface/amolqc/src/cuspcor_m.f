c
c module cuspopt handles the optimization of the cusp correction
c for GTOs.

c $Id: cuspcor_m.f,v 1.2 2008/02/07 16:51:25 luechow Exp $


c     09.09.1999	SM
c     Neues Modul zur Optimierung der Cusp-Korrekturfunktion (STO).
c     Optimiert wird die Abweichung der Funktion selbst (Varianz) von den
c     urspruenglichen Splinepunkten fuer mittelnahe Abstaende.
c     (alternativ auch Optimierung ueber Laplacian)
c     Zur Optimierung wird wird das Nelder-Mead-Verfahren benutzt.


c     --------------
      MODULE cuspopt
c     --------------

      use global
      use cubicspline
      use utilsmodule
      implicit none

      integer point1,point2             ! erster + letzter Punkt der in
                                        ! die Varianzberechnung eingeht
      real*8, allocatable :: cuspopt_xsav(:)    ! x-Koordinaten der Punkte
      real*8, allocatable :: cuspopt_ysav(:)    ! Funktionswerte
      real*8 nuccharge                  ! ??

      ! transition region is [x-transgap,x2], length=translen
      real*8, parameter :: transgap = 0.0015d0
      real*8, parameter :: translen = 0.0025d0

      CONTAINS

c     ---------------------------
      subroutine cuspopt_init(np)
c     ---------------------------

      integer np
      integer ierr

      allocate(cuspopt_xsav(np),cuspopt_ysav(np),stat=ierr)
      if (ierr /= 0) then
         call abortp('(csupopt_init): allocation failed')
      endif

      end subroutine cuspopt_init

c     --------------------------
      subroutine cuspopt_final()
c     --------------------------

      deallocate(cuspopt_xsav,cuspopt_ysav)

      end subroutine cuspopt_final

c========================================================

c     -----------------------------------------------------------
      SUBROUTINE cuspcorrect(ngto,cntrctn,bf,z,a,c,p1,p2,y,yfd,ysd,
     .                       np,k0)
c     -----------------------------------------------------------

c cuspcorrect attempts to correct the electron nuclear cusp.

c input parameter:
      integer, intent(in) :: ngto(:)
      real*8, intent(in)  :: cntrctn(:,:,:)
      integer np       ! actual # of points
      integer bf       ! idx of basis function
      real*8 z         ! desired STO orbital exponent at r=0
                       ! (== nuclear charge to satisfy cusp condition)
      real*8 a,c,b     ! params of function for small r (cusp correction)
                       ! note: a is also used as to flag:
                       ! a < 0.2 means: calculate optimal a and c.
                       ! 0 < a < 0.2: optimize 1S
                       ! a < 0 : optimize 2S
      real*8 p1,p2     ! rmin rmax for transition interval
c input/output parameter
      real*8 y(np)     ! in: original; out: modified values y_i = f(x_i)
      real*8 yfd(np)   ! in: original; out: modified values y_i = f'(x_i)
      real*8 ysd(np)   ! in: original; out: modified values y_i = f''(x_i)
c output parameter:
      integer k0       ! cusp correction from x_1 to x_k0


      real*8 ap,bp,cp,dp          ! coefficients of transition polynomial
      integer point1sav,point2sav,spoint
      integer k,l,j,ic
      real*8 y0,dev,x0(10),fmin
      real*8 yold(np),yfdold(np),ysdold(np),dx
      real*8 varloc,varlocmin,tmp,tmp1,tmp2
      real*8 x11,x12,x13,y1,y1fd,x21,x22,x23,y2,y2fd
      real*8 k1,k2,k3,k4,k5
!!!      external func1a
      interface
        real*8 function func1a(x,n)
          integer :: n
          real*8  :: x(n)
        end function func1a
      end interface

      save point1sav,point2sav

      if (np /= csplnpnt) then
         call abortp('(cuspcorrect): expecting np == csplnpnt')
      endif

      if (a < 0.2) then
         ! calculate optimal cusp correction function

         if (mytid==0 .and. logmode>1) then
            write(iul,*) ' constructing optimal cusp correction for bf=',bf
         endif

c        ! find first point for optimization of cusp correction function:
c        ! this is point where f'(x)/f(x) = -z (as holds for exp(-z*r))
         if (a > 0) then
            do k=1,np-1
               y0 = y(k)
               y1 = (y(k+1) - y(k)) / (csplx(k+1) - csplx(k))
               if (y1/y0 < -z) goto 10
            enddo
 10         continue
            if (k < np-1) then
               k0  = k
            else
               k0  = 0
            endif
         endif

         call cuspopt_init(csplnpnt)
         ! store x,y pairs for optimization of cusp function
         do k=1,csnpmax
            cuspopt_xsav(k)  = csplx(k)
            cuspopt_ysav(k)  = y(k)
         enddo


         ! first and last point in optimization
         if (a > 0d0) then
            point1 = k0
            point2 = (csplnpnt-1)*p2/(csalpha+p2)+1
            point1sav = point1  ! speichere Punkte fuer 2S-Funktion
            point2sav = point2
         else
            point1 = point1sav  ! lade Punkte der 1S-Funktion
            point2 = point2sav
         endif

         if (mytid==0 .and. logmode>1) then
            write(iul,'(1X,A,I5,1X,I5)')
     .           'first, last point in optimization:', point1,point2
         endif

         ! optimization of correction function
         if (a > 0d0) then
            x0(1) = 1d0*z-0.2
            x0(2) = z*sqrt(z/pi)
            x0(3) = 0d0
            nuccharge  = 1d0*z-0.1
         else
            x0(1) = 0.8d0*z
            x0(2) = -0.4d0*z*sqrt(z/pi)
            x0(3) = 0d0
            nuccharge  = 100d0
         endif

         call nelmead(x0,2,1000,0.8d0,1d-6,fmin,func1a,iul,logmode)
         x0(3) = 1d-3
         call nelmead(x0,3,1000,0.8d0,1d-6,fmin,func1a,iul,logmode)
         a     = x0(1)
         c     = x0(2)
         b     = x0(3)

         call cuspopt_final()    ! free space

         ! find optimal range for transition of 2nd derivative of function
         ! to 2nd derivative of STO
         spoint = (csplnpnt-1)*p1/(csalpha+p1)+1
         varlocmin = 1d100
         do l= 21,spoint
            varloc = 0d0
            do k=l-20,l+20
               tmp    = c*exp(-a*csplx(k))*(a**2)
               varloc = varloc + (tmp - ysd(k))**2
            enddo
            if(varloc <= varlocmin) then
               p1      = (csplx(l+2)+csplx(l-2))/2d0
               varlocmin = varloc
               if (logmode >= 4) write(iul,'(A5,I5,2(1X,F16.9))')
     .              'new',l,csplx(l),varlocmin
            endif
         enddo

         if (mytid==0 .and. logmode>1) then
            write(iul,'(1X,A)') 'parameters a,c,b,p1:'
            write(iul,'(4(G12.6,2X))') x0(1),x0(2),x0(3),p1
            write(iul,'(1X,A10,2X,1G12.6)') ' fmin :',fmin
         endif

         ! first and last point of transition region
         p1     = p1 - transgap
         p2     = p1 + translen
         point1 = k0
         point2 = (csplnpnt-1)*p2/(csalpha+p2)

         ! note: b is discarded!
         if (mytid==0 .and. logmode>1) then
            write(iul,'(1X,A)') 'parameters a,c,p1,p2:'
            write(iul,'(4(G12.6,2X))') a,c,p1,p2
         endif

      else

         if (mytid==0 .and. logmode==2) then
            write(iul,'(i5)',advance='no') bf
         endif

         if (mytid==0 .and. logmode>2) then
            write(iul,*) ' correcting cusp of basis function bf=',bf
            write(iul,'(1X,A)') 'parameters a,c,p1,p2:'
            write(iul,'(4(G12.6,2X))') a,c,p1,p2
         endif

         point1 = (csplnpnt-1)*p1/(csalpha+p1)
         point2 = (csplnpnt-1)*p2/(csalpha+p2)

         if (mytid == 0 .and. logmode>2) then
            write(iul,'(1X,A,1X,I5,1X,I5)')
     .           ' first, last point of transition interval :',
     .           point1,point2
         endif
      endif

      ! construction of transition polynomial for 2nd derivative
      ! cubic polynomial is used (with differentiable transitions)
      x11 = csplx(point1)
      x12 = csplx(point1)**2
      x13 = csplx(point1)**3
      y1   = c*exp(-a*csplx(point1))*a**2
      y1fd = c*exp(-a*csplx(point1))*(-a**3)
      x21 = csplx(point2)
      x22 = csplx(point2)**2
      x23 = csplx(point2)**3
      y2   = ysd(point2)
      y2fd = 0d0
      do ic=1,ngto(bf)
        tmp1 = cntrctn(1,ic,bf)
        tmp2 = cntrctn(2,ic,bf) * exp(-tmp1*x22)
        y2fd = y2fd + tmp2*(12d0*tmp1**2*x21-8d0*tmp1**3*x23)
      enddo
      k1 = (x13-x23)/(x11-x21) - 3d0*x22
      k2 = (x12-x22)/(x11-x21) - 2d0*x21
      k3 = (y1-y2)/(x11-x21) - y2fd
      k4 = 3d0*(x12-x22)/(2d0*(x11-x21))
      k5 = (y1fd-y2fd)/(2d0*(x11-x21))

      ! coefficients of cubic polynomial
      ap = (k5-k3/k2)/(k4-k1/k2)
      bp = k5 -k4*ap
      cp = y1fd-3d0*x12*ap-2d0*bp*x11
      dp = y1 -x13*ap-x12*bp-x11*cp

      if (mytid==0 .and. logmode>2) then
        write(iul,'(A,4(1X,D16.10))')
     . ' coeffs of transition polynomial:',ap,bp,cp,dp
      endif

      ! correction of spline points
      do k=csplnpnt-1,1,-1        ! wegen Uebergangsueberpruefung
         yold(k)   = y(k)
         yfdold(k) = yfd(k)
         ysdold(k) = ysd(k)
         if (k < point1) then
            ! STO interval
            ysd(k) = c*exp(-a*csplx(k))*a**2
            tmp1   = c*exp(-a*csplx(point1))*(-a)
     .           - ap/4d0*csplx(point1)**4-bp/3d0*csplx(point1)**3
     .           - cp/2d0*csplx(point1)**2-dp*csplx(point1)
     .           + ap/4d0*csplx(point2)**4+bp/3d0*csplx(point2)**3
     .           + cp/2d0*csplx(point2)**2+dp*csplx(point2) - yfd(point2)
            tmp2   = ap/4d0*csplx(point2)**4+bp/3d0*csplx(point2)**3
     .           + cp/2d0*csplx(point2)**2+dp*csplx(point2) - yfd(point2)

            yfd(k) = c*exp(-a*csplx(k))*(-a) - tmp1
            y(k)   = c*exp(-a*csplx(k))      - tmp1*csplx(k)
     .           - c*exp(-a*csplx(point1)) + tmp1*csplx(point1)
     .           + ap/20d0*csplx(point1)**5+bp/12d0*csplx(point1)**4
     .           + cp/6d0*csplx(point1)**3+dp/2d0*csplx(point1)**2
     .           - tmp2*csplx(point1)
     .           - ap/20d0*csplx(point2)**5-bp/12d0*csplx(point2)**4
     .           - cp/6d0*csplx(point2)**3-dp/2d0*csplx(point2)**2
     .           + tmp2*csplx(point2) + y(point2)
         else if (k < point2) then
            ! transition polynomial
            ysd(k) = ap*csplx(k)**3+bp*csplx(k)**2+cp*csplx(k)+dp
            tmp2   = ap/4d0*csplx(point2)**4+bp/3d0*csplx(point2)**3
     .           + cp/2d0*csplx(point2)**2+dp*csplx(point2) - yfd(point2)
            yfd(k) = ap/4d0*csplx(k)**4+bp/3d0*csplx(k)**3
     .           + cp/2d0*csplx(k)**2+dp*csplx(k) - tmp2
            y(k)   = ap/20d0*csplx(k)**5+bp/12d0*csplx(k)**4
     .           + cp/6d0*csplx(k)**3+dp/2d0*csplx(k)**2 - tmp2*csplx(k)
     .           - ap/20d0*csplx(point2)**5-bp/12d0*csplx(point2)**4
     .           - cp/6d0*csplx(point2)**3-dp/2d0*csplx(point2)**2
     .           + tmp2*csplx(point2) + y(point2)
         endif
         if (logmode > 5) write(iul,'(I5,2X,7(D16.8,2X))')
     .           k,csplx(k),yold(k),y(k),yfdold(k),yfd(k),ysdold(k),ysd(k)
      enddo

c     Ableitungen der Cubic-Splines zur Kontrolle
c
c      tmp1 = -a*c
c      tmp2 =  a*a*c
c
c      call cspline(3*bf-2,y,tmp1,tmp2)
c
c      do l= csplnpnt-1,1,-1
c        yold(l)   = y(l)
c        yfdold(l) = yfd(l)
c        ysdold(l) = ysd(l)
c        j = (csplnpnt-1)*csplx(l)/(csalpha+csplx(l))  + 1
c        k = 3*bf-2
c        dx = csplx(l) - csplx(j)
c        yfd(l)=  csplb(k,j)+2*dx*csplc(k,j) + 3*dx*dx*cspld(k,j)
c        ysd(l)=             2   *csplc(k,j) + 6*dx   *cspld(k,j)
c        if (logmode .ge. 3) write(iul,'(I5,2X,7(D16.8,2X))')
c     .          l,csplx(l),yold(l),y(l),yfdold(l),yfd(l),ysdold(l),ysd(l)
c      enddo
c
      if (mytid == 0) then
         call flush(iul)
      endif

      END SUBROUTINE cuspcorrect

      END MODULE cuspopt

c===================================================================

c     -----------------------------
      real*8 function func1a(x0,n0)
c     -----------------------------

c func1a calculates the sum of squared deviations for
c an exponential function (parameters x0) from the stored
c data points

      use cuspopt
      implicit none

      real*8 x0(10),a,b,c,func
      integer n0,k

      a = x0(1)
      c = x0(2)
      if (n0.eq.3) then
        b = x0(3)
      else
        b = 0d0
      endif

      func1a = 0d0
      do k  =  point1,point2
         func   = c*exp(-a*cuspopt_xsav(k)) + b
         func1a = func1a + (func-cuspopt_ysav(k))**2
      enddo

      end function func1a

