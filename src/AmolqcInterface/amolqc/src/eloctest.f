
      MODULE eloctest


      use rwSampleModule
      use wfdata
      use utilsmodule
      use elocaldata
      use elocal, only: eloc
      !!!use ecp, only: psit

      implicit none

      private
      public :: runEloctest

      CONTAINS

      !-------------------------------------!
      subroutine runEloctest(lines,nl,sample)
      !-------------------------------------!

      integer, intent(in)           :: nl
      character(len=120), intent(in) :: lines(nl)
      type(rwSample),intent(inout)  :: sample

      integer :: points=7          ! points for differentiation rule
      integer :: sampleSize=0      ! number of points from sample
      real*8  :: h=1.d-5           ! denominator for numerical differentiation
      logical :: numDerivs=.false. ! do numerical derivatives?
      logical :: update=.false.
      real*8 x(ne),y(ne),z(ne)
      integer iflag,i
      type(randomWalker), pointer  :: rwp

      call assert(ne>0,
     .   ' runEloctest: wave function must be initialized')
      call assert(getSampleSize(sample) > 0,
     .   ' runEloctest: no sample available')

      h=1.d-5
      call getdbla(lines,nl,'h=',h,iflag)
      points = 7
      call getinta(lines,nl,'rule=',points,iflag)
      sampleSize=1
      call getinta(lines,nl,'points=',sampleSize,iflag)
      numDerivs = .false.
      numDerivs = finda(lines,nl,'derivatives')
      update = .false.
      update = finda(lines,nl,'updatetest')

      if (logmode >= 2) then
         write(iul,'(/A/)') '  * * *  Eloctest  * * *'
         write(iul,*) ' calculating psi and E_local '
         write(iul,*) '  with contributions for sample points'
      endif

      rwp => getFirst(sample)
      do i=1,sampleSize
         call pos(rwp,x,y,z)
         call elocContribs(x,y,z)
         if (numDerivs) then
            call elocDerivTest(x,y,z,h,points)
         endif
         if (update) then
            call updatetest(x,y,z)
         endif
         if (.not.isNext(sample)) exit
         rwp => getNext(sample)
      enddo

      if (logmode >= 2) then
         write(iul,*) ' * * * end Eloctest * * *'
      endif

      end subroutine


      !----------------------------!
      subroutine elocContribs(x,y,z)
      !----------------------------!

      real*8, intent(in)  :: x(:),y(:),z(:)

      integer i
      real*8 psiv,vpot,psilapl
      real*8 jas,phi,jlapl,flapl,jgrad(3*ne),fgrad(3*ne),lapli(ne)
      real*8 expU
      type(eConfigArray) :: ec

      if (logmode >= 2) then
         write(iul,*) ' electron coordinates:'
         do i=1,ne
            write(iul,'(i5,3F8.3)') i,x(i),y(i),z(i)
         enddo
      endif

      call eConfigArray_new(ec,ne,1)
      call eConfigArray_set(ec,1,x,y,z)
      call eloc(0,ec,'none')

      write(iul,'(3(A8,G14.8))') ' Eloc = ',elEloc(1),
     .      ' Phi = ',elPhi(1),' U = ',elU(1)
      write(iul,'(3(A8,G14.8))') ' Psi = ',elPhi(1)*exp(elU(1)),
     .      ' Vpot = ',elVPot(1)+vpot0,' Vpp = ',elECPPot
      write(iul,'(A8,G14.8)') ' Vnuc = ',vpot0

      expU = exp(elU(1))
      vpot  = elVPot(1)
      jas   = elU(1)
      phi   = elPhi(1)
      psiv  = elPhi(1)*expU
      jlapl = elUlapl(1)
      flapl = elFlapl(1)
      do i=1,ne
        jgrad(3*i-2) = elUgrad(3*i-2,1)
        jgrad(3*i-1) = elUgrad(3*i-1,1)
        jgrad(3*i)   = elUgrad(3*i,1)
        fgrad(3*i-2) = elFgrad(3*i-2,1)
        fgrad(3*i-1) = elFgrad(3*i-1,1)
        fgrad(3*i)   = elFgrad(3*i,1)
      enddo
      psilapl = flapl*jas + phi*jlapl + 2*dot_product(fgrad,jgrad)

      write(iul,'(2(A8,G17.7))') ' jas =   ',elU(1),' jlapl = ',jlapl
      write(iul,'(3(A8,G17.7))') ' Ekin1 = ',elEloc(1)-elVPot(1)-vpot0,
     .   ' Ekin2 = ',-0.5d0*psilapl/psiv,' Ekin3 = ',sum(elEkini(1:ne,1))

      if (logmode >= 2) then
         write(iul,*) ' grad Psi / Psi:'
         do i=1,ne
            write(iul,*) i,elxDrift(i,1),elyDrift(i,1),elzDrift(i,1)
         enddo
         write(iul,*) ' grad Phi:'
         do i=1,ne
            write(iul,*) i,fgrad(3*i-2),fgrad(3*i-1),fgrad(3*i)
         enddo
         write(iul,*) ' grad exp(U):'
         do i=1,ne
            write(iul,*) i,jgrad(3*i-2),jgrad(3*i-1),jgrad(3*i)
         enddo
      endif

      call eConfigArray_destroy(ec)

      end subroutine



      !-----------------------------------------!
      subroutine elocDerivtest(xx,yy,zz,h,points)
      !-----------------------------------------!

      real*8, intent(inout)  :: xx(:),yy(:),zz(:)
      real*8, intent(in)     :: h
      integer, intent(in)    :: points
      integer i,k,pnts
      real*8, pointer        :: x(:),y(:),z(:)
      real*8 psiv,tmp,expU
      real*8 jas,phi,jlapl,flapl,jgrad(3*ne),fgrad(3*ne),lapli(ne)
      real*8 psi0(8,3*ne),fx(3*ne),fxx(3*ne),f1
      real*8 psif(8,3*ne),ff1
      real*8 psij(8,3*ne),fj1
      real*8 vpots0,vpot,grad(3*ne)
      real*8 gr0(3*ne),grj0(3*ne),grf0(3*ne)
      real*8 numEloc,eloc0,lapl0, flapl0, jlapl0
      real*8 maxError,maxErrorPhi,maxErrorU
      type(eConfigArray) :: ec

      call eConfigArray_new(ec,ne,1)
      call eConfigArray_set(ec,1,xx,yy,zz)
      call eConfigArray_getPtr(ec,1,x,y,z)   ! x,y,z pointer to ec positions
      call eloc(0,ec,'none')

      maxError = 0
      maxErrorPhi = 0
      maxErrorU = 0
      expU = exp(elU(1))
      vpot  = elVPot(1)
      jas   = elU(1)
      phi   = elPhi(1)
      psiv  = elPhi(1)*expU
      jlapl = elUlapl(1)
      flapl = elFlapl(1)
      do i=1,ne
        grad(3*i-2)  = elxDrift(i,1)
        grad(3*i-1)  = elyDrift(i,1)
        grad(3*i)    = elzDrift(i,1)
        jgrad(3*i-2) = elUgrad(3*i-2,1)
        jgrad(3*i-1) = elUgrad(3*i-1,1)
        jgrad(3*i)   = elUgrad(3*i,1)
        fgrad(3*i-2) = elFgrad(3*i-2,1)
        fgrad(3*i-1) = elFgrad(3*i-1,1)
        fgrad(3*i)   = elFgrad(3*i,1)
      enddo


c     // original values
      f1 = psiv
      fj1 = jas
      ff1 = phi
      eloc0 = elEloc(1)
      jlapl0 = jlapl
      flapl0 = flapl
      vpots0 = elVPot(1)+vpot0
      lapl0 = -2d0*(elEloc(1) - vpots0)*psiv
      do i=1,3*ne
         gr0(i) = grad(i)*psiv
         grj0(i) = jgrad(i)
         grf0(i) = fgrad(i)
      enddo

c     // Calculate points for numerical derivaties

       do k=1,ne
          do pnts=1,points/2
             x(k) = x(k) + pnts*h
             call eloc(0,ec,'none')
             psi0(2*pnts-1,3*k-2) = elPhi(1)*exp(elU(1))
             psij(2*pnts-1,3*k-2) = elU(1)
             psif(2*pnts-1,3*k-2) = elPhi(1)
             x(k) = x(k) - 2*pnts*h
             call eloc(0,ec,'none')
             psi0(2*pnts,3*k-2) = elPhi(1)*exp(elU(1))
             psij(2*pnts,3*k-2) = elU(1)
             psif(2*pnts,3*k-2) = elPhi(1)
             x(k) = x(k) + pnts*h

             y(k) = y(k) + pnts*h
             call eloc(0,ec,'none')
             psi0(2*pnts-1,3*k-1) = elPhi(1)*exp(elU(1))
             psij(2*pnts-1,3*k-1) = elU(1)
             psif(2*pnts-1,3*k-1) = elPhi(1)
             y(k) = y(k) - 2*pnts*h
             call eloc(0,ec,'none')
             psi0(2*pnts,3*k-1) = elPhi(1)*exp(elU(1))
             psij(2*pnts,3*k-1) = elU(1)
             psif(2*pnts,3*k-1) = elPhi(1)
             y(k) = y(k) + pnts*h

             z(k) = z(k) + pnts*h
             call eloc(0,ec,'none')
             psi0(2*pnts-1,3*k) = elPhi(1)*exp(elU(1))
             psij(2*pnts-1,3*k) = elU(1)
             psif(2*pnts-1,3*k) = elPhi(1)
             z(k) = z(k) - 2*pnts*h
             call eloc(0,ec,'none')
             psi0(2*pnts,3*k) = elPhi(1)*exp(elU(1))
             psij(2*pnts,3*k) = elU(1)
             psif(2*pnts,3*k) = elPhi(1)
             z(k) = z(k) + pnts*h
          enddo
       enddo

c-----calculate numerical derivatives
c      // second order formula for 1st derivative
c      // 3/5/7 point formula dependent on "points" for 2nd derivative

c     // total wavefunction
      flapl = 0d0
      do k=1,3*ne
          if (points .eq. 3) then
             fx(k)  = (psi0(1,k) - psi0(2,k))/(2d0*h)
             fxx(k) = ( -2d0*f1 + psi0(2,k) + psi0(1,k) ) / h**2
          else if  (points .eq. 5) then
             fx(k)  = ( 8d0*(psi0(1,k) - psi0(2,k)) -
     .                      (psi0(3,k) - psi0(4,k)) ) /(12d0*h)
             fxx(k) = ( -30d0*f1 + 16d0*(psi0(1,k) + psi0(2,k))
     .                 - (psi0(3,k) + psi0(4,k)) ) / (12d0*h**2)
          else if  (points .eq. 7) then
             fx(k)  = ( 45d0*(psi0(1,k) - psi0(2,k)) -
     .                   9d0*(psi0(3,k) - psi0(4,k)) +
     .                       (psi0(5,k) - psi0(6,k))) /(60d0*h)
             fxx(k) = ( -490d0*f1 + 270d0*(psi0(1,k) + psi0(2,k))
     .                 - 27d0*(psi0(3,k) + psi0(4,k))
     .                 + 2d0*(psi0(5,k) + psi0(6,k)) ) / (180d0*h**2)
          else if (points .eq. 9) then
             fx(k) = (672d0*(psi0(1,k) - psi0(2,k)) -
     .                168d0*(psi0(3,k) - psi0(4,k)) +
     .                 32d0*(psi0(5,k) - psi0(6,k)) -
     .                  3d0*(psi0(7,k) - psi0(8,k))) / (840d0*h)
             fxx(k) =(-71750d0*f1 + 40320d0*(psi0(1,k) + psi0(2,k))
     .                - 5040d0*(psi0(3,k) + psi0(4,k))
     .                +  640d0*(psi0(5,k) + psi0(6,k))
     .                -   45d0*(psi0(7,k) + psi0(8,k))) / (25200d0 * h**2)
          else
             call abortp('(eloctest): wrong input for points')
          endif
          flapl = flapl + fxx(k)
      enddo
      numEloc = -0.5d0*flapl/f1 + vpots0

      write(iul,*) ' '
      write(iul,*) ' numerical derivatives for Psi:'
      write(iul,*) ' computed, numerical, abs error, rel error'
      write(iul,'(A,4G17.7)') ' laplacian:',lapl0,flapl,lapl0-flapl,
     .    (lapl0-flapl)/lapl0
      write(iul,'(A,4G17.7)') ' Eloc:',eloc0,numEloc,eloc0-numEloc,
     .    (numEloc-eloc0)/eloc0
      if (logmode >= 3) then
         write(iul,*) ' gradient:'
         do i=1,3*ne
          write(iul,'(i5,4g17.7)') i,gr0(i),fx(i),(fx(i)-gr0(i)),
     .       (fx(i)-gr0(i))/gr0(i)
         enddo
      endif
      if (logmode >=2) then
         maxError = 0
         do i=1,3*ne
            maxError = max(maxError,abs((fx(i)-gr0(i))/gr0(i)))
         enddo
         write(iul,*) ' maximal relative gradient error:',maxError
      endif

c     // Jastrow part
      flapl = 0d0
      do k=1,3*ne
          if (points .eq. 3) then
             fx(k)  = (psij(1,k) - psij(2,k))/(2d0*h)
             fxx(k) = ( -2d0*fj1 + psij(2,k) + psij(1,k) ) / h**2
          else if  (points .eq. 5) then
             fx(k)  = ( 8d0*(psij(1,k) - psij(2,k)) -
     .                      (psij(3,k) - psij(4,k)) ) /(12d0*h)
             fxx(k) = ( -30d0*fj1 + 16d0*(psij(1,k) + psij(2,k))
     .                 - (psij(3,k) + psij(4,k)) ) / (12d0*h**2)
          else if  (points .eq. 7) then
             fx(k)  = ( 45d0*(psij(1,k) - psij(2,k)) -
     .                   9d0*(psij(3,k) - psij(4,k)) +
     .                       (psij(5,k) - psij(6,k))) /(60d0*h)
             fxx(k) = ( -490d0*fj1 + 270d0*(psij(1,k) + psij(2,k))
     .                 - 27d0*(psij(3,k) + psij(4,k))
     .                 + 2d0*(psij(5,k) + psij(6,k)) ) / (180d0*h**2)
          else if (points .eq. 9) then
             fx(k) = (672d0*(psij(1,k) - psij(2,k)) -
     .                168d0*(psij(3,k) - psij(4,k)) +
     .                 32d0*(psij(5,k) - psij(6,k)) -
     .                  3d0*(psij(7,k) - psij(8,k))) / (840d0*h)
             fxx(k) =(-71750d0*fj1 + 40320d0*(psij(1,k) + psij(2,k))
     .                - 5040d0*(psij(3,k) + psij(4,k))
     .                +  640d0*(psij(5,k) + psij(6,k))
     .                -   45d0*(psij(7,k) + psij(8,k))) / (25200d0 * h**2)
          else
             call abortp('(elocwf): wrong input for points')
          endif
          flapl = flapl + fxx(k)
      enddo
      write(iul,*) ' '
      write(iul,*) ' JASTROW part: numerical derivatives'
      write(iul,'(A,4G17.7)') ' laplacian:',jlapl0,flapl,jlapl0-flapl,
     .    (jlapl0-flapl)/jlapl0
      if (logmode >= 3) then
         write(iul,*) ' gradient:'
         do i=1,3*ne
          write(iul,'(I5,4G17.7)') i,grj0(i),fx(i),(fx(i)-grj0(i)),
     .       (fx(i)-grj0(i))/grj0(i)
         enddo
      endif
      if (logmode >=2) then
         maxError = 0
         do i=1,3*ne
            maxError = max(maxError,abs((fx(i)-grj0(i))/grj0(i)))
         enddo
         write(iul,*) ' maximal relative gradient error:',maxError
      endif

c     // Phi part
      flapl = 0d0
      do k=1,3*ne
          if (points .eq. 3) then
          fx(k)  = (psif(1,k) - psif(2,k))/(2d0*h)
             fxx(k) = ( -2d0*ff1 + psif(2,k) + psif(1,k) ) / h**2
          else if  (points .eq. 5) then
             fx(k)  = ( 8d0*(psif(1,k) - psif(2,k)) -
     .                      (psif(3,k) - psif(4,k)) ) /(12d0*h)
             fxx(k) = ( -30d0*ff1 + 16d0*(psif(1,k) + psif(2,k))
     .                 - (psif(3,k) + psif(4,k)) ) / (12d0*h**2)
          else if  (points .eq. 7) then
             fx(k)  = ( 45d0*(psif(1,k) - psif(2,k)) -
     .                   9d0*(psif(3,k) - psif(4,k)) +
     .                       (psif(5,k) - psif(6,k))) /(60d0*h)
             fxx(k) = ( -490d0*ff1 + 270d0*(psif(1,k) + psif(2,k))
     .                 - 27d0*(psif(3,k) + psif(4,k))
     .                 + 2d0*(psif(5,k) + psif(6,k)) ) / (180d0*h**2)
          else if (points .eq. 9) then
             fx(k) = (672d0*(psif(1,k) - psif(2,k)) -
     .                168d0*(psif(3,k) - psif(4,k)) +
     .                 32d0*(psif(5,k) - psif(6,k)) -
     .                  3d0*(psif(7,k) - psif(8,k))) / (840d0*h)
             fxx(k) =(-71750d0*ff1 + 40320d0*(psif(1,k) + psif(2,k))
     .                - 5040d0*(psif(3,k) + psif(4,k))
     .                +  640d0*(psif(5,k) + psif(6,k))
     .                -   45d0*(psif(7,k) + psif(8,k))) / (25200d0 * h**2)
          else
             call abortp('(eloctest): wrong input for points')
          endif
          flapl = flapl + fxx(k)
      enddo
      write(iul,*) ' '
      write(iul,*) ' DETERMINANT PART: numerical derivatives:'
      write(iul,'(A,4G17.7)') ' laplacian:',flapl0,flapl,flapl0-flapl,
     .    (flapl0-flapl)/flapl0
      if (logmode >= 3) then
         write(iul,*) ' gradient:'
         do i=1,3*ne
          write(iul,'(I5,4G17.7)') i,grf0(i),fx(i),(fx(i)-grf0(i)),
     .       (fx(i)-grf0(i))/grf0(i)
         enddo
      endif
      if (logmode >=2) then
         maxError = 0
         do i=1,3*ne
            maxError = max(maxError,abs((fx(i)-grf0(i))/grf0(i)))
         enddo
         write(iul,*) ' maximal relative gradient error:',maxError
      endif


      call eConfigArray_destroy(ec)

      end subroutine

c     -------------------------------
      subroutine updatetest(xx,yy,zz)
c     -------------------------------
      real*8,intent(inout) :: xx(:),yy(:),zz(:)
      real*8, parameter    :: dt = 0.01    ! time step for move
      real*8               :: x1(ne),y1(ne),z1(ne)
      real*8, pointer      :: x(:),y(:),z(:)
      real*8 phi,u
      integer ie
      type(eConfigArray) :: ec

      call eConfigArray_new(ec,ne,1)
      call eConfigArray_set(ec,1,xx,yy,zz)
      call eConfigArray_getPtr(ec,1,x,y,z)   ! x,y,z pointer to ec positions

      write(iul,*)
      write(iul,*) ' test of psi update with psit:: currently disabled'
      call abortp("UPDATE CURRENTLY NOT POSSIBLE")

      call eloc(0,ec,'none')
      write(iul,'(2(A8,G14.8))') ' Phi = ',elPhi(1),' U = ',elU(1)

      x1 = elxDrift(1:ne,1)*dt
      y1 = elyDrift(1:ne,1)*dt
      z1 = elzDrift(1:ne,1)*dt

      x(1) = x(1) + x1(1)
      y(1) = y(1) + y1(1)
      z(1) = y(1) + z1(1)

      ie = 1
      !!!!!!loc_method = "SD*J"

      !!!!!!call psit(.true.,ie,x,y,z,phi,u)
      
      write(iul,*) ' update with psit'
      write(iul,'(2(A8,G14.8))') ' Phi = ',phi,' U = ',u

      call eloc(0,ec,'none')
      write(iul,*) ' check with eloc'
      write(iul,'(2(A8,G14.8))') ' Phi = ',elPhi(1),' U = ',elU(1)
      call assertEqualRelative(phi,elPhi(1),1.d-5,
     .                         '(updatetest): phi not equal')
      call assertEqualRelative(u,elU(1),1.d-5,
     .                         '(updatetest): u not equal')

      do ie=2,ne
         x(ie) = x(ie) + x1(ie)
         y(ie) = y(ie) + y1(ie)
         z(ie) = y(ie) + z1(ie)
         !!!!!!!!!!!!!call psit(.false.,ie,x,y,z,phi,u)
         write(iul,*) ' update with psit ie=',ie
         write(iul,'(2(A8,G14.8))') ' Phi = ',phi,' U = ',u
      enddo

      call eloc(0,ec,'none')
      write(iul,*) ' check with eloc'
      write(iul,'(2(A8,G14.8))') ' Phi = ',elPhi(1),' U = ',elU(1)
      call assertEqualRelative(phi,elPhi(1),1.d-5,
     .                         '(updatetest): phi not equal')
      call assertEqualRelative(u,elU(1),1.d-5,
     .                         '(updatetest): u not equal')

      write(iul,*) ' PASSED: eloctest -- updatetest'

      call eConfigArray_destroy(ec)

      end subroutine


      END MODULE eloctest
