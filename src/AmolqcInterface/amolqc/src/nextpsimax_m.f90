!
! implementation of the L-BFGS algorithm for finding maxima
! core electrons are identified and put at the nuclei
! TODO: core identification with a threshold like H core identification
! (threshold dependend on Z)


module nextpsimaxModule
   use global
   use eConfigsModule
   use wfdata
   use findCoreModule
   use elocaldata, only: elxDrift,elyDrift,elzDrift,elPhi,elU,elEloc,elVPot,eloc_isAtSingularity
   use elocal, only: eloc
   use lbfgsb3

   implicit none

   private
   public :: nextpsimax


contains

   subroutine nextpsimax(x,y,z,maxiter,HThresh,eps,verb0,f,nCore,ierr,idx)
      real*8, intent(inout) :: x(:),y(:),z(:)  ! coords
      integer, intent(in)   :: maxiter         ! max # of iterations of minimizer
      real*8, intent(in)    :: HThresh         ! put electrons on H nuc below this distance
      real*8, intent(in)    :: eps(3)          ! epsilon (1:grad, 2:dist, 3:fctn)
      integer, intent(in)   :: verb0           ! verbose level, output only if >= 2
      real*8, intent(inout) :: f               ! function value at minimum
      integer, intent(out)  :: nCore           ! # of elecs at nuclei (on exit)
      integer, intent(out)  :: ierr            ! exit codes:
         ! 0: all crit met, 1: function value not converged, 2: grad not converged
         ! 3: distance not converged, 4: boundary hit, 5: some error, 6: singularity in eloc, 7: f==NAN
      integer, intent(inout), optional :: idx(:)   ! update permutation index due to core electron permutation
      integer             :: naCore, nbCore    ! alpha, beta elecs at nuclei
      integer             :: aCoreList(ncenter), bCoreList(ncenter)
      type(eConfigArray)  :: eca
      logical found
      real*8 deltaf, oldf
      real*8 lowerb, upperb
      integer, save :: sNCALLS = 0

!     LBFGS parameters
!     nmax  is the dimension of the largest problem to be solved.
!     mmax  is the maximum number of limited memory corrections.
!     lenwa is the corresponding real workspace required.
      integer, parameter  :: nmax=360
      integer, parameter  :: mmax = 17
      integer, parameter  :: lenwa = 2*mmax*nmax + 4*nmax + 11*mmax*mmax + 8*mmax
      character*60     task, csave
      logical          lsave(4)
      integer          i,n,m,iprint,restart,myiter,iulref,verb
      integer          nbd(nmax), iwa(3*nmax), isave(44)
      real*8           ff, factr, pgtol
      real*8           xx(nmax), l(nmax), u(nmax), gg(nmax), dsave(29)
      real*8           wa(lenwa)

      call assert(3*getNElec()<=nmax,'(nextpsimax): too many electrons - change nmax')

      ! discard point when electron hits lower or upper bound
      lowerb = -100.d0      ! should be determined from atom box!! could be parameters
      upperb =  100.d0

      sNCALLS = sNCALLS + 1

      verb = verb0 - 2    ! local verbosity: starting output from verb=1

      if (verb >= 3) then
         write(iull,*) " *** nextpsimax ",sNCALLS,"  *** "
         do i=1,getNElec()
            write(iull,'(i5,3F14.7)') i,x(i),y(i),z(i)
         enddo
         if (verb >= 4) then
            iulref = 300 + getMyTaskId()
            write(iulref,'(a,i6,a)') "LBFGS:",sNCALLS,"  1  F(Max):  0.0  found: 0 0 0 0 0"
            write(iulref,'(i5)') getNElec()
            do i=1,getNElec()
               write(iulref,'(3F14.7)') x(i),y(i),z(i)
            enddo
         endif
      endif

      call eConfigArray_new(eca,ne,1)
      call initCoreLists(naCore,aCoreList,nbCore,bCoreList)
      call findCoreElecs(naCore,aCoreList,nbCore,bCoreList,x,y,z,idx)
      call findHCoreElecs(naCore,aCoreList,nbCore,bCoreList,HThresh,x,y,z,found,idx)

      call eConfigArray_set(eca,1,x(1:ne),y(1:ne),z(1:ne))
      call eloc(0,eca,'none')

      call getxfg(x,y,z,naCore,nbCore,xx,ff,gg)

      if (verb >= 4) then
         write(iulref,'(a,i6,a,f13.3,a)') "LBFGS:",sNCALLS,"   2  F(Max): ",ff,"  found: 0 0 0 0 0"
         write(iulref,'(i5)') getNElec()
         do i=1,getNElec()
            write(iulref,'(3F14.7)') x(i),y(i),z(i)
         end do
         myiter = 2
         write(iull,'(3(a,i3),a,f13.6)') "myiter=",myiter," nac=",naCore," nbc=",nbCore," ff=",ff
      endif

      iprint = -1     ! We suppress the output for L-BFGS
      factr  = 0.0d0  ! use own stopping criteria
      pgtol  = 0.0d0
      n      = 3*(getNElec()-naCore-nbCore)  ! dimension of optimization problem
      m      =  5  ! number of limited memory corrections

      ! We now specify nbd which defines the bounds on the variables:
      ! l   specifies the lower bounds, u   specifies the upper bounds.
      do i = 1, n
         nbd(i) = 2
         l(i)   = lowerb-1.d-4   ! lower than lowerb
         u(i)   = upperb+1.d-4    ! higher than upperb
      enddo

      task = 'START'
      restart = 0
      oldf = ff
      deltaf = ff
      ierr = 5

      ! L-BFGS-B loop
      do
         if (n==0) then    ! all electrons are at cores including H. We are done
            if (verb >=3) write(iull,'(a)') ' -> L-BFGS not called. All electrons are at nuclei'
            ierr = 0
            exit
         end if
         call setulb(n,m,xx,l,u,nbd,ff,gg,factr,pgtol,wa,iwa,task,iprint, &
                     csave,lsave,isave,dsave)
         if (.not.all(abs(xx(1:n))<huge(1.d0))) then
            if (verb>=3) write(iull,'(a,i5,a)') ' -> L-BFGS produced NaN or Inf at iter=',isave(34), &
              ': leaving L-BFGS and discarding point'
            ff = 99999.99999d0
            ierr = 5
            exit
         else if (any(xx(1:n)<=lowerb) .or. any(xx(1:n)>=upperb)) then
            if (verb>=3) write(iull,'(a,i5,a)') ' -> boundary hit at iter=',isave(34),': leaving L-BFGS and discarding point'
            ff = 88888.88888d0
            ierr = 4
            exit
         endif
         if (task(1:2) .eq. 'FG') then
            !!write(iull,*) "DBG:"
            !!write(iull,'(5g12.5)') xx(1:n)

            call putx(xx,naCore,nbCore,x,y,z)

            !!do i=1,getNElec()
            !!   write(iull,'(i5,3F14.7)') i,x(i),y(i),z(i)
            !!end do

            call findHCoreElecs(naCore,aCoreList,nbCore,bCoreList,HThresh,x,y,z,found,idx)

            !!do i=1,getNElec()
            !!   write(iull,'(i5,3F14.7)') i,x(i),y(i),z(i)
            !!end do

            if (found) then
               task = 'START'   ! restart with less electrons
               n =  3*(getNElec()-naCore-nbCore)
               factr  = 0.0d0  ! own stopping criteria
               pgtol  = 0.0d0
               m = 5
               restart = restart + 1
               if (verb >= 4) write(iull,*) ' -> restarting L-BFGS with n=',n,' restart=',restart
            endif
            call eConfigArray_set(eca,1,x(1:ne),y(1:ne),z(1:ne))
            call eloc(0,eca,'none')
            if (eloc_isAtSingularity()) then
               if (verb >= 3) write(iull,*) ' -> hitting singularity'
               ff = 99999.99999d0
               ierr = 6
               exit
            endif
            call getxfg(x,y,z,naCore,nbCore,xx,ff,gg)

            if (verb >= 4) then
               myiter = myiter + 1
               write(iulref,'(a,i6,i4,a,f13.3,a)') "LBFGS:",sNCALLS,myiter," F(Max): ",ff,"  found: 0 0 0 0 0"
               write(iulref,'(i5)') getNElec()
               do i=1,getNElec()
                  write(iulref,'(3F14.7)') x(i),y(i),z(i)
               enddo
               write(iull,'(3(a,i3),a,f13.6)') "myiter=",myiter," nac=",naCore," nbc=",nbCore," ff=",ff
            endif

            deltaf = abs(ff-oldf)
            oldf = ff
         else if (task(1:5) .eq. 'NEW_X') then
            if (dsave(13) <= eps(1)*(1.0d0 + abs(ff)) .and.  &     ! grad crit
                dsave(4)  <= eps(2)                   .and.  &     ! dist crit
                deltaf    <= eps(3)*(1.d0 + abs(ff))) then         ! delta f crit
               ierr = 0
               task='STOP: ALL CONVERGENCE CRITERIA MET'
            else if (isave(34) >= maxiter) then
               if (deltaf > eps(3)*(1.d0 + abs(ff))) then
                  ierr = 1
               else if (dsave(13) > eps(1)*(1.0d0 + abs(ff))) then
                  ierr = 2
               else
                  ierr = 3
               end if
               task='STOP: TOTAL NO. OF FUNCTION EVALUATIONS EXCEEDS LIMIT'
            endif
            ! the current iteration number, isave(30),
            ! the total number of f and g evaluations, isave(34),
            ! the value of the objective function f,
            ! the norm of the projected gradient,  dsave(13)
            if (verb >= 4) write(iull,'(2(a,i6),3(a,g15.5))') ' Iter=',isave(30),' nfg=',isave(34), &
               ' f=',ff, 'max grad =',dsave(13),' dist =',dsave(4)
         else
            if (task(1:4) == 'STOP') then
               if (verb >= 4) write(iull,*) 'leaving L-BFGS with ',trim(task)
            else
               ierr = 5
               if (verb >= 4) write(iull,*) 'leaving L-BFGS with ',trim(task)
            endif
            exit
         endif
      enddo

      call putx(xx,naCore,nbCore,x,y,z)
      f = ff
      nCore = naCore + nbCore

      if (isnan(f)) ierr = 7

      if (verb >= 1) then
         write(iull,'(a,f12.5,4(a,i4))') ' nextpsimax result: final f=',f,' iter=',isave(34),' ierr=',ierr, &
                                         ' nCore=',nCore,' restarts=',restart
         if (verb >= 2) then
            write(iull,'(3(a,g15.5))') ' with criteria: grad=',dsave(13),' dist=',dsave(4),' delta f=',deltaf
            if (verb >= 3) then
               write(iull,*) "with x,y,z:"
               do i=1,getNElec()
                  write(iull,'(i5,3f12.5)') i,x(i),y(i),z(i)
               enddo
               if (verb >= 4) then
                  write(iull,*) "with xx:"
                  do i=1,getNElec()-naCore-nbCore
                     write(iull,'(i5,3f12.5)') i,xx(3*i-2:3*i)
                  enddo
                  write(iull,*) "with gg:"
                  do i=1,getNElec()-naCore-nbCore
                     write(iull,'(i5,3f12.5)') i,gg(3*i-2:3*i)
                  enddo
               endif
               write(iull,*)
            endif
         endif
      endif

      call eConfigArray_destroy(eca)

   end subroutine nextpsimax

   subroutine getxfg(x,y,z,naCore,nbCore,xx,ff,gg)
      real*8, intent(in)    :: x(:),y(:),z(:)
      integer, intent(in)   :: naCore,nbCore
      real*8, intent(inout) :: xx(:),ff,gg(:)
      integer ia,ib,ib0,ib1,ib2,na
      na = getNAlpha()
      ia = na - naCore
      ib = getNBeta() - nbCore
      ib0 = 3*ia
      ib1 = na + nbCore + 1
      ib2 = getNElec()
      xx(1:ia)                = x(naCore+1:na)
      xx(ia+1:2*ia)           = y(naCore+1:na)
      xx(2*ia+1:3*ia)         = z(naCore+1:na)
      xx(ib0+1:ib0+ib)        = x(ib1:ib2)
      xx(ib0+ib+1:ib0+2*ib)   = y(ib1:ib2)
      xx(ib0+2*ib+1:ib0+3*ib) = z(ib1:ib2)
      ff = - 2.d0*(log(abs(elPhi(1))) + elU(1))
      gg(1:ia)                = -elxDrift(naCore+1:na,1)
      gg(ia+1:2*ia)           = -elyDrift(naCore+1:na,1)
      gg(2*ia+1:3*ia)         = -elzDrift(naCore+1:na,1)
      gg(ib0+1:ib0+ib)        = -elxDrift(ib1:ib2,1)
      gg(ib0+ib+1:ib0+2*ib)   = -elyDrift(ib1:ib2,1)
      gg(ib0+2*ib+1:ib0+3*ib) = -elzDrift(ib1:ib2,1)
      gg = 2.d0*gg
   end subroutine getxfg

   subroutine putx(xx,naCore,nbCore,x,y,z)
      real*8, intent(in)     :: xx(:)
      integer, intent(in)    :: naCore,nbCore
      real*8, intent(inout)  :: x(:),y(:),z(:)
      integer ia,ib,ib0,ib1,ib2,na
      na = getNAlpha()
      ia = na - naCore
      ib = getNBeta() - nbCore
      ib0 = 3*ia
      ib1 = na + nbCore +1
      ib2 = getNElec()
      x(naCore+1:na) = xx(1:ia)
      y(naCore+1:na) = xx(ia+1:2*ia)
      z(naCore+1:na) = xx(2*ia+1:3*ia)
      x(ib1:ib2) = xx(ib0+1:ib0+ib)
      y(ib1:ib2) = xx(ib0+ib+1:ib0+2*ib)
      z(ib1:ib2) = xx(ib0+2*ib+1:ib0+3*ib)
   end subroutine putx

end module nextpsimaxModule
