!
! module to perform the maximization of psi**2
! implemented as minimization of -ln(psi**2)
!
! currently uses:
!   swarm optimization RPSO
!   L-BFGS, BFGS
!   Nelder Mead (NELM)
!
! 2015/3: split into psimax_m and MAXIMA_m
!   psimax does only the actual minimization of -ln|psi^2|
!   psimax depends on maxanalysis and maxbasin which it calls
!   psimax is the interface to maxanalysis and maxbasin, including epart
!   minimization is done fully in parallel
!
! 2013: there are basically two codes. Both use:
!   psimax_init: initialization from .in file
!   psimax_opt:  local maximization for one walker called from (e.g.) qmc
!     old code (RP) RPSO/BFGS/NELM using psimax_calc / psimax_add_to_list
!     new code (AL) L-BFGS using nextpsimax / psimax_addToList
!        new code uses sophisticated "referenceContainer" object for storing/sorting references
!
! parallel version:
!   both old and new codes run optimizer locally and collect local maxima on MASTER (in add_to_list/addTOList)
!   ONLY the MASTER keeps/stores/sorts references
!   new code: ONLY MASTER has/needs referenceContainer object allocated (mRC_p)

module psimax

   use minimization   ! old code: BFGS
   use nextpsimaxModule, only: nextpsimax  ! new code: L-BFGS
   use findcores, only: set_core_rule, findcores_find, findcores_sort   ! old code only
   use RandomWalkerModule
   use wfdata, only: atoms
   use maxanalysisModule, only: maxana_isInitialized, maxana_destroy, maxana_add, maxana_writeResults, &
       maxana_getFirstF, maxana_getDiffMax
   use maxbasinsModule, only: maxbas_isInitialized, maxbas_destroy, maxbas_add, maxbas_writeResults

   implicit none
   private
   public :: psimax_init, psimax_destroy, psimax_opt, psimax_writeResults, psimax_writeParams, psimax_doMaxAnalysis, &
             psimax_doMaxSearch, psimax_getDiffMax, psimax_getFirstF, rebuildvec, writeCoordsToVec

   integer, parameter  :: RPSO=1,BFGS=2,NELM=3,LBFGS=4
   integer, parameter  :: MAXIMA=1,BASINS=2
   integer             :: mOptMode = LBFGS  ! 2=bfgs, 3=nelder/mead (nelm), (1=rpso disabled)
   integer             :: mVerbose=0
   integer             :: mAnalyseMode=0   ! allowed values: MAXIMA|BASINS

   ! parameter for L-BFGS
   integer             :: mMaxIter = 150
   real*8              :: mHDist = -1.d0   ! if r_ai is lower (for H atom) electron is put at nucleus (<0 means turned off)
   real*8              :: mEPSGrad = 1.d-6    ! grad crit (relative)
   real*8              :: mEPSDist = 1.d-5    ! dist crit (absolute)
   real*8              :: mEPSFctn = 1.d-6    ! delta f crit (relative)

   real*8              :: mFtol = 1.d-5
   real*8              :: mTolFac = 10   !
   real*8              :: mTol         ! Tolerance=mFtol/mTolFac for the optimizer (factor 10 is empirical)
   real*8              :: mGeomTol = 1d-5 ! Crit. for geometric similarity (<mGeomTol --> equal)

   logical             :: mAllowSingularities = .false.

   ! Nelder Mead
   real*8         :: mDelta   = 0.2d0
   logical        :: mNelmIt  = .false.

   ! BFGS and Nelder Mead
   ! Separate Cores, i.e. put two electrons at nucleus for Z>2
   logical       :: mConstraint = .true.
   integer       :: mNocel = 0
   type(coord)   :: mCoreCoords
   integer       :: mCoreRule = 1  ! 0=no core sep, 1=He core, 2=Ne core???
   ! Bound Constraint
   real*8        :: mBound = 20d0


   ! counters variables:
   real*8       :: sCounter(0:9) = 0.d0  ! counters: exit codes 0:7, 8 all, 9 successful minimizer calls

contains


   subroutine psimax_init(lines,nl)
      character(len=120), intent(in) :: lines(:)
      integer, intent(in)            :: nl
      character(len=10)              :: s
      integer                        :: iflag,maxFull,iflag1,iflag2,iflag3,idx,nre,nrefs,i,j
      integer, pointer               :: ire(:,:) => null()
      real*8                         :: st,distThr,tolDist2,tolMeanDist,tolSim,tolSame
      real*8                         :: nucThresh,coreThresh,bondThresh
      logical                        :: finda,doSortFreq,l1,l2
      character(len=3)               :: maxMode
      character(len=40)              :: reffile,exclfile

      mVerbose = logmode
      call getinta(lines,nl,'verbose=',mVerbose,iflag)

      call getinta(lines,nl,'max_iter=',mMaxIter,iflag)
      if (iflag==0) then
         call minimization_setIterMax(mMaxIter)
      end if

      ! Nelder Mead and BFGS only
      call getdbla(lines,nl,'step_max=',st,iflag)
      if (iflag==0) call minimization_setStepMax(st)

      ! LBFGS options
      mHDist = -1.d0
      call getdbla(lines,nl,'H_dist=',mHDist,iflag)
      mHDist = mHDist / bohr2angs
      mEPSGrad = 1.d-5
      call getdbla(lines,nl,'max_grad=',mEPSGrad,iflag)
      mEPSDist = 1.d-4
      call getdbla(lines,nl,'max_dist=',mEPSDist,iflag)
      mEPSFctn = 1.d-5
      call getdbla(lines,nl,'max_f=',mEPSFctn,iflag)
      mAllowSingularities = .false.
      if (finda(lines,nl,'no_allow_sing')) then
         mAllowSingularities = .false.
      else if (finda(lines,nl,'allow_sing')) then
         mAllowSingularities = .true.
      end if

      mOptMode = LBFGS
      call getstra(lines,nl,'method=',s,iflag)
      if(iflag == 0)then
         !if (trim(s)=='rpso') then
         !    !mOptMode = 1
         !    write(iul,*) "RPSO currently inactive due to stability/testing issues"
         !    stop
         if (s=='lbfgs')then
            mOptMode = LBFGS
         else if (s=='nelm') then
            mOptMode = NELM
         else if (s=='bfgs') then
            mOptMode = BFGS
         else
            call abortp('$init_maxima_search: unknown mode given')
         endif
      endif


      ! Nelder Mead and BFGS only
      if(finda(lines,nl,'numgrad'))then
          call minimization_setGrad(.true.)
      endif

      ! Nelder Mead and BFGS only
      if (finda(lines,nl,'no_core_sep')) then
         mConstraint = .false.
         mCoreRule = 0
         if (mOptMode==LBFGS) call abortp("init_maxima_search: no_core_sep no implemented for LBFGS")
      else if (finda(lines,nl,'core_sep')) then
         mConstraint = .true.
         mCoreRule = 1
      end if
      if (mConstraint) then
          call set_core_rule(mCoreRule)
          call findcores_find(mNocel)
          call psimax_set_core_el_coords()
      endif

      ! Nelder Mead
      call getdbla(lines,nl,'border=',mBound,iflag)
      call getdbla(lines,nl,'delta=',mDelta,iflag)

      ! BFGS code
      call getdbla(lines,nl,'ftol=',mFtol,iflag)
      call getdbla(lines,nl,'tolfac=',mTolFac,iflag)
      mTol = mFtol / mTolFac

      sCounter = 0.d0;

      l1 = maxana_isInitialized()
      l2 = maxbas_isInitialized()
      if (.not.(l1.or.l2) .or. (l1.and.l2)) then
         call abortp(' $init_maxima_search requires previous maxima or basin analysis initialization')
      endif
      if (l1) then
         mAnalyseMode = MAXIMA
      else if (l2) then
         mAnalyseMode = BASINS
         if (mOptmode /= LBFGS) then
            call abortp(" $init_maxima_search: basin analysis requires currently method=lbfgs")
         endif
      endif
   end subroutine psimax_init


   subroutine psimax_destroy()
      sCounter = 0.d0
      select case(mAnalyseMode)
      case (MAXIMA)
         call maxana_destroy()
      case (BASINS)
         call maxbas_destroy()
      end select
   end subroutine psimax_destroy


   function psimax_doMaxSearch() result(res)
      logical :: res
      res = (mAnalyseMode > 0)
   end function psimax_doMaxSearch


   function psimax_doMaxAnalysis() result(res)
      logical :: res
      res = (mAnalyseMode == MAXIMA)
   end function psimax_doMaxAnalysis


   function psimax_getFirstF() result(res)
      real*8 :: res
      res = 0.d0
      if (mAnalyseMode==MAXIMA) then
         res = maxana_getFirstF()
      endif
   end function psimax_getFirstF


   function psimax_getDiffMax() result(res)
      integer :: res
      res = 0
      if (mAnalyseMode==MAXIMA) then
         res = maxana_getDiffMax()
      endif
   end function psimax_getDiffMax


   subroutine psimax_writeParams(iu)
      integer, intent(in) :: iu
      character(len=5)    :: s
      character(len=10)   :: s1
      character(len=120)   :: s2

      write(iu,'(A/)') '    maxima search parameters:'

      select case (mOptMode)
      case (RPSO)
         s = 'rpso'
      case (BFGS)
         s = 'bfgs'
      case (NELM)
         s = 'nelm'
      case (LBFGS)
         s = 'lbfgs'
      end select

      if (mOptmode==LBFGS) then

         select case(mAnalyseMode)
         case (MAXIMA)
            write(iu,'(1X,A34)') 'L-BFGS search with maxima analysis'
         case (BASINS)
            write(iu,'(1X,A34)') 'L-BFGS search with basin analysis'
         case default
            write(iu,'(1X,A34)') 'L-BFGS search with no analysis'
         end select

         if (mConstraint) then
            write(iu,'(4X,A)') 'maximization with 1s core electrons fixed at nucleus'
         else
            write(iu,'(4X,A)') 'maximization for all electrons'
         end if
         write(iu,'(1X,A21,I5)') 'max_iter =',mMaxIter
         write(iu,'(2(1X,A21,E12.2,2X))') " conv. grad tol =",mEPSGrad," conv. dist tol =",mEPSDist
         write(iu,'(1X,A21,G12.2,3X,A21,F12.4)') " conv. func tol =",mEPSFctn," H_dist (A) =",mHDist*bohr2angs
         write(iu,*)
      else
         write(iu,'(1X,A21,1X,A5,9X,A21,I5)') 'mode =',s
         write(iu,'(1X,A21,I6,9X,A21,F10.4)') 'max_iter =',minimization_getIterMax(), &
         ' step_max =',minimization_getStepMax()
         if (mConstraint) then
            write(iu,'(1X,A)') 'maximization with 1s core electrons fixed at nucleus'
         else
            write(iu,'(1X,A)') 'maximization for all electrons'
         end if

         select case (mOptMode)
         case (BFGS)
            write(iu,'(2(1X,A21,G12.2,3X))') " conv. func tol =",mTol
         case (NELM)
            write(iu,'(2(1X,A21,G12.2,3X))') " conv. func tol =",mTol," initial delta =",mDelta
         case default
            call abortp('psimax_writeParams: internal error, illegal value')
         end select

         write(iu,*)
      endif
   end subroutine psimax_writeParams


   subroutine psimax_opt(rwap)
      type(RandomWalker), pointer :: rwap(:)
      type(RandomWalker), pointer :: rwp
      real*8 :: F,F0
      real*8 :: x(ne),y(ne),z(ne)
      integer:: w,i,iflag,v,nc
      integer:: idx(ne)
      real*8 :: eps(3)
      logical lold
      type(Coord) :: nall

      v = mVerbose

      if (mAllowSingularities) then
         lold = eloc_getStopAtSingularity()
         call eloc_setStopAtSingularity(.false.)
      end if

      do w=1,size(rwap)

         idx = (/ (i,i=1,ne) /)

         select case (mOptMode)
         case (LBFGS)
            eps(1) = mEPSGrad; eps(2) = mEPSDist; eps(3) = mEPSFctn
            call pos(rwap(w),x,y,z)

            if (mVerbose>=4) then
               write(iull,*) 'before nextpsimax:'
               do i=1,ne
                  write(iull,'(i5,3f15.6)') i,x(i),y(i),z(i)
               end do
            end if

            call nextpsimax(x,y,z,mMaxIter,mHDist,eps,v,F,nc,iflag,idx)

            if (mVerbose>=4) then
               write(iull,*) 'after nextpsimax:'
               write(iull,'(50i3)') (/ (i,i=1,ne) /)
               write(iull,'(50i3)') idx(1:ne)
               do i=1,ne
                  write(iull,'(i5,3f15.6)') i,x(i),y(i),z(i)
               end do
               write(iull,*) ' ---------------'
            end if
         case (BFGS)
            allocate(nall%x(ne),nall%y(ne),nall%z(ne))
            call pos(rwap(w),x,y,z)
            if (mConstraint) call findcores_sort(nall)
            call psimax_calc(x,y,z,F,iflag)
            deallocate(nall%x,nall%y,nall%z)
         case (NELM)
            allocate(nall%x(ne),nall%y(ne),nall%z(ne))
            call pos(rwap(w),x,y,z)
            if (mConstraint) call findcores_sort(nall)
            call psimax_calc(x,y,z,F,iflag)
            deallocate(nall%x,nall%y,nall%z)
         case default
            call abortp('(psimax_opt): internal error: illegal mOptMode')
         end select

         sCounter(8) = sCounter(8) + 1.d0
         sCounter(iflag) = sCounter(iflag) + 1.d0
         if (iflag == 0) sCounter(9) = sCounter(9) + 1.d0

         select case (mAnalyseMode)
         case (MAXIMA)
            call maxana_add(x,y,z,F,iflag)
         case (BASINS)
            rwp => rwap(w)
            call maxbas_add(x,y,z,F,rwp,idx,iflag)
         end select

      enddo

      if (mAllowSingularities) then
         call eloc_setStopAtSingularity(lold)
      endif
   end subroutine psimax_opt



   subroutine psimax_writeResults()
      integer, parameter :: iu=12
      real*8          :: min_save,rcvCount(0:9)
      integer         :: i,k
      real*8          :: F

      call myMPIReduceSumDouble(sCounter,rcvCount,10)

      if (MASTER) then

         write(iul,*)
         write(iul,*) "Summary for maxima search:"
         write(iul,*)


         write(iul,'(a,f15.0)') "  # minimizer calls:",rcvCount(8)
         write(iul,'(a,f15.0)') "  # maxima analyzed:",rcvCount(9)
         write(iul,'(a)') "  exit code percentages of minimizer:"
         select case (mOptMode)
         case (LBFGS)
            write(iul,'(a30,f10.2)') " fully converged: ",100*rcvCount(0)/rcvCount(8)
            write(iul,'(a30,f10.2)') " fctn value not converged: ",100*rcvCount(1)/rcvCount(8)
            write(iul,'(a30,f10.2)') " gradient not converged: ",100*rcvCount(2)/rcvCount(8)
            write(iul,'(a30,f10.2)') " distance not converged: ",100*rcvCount(3)/rcvCount(8)
            write(iul,'(a30,f10.2)') " boundary hit: ",100*rcvCount(4)/rcvCount(8)
            write(iul,'(a30,f10.2)') " error occurred: ",100*rcvCount(5)/rcvCount(8)
            write(iul,'(a30,f10.2)') " eloc singularity: ",100*rcvCount(6)/rcvCount(8)
            write(iul,'(a30,f10.2)') " fctn value is nan: ",100*rcvCount(7)/rcvCount(8)
         case (BFGS)
            write(iul,'(a30,f10.2)') " fctn value converged: ",100*rcvCount(0)/rcvCount(8)
            write(iul,'(a30,f10.2)') " gradient converged: ",100*rcvCount(1)/rcvCount(8)
            write(iul,'(a30,f10.2)') " fctn value not a number: ",100*rcvCount(2)/rcvCount(8)
            write(iul,'(a30,f10.2)') " not converged: ",100*rcvCount(3)/rcvCount(8)
         case (NELM)
            write(iul,'(a30,f10.2)') " fctn value converged: ",100*rcvCount(0)/rcvCount(8)
            write(iul,'(a30,f10.2)') " not converged: ",100*rcvCount(1)/rcvCount(8)
         end select
         write(iul,*)
      endif

      select case (mAnalyseMode)
      case (MAXIMA)
         call maxana_writeResults()
      case (BASINS)
         call maxbas_writeResults()
      end select

   end subroutine psimax_writeResults


   !!! old code  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine psimax_calc(x,y,z,f,iflag)
       ! performs a local minimiation of -ln(Psi**2), Psi unnormalized
       real*8,intent(inout) :: x(ne),y(ne),z(ne)   ! in: initial coords
                                                   ! out: optimized coords
       real*8,intent(out)   :: f                   ! resulting function value
       integer              :: iflag               ! convergence information
       real*8               :: vec(3*ne-6*mNocel),g(3*ne-6*mNocel)
       real*8               :: func_max
       external func_max

       vec  = 0d0

       call writeCoordsToVec(vec,x,y,z)
       select case (mOptMode)
          !case (RPSO)
          !  call initSwarm(mBound,1,n,mTol,1)
          !  call swarmOpt(vec)
          case (BFGS)
             call dfpmin(vec,3*ne-6*mNocel,mTol,func_max,mVerbose,iflag)
             f=func_max(vec,3*ne-6*mNocel)
          case (NELM)
             call psimax_nelm_driver(vec,f,iflag)
          case default
             call abortp("psimax_calc: illegal mOptmode")
       end select

       call readCoordsFromVec(vec,x,y,z)
   end subroutine psimax_calc


   subroutine writeCoordsToVec(vec,x,y,z)
       real*8,intent(in)  :: x(:),y(:),z(:)
       real*8,intent(out) :: vec(:)
       integer            :: istat,n
       integer            :: i,j
       n = size(x)

       vec = 0.d0

       i = 1
       j = nalpha-mNocel
       vec(1:j)  = x(1:nalpha-mNocel)

       i = nalpha-mNocel+1
       j = ne - 2*mNocel
       vec(i:j) = x(nalpha+1:ne-mNocel)

       i =  ne-2*mNocel+1
       j =  ne-2*mNocel+nalpha-mNocel
       vec(i:j) = y(1:nalpha-mNocel)

       i= ne-2*mNocel+nalpha-mNocel+1
       j= 2*(ne-2*mNocel)
       vec(i:j) = y(nalpha+1:ne-mNocel)

       i = 2*(ne-2*mNocel)+1
       j = 2*(ne-2*mNocel)+nalpha-mNocel
       vec(i:j) = z(1:nalpha-mNocel)

       i = 2*(ne-2*mNocel)+nalpha-mNocel+1
       j = 3*(ne-2*mNocel)
       vec(i:j) = z(nalpha+1:ne-mNocel)


   end subroutine writeCoordsToVec


   subroutine readCoordsFromVec(vec,x,y,z)
       real*8,intent(in)               :: vec(:)
       real*8,intent(out)              :: x(:),y(:),z(:)
       integer                         :: istat,i,j,n

       n = size(vec) / 3

       x = 0.d0; y=0.d0; z=0.d0

       i = 1
       j = nalpha-mNocel
       x(1:nalpha-mNocel) = vec(i:j)

       i = nalpha-mNocel+1
       j = ne-2*mNocel
       x(nalpha+1:ne-mNocel) = vec(i:j)

       i = ne-2*mNocel+1
       j = ne-2*mNocel+nalpha-mNocel
       y(1:nalpha-mNocel) = vec(i:j)

       i = ne-2*mNocel+nalpha-mNocel+1
       j = 2*(ne-2*mNocel)
       y(nalpha+1:ne-mNocel) = vec(i:j)

       i = 2*(ne-2*mNocel)+1
       j = 2*(ne-2*mNocel)+nalpha-mNocel
       z(1:nalpha-mNocel) = vec(i:j)

       i = 2*(ne-2*mNocel)+nalpha-mNocel+1
       j = 3*(ne-2*mNocel)
       z(nalpha+1:ne-mNocel) = vec(i:j)

       if(mConstraint) then
           x(nalpha-mNocel+1:nalpha) = mCoreCoords%x(1:mNocel)
           y(nalpha-mNocel+1:nalpha) = mCoreCoords%y(1:mNocel)
           z(nalpha-mNocel+1:nalpha) = mCoreCoords%z(1:mNocel)

           x(ne-mNocel+1:ne) = mCoreCoords%x(1:mNocel)
           y(ne-mNocel+1:ne) = mCoreCoords%y(1:mNocel)
           z(ne-mNocel+1:ne) = mCoreCoords%z(1:mNocel)
       endif
   end subroutine readCoordsFromVec

   subroutine rebuildvec(vec_core,vec)
       real*8,intent(in) :: vec(3*ne-mNocel)
       real*8,intent(out) :: vec_core(3*ne)
       integer           :: i,j

       vec_core(1:nalpha-mNocel) = vec(1:nalpha-mNocel)
       vec_core(nalpha+1:ne-mNocel) = vec(nalpha-mNocel+1:ne-2*mNocel)
       vec_core(ne+1:ne+nalpha-mNocel) = vec(ne-2*mNocel+1:ne-2*mNocel+nalpha-mNocel)
       vec_core(ne+nalpha+1:2*ne-mNocel) = vec(ne-2*mNocel+nalpha-mNocel+1:2*(ne-2*mNocel))
       vec_core(2*ne+1:2*ne+nalpha-mNocel) =  vec(2*(ne-2*mNocel)+1:2*(ne-2*mNocel)+nalpha-mNocel)
       vec_core(2*ne+nalpha+1:3*ne-mNocel) = vec(2*(ne-2*mNocel)+nalpha-mNocel+1:3*(ne-2*mNocel))
       if(mConstraint)then
           vec_core(nalpha-mNocel+1:nalpha) = mCoreCoords%x(1:mNocel)
           vec_core(ne-mNocel+1:ne) = mCoreCoords%x(1:mNocel)
           vec_core(ne+nalpha-mNocel+1:ne+nalpha) = mCoreCoords%y(1:mNocel)
           vec_core(2*ne-mNocel+1:2*ne) =  mCoreCoords%y(1:mNocel)
           vec_core(2*ne+nalpha-mNocel+1:2*ne+nalpha) = mCoreCoords%z(1:mNocel)
           vec_core(3*ne-mNocel+1:3*ne) =  mCoreCoords%z(1:mNocel)
       endif
   end subroutine rebuildvec

   subroutine psimax_set_core_el_coords()
       integer     :: c,i


       c = 0
       if(mCoreRule == 1) then
           do i=1,ncenter
               if(atoms(i)%elemIdx > 2) then
                   c = c +1
               endif
           enddo
       elseif(mCoreRule == 3) then
           do i=1,ncenter
               c = c +1
           enddo
       endif

       allocate(mCoreCoords%x(c),mCoreCoords%z(c),mCoreCoords%y(c))

       c = 0
       if(mCoreRule == 1) then
           do i=1,ncenter
               if(atoms(i)%elemIdx > 2) then
                   c = c +1
                   mCoreCoords%x(c) = atoms(i)%cx
                   mCoreCoords%y(c) = atoms(i)%cy
                   mCoreCoords%z(c) = atoms(i)%cz
               endif
           enddo
       elseif(mCoreRule ==3) then
           do i=1,ncenter
               c = c +1
               mCoreCoords%x(c) = atoms(i)%cx
               mCoreCoords%y(c) = atoms(i)%cy
               mCoreCoords%z(c) = atoms(i)%cz
           enddo

       endif

   end subroutine psimax_set_core_el_coords


   subroutine psimax_nelm_driver(vec,F,iflag)
       real*8,intent(inout) :: vec(:)
       real*8,intent(out)   :: F
       integer, intent(out) :: iflag   ! convergence information (0=converged, 1=unconverged)
       integer              :: k,N,iter
       real*8               :: func_max,Ftmp,delta
       external func_max
       delta = mDelta
       N =  3*ne-6*mNocel
       F=func_max(vec,N)
       Ftmp =1d30
       iflag = 0
       if(mNelmIt)then
           do k=1,100
               iter = mMaxIter
               call nelmead(vec,N,iter,mDelta,mTol,F,func_max,iul,logmode)
               delta = delta /(2d0*dble(k))
               if(DABS(F-Ftmp)<=mTol) then
                   F = Ftmp
                   exit
               endif
               Ftmp = F
               if(k == 100) iflag = 1
           enddo
       else
           iter = mMaxIter
           call nelmead(vec,N,iter,mDelta,mTol,F,func_max,iul,logmode)
           if (iter>=mMaxIter) iflag = 1
       endif

   end subroutine psimax_nelm_driver



end module psimax


real*8 function func_max(vec,n)
   use psimax, only: rebuildvec
   use elocaldata
   use elocal, only: eloc
   use eConfigsModule
   implicit none
   integer,intent(in)   :: n
   real*8,intent(in)    :: vec(n)
   real*8               :: vec_core(3*ne)
   integer              :: i
   type(eConfigArray)   :: ec

   call rebuildvec(vec_core,vec)
   call eConfigArray_new(ec,ne,1)
   call eConfigArray_set(ec,1,vec_core(1:ne),vec_core(ne+1:2*ne),vec_core(2*ne+1:3*ne))
   call eloc(0,ec,'none')
   call eConfigArray_destroy(ec)
   func_max = -2d0*(log(elPhi(1))+elU(1))
end function func_max

subroutine grad_max(vec,grad,n)
   use psimax, only: rebuildvec, writeCoordsToVec
   use elocaldata
   use elocal, only: eloc
   use eConfigsModule
   implicit none
   integer,intent(in)   :: n
   real*8,intent(in)    :: vec(n)
   real*8,intent(out)   :: grad(n)
   real*8               :: vec_core(3*ne)
   type(eConfigArray)   :: ec
   call rebuildvec(vec_core,vec)
   call eConfigArray_new(ec,ne,1)
   call eConfigArray_set(ec,1,vec_core(1:ne),vec_core(ne+1:2*ne),vec_core(2*ne+1:3*ne))
   call eloc(0,ec,'none')

   call eConfigArray_destroy(ec)
   call writeCoordsToVec(grad(1:n),elxDrift(1:ne,1),elyDrift(1:ne,1),elzDrift(1:ne,1))

   grad(1:n) = -2.d0 * grad(1:n)
end subroutine grad_max





