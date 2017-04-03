module assignment
    use findcores
    use global
    use RandomWalkerModule
    use psimax
    use wfdata
    use statistics
    use hungarian
    use Utils, only: tokenize,intToStr
    implicit none

    private
    public:: assignment_init,assignment_get,assignment_set_reference,assignment_read_reference,assignment_get_reference
    public:: assignment_assignstat_write


    integer,allocatable      ::    mAssign(:)          ! Contains assignment vector + one last element for mCoreRule
    type(coord)              ::    nall                ! nall: contains electron configuration

    !Read Reference
    type(coord),allocatable  ::    mnall(:)            ! mae: alpha electrons of maxima, mbe: beta elect.
    integer,allocatable      ::    mRef(:,:)           ! Array of reference indices
    integer                  ::    mNumRefs            ! Number of references found
    integer, parameter       ::    mMaxRefs = 100      ! Maximum number of references
    integer,allocatable      ::    mUsedRefs(:)        ! Used reference for electron n
    integer                  ::    mnoRef              ! number of references to read resp. number of references read
    logical                  ::    mReOpt=.false.

    !Separation of Core Electrons
    integer                  ::    mCoreEl   = 0       ! HALF the number of core electrons

    !Distance Criterion for n references
    real,allocatable         ::    mDmatA(:,:,:)         !Distance Matrix for alpha electrons
    real,allocatable         ::    mDmatB(:,:,:)         !Distance Matrix for beta  electrons


    !Parameters for Munkres
    logical                  ::    mMunkres=.true.    ! Use Munkres for assignment
    integer                  ::    mMode=1            ! Mode for Munkres algorithm, s. subroutine munkres

    !Greedy
    logical                  ::    mGreedy=.false.    ! Do assignment with simple greedy algorithm

    !Assignment Statistics
    logical                      ::    mAssignStat = .false. ! Do statistic over the distance between reference point and assigned electron
    type(simpleStat),allocatable ::    mAsgnStat(:)

    ! Assignment cost matrix type
    integer                      :: mCostMatrix = 1   ! 1: squared distance, 2: distance, 3: meanDist (Rene)

contains

   !-----------------------------------
   subroutine assignment_init(lines,nl)
   !-----------------------------------
      integer, intent(in)            :: nl
      character(len=120), intent(in) :: lines(:)
      integer                        :: iflag,alstat
      logical                        :: finda
      integer                        :: i
      character(len=7)               :: s         !longest value for s is 'squared' -> len=7
      integer                        :: tmp,tmp2  !need for sorting references

      if (.not.finda(lines,nl,'nocoresep'))then
        call findcores_find(mCoreEl)
      endif

      call getstra(lines,nl,'assign=',s,iflag)
      if (iflag == 0) then
         if (s=='squared') then
            mCostMatrix = 1
         else if (s=='simple') then
            mCostMatrix = 2
         else if (s=='mean') then
            mCostMatrix = 3
            mAssignStat = .true.
         else
            call abortp("illegal value given for assign=")
         end if
      else
         mCostMatrix = 1
      end if

      ! Reoptimize reference read from file
      if (finda(lines,nl,'reopt')) mReOpt = .true.

      ! search for add. references
      do i = 2, mMaxRefs
         if(.not. finda(lines,nl,'ref'// trim(intToStr(i)) //'_nr=')) then
            mNumRefs = i - 1
            exit
         endif
      enddo

      ! index of reference to read from file - reading n ref numbers
      allocate(mRef(mNumRefs,2),stat=alstat)
      call assert(alstat==0,"assignment::init: allocate 0 failed")

      ! no reference given -> default reference 1 1 to be used
      mRef(1,1) = 1
      mRef(1,2) = 1

      call getinta(lines,nl,'ref_nr=',mRef(1,1),iflag)
      call getinta(lines,nl,'ref_nr2=',mRef(1,2),iflag)

      ! read references
      do i = 2, mNumRefs
         ! default value for index
         tmp2 = 1
         call getinta(lines,nl,'ref'// trim(intToStr(i)) //'_nr=',mRef(i,1),iflag)
         call getinta(lines,nl,'ref'// trim(intToStr(i)) //'_nr2=',mRef(i,2),iflag)
      enddo

      ! sort rerences - first one to find: mRef(1)
      if (mNumRefs > 1) call sort_reference ()

      ! munkres is default
      if (finda(lines,nl,'greedy')) then
         mGreedy = .true.
         mMunkres = .false.
      else
         mGreedy = .false.
         mMunkres = .true.
      endif
      ! allocate and reset everything
      if (finda(lines,nl,'asgnstat')) mAssignStat = .true.
      allocate(mAssign(ne+1),mUsedRefs(ne),stat=alstat)
      call assert(alstat==0,"assignment::init: allocate 1 failed")
      allocate(nall%x(ne),nall%y(ne),nall%z(ne),stat=alstat)
      call assert(alstat==0,"assignment::init: allocate 2 failed")

      allocate(mDmatA(mNumRefs,nalpha-mCoreEl,nalpha-mCoreEl), &
               mDmatB(mNumRefs,nbeta-mCoreEl,nbeta-mCoreEl),stat=alstat)
      call assert(alstat==0,"assignment::init: allocate 3 failed")

      allocate(mnall(mNumRefs), stat=alstat)
      call assert(alstat==0,"assignment::init: allocate 4 failed")
      do i = 1, mNumRefs
         allocate(mnall(i)%x(ne),mnall(i)%y(ne),mnall(i)%z(ne),stat=alstat)
         call assert(alstat==0,"assignment::init: allocate 5 failed")
      enddo

      if (mAssignStat) then
         allocate(mAsgnStat(ne))
         do i=1,ne
             call reset(mAsgnStat(i))
         enddo
      endif


      contains

         subroutine sort_reference()
           ! subroutine sorting mRef w.r.t. reference indices
           ! algorithm: bubble sort
           integer   :: n,m,tmp
           logical   :: swapped

           n = mNumRefs
           do
              swapped = .false.
              do m=1, (n-1)
                 if ( mRef(m,1) > mRef(m+1,1) .or. &
                    ( mRef(m,1) == mRef(m+1,1) .and. &
                      mRef(m,2) > mRef(m+1,2) ) ) then
                    tmp = mRef(m,1)
                    mRef(m,1) = mRef(m+1,1)
                    mRef(m+1,1)=tmp
                    tmp = mRef(m,2)
                    mRef(m,2) = mRef(m+1,2)
                    mRef(m+1,2)=tmp
                    swapped = .true.
                 endif
              end do
              n=n-1
              if (.not. swapped) exit
           end do

         end  subroutine


   end subroutine assignment_init

    !---------------------------------------------
    subroutine assignment_read_reference(lines,nl)
    !---------------------------------------------
        !Reads Reference from ref-File
        integer, intent(in)           :: nl
        character(len=120), intent(in) :: lines(:)

        integer, parameter            :: iomax=45
        integer                       :: i,j,r,iflag,io,n,alstat
        real*8                        :: x_tmp(ne),y_tmp(ne),z_tmp(ne)
        real*8                        :: F(mNumRefs),Ftmp
        integer                       :: nomax,noel
        character(len=120)            :: line
        character(len=20)             :: words(12)
        integer                       :: refNo,refNo2,foundtmp,found(mNumRefs),foundrefs
        character(len=80)             :: line2
        character(len=40)             :: fname
        integer                       :: asgn_va(nalpha), asgn_vb(nbeta)
        real                          :: temp
        integer                       :: tindx1, tindx2             !temp. indices used for Dmat - calculation
        real,allocatable              ::    mDistA(:,:)             !Distance Matrix for alpha references
        real,allocatable              ::    mDistB(:,:)             !Distance Matrix for beta references

        ! allocate distance-matrices for calculation
        allocate(mDistA(nalpha,nalpha),mDistB(nbeta,nbeta),stat=alstat)
        call assert(alstat==0,"assignment::read_reference: allocate 1 failed")

        call getstra(lines,nl,'ref_file=',fname,iflag)
        if (iflag /= 0) then
            fname = trim(baseName)//'.ref'
        end if
        open(iomax,file=fname,status='old',iostat=io)
        call assert(io==0,' ref file '//trim(fname)//' does not exist')

        read(iomax,*) n   ! ncenter
        do i=1,n
            read(iomax,*) j  ! index,elem,x,y,z
        end do

        foundrefs = 0
        read(iomax,*) nomax
        do  !FIXME loop does not end, if worng refernece is given!
            read(iomax,'(A)',iostat=io) line
            call assert(io==0,' error while reading reference file')
            call tokenize(line,words,n)
            read(words(2),*) refNo
            read(words(3),*) refNo2
            read(words(5),*) Ftmp
            read(words(7),*) foundtmp
            read(iomax,*) noel
            do i=1,noel
                read(iomax,*) x_tmp(i),y_tmp(i),z_tmp(i)
            enddo

            do i = foundrefs + 1, mNumRefs ! XXX sort necessary? i = 1, ...
               if(mRef(i,1) == refNo .and. mRef(i,2) == refNo2) then
                  mnall(i)%x = x_tmp
                  mnall(i)%y = y_tmp
                  mnall(i)%z = z_tmp
                  F(i) = Ftmp
                  found(i) = foundtmp
                  foundrefs = foundrefs + 1
                  exit
               endif
            enddo
            if(foundrefs == mNumRefs) exit
        enddo

        if (mReOpt) then
           call psimax_init(lines,nl)
           do i = 1, mNumRefs
              call psimax_calc(mnall(i)%x(1:ne),mnall(i)%y(1:ne),mnall(i)%z(1:ne),F(i),iflag)
              if (MASTER .and. logmode >= 2) write(iul,'(A,G13.6)') "Reoptimized F = ",F(i)
           enddo
        endif

        if (MASTER .and. logmode >=2) then
            write(iul,'(3A)') "reference file ",trim(fname)," read!"
            do i = 1, mNumRefs
                write(iul,'(A,I3,A,2I3,A,F13.5,A,I8)') "reference ", i, " indices: ", &
                  mRef(i,1),mRef(i,2)," with function value: ",F(i)," found:",found(i)
                do j=1,ne
                   write(iul,'(3F12.5)') mnall(i)%x(j), mnall(i)%y(j), mnall(i)%z(j)
                enddo
            enddo
        endif
        close(iomax)

        do r = 2, mNumRefs
          !call munkres algorithm for assignment of ref i -> ref 1
          !calculate distance matrix: ref 1 - ref i
          !alpha
          do i=1,nalpha
             do j=1,nalpha
                mDistA(i,j) = sqrt( (mnall(1)%x(i)-mnall(r)%x(j))**2 + &
                                    (mnall(1)%y(i)-mnall(r)%y(j))**2 + &
                                    (mnall(1)%z(i)-mnall(r)%z(j))**2 )
             enddo
          enddo
          !beta
          do i=(nalpha+1),ne
             do j=(nalpha+1),ne
                tindx1 = i-nalpha
                tindx2 = j-nalpha
                mDistB(tindx1,tindx2) = sqrt( (mnall(1)%x(i)-mnall(r)%x(j))**2 + &
                                              (mnall(1)%y(i)-mnall(r)%y(j))**2 + &
                                              (mnall(1)%z(i)-mnall(r)%z(j))**2 )
             enddo
          enddo

          !call munkres to get correct assignment and permute

          !alpha
          call munkres(mMode, mDistA, nalpha, nalpha, asgn_va, temp)
          !permute reference
          ! save reference r (which will be sorted) to temp var - only alpha maxima
          do i=1,nalpha
             x_tmp(i) = mnall(r)%x(i)
             y_tmp(i) = mnall(r)%y(i)
             z_tmp(i) = mnall(r)%z(i)
          enddo
          do i=1,nalpha
             mnall(r)%x(i) = x_tmp(asgn_va(i))
             mnall(r)%y(i) = y_tmp(asgn_va(i))
             mnall(r)%z(i) = z_tmp(asgn_va(i))
          enddo

          !beta
          call munkres(mMode, mDistB, nbeta, nbeta, asgn_vb, temp)
          !permute reference
          ! save reference r - only beta maxima
          do i=nalpha+1,ne
             x_tmp(i) = mnall(r)%x(i)
             y_tmp(i) = mnall(r)%y(i)
             z_tmp(i) = mnall(r)%z(i)
          enddo
          do i=nalpha+1,ne
             mnall(r)%x(i) = x_tmp(asgn_vb(i-nalpha)+nalpha)
             mnall(r)%y(i) = y_tmp(asgn_vb(i-nalpha)+nalpha)
             mnall(r)%z(i) = z_tmp(asgn_vb(i-nalpha)+nalpha)
          enddo
        enddo

        !print sorted references
        if (MASTER .and. logmode >=2 .and. mNumRefs > 1) then
          write(iul,'(A)') "reference sorted using munkres-algorithm"
          do r = 1, mNumRefs
             write(iul,'(A,I3,A,2I3)') "reference", r, " indices: ",mRef(r,1),mRef(r,2)
             do i=1,ne
                 write(iul,'(3F12.5)') mnall(r)%x(i),mnall(r)%y(i),mnall(r)%z(i)
             enddo
             write(iul,*)
          enddo
        endif

        if(mCoreEl>0) then
          do r = 1, mNumRefs
             call findcores_sort(mnall(r))
          enddo
        endif

        !dellocate Distance matrices
        deallocate(mDistA,mDistB)

    end subroutine assignment_read_reference

    subroutine assignment_get_reference(ref,x,y,z)
        integer,intent(in) :: ref
        real*8,intent(out) :: x(ne),y(ne),z(ne)

        x = mnall(ref)%x
        y = mnall(ref)%y
        z = mnall(ref)%z
    end subroutine assignment_get_reference

    !-----------------------------------------------
    subroutine assignment_set_reference(ref,x,y,z,print)
        !-------------------------------------------
        !Set Reference from previous max/coc calculation
        integer,intent(in) :: ref
        real*8,intent(in)  :: x(:),y(:),z(:)
        integer :: i,idx
        logical,intent(in) :: print

        mnall(ref)%x = x
        mnall(ref)%y = y
        mnall(ref)%z = z

        if (mCoreEl>0) call findcores_sort(mnall(1))
        if (MASTER .and. logmode >= 2)then
            write(iul,'(A)') " new assignment reference:"
            do i=1,ne
                write(iul,'(3F12.5)') mnall(1)%x(i),mnall(1)%y(i),mnall(1)%z(i)
            end do
            write(iul,*)
        end if
    end subroutine assignment_set_reference

    !------------------------------------
    subroutine assignment_get(rwp,asgn,bl)
    !------------------------------------
         !Get the asignment
         type(RandomWalker),pointer,intent(in) :: rwp(:)
         integer,intent(inout)                 :: asgn(:,:)
         integer                               :: aidx_a(nalpha-mCoreEl),aidx_b(nbeta-mCoreEl)
         integer                               :: aidx_a2(nalpha-mCoreEl),aidx_b2(nbeta-mCoreEl)
         integer                               :: aidx_save(1,ne)
         real                                  :: dist1,dist2
         integer                               :: i,k,w
         integer                               :: no
         integer                               :: atmp(ne+1)
         integer                               :: usedRefA, usedRefB
         logical                               :: bl
         real*8                                :: dummy

         call assert(associated(rwp),'assignment_get:: rw pointer not associated')
         do w = 1,size(rwp)
            usedRefA = 1
            usedRefB = 1
            !Get Electron Configuration from random walker
            call pos(rwp(w),nall%x,nall%y,nall%z)
            !Push the alpha/beta core electrons to the end of the alpha, beta lists
            if (mCoreEl>0) then
                call findcores_sort(nall)
            end if
            ! get the dist. matrix for alpha and beta
            select case (mCostMatrix)
               case (1)
                  call assignment_getSimpleDistMat()
               case (2)
                  call assignment_getAbsDistMat()
               case (3)
                  if(bl) then
                     call assignment_getMeanDistMat()
                  else
                     call assignment_getSimpleDistMat()
                  endif
               case default
                  call abortp("assignment_get: illegal value for mCostMatrix")
            end select

            if(mMunkres)then
               !assign alpha
               no   = nalpha-mCoreEl
               call munkres(mMode, mDmatA(1,:,:), no, no, aidx_a, dist1)
               do i = 2, mNumRefs
                  call munkres(mMode, mDmatA(i,:,:), no, no, aidx_a2, dist2)
                  if (dist1>dist2) then  !reference i: lower distance
                     usedRefA = i
                     dist1 = dist2
                     aidx_a=aidx_a2
                  endif
               enddo

               !assign beta
               no   = nbeta-mCoreEl
               call munkres(mMode, mDmatB(1,:,:), no, no, aidx_b, dist1)
               do i = 2, mNumRefs
                  call munkres(mMode, mDmatB(i,:,:), no, no, aidx_b2, dist2)
                  if (dist1>dist2) then  !reference i: lower distance
                     usedRefB = i
                     dist1 = dist2
                     aidx_b=aidx_b2
                  endif
               enddo

               !write the assignment vector for the whole electron configuration
               mAssign(1:nalpha-mCoreEl) = aidx_a
               mAssign(nalpha+1:ne-mCoreEl) = aidx_b + nalpha

               if(mCoreEl>0)then
                 do i=nalpha-mCoreEl+1,nalpha
                     mAssign(i) = i
                 enddo
                 do i=ne-mCoreEl+1,ne
                     mAssign(i) = i
                 enddo
               endif
               atmp(1:ne) = mAssign(1:ne)
            endif
            if(mGreedy)then
                call greedy(dummy)
            endif

            !Last Entry is the number of Core Electrons, nec. for other routines
            mAssign(size(mAssign))  = mCoreEl
            if(mAssignStat .and. mMunkres) then
               call assignment_assignstat_add(usedRefA,usedRefB)
            endif
            asgn(:,w) = mAssign
        enddo
    end subroutine assignment_get

    subroutine assignment_assignstat_add(refA, refB)
        integer,intent(in) :: refA, refB
        integer :: i,idx,r

        do i=1,nalpha-mCoreEl
            idx=mAssign(i)
            call addData(mAsgnStat(i),dble(mDmatA(refA,idx,i)))
        enddo

        do i=1,nbeta-mCoreEl
            idx=mAssign(nalpha+i)
            call addData(mAsgnStat(nalpha+i),dble(mDmatB(refB,idx-nalpha,i)))
        enddo
    end subroutine assignment_assignstat_add

    subroutine assignment_assignstat_write()
        integer :: i

        if(mMunkres .and. mAssignStat)then

            write(iul,'(A)') ''
            write(iul,'(A)') 'Assignment Statistics for Munkres Algorithm'
            write(iul,'(A)') 'Gives the mean distance for electrons assigned to a reference'
            write(iul,'(A)') 'IDX    Mean Distance    StdDev '
            do i=1,ne
                write(iul,'(I3,X,F13.5,X,F13.5)') i,meanAllNodes(mAsgnStat(i)),stdDevMeanAllNodes(mAsgnStat(i))
            enddo
        endif
    end subroutine assignment_assignstat_write


    ! Routines for Distance Matrix

    !---------------------------------------
    subroutine assignment_getSimpleDistMat()
        !---------------------------------------
        !Calculate  simple quadratic euclidian distance
        integer                :: i,j,r
        integer                :: limit,limit2,idx,idx2

        limit = nalpha-mCoreEl
        do i=1,limit
            do j=1,limit
               do r = 1, mNumRefs
                  mDmatA(r,i,j)=(mnall(r)%x(i)-nall%x(j))**2 + &
                                (mnall(r)%y(i)-nall%y(j))**2 + &
                                (mnall(r)%z(i)-nall%z(j))**2
               enddo
            enddo
        enddo
        limit = nalpha+1
        limit2 = ne-mCoreEl
        do i=limit,limit2
            do j=limit,limit2
               idx = i - nalpha
               idx2 = j - nalpha

               do r = 1, mNumRefs
                  mDmatB(r,idx,idx2)=(mnall(r)%x(i)-nall%x(j))**2 + &
                                     (mnall(r)%y(i)-nall%y(j))**2 + &
                                     (mnall(r)%z(i)-nall%z(j))**2
               enddo
            enddo
        enddo
    end subroutine assignment_getSimpleDistMat

    !------------------------------------
    subroutine assignment_getAbsDistMat()
        !--------------------------------
        !Calculate simple distance matrix
        integer                :: i,j,r
        integer                :: limit,limit2,idx,idx2

        limit = nalpha-mCoreEl
        do i=1,limit
            do j=1,limit
               do r = 1, mNumRefs
                  mDmatA(r,i,j)=sqrt((mnall(r)%x(i)-nall%x(j))**2 + &
                                     (mnall(r)%y(i)-nall%y(j))**2 + &
                                     (mnall(r)%z(i)-nall%z(j))**2)
               enddo
            enddo
        enddo
        limit = nalpha+1
        limit2 = ne-mCoreEl
        do i=limit,limit2
            do j=limit,limit2
               idx = i - nalpha
               idx2 = j - nalpha
               do r = 1, mNumRefs
                  mDmatB(r,idx,idx2)=sqrt((mnall(r)%x(i)-nall%x(j))**2 + &
                                          (mnall(r)%y(i)-nall%y(j))**2 + &
                                          (mnall(r)%z(i)-nall%z(j))**2)
               enddo
            enddo
        enddo
    end subroutine assignment_getAbsDistMat

        !---------------------------------------
    subroutine assignment_getMeanDistMat()
        !---------------------------------------
        !weight euclidian dist by meandist

        ! not sure if correct with more than one reference
        integer                :: i,j,r
        integer                :: limit,limit2,idx,idx2

        limit = nalpha-mCoreEl
        do i=1,limit
            do j=1,limit
               do r = 1, mNumRefs
                  mDmatA(r,i,j)=(mnall(r)%x(i) - nall%x(j))**2 + &
                                (mnall(r)%y(i) - nall%y(j))**2 + &
                                (mnall(r)%z(i) - nall%z(j))**2
                  mDmatA(r,i,j)=mDmatA(r,i,j)/meanAllNodes(mAsgnStat(j))
               enddo
            enddo
        enddo
        limit = nalpha+1
        limit2 = ne-mCoreEl
        do i=limit,limit2
            do j=limit,limit2
                idx = i - nalpha
                idx2 = j - nalpha
                do r = 1, mNumRefs
                   mDmatB(r,idx,idx2)=(mnall(r)%x(i) - nall%x(j))**2 + &
                                      (mnall(r)%y(i) - nall%y(j))**2 + &
                                      (mnall(r)%z(i) - nall%z(j))**2
                   mDmatB(r,i,j)=mDmatB(r,i,j)/meanAllNodes(mAsgnStat(j))
                enddo
            enddo
        enddo
    end subroutine assignment_getMeanDistMat

    ! Routines for Assignment

    subroutine greedy(dist)
        real*8,intent(out) :: dist
        real*8             :: minDist
        integer            :: idx1(2)
        integer            :: asgnA(mNumRefs,nalpha-mCoreEl)
        integer            :: asgnB(mNumRefs,nbeta-mCoreEl)
        integer            :: i,minRef,r
        integer            :: asgnA_test(nalpha-mCoreEl)
        integer            :: asgnB_test(nbeta-mCoreEl)
        logical            :: comp

        do r = 1, mNumRefs
           dist = 0.0
           do i=1,size(mDmatA(r,:,1))
               idx1 = minloc(mDmatA(r,:,:),MASK = mDmatA(r,:,:) .ge. 0)
               dist = dist + mDmatA(r,idx1(1),idx1(2))
               mDmatA(r,:,idx1(2)) = -1
               mDmatA(r,idx1(1),:) = -1
               asgnA(r,idx1(2)) = idx1(1)
           enddo

           do i=1,size(mDmatB(r,:,1))
               idx1 = minloc(mDmatB(r,:,:),MASK = mDmatB(r,:,:) .ge. 0)
               dist = dist + mDmatB(r,idx1(1),idx1(2))
               mDmatB(r,:,idx1(2)) = -1
               mDmatB(r,idx1(1),:) = -1
               asgnB(r,idx1(2)) = idx1(1)
           enddo
           if(dist < minDist) minRef = r
        enddo

        mAssign(1:nalpha-mCoreEl) = asgnA(minRef,:)
        mAssign(nalpha+1:ne-mCoreEl) = asgnB(minRef,:) + nalpha
        if(mCoreEl>0)then
            do i=nalpha-mCoreEl+1,nalpha
                mAssign(i) = i
            enddo
            do i=ne-mCoreEl+1,ne
                mAssign(i) = i
            enddo
        endif

    end subroutine greedy

end module assignment
