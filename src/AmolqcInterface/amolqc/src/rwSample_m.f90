
MODULE RWSampleModule

! $Id: rwSample_m.f90,v 1.1.1.1 2007/04/25 13:42:20 luechow Exp $

! $Log: rwSample_m.f90,v $
! Revision 1.1.1.1  2007/04/25 13:42:20  luechow
! QMC program amolqc. rewritten in 2006. AL
!

  use global
  use RandomWalkerModule
  use utilsModule, only: myran
  use statistics
  use double_histogram_statistics
  use indexed_histogram_statistics
  use intListModule, ilm_isEmpty => isEmpty
  use sortingModule, only: entry, quickSortEntry, entryToVector, entryFromVector

  ! This implementation: array
  ! same order of elements as old static implementation
  ! This implementation uses copying of elements rather
  ! than copying of pointers

  implicit none

  type RWSample
     private
     integer                         :: sampleSize=0
     integer                         :: current=0
     type(RandomWalker), pointer     :: rwArray(:)=>null()
  end type RWSample

  integer :: mHistogramCounter=0

  private :: mHistogramCounter

  interface assignment(=)
     module procedure sample_assign
  end interface

CONTAINS

  !------------------------!
  subroutine initSample(s,n)
  !------------------------!
    type(RWSample), intent(inout) :: s
    integer, intent(in)           :: n        ! max size of rw array
    integer al

    call assert(nElecs()>0,"RWSample:init: number of electrons in RW module not set")
    call assert(.not.associated(s%rwArray),"RWSample:init: already initialized")
    call assert(n>0,"RWSample:init: initialization with zero size")
    allocate(s%rwArray(n),stat=al)
    if (al /= 0) call error("RWSample:init: allocation failed")
    s%sampleSize = 0
    s%current    = 0
  end subroutine initSample

  !-------------------------!
  subroutine destroySample(s)
  !-------------------------!
    type(RWSample), intent(inout) :: s
    integer i,al

    do i=s%sampleSize,1,-1
       call rw_destroy(s%rwArray(i))
    enddo
    deallocate(s%rwArray,stat=al)
    nullify(s%rwArray)
    if (al/=0) call error("RWSample:destroy: deallocation failed")
  end subroutine destroySample

  !------------------------------!
  elemental subroutine sample_assign(s2,s1)
  !------------------------------!
    type(RWSample), intent(inout) :: s2
    type(RWSample), intent(in)    :: s1
    integer                       :: i

    if(.not.associated(s1%rwArray)) then
      return
    endif

    if(.not.associated(s2%rwArray))then
        allocate(s2%rwArray(s1%sampleSize))
    endif
    do i=1,s1%sampleSize
     s2%rwArray(i) = s1%rwArray(i)
    enddo
    s2%sampleSize = s1%sampleSize
    s2%current = 0

  end subroutine sample_assign

!--------------- elementary functions with access to data structure -------

  !------------------!
  function getFirst(s)
  !------------------!
    type(RandomWalker), pointer        :: getFirst
    type(RWSample), intent(inout)      :: s
    s%current = 1
    getFirst => s%rwArray(1)
  end function getFirst


  !--------------------------------------!
  subroutine getFirstWalkerBlock(s,bs,fbp)
  !--------------------------------------!
    type(RWSample), intent(inout)      :: s
    integer, intent(inout)             :: bs ! block size
    type(RandomWalker), pointer        :: fbp(:)
    bs = min(bs,s%sampleSize)
    fbp => s%rwArray(1:bs)
    s%current = bs
  end subroutine getFirstWalkerBlock


  !---------------------------------------!
  subroutine getFirstWalker(s,idx,dataPtr)
  !---------------------------------------!
    type(RWSample), intent(inout)  :: s
    integer, intent(out)           :: idx
    type(RandomWalker), pointer    :: dataPtr
    s%current = 1
    idx = 1
    dataPtr => s%rwArray(1)
  end subroutine getFirstWalker


  !-----------------!
  function getNext(s)
  !-----------------!
    type(RandomWalker), pointer        :: getNext
    type(RWSample), intent(inout)      :: s
    s%current = s%current + 1
    getNext => s%rwArray(s%current)
  end function getNext

  !-------------------------------------!
  subroutine getNextWalkerBlock(s,bs,nbp)
  !-------------------------------------!
    type(RWSample), intent(inout)      :: s
    integer, intent(inout)             :: bs  ! block size
    type(RandomWalker), pointer        :: nbp(:)
    bs = min(bs,s%sampleSize-s%current)
    if (bs > 0) then
       nbp => s%rwArray(s%current+1:s%current+bs)
    else
       nbp => null()
    endif
    s%current = s%current + bs
  end subroutine getNextWalkerBlock



  !--------------------------------------!
  subroutine getNextWalker(s,idx,dataPtr)
  !--------------------------------------!
    type(RWSample), intent(inout)  :: s
    integer, intent(out)           :: idx
    type(RandomWalker), pointer    :: dataPtr

    s%current = s%current + 1
    idx = s%current
    dataPtr => s%rwArray(s%current)
  end subroutine getNextWalker


  !-----------------!
  function getLast(s)
  !-----------------!
    type(RandomWalker), pointer        :: getLast
    type(RWSample), intent(inout)      :: s
    getLast => s%rwArray(s%sampleSize)
  end function getLast


  !--------------------------------------!
  subroutine getLastWalker(s,idx,dataPtr)
  !--------------------------------------!
    type(RWSample), intent(inout) :: s
    integer, intent(out)          :: idx
    type(RandomWalker), pointer   :: dataPtr
    idx = s%sampleSize
    dataPtr => s%rwArray(s%sampleSize)
  end subroutine getLastWalker

  !-----------------------!
  function getWalker(s,idx)
  !-----------------------!
    type(RandomWalker), pointer        :: getWalker
    type(RWSample), intent(inout)      :: s
    integer,intent(in)                  :: idx

    getWalker => s%rwArray(idx)

  end function getWalker


  !------------------------!
  logical function isNext(s)
  !------------------------!
    type(RWSample), intent(in)         :: s
    isNext = ( s%sampleSize > s%current )
  end function isNext


  !------------------------!
  logical function isLast(s)
  !------------------------!
    type(RWSample), intent(in)         :: s
    isLast = ( s%sampleSize == s%current )
  end function isLast


  !-------------------------!
  logical function isValid(s)
  !-------------------------!
    type(RWSample), intent(in)         :: s
    isValid = ( s%current <= s%sampleSize )
  end function isValid


  !-------------------------!
  logical function isEmpty(s)
  !-------------------------!
    type(RWSample), intent(in)      :: s
    isEmpty = (s%sampleSize == 0)
  end function isEmpty


  !-------------------------------!
  integer function getSampleSize(s)
  !-------------------------------!
    type(RWSample), intent(in)      :: s
    getSampleSize = s%sampleSize
  end function getSampleSize


  !----------------------------------!
  integer function getMaxSampleSize(s)
  !----------------------------------!
    type(RWSample), intent(in)      :: s
    getMaxSampleSize = size(s%rwArray)
  end function getMaxSampleSize

  !---------------------------------------!
  subroutine getSampleSizes(s, sizes)
  !---------------------------------------!
    type(RWSample), intent(in) :: s
    integer, intent(inout) :: sizes(:)
    integer :: ierr
    call assert(size(sizes) == nproc, "Incorrectly allocated array in getSampleSizes")
    call myMPIGatherInteger(getSampleSize(s),1,sizes,ierr)
  end subroutine getSampleSizes

  !---------------------------------------!
  integer function getSampleSizeAllNodes(s)
  !---------------------------------------!
    type(RWSample), intent(in)      :: s
    integer sum,allSum
    sum = s%sampleSize
    call myMPIAllReduceSumInteger(sum,allSum,1)
    getSampleSizeAllNodes = allSum
  end function getSampleSizeAllNodes


  !--------------------------------!
  real*8 function getSampleWeight(s)
  !--------------------------------!
    type(RWSample), intent(in)      :: s
    real*8                          :: sum
    integer                         :: i
    sum = 0
    do i=1,s%sampleSize
       sum = sum + wgt(s%rwArray(i))
    enddo
    getSampleWeight = sum
  end function getSampleWeight


  !----------------------------------------!
  real*8 function getSampleWeightAllNodes(s)
  !----------------------------------------!
    type(RWSample), intent(in)      :: s
    real*8                          :: sum,allSum
    integer                         :: i
    sum = 0
    do i=1,s%sampleSize
       sum = sum + wgt(s%rwArray(i))
    enddo

    call myMPIAllReduceSumDouble(sum,allSum,1)
    getSampleWeightAllNodes = allSum
  end function getSampleWeightAllNodes


  !---------------------------------------------------!
  function getMeanWeight(s)
  !---------------------------------------------------!
    type(RWSample), intent(in)      :: s
    real*8                          :: getMeanWeight
       real*8                       :: var,stdDev

    integer                         :: i
    type(SimpleStat)                :: stat

    call reset(stat)
    do i=1,s%sampleSize
       call addData(stat,wgt(s%rwArray(i)))
    enddo

    getMeanWeight = mean(stat)


  end function getMeanWeight


  !----------------------------------------!
  real*8 function getMeanWeightAllNodes(s)
  !----------------------------------------!
    type(RWSample), intent(in)      :: s
    real*8                          :: sum,allSum
    integer                         :: i
    sum = 0
    do i=1,s%sampleSize
       sum = sum + wgt(s%rwArray(i))
    enddo

    call myMPIAllReduceSumDouble(sum,allSum,1)
    getMeanWeightAllNodes = allSum/getSampleSizeAllNodes(s)

  end function getMeanWeightAllNodes



  !---------------------------------------------------!
  subroutine getSampleEnergyAndVariance(s,E,var,stdDev)
  !---------------------------------------------------!
    type(RWSample), intent(in)      :: s
    real*8,intent(out)              :: E,var,stdDev

    integer                         :: i
    type(weightStat)                :: stat

    call reset(stat)
    do i=1,s%sampleSize
       call addData(stat,E_local(s%rwArray(i)),wgt(s%rwArray(i)))
    enddo

    E = mean(stat)
    var = variance(stat)
    stdDev = stdDevMean(stat)
  end subroutine getSampleEnergyAndVariance


  !-----------------------------------------------------------!
  subroutine getSampleEnergyAndVarianceAllNodes(s,E,var,stdDev)
  !-----------------------------------------------------------!
    type(RWSample), intent(in)      :: s
    real*8,intent(out)              :: E,var,stdDev

    integer                         :: i
    type(weightStat)                :: stat

    call reset(stat)
    do i=1,s%sampleSize
       call addData(stat,E_local(s%rwArray(i)),wgt(s%rwArray(i)))
    enddo

    E = meanAllNodes(stat)
    var = varianceAllNodes(stat)
    stdDev = stdDevMeanAllNodes(stat)
  end subroutine getSampleEnergyAndVarianceAllNodes


  !----------------------------------------------------!
  subroutine writeSampleHistogram(s,emean,var,fac,nbins)
  !----------------------------------------------------!
    type(RWSample), intent(in)      :: s
    real*8, intent(in)              :: emean          ! mean energy
    real*8, intent(in)              :: var            ! mean energy
    real*8, intent(in)              :: fac            ! histogram interval emean +/- fac*sqrt(var)
    integer, intent(in)             :: nbins

    integer                         :: i
    real*8                          :: emin,emax
    type(DoubleHistogram)           :: stat
    character(len=4)                :: tid
    character(len=2)                :: countstr

    emin = emean - fac*sqrt(var)
    emax = emean + fac*sqrt(var)
    mHistogramCounter = mHistogramCounter + 1
    write(countstr,'(i2.2)') mHistogramCounter

    call dh_createHistogram(stat,emin,emax,nbins)

    do i=1,s%sampleSize
       call dh_addData(stat,E_local(s%rwArray(i)),wgt(s%rwArray(i)))
    enddo

    if (PARALLEL_RUN .and. logmode >= 3) then
       call getMyTaskIdChar(tid)
       open(31,file=trim(basename)//trim(tid)//'-hist'//countstr//'.txt')
       call dh_writeMeanHistogram(stat,31)
       close(31)
    else if (MASTER .and. logmode==2) then
       open(31,file=trim(basename)//'-hist'//countstr//'.txt')
       call dh_writeMeanHistogram(stat,31)
       close(31)
    endif

  end subroutine writeSampleHistogram


  !---------------------------------------------------------------!
  subroutine removeSampleOutliers(s,emean,var,fac,rfac,riter,nbins)
  !---------------------------------------------------------------!

    type(RWSample), intent(inout)   :: s              ! sample
    real*8, intent(inout)           :: emean          ! mean energy, updated after removal of outliers
    real*8, intent(inout)           :: var            ! mean energy, updated after removal of outliers
    real*8, intent(in)              :: fac            ! histogram interval emean +/- fac*sqrt(var)
    real*8, intent(in)              :: rfac           ! allow rfac more in percentile than normally expected
    integer, intent(in)             :: riter
    integer, intent(in)             :: nbins
    real*8                          :: emin,emax
    type(IntList)                   :: leftList,rightList
    type(IntNode), pointer          :: p
    real*8                          :: leftProb,rightProb,sigma
    integer                         :: i,j,k,it,maxData,idx,keepLeft,keepRight,removeLeft,removeRight,start
    type(IndexedHistogram)          :: stat

    maxData = s%sampleSize
    emin = emean - fac*sqrt(var)
    emax = emean + fac*sqrt(var)
    call idh_createHistogram(stat,emin,emax,nbins,maxData)

    do it=1,riter
       do i=1,s%sampleSize
          call idh_addData(stat,E_local(s%rwArray(i)),wgt(s%rwArray(i)))
       enddo

       call idh_getOutliers(stat,leftList,leftProb,rightList,rightProb)
       emean = idh_mean(stat); var = idh_variance(stat)

       if (logmode>=3) then
          write(iull,*) ' Outlier analysis (this process only)'
          write(iull,*) ' E_mean = ',emean,' var = ',var
          write(iull,*) ' left outliers:'
          write(iull,'(a,i5,g11.3)') ' count / expected count in normal distribution: ',getSize(leftList),leftProb
          write(iull,*) ' right outliers:'
          write(iull,'(a,i5,g11.3)') ' count / expected count in normal distribution: ',getSize(rightList),rightProb
       endif

       keepLeft = int(leftProb * rfac)
       removeLeft = getSize(leftList) - keepLeft
       keepRight = int(rightProb * rfac)
       removeRight = getSize(rightList) - keepRight

    if (removeLeft<=0 .and. removeRight<=0) exit
       ! erasing randomly "keepLeft" elements from left outlier list
       if (keepLeft >= getSize(leftList)) then
          call clear(leftList)
       else
          do i=1,keepLeft
             k = int(myran()*getSize(leftList)) + 1
             p => beginPtr(leftList)
             do j=1,k-1
                call nextPtr(p)
             enddo
             p => erase(leftList,p)
          enddo
       endif
       if (keepRight >= getSize(rightList)) then
          call clear(rightList)
       else
          do i=1,keepRight
             k = int(myran()*getSize(rightList)) + 1
             p => beginPtr(rightList)
             do j=1,k-1
                call nextPtr(p)
             enddo
             p => erase(rightList,p)
          enddo
       endif

       call mergeList(leftList,rightList)
       call sort(leftList)
       p => endPtr(leftList)
       do
           if (isBegin(leftList,p)) exit
           call previousPtr(p)
           call deleteWalker(s,dataPtr(p))
       enddo

       if (logmode>=3) then
          write(iull,*) ' removed ',max(removeLeft,0),' left outliers and ',max(removeRight,0),' right outliers '
          write(iull,*) ' new sample size = ',getSampleSize(s)
       endif

       call clear(leftList)
       call clear(rightList)
       emin = emean - fac*sqrt(var)
       emax = emean + fac*sqrt(var)
       call idh_reset(stat,emin,emax)
    enddo

    call destroyList(leftList)
    call destroyList(rightList)

    call getSampleEnergyAndVarianceAllNodes(s,emean,var,sigma)
    if (MASTER .and. logmode >=2) then
       write(iul,'(/a)') ' Outlier removal:'
       write(iul,'(a,i8,a,f15.5,a,f10.5,a,f15.2)') ' final size=',getSampleSize(s),' E_mean=',emean,' +/- ',sigma,' var=',var
    endif
  end subroutine removeSampleOutliers


  !----------------------------------------------------!
  subroutine removeSampleOutliers2(s,tol,global,replace)
  !----------------------------------------------------!

#ifdef PARALLEL
     include 'mpif.h'
#endif
     type(RWSample), intent(inout) :: s
     real*8, intent(in)            :: tol              ! tolerance for outlier removal (factor of sigma)
     logical, intent(in)           :: global           ! treat distributed sample as one sample or determine outliers only locally
     logical, intent(in)           :: replace          ! replace removed walkers by duplicating random walkers
     type(entry), pointer      :: walkerList(:)
     real*8  x,p,loga,frac,inverse_erf,quantile,b
     real*8  mean,emean,var,sigma
     integer alstat,i,n,listSize,tid,localIdx,ssize,oldSize
     integer, allocatable          :: killList(:),localList(:),killIdx(:),countOverflow(:)
     real*8, parameter  :: a=0.147            ! a const from inv_erf approx
     integer, parameter :: maxKill = 15       ! percent to be removed at most
     ! variables for parallel merge
     integer            :: stride    ! current stride for joining tid: 2**iter
     integer            :: nap       ! number of "active" processes participating in join
     integer            :: idlemax   ! the highest "active" tid may be "idle"
     integer            :: tag1, nmax, nmax0, count, nSize
#ifdef PARALLEL
     integer            :: ierr, status(MPI_STATUS_SIZE)
#endif
     real               :: starttime, endtime
     real*8, pointer    :: vec(:),vecrecv(:),vecjoin(:)
     integer            :: thissize(1),sizes(maxMPINodes),maxsize

     call assert(s%sampleSize * nproc >= 10,"remove_outliers with trimmed mean requires at least 10 walkers")

     allocate(walkerList(s%sampleSize),stat=alstat)
     if (alstat /= 0) call abortp("rwSample:removeSampleOutliers2:allocate failed")

     !!print*, "DBG: in ro2",tol,global,replace
     do i=1,s%sampleSize
        walkerList(i) = entry(idx=i,tid=mytid,data=E_local(s%rwArray(i)))
     enddo

     call quickSortEntry(walkerList)

     !! parallel (global) version: merge locally sorted lists hierarchically
     !! 0-1, 2-3, 4-5, 6-7 ...
     !! 0-1-2-3, 4-5-6-7 ...
     !! until master has fully sorted list. careful about the size of the lists!!
     !! "idlemax" indicates the last process that has no partner (occurs when nproc is not
     !! 2**n)

     if (PARALLEL_RUN .and. global) then

        n = s%sampleSize
        thissize(1) = n
#ifdef PARALLEL
        call MPI_GATHER(thissize,1,MPI_INTEGER,sizes,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
        if (MASTER) then
           maxsize = maxval(sizes(1:nproc))
           thissize(1) = maxsize
        endif
        call MPI_BCAST(thissize,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
#endif
        nmax0 = thissize(1)
        nmax = nmax0

        allocate(vec(2*nmax))
        call entryToVector(walkerList,vec,n)
        nullify(vecrecv)
        nullify(vecjoin)
        deallocate(walkerList)
        ! from here use size of linear vec = 2*size of entry vector
        n = 2*s%sampleSize
        nmax = 2*nmax

        if (MASTER) call cpu_time(starttime)
        stride = 1     ! for joining mytid. Is 2**iter
        nap = nproc
        do
        if (stride >= nproc) exit
           idlemax = stride * ((nap/2)*2)
           !!if (MASTER) print*, "nap,idlemax=",nap,idlemax
           if (mod(mytid,stride)==0) then
              ! active
              if (mytid/=idlemax) then
                 if (mod(mytid,2*stride)==0) then
                    !!print*, mytid, " <- ",mytid+stride
                    tag1 = 2*stride
                    !! get size of vector to receive
                    allocate(vecrecv(nmax))
#ifdef PARALLEL
                    call MPI_RECV(vecrecv,nmax,MPI_DOUBLE_PRECISION,mytid+stride,tag1,MPI_COMM_WORLD,status,ierr)
                    call MPI_GET_COUNT(status,MPI_DOUBLE_PRECISION,count,ierr)
                    !!print*, "recv:",mytid,count,status(MPI_SOURCE),status(MPI_TAG)
#endif
                    allocate(vecjoin(2*nmax))
                    call mergeArrays(vec,n,vecrecv,count,vecjoin)
                    !!print*, "merge:",mytid,n,count,(vecjoin(i),i=1,n+count)
                    deallocate(vec,vecrecv)
                    vec => vecjoin
                    nullify(vecjoin)
                    n = n + count
                    nmax = 2*nmax
                 else
                    !!print*, mytid, " -> ",mytid-stride
                    tag1 = 2*stride
#ifdef PARALLEL
                    call MPI_SEND(vec,n,MPI_DOUBLE_PRECISION,mytid-stride,tag1,MPI_COMM_WORLD,ierr)
#endif
                    !!print*, "sent:", mytid, n," elements with tag ",tag1
                 end if
              else
                 !!print*, mytid, "is idle"
              end if
           end if
           nap = (nap+1)/2
           stride = stride * 2
        end do


        if (MASTER .and. logmode>=2) then
           call cpu_time(endtime)
           write(iul,*) " parallel merge required: ",endtime-starttime," sec."
        end if

        if (MASTER) then
           ! back to entry size
           n = n/2
           allocate(walkerList(n))
           call entryFromVector(walkerList,vec,n)
        end if

        deallocate(vec)

     end if  ! PARALLEL_RUN

     if (MASTER) then
        ! calculate 10% trimmed mean and variance
        mean = 0; var = 0
        ssize = size(walkerList)        ! in parallel case /= s%sampleSize!
        do i=ssize/10,(9*ssize)/10
           mean = mean + walkerList(i)%data
           var = var + walkerList(i)%data**2
        end do
        n = (9*ssize)/10 - ssize/10 + 1
        mean = mean / n
        var = var / n - mean**2
        sigma = sqrt(var)
        if (logmode >= 2) then
           write(iul,'(2(a,f12.3))') "  trimmed mean = ",mean,"  sigma = ",sigma
        end if

        !! now remove outliers from left and right

        listSize = (ssize/nproc * maxKill) / 100
        !!print*,"DBG: listsize=",listSize,nproc,ssize
        allocate(killList(nproc*listSize), killIdx(0:nproc-1), countOverflow(0:nproc-1))
        killIdx = 0

        ! with sample size N, one walker can be expected in the 1/N quantil
        ! zero walkers in the 1/(10*N) quantil (corresponds to 1/(5*N) probability
        ! if right side is also taken into account. See elaboration
        p = 0.1d0 / ssize
        x = p - 1.d0
        frac = 2d0/(Pi*a)
        loga = log(1d0-x*x)
        inverse_erf = sign(1d0,x)*sqrt( sqrt( (frac + 0.5d0*loga)**2 - loga/a ) &
                      - (frac + 0.5d0*loga) )
        quantile = -sqrt(2d0)*inverse_erf
        b = tol*quantile
        !!print*,"DBG:",p,quantile,b

        do i=1,ssize/2
        if (walkerList(i)%data > mean - b*sigma) exit
           tid = walkerList(i)%tid
           killIdx(tid) = killIdx(tid) + 1
           if (killIdx(tid) <= listSize) then
              killList(tid*listSize+killIdx(tid)) = walkerList(i)%idx
           end if
        end do
        do i=ssize,ssize/2,-1
        if (walkerList(i)%data < mean + b*sigma) exit
           tid = walkerList(i)%tid
           killIdx(tid) = killIdx(tid) + 1
           if (killIdx(tid) <= listSize) then
              killList(tid*listSize+killIdx(tid)) = walkerList(i)%idx
           end if
        end do
        !!print*,"DBG: killIdx=",killIdx
        !!print*,"DBG: # to delete:",sum(killIdx)
        countOverflow = 0
        where (killIdx > listSize)
           countOverflow = 1
        end where
        !!print*,"DBG: processes with overflow:",sum(countOverflow)
        ! cut killIdx to listSize
        killIdx = min(killIdx,listSize)
        if (logmode>=2) then
           if (replace) then
              write(iul,'(i7,a)') sum(killIdx)," walkers will be deleted and replaced"
           else
              write(iul,'(i7,a)') sum(killIdx)," walkers will be deleted and not replaced"
           end if
           write(iul,'(i7,a,i7,a)') sum(countOverflow)," processes attempted to exceed the limit of ",listSize," deletions"
        end if
        deallocate(countOverflow)
     else
        allocate(killList(1),killIdx(0:0))
     end if

#ifdef PARALLEL
     call MPI_BCAST(listSize,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
#endif
     allocate(localList(listSize))
#ifdef PARALLEL
     call MPI_SCATTER(killList,listSize,MPI_INTEGER,localList,listSize,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
     call MPI_SCATTER(killIdx,1,MPI_INTEGER,localIdx,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
#else
     localList = killList
     localIdx = killIdx(0)
#endif
     oldSize = s%sampleSize
     if (localIdx > 0) then
        call deleteWalkerList(s,localList,localIdx)
     end if

     if (replace .and. localIdx>0) then
        call changeSampleSizeTo(s,oldSize)
     end if

     deallocate(killList,localList,killIdx)

#ifdef PARALLEL
     call MPI_BARRIER(MPI_COMM_WORLD,ierr)
#endif
     if (logmode >= 3) then
        write(iull,'(/a)') ' local sample after outlier removal:'
        call getSampleEnergyAndVariance(s,emean,var,sigma)
        write(iull,'(a,i6,a,i8,a,f15.5,a,f10.5,a,f15.2)') ' tid=',mytid,': final size=', &
           getSampleSize(s),' E_mean=',emean,' +/- ',sigma,' var=',var
     end if
     call getSampleEnergyAndVarianceAllNodes(s,emean,var,sigma)
     nSize = getSampleSizeAllNodes(s)
     if (MASTER .and. logmode >=2) then
        write(iul,'(/a)') '  after outlier removal:'
        write(iul,'(a,i8,a,f15.5,a,f10.5,a,f15.2)') '  final total walker size=',nSize, &
           ' E_mean=',emean,' +/- ',sigma,' var=',var
     end if
     call setCurrentResult(emean,sigma,var)


   contains

   subroutine mergeArrays(vec,n,vecrecv,nr,vecjoin)
      real*8, pointer :: vec(:),vecrecv(:),vecjoin(:)
      integer, intent(in) :: n, nr   ! actual sizes of vectors
      integer i,i1,i2,sizejoin
      sizejoin = n + nr
      i1 = 2; i2 = 2
      do i=2,sizejoin,2
         if (i1 > n) then
            vecjoin(i-1:i) = vecrecv(i2-1:i2)
            i2 = i2+2
         else if (i2 > nr) then
            vecjoin(i-1:i) = vec(i1-1:i1)
            i1 = i1+2
         else if (vec(i1) < vecrecv(i2)) then
            vecjoin(i-1:i) = vec(i1-1:i1)
            i1 = i1+2
         else
            vecjoin(i-1:i) = vecrecv(i2-1:i2)
            i2 = i2+2
         endif
      end do
   end subroutine mergeArrays


  end subroutine removeSampleOutliers2

!---------- elementary functions that change sample ------------------


  !-----------------------------!
  subroutine deleteWalker(s,idx)
  !-----------------------------!
    type(RWSample), intent(inout)               :: s
    integer, intent(in)                         :: idx

    call assert(idx>0.and.idx<=s%sampleSize,   &
         "RWSample: deleting illegal element")
    call rw_destroy(s%rwArray(idx))
    if (idx < s%sampleSize) s%rwArray(idx) = s%rwArray(s%sampleSize)
    s%sampleSize = s%sampleSize - 1
  end subroutine deleteWalker


  !-------------------------------------!
  subroutine deleteWalkerList(s,list,idx)
  !-------------------------------------!
    type(RWSample), intent(inout)               :: s
    integer, intent(inout)                      :: list(:)
    integer, intent(in)                         :: idx
    integer i,j

    call assert(idx>=0.and.size(list)>=idx,   &
         "RWSample:deletaWalkerList: illegal arguments")
    do i=1,idx
       call rw_destroy(s%rwArray(list(i)))
       if (list(i) < s%sampleSize) then
          s%rwArray(list(i)) = s%rwArray(s%sampleSize)
          do j=i+1,idx
             if (list(j)==s%sampleSize) then
                list(j) = list(i)
             end if
          end do
       end if
       s%sampleSize = s%sampleSize - 1
    end do
  end subroutine deleteWalkerList


  !---------------------------!
  subroutine appendWalker(s,rw)
  !---------------------------!
    ! note appends _copy_ of walker rw
    type(RWSample), intent(inout)          :: s
    type(RandomWalker), intent(in)         :: rw

    if (s%sampleSize >= size(s%rwArray)) call reallocateSample(s,s%sampleSize + 1)
    s%sampleSize = s%sampleSize + 1
    s%rwArray(s%sampleSize) = rw
  end subroutine appendWalker

  !---------------------------!
  subroutine reallocateSample(s,n)
  !---------------------------!
    type(RWSample), intent(inout) :: s
    integer, intent(in)           :: n
    type(RandomWalker), pointer :: rwArray(:)
    integer :: i, ierr, walkersKept

    call assert(associated(s%rwArray), "Trying to reallocate uninitialized sample")

    !!!write(999,*) "WARNING:reallocating Sample:",getMyTaskId(),s%sampleSize,size(s%rwArray),n

    ! allocate to target size, not walkersKept. this way, the sample array can
    ! be over-allocated, which may prevent having to resize the array multiple
    ! times
    allocate(rwArray(n), stat=ierr)
    call assert(ierr == 0, "allocation failed in reallocateSample")

    ! makes *copies* of the walkers. can't use pointer because we have to
    ! re-allocate s%rwArray, but the only function capable of effectively
    ! re-allocating (move_alloc) is only allowed on allocatable arrays, not
    ! on pointer arrays. so to re-allocate we have to effectively deallocate
    ! and then allocate again, which makes using pointers impossible
    walkersKept = min(s%sampleSize, n)

    ! DO NOT TURN THIS INTO A LIST ASSIGNMENT
    ! BOTH GFORTRAN AND IFORT GENERATE BROKEN CODE FOR THAT
    do i = 1, walkersKept
      rwArray(i) = s%rwArray(i)
    enddo

    ! deallocate initial walkers
    do i = s%sampleSize, 1, -1
      call rw_destroy(s%rwArray(i))
    enddo
    deallocate(s%rwArray)

    s%rwArray => rwArray
    s%sampleSize = walkersKept
  end subroutine reallocateSample

  !-----------------------------------!
  subroutine replaceWalker(s,idx1,idx2)
  !-----------------------------------!

    type(RWSample), intent(inout)          :: s
    integer, intent(in)                    :: idx1,idx2

    call assert(idx1>0.and.idx1<=s%sampleSize,"RWSample: copying illegal element idx1")
    call assert(idx2>0.and.idx2<=s%sampleSize,"RWSample: copying illegal element idx2")
    call rw_destroy(s%rwArray(idx1))
    s%rwArray(idx1) = s%rwArray(idx2)

  end subroutine replaceWalker

  !-----------------------------------!
  subroutine copyCurrentWalkerTo(s,idx)
  !-----------------------------------!

    type(RWSample), intent(inout)          :: s
    integer, intent(in)                    :: idx

    call assert(idx>0.and.idx<=s%sampleSize,"RWSample: copying illegal element")
    s%rwArray(idx) = s%rwArray(s%current)
  end subroutine copyCurrentWalkerTo

!------------ compound functions using elementary functions -------------


  !-------------------------------!
  subroutine deleteCurrentWalker(s)
  !-------------------------------!
    type(RWSample), intent(inout)      :: s

    call deleteWalker(s,s%current)
  end subroutine deleteCurrentWalker



  !------------------------!
  subroutine resetWeights(s)
  !------------------------!
    type(RWSample), intent(inout)      :: s
    integer i

    do i=1,s%sampleSize
       call setWeight(s%rwArray(i),1.d0)
    enddo
  end subroutine resetWeights

  !------------------------------!
  subroutine resetPersistencies(s)
  !------------------------------!
    type(RWSample), intent(inout)      :: s
    integer i

    do i=1,s%sampleSize
       call resetPersistency(s%rwArray(i))
    enddo
  end subroutine resetPersistencies


  !------------------------!
  subroutine resettoMean(s)
  !------------------------!
    type(RWSample), intent(inout)      :: s
    integer i
    real*8 :: w

    w = getMeanWeight(s)


    do i=1,s%sampleSize
       call setWeight(s%rwArray(i),w)
    enddo
  end subroutine resettoMean

  !-----------------------------!
  subroutine recalculateSample(s)
  !-----------------------------!
    type(RWSample), intent(inout)      :: s
    integer i
    do i=1,s%sampleSize
       call recalculateEloc(s%rwArray(i))
    enddo
  end subroutine recalculateSample


  !------------------------!
  subroutine deleteSample(s)
  !------------------------!
    type(RWSample), intent(inout)      :: s
    integer n

    do n=s%sampleSize,1,-1
       call deleteWalker(s,n)
    enddo
  end subroutine deleteSample


  !-------------------------------!
  subroutine multiplyWalker(s,rw,n)
  !-------------------------------!
    ! append (n-1) copies of walker rw
    type(RWSample), intent(inout)          :: s
    type(RandomWalker), target, intent(in) :: rw
    integer, intent(in)                    :: n
    type(RandomWalker), pointer            :: newRW
    integer i

    allocate(newRW)
    newRW = rw          ! deep copy
    do i=2,n
       call appendWalker(s,newRW)
    end do
    deallocate(newRW)
  end subroutine multiplyWalker


  !---------------------------------!
  integer function countSampleSize(s)
  !---------------------------------!
    type(RWSample), intent(inout)      :: s
    countSampleSize = s%sampleSize
  end function countSampleSize

  !--------------------------------!
  subroutine changeSampleSizeTo(s,n)
  !--------------------------------!
    type(RWSample), intent(inout)  :: s
    integer, intent(in)            :: n
    integer i,k,mNew,nOld,sSize
    real*8 xi
    type(RandomWalker), pointer    :: newRW
    real*8, parameter              :: REALLOC_FACTOR = 1.3d0

    call assert(s%sampleSize>=1,"changeSampleSizeTo: impossible for empty sample")

    if (n < s%sampleSize) then
       ! reduce sample
       nOld = s%sampleSize
       do i=s%sampleSize,n+1,-1
          call deleteWalker(s,i)
       enddo
       if (logmode >=3) then
          write(iul,'(a,i7,a,i7,a)') ' sample size has been reduced by ',nOld-n,' to ', n,' walkers per node'
       endif
    else
       ! multiply random walkers
       sSize = s%sampleSize
       mNew = n - sSize
       if (n > size(s%rwArray)) call reallocateSample(s, int(REALLOC_FACTOR*n))
       do i=1,mNew
          xi = myran()
          k = floor(sSize*xi)+1
          allocate(newRW)
          newRW = s%rwArray(k)
          call appendWalker(s,newRW)
       enddo

       if (logmode >= 3 .and. mNew > 0) then
          write(iul,'(a,i7,a,i7,a)') ' sample size has been increased by ',mNew,' (random duplication) to ', n,' walkers per node'
       endif

    endif
    call assert(n==s%sampleSize,"changeSampleSizeTo: failed")

  end subroutine changeSampleSizeTo

  !--------------------------!
  subroutine reduceToLast(s,n)
  !--------------------------!
    type(RWSample), intent(inout) :: s
    integer, intent(in) :: n
    type(RandomWalker), pointer :: rw(:)
    integer :: i, ierr

    call assert(n <= s%sampleSize, "reduceToLast: cannot reduce to larger sample")
    if (s%sampleSize == n) return

    ! copy list of walkers to keep
    ! last walkers, in reverse order
    allocate(rw(n), stat=ierr)
    call assert(ierr == 0, "allocation failed in reduceToLast")
    do i = 1, n
      rw(i) = s%rwArray(s%sampleSize - n + i)
    enddo

    ! deallocate initial walkers
    do i = 1, s%sampleSize
      call rw_destroy(s%rwArray(i))
    enddo
    deallocate(s%rwArray)

    s%rwArray => rw
    s%sampleSize = n
  end subroutine reduceToLast

  !--------------------------------!
  subroutine sendWalkersTo(s,n,node)
  !--------------------------------!
    ! send n walkers (the last n) to 'node'
    type(RWSample), intent(inout)      :: s
    integer, intent(in)                :: n, node
    integer i,idx

    do i=1,n
       idx = s%sampleSize
       call sendRandomWalkerTo(s%rwArray(idx),node,i)
       call deleteWalker(s,idx)
    enddo

  end subroutine sendWalkersTo


  !-------------------------------------!
  subroutine receiveWalkersFrom(s,n,node)
  !-------------------------------------!
    ! receive n walkers from 'node' and append to sample
    type(RWSample), intent(inout)      :: s
    integer, intent(in)                :: n, node
    integer i,idx
    type(RandomWalker), pointer            :: newRW

    allocate(newRW)
    do i=1,n
       call receiveRandomWalkerFrom(newRW,node,i)
       call appendWalker(s,newRW)
    enddo
    call rw_destroy(newRW)
    deallocate(newRW)

  end subroutine receiveWalkersFrom


  !------------------------!
  subroutine SaveFathers(s,i)
  !------------------------!

!     set new father-walker and save father-weights and position
      implicit none
!     input parameter:
      integer, intent(in)           :: i
      type(RWSample),intent(inout)  :: s
!     internal parameters:
      integer             :: n,f

      call assert(i<= nDWGen(),"(SaveFathers): invalid index i")

      do n = 1,s%samplesize
         call setfather(s%rwArray(n),i,n)
      enddo
      end subroutine SaveFathers


  !------------------------------------!
  subroutine GetSonWeight(i,s,wsons)
  !------------------------------------!

!     get position vectors and weights of father-walkers and w.sums of suns
      implicit none
!     input parameter:
      integer,intent(in)         :: i             ! father tag
      type(RWSample),intent(in)  :: s
!     in/output parameter:
      real*8, intent(inout) :: wsons(:)   ! w.sum of sons

      integer :: n,f

      call assert(i<= nDWGen(),"(SaveFathers): invalid index i")

     wsons(:) = 0d0

    do n = 1,s%samplesize
      wsons(father(s%rwArray(n),i))=wsons(father(s%rwArray(n),i))+wgt(s%rwArray(n))
    enddo

      end subroutine GetSonWeight

  !------------------------------!
  subroutine loadBalanceSamples(s)
  !------------------------------!
    type(RWSample), intent(inout)      :: s
    integer                            :: sizes(nproc)
    integer                            :: senders(2*nproc)
    integer                            :: receivers(2*nproc)
    integer                            :: idx(nproc)
    integer i,ierr,sender(2),receiver(2),n,targetSize,sdr,rcv
    integer walkers,walkersToSend,walkersToReceive,walkersMax
    real*8                             :: sizeMean,sizeStdDev
    type(simpleStat)                   :: sizeStat

    if (nproc==1) return

    call myMPIGatherInteger(s%samplesize,1,sizes,ierr)

    senders = 0
    receivers = 0
    call reset(sizeStat)
    if (mytid==0) then
       ! only master has sizes
       ! master calculates for each node (sdr) how many walkers to send
       ! to which node (rcv). Note: each node sends walkers to only
       ! *one* node. Each node receives from only *one* node.
       ! This yields minimal communication but not optimal load balancing
       do i=1,nproc
          call addData(sizeStat,dble(sizes(i)))
          idx(i) = i
       enddo
       sizeMean = mean(sizeStat)
       targetSize = nint(sizeMean)
       sizeStdDev  = sqrt(variance(sizeStat))

       call mySort(idx,nproc)   ! order sizes using index 'idx' decreasing

       ! send walkers from  1->n, 2->n-1 ...
       do i=1,nproc/2
          sdr = idx(i)                                      ! sender index
          walkersToSend = sizes(sdr) - targetSize
          rcv = idx(nproc+1-i)                              ! receiver index
          walkersToReceive = targetSize - sizes(rcv)
          walkers = (walkersToSend + walkersToReceive) / 2  ! send mean
          if (i==1) walkersMax = walkers
          if (walkers <= 0) exit
          if (sizes(rcv)+walkers > 0.9*getMaxSampleSize(s)) write(*,*) 'DBG:',rcv,sizes(rcv),walkers
          senders(2*sdr-1)   = walkers
          senders(2*sdr)     = rcv - 1   ! receiver task id
          receivers(2*rcv-1) = walkers
          receivers(2*rcv)   = sdr - 1   ! sender task id
       enddo
       if (logmode>=3) then
          write(iul,'(A,I4,A,I4)') ' loadBalance: ',i-1,' nodes sending walkers. max=',walkersMax
       endif
    endif

    call myMPIScatterInteger(senders,2,sender,ierr)
    call myMPIScatterInteger(receivers,2,receiver,ierr)

    if (sender(1)>0) then
       call sendWalkersTo(s,sender(1),sender(2))
    endif
    if (receiver(1)>0) then
       call receiveWalkersFrom(s,receiver(1),receiver(2))
    endif

  contains

    subroutine mySort(idx,n)

    ! internal function to sort sizes(idx(i)) with decreasing order
    ! uses bubble sort -> replace by heap/quick sort for large n

    integer n
    integer idx(n)
    integer i,j,tmp

    do i=2,n
       do j=n,i,-1
          if (sizes(idx(j)) > sizes(idx(j-1))) then
             tmp = idx(j); idx(j)=idx(j-1); idx(j-1)=tmp
          endif
       enddo
    enddo

    end subroutine mySort


  end subroutine loadBalanceSamples


  !----------------------------!
  subroutine displaySample(s,iu)
  !----------------------------!
    type(RWSample), intent(inout)  :: s
    integer, optional              :: iu
    integer i

    if (.not.MASTER) return    ! to be replaced by parallel version

    if (present(iu)) then
       write(iu,*) 'displaying sample:'
       do i=1,s%sampleSize
          write(iu,*) i
          call display(s%rwArray(i),2,iu)
       enddo
    else
       write(*,*) 'displaying sample:'
       do i=1,s%sampleSize
          write(*,*) i
          call display(s%rwArray(i),2)
       enddo
    endif
  end subroutine displaySample



  subroutine writeSamplePositions(s,iu,tsize)
  !-----------------------------------------!
     ! write electron positions of sample unformatted to open data stream iu
     ! walker positions are written individually as 3n-dim vectors.
     ! The sample is written in the order node0,node1,..nodeN.
     ! It would be more efficient to transfer walker position blocks via MPI
     ! rather than individual walkers.
#ifdef PARALLEL
     include 'mpif.h'
#endif
     type(RWSample), intent(inout)   :: s
     integer, intent(in)             :: iu
     integer, intent(out)            :: tsize
     real*8                          :: x(nElecs()),y(nElecs()),z(nElecs()),xall(3*nElecs())
     integer i,n,tag1,tag2,proc,count,ssize
#ifdef PARALLEL
     integer            :: ierr, status(MPI_STATUS_SIZE)
#endif
     n = nElecs()
     call assert(n>0,"RWSample:writeSample: number of electrons in RW module not set")
     ssize = 0
     tag1 = 0; tag2 = 0

     tsize = getSampleSizeAllNodes(s)

     if (MASTER) then

        write(iu) n,tsize
        ! first write local sample
        do i=1,s%sampleSize
           call pos(s%rwArray(i),x,y,z)
           xall(1:n) = x
           xall(n+1:2*n) = y
           xall(2*n+1:3*n) = z
           write(iu) xall
        end do


#ifdef PARALLEL
        ! now read samples from NODES
        do proc=1,nproc-1
           call MPI_RECV(ssize,1,MPI_INTEGER,proc,tag1,MPI_COMM_WORLD,status,ierr)
           do i=1,ssize
              tag2 = i
              call MPI_RECV(xall,3*n,MPI_DOUBLE_PRECISION,proc,tag2,MPI_COMM_WORLD,status,ierr)
              call MPI_GET_COUNT(status,MPI_DOUBLE_PRECISION,count,ierr)
              call assert(count==3*n,'writeSamplePositions: error MPI communication (i)')
              write(iu) xall
           end do
        end do
#endif

     else ! SLAVE

#ifdef PARALLEL
        call MPI_SEND(s%sampleSize,1,MPI_INTEGER,0,tag1,MPI_COMM_WORLD,ierr)
        do i=1,s%sampleSize
           call pos(s%rwArray(i),x,y,z)
           xall(1:n) = x
           xall(n+1:2*n) = y
           xall(2*n+1:3*n) = z
           tag2 = i
           call MPI_SEND(xall,3*n,MPI_DOUBLE_PRECISION,0,tag2,MPI_COMM_WORLD,ierr)
        end do
#endif

     end if ! (MASTER)

  end subroutine writeSamplePositions


  subroutine readSamplePositions(s,iu,nsize)
  !----------------------------------------!
     ! read electron positions of sample unformatted to open data stream iu
     ! parallels writeSamplePositions. See format description there
#ifdef PARALLEL
     include 'mpif.h'
#endif
     type(RWSample), intent(inout)   :: s
     integer, intent(in)             :: iu
     integer, intent(inout)          :: nsize  ! total size per node to be read from open file, default=0: all positions
                                               ! on output: initial size on node
     real*8                          :: x(nElecs()),y(nElecs()),z(nElecs()),xall(3*nElecs())
     integer i,n,nn,tag1,tag2,proc,cnt,count,addsize,ssize,tsize
#ifdef PARALLEL
     integer            :: ierr, status(MPI_STATUS_SIZE)
#endif

     tag1 = 0; tag2 = 0
     n = nElecs()
     ssize = 0
     cnt = 0
     call assert(n>0,"RWSample:readSamplePositions: number of electrons in RW module not set")

     if (MASTER) then

        read(iu) nn,tsize
        call assert(nn==n,"RWSample:readSamplePositions: number of electrons differs from pos file")
        if (nsize>0) then
           call assert(nsize*nproc <= tsize,"RWSample:readSamplePositions: not enough walkers in pos file")
           tsize = nsize*nproc
        end if
        ssize = tsize / nproc            ! local sample size
        addsize = mod(tsize,nproc)       ! the last addsize procs get an additional walker

        call assert(ssize <= size(s%rwArray),"RWSample:readSamplePositions: size not large enough")

        ! destroy current sample
        do i=s%sampleSize,1,-1
           call rw_destroy(s%rwArray(i))
        enddo

        ! first read local sample
        do i=1,ssize
           read(iu) xall
           cnt = cnt + 1
           x = xall(1:n)
           y = xall(n+1:2*n)
           z = xall(2*n+1:3*n)
           call rw_new(s%rwArray(i))
           call resetTo(s%rwArray(i),x,y,z)
        end do
        s%sampleSize = ssize

        ! now read further walkers and distribute to NODES
#ifdef PARALLEL
        do proc=1,nproc-1
           if (proc >= nproc-addsize) ssize = ssize + 1
           call MPI_SEND(ssize,1,MPI_INTEGER,proc,tag1,MPI_COMM_WORLD,ierr)
           do i=1,ssize
              read(iu) xall
              cnt = cnt + 1
              tag2 = i
              call MPI_SEND(xall,3*n,MPI_DOUBLE_PRECISION,proc,tag2,MPI_COMM_WORLD,ierr)
           end do
        end do
#endif

        call assert(cnt==tsize,"RWSample:readSamplePositions: internal error")

     else ! SLAVE

#ifdef PARALLEL
        ! destroy current sample
        do i=s%sampleSize,1,-1
           call rw_destroy(s%rwArray(i))
        enddo

        call MPI_RECV(ssize,1,MPI_INTEGER,0,tag1,MPI_COMM_WORLD,status,ierr)
        call MPI_GET_COUNT(status,MPI_INTEGER,count,ierr)
        call assert(count==1,'writeSamplePositions: error MPI communication (i)')
        call assert(ssize <= size(s%rwArray),"RWSample:readSamplePositions: size not large enough")

        ! receive sample from master
        do i=1,ssize
           tag2 = i
           call MPI_RECV(xall,3*n,MPI_DOUBLE_PRECISION,0,tag2,MPI_COMM_WORLD,status,ierr)
           call MPI_GET_COUNT(status,MPI_DOUBLE_PRECISION,count,ierr)
           call assert(count==3*n,'writeSamplePositions: error MPI communication (i)')
           x = xall(1:n)
           y = xall(n+1:2*n)
           z = xall(2*n+1:3*n)
           call rw_new(s%rwArray(i))
           call resetTo(s%rwArray(i),x,y,z)
        end do
        s%sampleSize = ssize

#endif

     end if ! (MASTER)

     nsize = ssize

  end subroutine readSamplePositions


  subroutine writeSampleCommand(lines,nl,s)
  !---------------------------------------!

    integer, intent(in)           :: nl
    character(len=120), intent(in) :: lines(nl)
    type(RWSample), intent(inout) :: s
    character(len=70)             :: posFile

    integer, parameter            :: iu=19
    integer iflag,ssize

    call assert(.not.isEmpty(s),'writeSampleCommand: empty sample')

    if (MASTER) then
       call getstra(lines,nl,'file=',posFile,iflag)
       if (iflag /= 0) then
          if (baseName=="") call abortp('writeSampleCommand: baseName not set')
          posFile = trim(baseName)//'.pos'
       endif
       open(iu,file=posFile,form='unformatted')
    end if

    call writeSamplePositions(s,iu,ssize)

    close(iu)

    if (logmode >= 2) then
       write(iul,'(a,i10,2a)') ' positions of current sample (size =',ssize,') written to file ',trim(posFile)
    endif

  end subroutine writeSampleCommand


end MODULE RWSampleModule

!-----------------------!
subroutine testBalance(s)
!-----------------------!
  use RWSampleModule
  implicit none
  type(RWSample), intent(inout)   :: s
  real*8  E,EAll,stdDev,stdDevAll,var,varAll
  real*8, allocatable :: v(:)
  integer n,nAll,sender,receiver,ierr,sizeRW
  type(RandomWalker), pointer :: rwp,newRW

  nAll = getSampleSizeAllNodes(s)
  call getSampleEnergyAndVarianceAllNodes(s,EAll,varAll,stdDevAll)

  call getSampleEnergyAndVariance(s,E,var,stdDev)
  n = getSampleSize(s)
  if (mytid==0) then
     write(*,*) ' testBalance sample'
     write(*,*)
     write(*,'(A,I5,2F13.4)') ' all nodes:',nAll,EAll,stdDevAll

     write(*,*) ' node    size     E    sigma'
     write(*,*) ' ---------------------------'
  endif
  write(*,'(2I5,2F13.4)') mytid,n,E,stdDev

  call myMPIBarrier(ierr)

  if (nproc < 3) then
     if (mytid==0) write(*,*) ' balance test requires at least 3 nodes required'
     return
  endif

!  sizeRW = sizeOfRandomWalker()
!  allocate(v(sizeRW))
!  rwp => getFirst(s)
!  allocate(newRW)
!
!  call writeRandomWalkerToVector(rwp,v)
!  call readRandomWalkerFromVector(newRW,v)
!
!  write(*,*) mytid,E_local(rwp),E_local(newRW)
!
!  call appendWalker(s,newRW)
!
!  call rw_destroy(newRW)
!  deallocate(newRW)


  sender = 1
  receiver = 2
  if (mytid==sender) call sendWalkersTo(s,2,receiver)
  if (mytid==receiver) call receiveWalkersFrom(s,2,sender)

  call myMPIBarrier(ierr)

  nAll = getSampleSizeAllNodes(s)
  call getSampleEnergyAndVarianceAllNodes(s,EAll,varAll,stdDevAll)

  call getSampleEnergyAndVariance(s,E,var,stdDev)
  n = getSampleSize(s)
  if (mytid==0) then
     write(*,*) ' after moving 2 walkers from 1->2:'
     write(*,*)
     write(*,'(A,I5,2F13.4)') ' all nodes:',nAll,EAll,stdDevAll

     write(*,*) ' node    size     E    sigma'
     write(*,*) ' ---------------------------'
  endif
  write(*,'(2I5,2F13.4)') mytid,n,E,stdDev

  call myMPIBarrier(ierr)

  call loadBalanceSamples(s)

  nAll = getSampleSizeAllNodes(s)
  call getSampleEnergyAndVarianceAllNodes(s,EAll,varAll,stdDevAll)

  call getSampleEnergyAndVariance(s,E,var,stdDev)
  n = getSampleSize(s)
  if (mytid==0) then
     write(*,*) ' after load balance call:'
     write(*,*)
     write(*,'(A,I5,2F13.4)') ' all nodes:',nAll,EAll,stdDevAll

     write(*,*) ' node    size     E    sigma'
     write(*,*) ' ---------------------------'
  endif
  write(*,'(2I5,2F13.4)') mytid,n,E,stdDev



end subroutine testBalance
