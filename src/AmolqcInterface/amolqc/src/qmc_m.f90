
MODULE qmc

! containing qmc parameters and qmc run routines
! correction in Branching -> NSplit ps=0
! correction of bgte in qmc_m.f90 and propagator_m.f90 included

  use omp_lib
  use global
  use wfdata, only: do_epart
  use properties
  use statistics
  use newstatistics, only: stat,vectorstat
  use RandomWalkerModule
  use RWSampleModule
  use Propagator
  use RWStatistics
  use utilsmodule
  use Utils, only: getTimeString
  use elocaldata, only: setElocCutoff, getElocCutoff, eloc_initialize, &
         mElocCut, mElocCutCount, mDriftCut, mDriftCutCount
  use elocal
  use aosdata, only: aos_initialize
  use mos, only: mos_initialize
  use wfmodule, only: trajectory_head
  use reconfg
  use psimax
  use sed
!  use epart
  use epartModule
!  use assignment
  use assign
  use coc, only: coc_isConverged, coc_add
  use QmcSample, only: currentSampleSize

!  use ce_new_m
  implicit none

  private
  public  :: qmc_init, initWalkerStat, initTrajectory, &
    qmc_writeParams, qmc_run, qmc_Energy, qmc_StdDev, qmc_Variance, qmc_cocIsConverged

  integer, parameter :: WEIGHT_NONE=0, WEIGHT_REY=1, WEIGHT_UMR=2, WEIGHT_ACCEPT=3
  integer, parameter :: POP_NONE=0, POP_GLOBAL=1, POP_LOCAL=2

  ! parameters (constant during run)
  real*8      :: mWFac=1.d0                ! scale factor in Umrigar pop control
  real*8      :: mLThresh=0.5d0, mUThresh=2.0d0  ! thresholds for branching
  real*8      :: mTargetAR=0.5d0           ! target acceptance ratio (AR) for AR adaptation
  real*8      :: mEpsAR=0.01d0             ! required accuracy for AR
  real*8      :: mStdDev=0.05d0            ! target std deviation for result
  integer*8   :: mSteps=1000               !
  integer*8   :: mBlockLen=100             ! blocking of steps
  integer*8   :: mStepsDiscard=0           ! steps to discard initially
  integer*8   :: mStepStride=10            ! stride for data collection
  integer     :: mTargetSampleSize=0             ! actual size for VMC, target size for DMC, this node
  integer     :: mTargetSampleSizeAllNodes=0     ! same for all nodes
  integer     :: mPerMax=10                ! max allowed persistency
  integer     :: mWeight=WEIGHT_NONE       ! weighting strategy (none=0,Rey=1,Umr,Accept)
  integer     :: mPopctrl=POP_GLOBAL       ! 0: noPopctrl,=1 : Global Popctrl in MPIJob, =2 : separated/local for each node
  integer     :: mRcf=1                    ! 1: Our Method 2: Implementation following Caffarel
  integer     :: mWalkerBlock=1            ! # of walkers to be moved and calculated "simultaneously"
  integer     :: mAutoCorrMax=500          ! calc autocorrelation up to mAutoCorrMax steps
  integer     :: mShowDetails=0            ! details in output
  logical     :: mCheckStdDev=.false.      ! check if target std dev is reached
  logical     :: mBranch=.false.
  logical     :: mReconf=.false.
  logical     :: mJoin=.true.
  logical     :: mKill=.false.
  logical     :: mShowSteps=.false.
  logical     :: mSplit=.true.
  logical     :: mKillPersist=.false.
  logical     :: mWalkerStatistics = .false.  ! do topological analysis for walkers?
  logical     :: mTrajectory = .false.        ! record trajectories of walkers?
  logical     :: mLoadBalance = .false.       ! balance sample sizes in parallel runs
  logical     :: mAdaptTau = .false.          ! adapt tau at beginning of run to some AR
  logical     :: mAutoCorrelation = .false.   ! calc autocorrelation of E_L data (only VMC)
  character*1 :: mStatType='c'
  character*3 :: mMethod='VMC'        ! QMC method

  !Future Walking
  logical    :: mFuture = .false.            ! switch for Future-Walking
  integer    :: mDprops = 0                  ! number of Generations (saved fathers managed)
  integer    :: mStprops = 0                 ! Offset between fathers taken (tau(corr) should be taken into account)
                                              ! stprops*dprops should equal Block length
  integer    :: mTau2   = 0                 ! Steps between output pro Generation (tau itself is like die usual tau)
  integer    :: mNotau2 = 0                  ! no. of outputs pro stprops*dprops

  ! status parameters
  integer :: sBlock=0
  integer :: sSampleSize=0
  integer :: sOldLogmode=0
  real*8  :: sERef=0
  integer :: sMode=0            ! propagator mode
  real*8  :: sEMean=0           ! QMC result: Energy
  real*8  :: sEMeanStdDev = 0   ! std deviation (of mean)
  real*8  :: sVar=0             ! variance (of local energy)
  type(simpleStat),save :: sBlockStat
  type(weightStat),save :: sTotalStat
  !Assignment parameters
  logical             :: mAssignment=.false.
  ! SED parameters
  logical             :: mSED = .false.
  ! Energy partitioning parameter
  logical             :: mEpartSED = .false.
  ! Center of Charge parameter
  logical             :: mCOC = .false.
  logical             :: sCOC_isConverged = .false.
  ! Parameters for Maxima Search
  logical            :: mMaxSearch = .false.
  logical            :: mMaxAnalysis = .false.
  integer            :: mBlockOffset = 1
  real*8,allocatable :: x_max(:),y_max(:),z_max(:)
  real*8             :: best
  ! accumulate seen samples into a history
  logical            :: mAccumulate = .false.
  integer            :: mAccumulateSampleSize = 0

contains


  subroutine qmc_init(lines,nl,sample)
  !----------------------------------!

    type(RWSample), intent(inout) :: sample

    character(len=120), intent(in) :: lines(:)
    integer, intent(in)           :: nl

    integer iflag,i
    character(len=3)           :: s,mt
    character(len=12)          :: s1
    logical                    :: found,found1,found2
    real*8 tau,ds,cf
    integer wgt,bgte,move,tmove,tauFlag,iflag1,vb,blockdiscard,tmoverej, tmovew, tmovecross
    integer*8 accSize
    logical join,kill,split,ar,rc,eLocalCutOff,changed
    integer :: int_rcv(nproc), ierr

    mMaxSearch = psimax_doMaxSearch()
    mMaxAnalysis = psimax_doMaxAnalysis()
    mSED = finda(lines,nl,'$sed')
    mEpartSED = .false.
    mCOC = finda(lines,nl,'$iterate_coc')

    ! set default values for standard methods
    found1 = finda(lines,nl,'vmc')
    found2 = finda(lines,nl,'VMC')
    if (found1 .or. found2 .or. mSED .or. mCOC) then
       ! set VMC default values)
       mMethod='VMC';ar=.true.;rc=.false.; move=2; mWeight=WEIGHT_NONE; mBranch=.false.; mKillPersist=.false.
       mPerMax=0; mLoadBalance=.false.;eLocalCutOff=.false.
       mBlockLen=200; mStepsDiscard=1000; mAutoCorrelation=.true.
       mTargetAR=0.5d0
    endif

    found1 = finda(lines,nl,'dmc')
    found2 = finda(lines,nl,'DMC')
    mFuture = finda(lines,nl,'FWDMC')
    if (found1 .or. found2 .or. mFuture) then
       ! DMC calculation
       mMethod='DMC';ar=.true.;rc=.true.; move=1; mWeight=WEIGHT_REY; mBranch=.true.; mJoin=.true.; mKill=.false.
       mSplit=.true.;mKillPersist=.true.;mPerMax=10;eLocalCutOff=.true.
       mLoadBalance=.false.; if (nproc>1) mLoadBalance=.true.
       mBlockLen=1000; mStepsDiscard=3000; mAutoCorrelation=.false.
       mTargetAR=0.90d0
    endif

    ! QMC parameters
    mTargetSampleSize = 0   ! default value (0) is: "keep current sample size"
    call getinta(lines,nl,'walker=',mTargetSampleSize,iflag)
    mTargetSampleSizeAllNodes = mTargetSampleSize * nproc

    ! propagate walkers (parallel) in small blocks
    call getinta(lines,nl,'walker_block=',mWalkerBlock,iflag)
    if (iflag == 0) then
       call eloc_initialize(mWalkerBlock)
       call aos_initialize(mWalkerBlock)
       call mos_initialize(mWalkerBlock)
    endif

    ! total # of QMC steps
    call getint8a(lines,nl,'steps=',mSteps,iflag)
    mShowSteps = finda(lines,nl,'show_steps')

    ! block length for statistics of correlated data
    call getint8a(lines,nl,'block_len=',mBlockLen,iflag)

    ! allow finishing a run if a given std dev is reached
    call getdbla(lines,nl,'std_dev=',mStdDev,iflag)
    if (iflag == 0) then
       mCheckStdDev = .true.
    else
       mCheckStdDev = .false.
    endif

    ! step_stride for all calls using current walkers
    call getint8a(lines,nl,'step_stride=',mStepStride,iflag)

    ! accumulation of large samples with size acc_size for optimization
    mAccumulate = finda(lines,nl,'accumulate')
    if (mAccumulate .and. mMethod /= 'VMC') call abortp('$qmc: accumulate requires a VMC run')
    call getint8a(lines,nl,'acc_size=',accSize,iflag)
    if (iflag==0) mSteps = (accSize - 1) * mStepStride / getSampleSize(sample) + mStepsDiscard

    ! discard first steps in all statistics
    call getint8a(lines,nl,'discard=',mStepsDiscard,iflag)
    found = finda(lines,nl,'discard_all')
    if (found) mStepsDiscard = mSteps

    ! allow to kill persistent walkers after mPerMax steps
    call getinta(lines,nl,'persist=',mPerMax,iflag)
    mKillPersist = .false.
    if (mPerMax > 0)  mKillPersist = .true.

    ! DMC branching
    call getloga(lines,nl,'branch=',mBranch,iflag)
    call getdbla(lines,nl,'lower_thresh=',mLThresh,iflag)
    call getdbla(lines,nl,'upper_thresh=',mUThresh,iflag)
    call getloga(lines,nl,'join=',mJoin,iflag)
    call getloga(lines,nl,'kill=',mKill,iflag)
    if (mKill .and. mJoin) call abortp('initQMCParams: join and kill exclude each other')
    call getloga(lines,nl,'split=',mSplit,iflag)

    ! DMC load balancing
    found = finda(lines,nl,'no_load_balance')
    if (found) then
       mLoadBalance = .false.
    else
       found = finda(lines,nl,'load_balance')
       if (found) mLoadBalance = .true.
    endif

    ! propagator parameters
    call getstra(lines,nl,'move=',s,iflag)
    if (iflag == 0) then
       if (s=='Rey' .or. s=='rey') move = 1
       if (s=='Umr' .or. s=='umr') move = 2
       if (s=='Two' .or. s=='two') move = 3
       if (s=='Gss' .or. s=='gss') move = 4
    endif
    if (move==2) then
       found = finda(lines,nl,'no_exp')
       if (found) call setNoExp()
    endif
    call getstra(lines,nl,'weight=',s,iflag)   ! weight overrides default values
    if (iflag == 0) then
       if (s=='none') then
          mWeight = WEIGHT_NONE
       else if (s=='Rey' .or. s=='rey') then
          mWeight = WEIGHT_REY
       else if (s=='Umr' .or. s=='umr') then
          mWeight = WEIGHT_UMR
       else if (s=='Acc' .or. s=='acc') then
          mWeight = WEIGHT_ACCEPT
       else
          call abortp('qmc_init: unknown weight given')
       endif
    endif
    if (move==2 .and. .not.(mWeight==WEIGHT_UMR .or. mWeight==WEIGHT_NONE) .or. mWeight==WEIGHT_UMR .and. .not.move==2) then
       call abortp("$qmc: Umrigar move only with Umrigar weight allowed")
    end if

    ! T moves (only with ecps)
    tmove = 0
    call getstra(lines,nl,'T_moves=',s1,iflag)
    if (iflag == 0) then
       if (s1=='simple') tmove = 1
       if (s1=='sc') tmove = 2
    end if

    if (tmove > 0) then
       tmoverej = 0
       call getinta(lines,nl,'T_move_reject=',tmoverej,iflag)
       select case (tmoverej)
          case (0)
            tmove = tmove
          case (1)
            tmove = tmove + 2
          case (2)
            tmove = tmove + 4
          case default
            call abortp("illegal value for T_move_reject given")
       end select
       tmovew = 0
       call getinta(lines,nl,'T_move_wgt=',tmovew,iflag)
       if (iflag == 0 .and. tmoverej > 0) call abortp("$qmc: T_move_weight only for T_move_reject=0")
       if (iflag == 0) then
          if (tmovew == 1) tmove = tmove + 6
       end if
       call getinta(lines,nl,'T_move_cross=',tmovecross,iflag)
       if (iflag == 0) then
          if (tmovecross == 1) tmove = tmove + 8
          if (tmovecross == 2) tmove = tmove + 10
       end if
    end if

    call getdbla(lines,nl,'wfac=',mWFac,iflag)
    E_trial = 0.d0
    call getdbla(lines,nl,'eref=',E_trial,iflag)
    call getdbla(lines,nl,'E_ref=',E_trial,iflag1)
    if (iflag /= 0 .and. iflag1 /= 0 .and. mWeight /= WEIGHT_NONE) then
       E_trial = qmc_Energy()
       if (E_trial==0) call abortp('(qmc_init): E_ref or previous vmc calc for dmc required')
    end if
    mt = "all"
    call getstra(lines,nl,'move_typ=',mt,iflag)
    call getinta(lines,nl,'blocks_tau_eff=',bgte,iflag)
    if (iflag /= 0) bgte = -1

    ! set or guess time step tau
    mAdaptTau = .true.
    call getdbla(lines,nl,'initial_tau=',tau,iflag)
    if (iflag /= 0) then
       call getdbla(lines,nl,'tau=',tau,tauFlag)
       if (tauFlag == 0) then
          mAdaptTau = .false.
        else
          tau = propagator_timeStep()
       endif
    end if
    call getdbla(lines,nl,'accept_ratio=',mTargetAR,iflag)

    ds = 1.d0
    call getdbla(lines,nl,'drift_scal=',ds,iflag)
    found = finda(lines,nl,'no_accept_step')
    if (found) then
       ar = .false.
    else
       found = finda(lines,nl,'accept_step')
       if (found) ar = .true.
    endif
    found = finda(lines,nl,'allow_cross')
    if (found) rc = .false.
    found = finda(lines,nl,'reject_cross')
    if (found) rc = .true.

    ! Rothstein/Vrbik E_local/Drift cutoff
    if (finda(lines,nl,'no_elocal_cutoff')) eLocalCutOff=.false.
    cf = 1.d0
    call getdbla(lines,nl,'elocal_cutoff=',cf,iflag)
    call setElocCutOff(eLocalCutOff,tau,cf)

    ! turn on/off autocorrelation calculation
    found = finda(lines,nl,'no_auto_corr')
    if (found) then
       mAutoCorrelation = .false.
    else
       found = finda(lines,nl,'auto_corr')
       if (found) mAutoCorrelation = .true.
    endif

    call getinta(lines,nl,'auto_corr_max=',mAutoCorrMax,iflag)
    if (iflag == 0) mAutoCorrelation = .true.

    ! Initialization of Future Walking

    if(mFuture) then
     call getinta(lines,nl,'tau2=',mTau2,iflag)
     if(iflag /=0) call abortp("qmc_future: tau2 required")
     call getinta(lines,nl,'notau2=',mNotau2,iflag)
     if(iflag /=0) call abortp("qmc_future: notau2 required")
     call getinta(lines,nl,'stprops=',mStprops,iflag)
     if(iflag /=0) call abortp("qmc_future: stprops required")
    endif

    ! Popcontrol and Init Parallel. Old Scheme

    call getinta(lines,nl,'popctrl=',mPopctrl,iflag)
    if (iflag /=0) mPopctrl = POP_GLOBAL

    if (finda(lines,nl,'old_par')) then
      mLoadBalance = .false.
      mPopctrl = POP_LOCAL
    endif

    ! Init for Stochastic Reconfiguration

    if (finda(lines,nl,'rcf')) then
      call getinta(lines,nl,'rcf=',mRcf,iflag)
      mReconf = .true.
      mBranch = .false.
      !!!mKillPersist = .false.
      !!!mPopctrl = 0
      mLoadBalance = .false.
    endif

    sOldLogmode = logmode
    call getinta(lines,nl,'verbose=',vb,iflag)
    if (iflag==0) then
       logmode = vb
    endif

    found = finda(lines,nl,'show_details')
    if (found) mShowDetails = 1

    ! energy partitioning within SED run
    if (mSED) then
        call sed_init(lines,nl)
        if (do_epart) then
           mEpartSED=.true.
           call epart_init(lines,nl)
        endif
    endif

    if (mSED .or. mCOC) then
        mAssignment = .true.
        call assign_init(lines,nl)
        if (finda(lines,nl,'readref') .or. finda(lines,nl,'ref_file')) then
            call assign_readRef(lines,nl,changed)
            if (changed.and.logmode>=2) then
               write(iul,*) "Warning: core electrons from reference have been permuted"
            end if
        else if (allocated(x_max)) then
            call assign_setRef(x_max,y_max,z_max,changed)
            if (changed.and.logmode>=2) then
               write(iul,*) "Warning set reference: core electrons from reference have been permuted"
            end if
        else
            call abortp('(qmc_init): no reference coordinates available ')
        endif
    endif

    call propagator_reset()
    blockdiscard = mStepsDiscard / mBlockLen
    call propagator_init(mWeight,move,tmove,blockdiscard,bgte,tau,ds,ar,rc,mt,mWalkerBlock)

    !Max Analysis - Check for consistent sample size on each core
    !happens if outliers are removed and not replaced on any core
    if (mMaxAnalysis) then
      call myMPIGatherInteger(getSampleSize(sample),1,int_rcv,ierr)
      if (MASTER) then
        do i=1,nproc
          if (int_rcv(i) /= getSampleSize(sample) ) then
            call abortp("(qmc_init): Different sample size on each core for maxima calculation.")
          endif
        enddo
      endif
      call myMPIBarrier(ierr)
    endif

  end subroutine qmc_init


  ! access functions
  real*8 function qmc_Energy()
     qmc_Energy = sEMean
  end function qmc_Energy

  real*8 function qmc_StdDev()
     qmc_StdDev = sEMeanStdDev
  end function qmc_StdDev

  real*8 function qmc_Variance()
     qmc_Variance = sVar
  end function qmc_Variance

  logical function qmc_cocIsConverged()
     qmc_cocIsConverged = sCOC_isConverged
  end function qmc_cocIsConverged



  subroutine qmc_writeParams(iu)
  !----------------------------!

    integer, intent(in) :: iu
    character(len=10)   :: weight
    character(len=20)   :: s
    real*8 cf,tau
    logical yn

    select case(mWeight)
    case(WEIGHT_NONE); weight="none"
    case(WEIGHT_REY); weight="Reynolds"
    case(WEIGHT_UMR); weight="Umrigar"
    case(WEIGHT_ACCEPT); weight="Accept"
    case default; call abortp("qmc_writeParams: illegal weight")
    end select

    if (mSed .and. .not. mEpartSED) then
       s = 'SED/'//mMethod
    else if (mSed .and. mEpartSED) then
       s = 'SED Epart/'//mMethod
    else
       s = mMethod
    end if
    write(iu,'(/3A/)') '   * * *  ',trim(s),' calculation  * * *'
    write(iu,'(A/)') '    QMC parameters:'

    write(iu,'(1X,A21,F12.5,3X,A21,L12)') 'tau =',propagator_timeStep(),' adapt tau =',mAdaptTau
    write(iu,'(1X,2(A21,I12,3X))') ' total walker =',mTargetSampleSizeAllNodes,    &
         ' local walker =',mTargetSampleSize
    write(iu,'(1X,A21,I12,3X,A21,I12)') ' steps =',mSteps,' discard =',mStepsDiscard
    write(iu,'(1X,2(A21,I12,3X))') ' block_len =',mBlockLen,' walker_block =',mWalkerBlock
    write(iu,'(1X,2(A21,I12,3X))') ' step_stride =',mStepStride
    if (mCheckStdDev) write(iu,'(1X,A21,F12.5)') ' target std dev =',mStdDev
    if (mAdaptTau) write(iu,'(1X,A21,F12.5)') 'target accept ratio =',mTargetAR
    write(iu,'(1X,2(A21,F12.5,3X))') 'E_ref =',E_trial,'wfac =',mWFac
    call getElocCutOff(yn,tau,cf)
    write(iu,'(1X,A21,L12,3X,A21,F12.5,3X)') 'E_loc_cutoff =',yn,'factor =',cf
    write(iu,'(1X,A21,L12,3X,A21,I12)') ' kill_persist =',mKillPersist,' max_persist =',mPerMax
    write(iu,'(2(1X,A21,L12,2X))') ' load balance =',mLoadBalance,' branch =',mBranch
    select case (mPopctrl)
    case (POP_NONE)
       write(iu,'(1X,A21,L12,3X,A21,I12)') ' future walking =',mFuture,' pop ctrl = none'
    case (POP_GLOBAL)
       write(iu,'(1X,A21,L12,3X,A21,I12)') ' future walking =',mFuture,' pop ctrl = global'
    case (POP_LOCAL)
       write(iu,'(1X,A21,L12,3X,A21,I12)') ' future walking =',mFuture,' pop ctrl = local'
    case default
       call abortp("qmc_writeParams: illegal mPopctrl value")
    end select
    write(iu,'(1X,A21,L12,3X,A21,I12)') ' Reconf =',mReconf,' RcfMethod =',mRcf
    write(iu,'(1X,A21,L12,3X,A21,I12)') 'accumulate =',mAccumulate
    write(iu,*)

    if (mFuture) call propOutput(iu)

    call propagator_writeParams(iu)

    if (mMaxSearch) call psimax_writeParams(iu)

    if (mSed) call sed_write_init(iu)


  end subroutine qmc_writeParams



  subroutine qmc_run(sample)
  !------------------------!

    type(RWSample), intent(inout) :: sample

    integer                     :: block,ps,bs,i,j,k,idx,tauFoundStep
    integer*8                   :: n,st
    real*8                      :: elocCounter,ERef,varAllNodes
    real*8                      :: EMeanAllNodes,sampleSizeAllNodes,eTotal,vTotal
    real*8                      :: ElocalStat,Emeanlocal,E,var,stddev,stdDevAllNodes
    real*8                      :: sampleWeightAllNodes,sampleWeight
    real*8                      :: Samplesize, EL,EECPL,EECPNL
    real*8                      :: tcorr,ACvar,ACNcorr
    real*8                      :: accepted ! 0.0 for (all) rejected, 1.0 for (all) accepted
    real*8                      :: wtimer1,wtimer2,ewtime,twtime
    real*8, allocatable         :: Estep(:),AC(:,:),ACRingBuffer(:,:),ACresult(:)  ! autocorrelation date
    logical                     :: first,ACavail,tauFound
    type(RandomWalker), pointer :: rwp
    type(RandomWalker), pointer :: rwbp(:)   ! pointer to walker block==array
    type(weightStat) :: stepStat          ! < E_local > (weighted)
    type(weightStat) :: testStat
!!!    type(simpleStat) :: totalAccStat      ! < acceptance ratio >
    type(simpleStat) :: blockAccStat      ! < acceptance ratio > of block
    type(simpleStat) :: adaptAccStat      ! < acceptance ratio > of some steps for adapting tau
    type(simpleStat) :: eRefStat          ! < E_ref >
    type(simpleStat) :: wgtStat           ! < total-Weight >
    type(simpleStat) :: blockStatlocal    ! for local (Procs) BlockStat
    type(stat)       :: autocorrEStat
    type(vectorstat) :: autocorrStat
    logical          :: converged
    integer,allocatable          ::  asgn(:,:)
    logical          :: t, cutoff
    type(RWSample) :: history ! used to accumulate samples
    integer :: histSize, tmove
    integer*8 :: dataAllNodes,tStep
    character(len=17)  :: s


#ifdef WTIMER
    if (MASTER) then
      wtimer1 = omp_get_wtime()
      wtimer = 0.d0
    endif
#endif

    if (mAccumulate) call internal_accumulateInit()

    if (mWalkerStatistics) call rwStatInit(mStatType)

    if (mReconf .and. mRcf>2) call qmc_reconfInit(sample,mRcf)

    if (mFuture) call future_init(sample,mStprops,mTau2,mNotau2)

    if (mTargetSampleSize==0) mTargetSampleSize = getSampleSize(sample)
    if (mTargetSampleSizeAllNodes==0) mTargetSampleSizeAllNodes = getSampleSizeAllNodes(sample)
    if (mMethod == 'VMC' .and. mTargetSampleSize /= getSampleSize(sample))  &
         call abortp('qmc_run: inconsistent walker numbers in VMC')

    if (mWeight==WEIGHT_NONE) then
       call resetWeights(sample)
       call resetPersistencies(sample)
    endif

    if (mAutoCorrelation) then
       n = getSampleSize(sample)
       allocate(Estep(n),AC(0:mAutoCorrMax,n),ACRingBuffer(mAutoCorrMax,n),ACresult(0:mAutoCorrMax))
       call autocorrStat%create(mAutoCorrMax+1)
       call autocorrEStat%create()
    end if

    call reset(sBlockStat)
    call reset(sTotalStat)
    call reset(blockStatlocal)
    call reset(eRefStat)
    call reset(wgtStat)
    call reset(adaptAccStat)
    call reset(testStat)
    mElocCut = 0; mElocCutCount = 0
    mDriftCut = 0; mDriftCutCount = 0

    if (mEPartSED .or. mSED .or. mCOC) allocate(asgn(ne,mWalkerBlock))

    sERef = E_trial
    sMode = 0           !
    tStep = 0           ! counter for total steps (no blocking) after discard (for autocorrelation)
    block = 0


    if (MASTER .and. logmode >= 2) then
       call qmc_writeParams(iul)
       if (mMaxAnalysis) then
          write(iul,'(/A)') '             step       size              <E>                   <V>     AR   #max  fmin'
          write(iul,'(A)')  ' --------------------------------------------------------------------------------------'
       else
          write(iul,'(/A)') '             step       size              <E>                   <V>     AR '
          write(iul,'(A)')  ' --------------------------------------------------------------------------'
       end if
    endif

    tauFoundStep = 0
    elocCounter = 0
    ERef = sERef
    call reset(stepStat)
    call reset(blockAccStat)

    STEPS: do st=1,mSteps

       ! loop over walker (in blocks of bs)
       idx = 0
       bs = mWalkerBlock  ! # of walker for simultaneous propagation
       call getFirstWalkerBlock(sample,bs,rwbp)
       first = .true.
       WALKERLOOP: do

          call propagateAndWeight(rwbp,ERef,sMode,accepted)
          elocCounter = elocCounter + bs

          do i=1,bs
             call addData(stepStat,E_local(rwbp(i)),wgt(rwbp(i)))
!              if (logmode >= 4) then
!                 if (wgt(rwbp(i)) > 3.d0) then
!                    call eloc_eloc_getCurrentElocData1(tphi,tu,teloc,tvpot,txdrift,tydrift,tzdrift)
!                    call eloc_eloc_getCurrentElocEpart1(tkin,tvee,tvne)
!                    call eloc_getECPContribs(tEECPL,tEECPNL)
!                    write(iull,'(2i7,g10.2,3f12.6)') st, idx+i, wgt(rwbp(i)), EL, EECPL, EECPNL
!                 end if
!              end if
             if (mAutoCorrelation) Estep(idx+i) = E_local(rwbp(i))     ! save E_L for autocorrelation
          enddo
          idx = idx + bs
          call addData(blockAccStat,accepted)
          call addData(adaptAccStat,accepted)

          if (mod(st,mStepStride)==0 .and. st > mStepsDiscard) then
             call internal_doStuff()
          endif

       if (.not. isNext(sample)) exit

          bs = mWalkerBlock
          call getNextWalkerBlock(sample,bs,rwbp)
          first = .false.

       enddo WALKERLOOP

       if (mWeight /= WEIGHT_NONE .and. mPopctrl /= POP_NONE) then
          call pop_control(sample,sampleWeightAllNodes,ERef,st)
          if (st > mStepsDiscard) then
             call addData(eRefStat,ERef)
             call addData(wgtStat,sampleWeightAllNodes)
          endif
       endif

       if (mKillPersist .and. .not. mReconf) then
          call kill_persist(sample)
          if (mMethod == 'VMC') call changeSampleSizeTo(sample,mTargetSampleSize)
       endif

       if (mShowSteps .or. logmode>=5) then
          call getSampleEnergyAndVarianceAllNodes(sample,E,var,stddev)
          if (MASTER) write(iul,'(a,2i6,g20.10,3g14.4,i8)') 'STEP:',block+1,st,E,var,stddev
       endif

       if (mBranch) call qmc_branch(sample)

       if (mReconf) call qmc_rcf(sample)

       if (mLoadBalance) call loadBalanceSamples(sample)

       if (mAdaptTau) then
          call qmc_adaptTau(adaptAccStat,tauFound)
          if (tauFound) tauFoundStep = st
       endif

       if (st > mStepsDiscard) then
          if (mFuture) call futurewlk(sample)
          tStep = tStep + 1
          if (mAutoCorrelation) call autocorrAdd()
       endif

       if (mod(st,mBlockLen)==0) then ! end of block

          EMeanlocal = mean(stepStat)
          EMeanAllNodes = meanAllNodes(stepStat)
          stdDevAllNodes = stdDevMeanAllNodes(stepStat)
          varAllNodes = varianceAllNodes(stepStat)

          if (st <= mStepsDiscard) then
             if(mPopctrl == POP_GLOBAL) then
                sERef = EMeanAllNodes
             else if (mPopctrl == POP_LOCAL) then
                sERef = EMeanlocal
             endif
          else
             sTotalStat = sTotalStat + stepStat     ! local statistics
             call addData(sBlockStat,EMeanAllNodes)
             call addData(blockStatlocal,Emeanlocal)
             if (mPopctrl == POP_GLOBAL) then
                sERef = mean(sBlockStat)
             else if (mPopctrl == POP_LOCAL) then
                sERef = mean(blockStatlocal)
             endif
          endif

          if (mWeight == WEIGHT_UMR) call propagator_setERef0(sERef)  ! for Umrigar weighting

          sampleSizeAllNodes = getSampleSizeAllNodes(sample)
          if (logmode>=2) then
             if (mMaxAnalysis .and. st>mStepsDiscard) then
                write(iul,'(i15,i10,f16.5,a4,f10.5,f10.3,f8.3,i4,f13.5)') st,nint(sampleSizeAllNodes), &
                     EMeanAllNodes,' +/-',stdDevAllNodes,varAllNodes,mean(blockAccStat), &
                     psimax_getDiffMax(),psimax_getFirstF()
             else
                write(iul,'(i15,i10,f16.5,a4,f10.5,f10.3,f8.3)') st,nint(sampleSizeAllNodes), &
                     EMeanAllNodes,' +/-',stdDevAllNodes,varAllNodes,mean(blockAccStat)
             end if
          end if

          block = st/mBlockLen
          call propagator_endOfBlock(block)  ! allow propagator to do stuff at end of block

          if (mWalkerStatistics) then
             call rwStatPrint(42,st)
             call rwStatReset()
          endif


          if (mCheckStdDev .and. st > mStepsDiscard+10*mBlockLen) then
             if (stdDevMean(sBlockStat) < mStdDev) exit STEPS
          endif

          if (mCOC .and. st > mStepsDiscard+2*mBlockLen) then
              sCOC_isConverged = coc_isConverged()
              if (sCOC_isConverged) exit STEPS
          endif

          elocCounter = 0
          ERef = sERef
          call reset(stepStat)
          call reset(blockAccStat)

       endif ! end of block

    enddo STEPS

    eTotal = meanAllNodes(sTotalStat)
    vTotal = varianceAllNodes(sTotalStat)
    dataAllNodes = dataCountAllNodes(sTotalStat)

    if (st > mStepsDiscard + mBlockLen) then

      if (mAutoCorrelation) then
         ACavail = autocorrStat%count() > 0
         if (ACavail) ACresult = autocorrStat%mean() - (autocorrEStat%mean())**2
      end if

      if (MASTER) then

        if (mAutoCorrelation .and. ACavail) then
           ACvar = ACresult(0)
           ACresult = ACresult / ACresult(0)
           do i=1,mAutoCorrMax
              if (ACresult(i) < 0.05) exit
           end do
           if (i==1) then
              ACNcorr = 1
           else if (i==mAutoCorrMax+1) then
              ACNcorr = 0     ! denotes larger than mAutoCorrMax
           else
              ! linear interpolation between R(i-1) and R(i)
              ACNcorr = i-1 + (0.05-ACresult(i-1))*1.d0/(ACresult(i)-ACresult(i-1))
           end if
           !! write preliminarily full autocorrelation data
           if (logmode >= 3) then
              write(987,*) ' full autocorrelation data:'
              write(987,*) i
              do i=0,size(ACresult)-1
                 write(987,'(i5,f12.5)') i,ACresult(i)
              end do
           end if
        end if

        ! calculation of correlation length
        ! note that the block statistic contains locally the block means of ALL nodes (see above)
        tcorr = variance(sBlockStat)/vTotal * &
                dble(dataAllNodes)/dble(datacount(sBlockStat))

        ! write final output
        write(iul,'(//2A/)') '  FINAL RESULT:'
        write(iul,'(a,f13.5,a,f8.5,a)') ' total energy                 = ',eTotal,' +/-',  &
           stdDevMean(sBlockStat),' E_h'
        write(iul,'(a,f13.5,A)') ' block average energy         = ',mean(sBlockStat),' E_h'
        write(iul,'(a,f13.5,A)') ' variance (of wave function)  = ',vTotal,' E_h^2'
        if (mWeight/=WEIGHT_NONE .and. mPopctrl /= POP_NONE) then
           write(iul,'(a,f13.5,a,f8.5,a)')       ' mean E_ref (sigma_i)         = ',mean(eRefStat),' +/-', &
                sqrt(variance(eRefStat)),' E_h'
           write(iul,'(a,f13.2,a,f8.2)')         ' mean weight (sigma_i)        = ',mean(wgtStat),' +/-', &
                sqrt(variance(wgtStat))
           write(iul,'(a,f13.2,a,f13.2)')        ' minimum weight               = ',minValue(wgtStat), &
                ' maximum weight = ',maxValue(wgtStat)
           write(iul,'(a,f13.4)')                ' tau_acc                      = ',propagator_timeStepAcc()
        endif
        if (tauFoundStep > 0) then
           write(iul,'(a,f13.4,a,i12)')          ' tau (adapted)                = ',propagator_timeStep(), &
                                                 ' fixed at step ',tauFoundStep
        endif
        if (mAutoCorrelation .and. ACavail) then
           if (ACNcorr > 0) then
              write(iul,'(a,f9.1)') ' N_corr (<5%)                 = ',ACNcorr
           else
              write(iul,'(a,i9)') ' N_corr (<5%)                 > ',mAutoCorrMax
           end if
        end if
        write(iul,'(a,f9.1)')  ' N_corr (global)              = ',tcorr
        if (mShowDetails >= 1) then
           if (propagator_getTMove() > 0) then
              write(iul,'(a,f9.5)') ' T move ratio                 = ', propagator_getTMoveRatio()
           end if
           if (mElocCutCount > 0) then
              write(iul,'(a,f9.5)') ' E_local cutoff ratio         = ', dble(mElocCut) / dble(mElocCutCount)
              write(iul,'(a,f9.5)') ' drift cutoff ratio           = ', dble(mDriftCut) / dble(mDriftCutCount)
           end if
        end if
        if (datacount(sBlockStat) < 20) then
           write(iul,'(/a/a)') ' WARNING: stddev and global N_corr may be unreliable', &
             '  (number of blocks not discarded < 20)'
        end if
        if (20*tcorr > mBlockLen) then
           write(iul,'(/a/a)') ' WARNING: stddev and global N_corr may be unreliable', &
             '  (block_len < 20*N_corr)'
        end if

      end if

    else

      if (MASTER) then
        write(iul,'(/A,F15.5,A,F15.5)') ' qmc: Emean = ',EMeanAllNodes,' var = ',varAllNodes
      endif

    endif

    if (mFuture)then
       call propCalculate()
       if (MASTER) call propPrint()
       call deallocFWArrays()
    endif

    ! save result for access with qmc_Energy etc.
    if (st > mStepsDiscard + mBlockLen) then
       sEMean = eTotal
       sEMeanStdDev = stdDevMean(sBlockStat)
       sVar = variance(sTotalStat)
    else
       sEMean = EMeanAllNodes
       sVar = varAllNodes
       sEMeanStdDev = 0
    endif
    ! save result also in global module
    call setCurrentResult(sEMean,sEMeanStdDev,sVar)

    if (mMaxSearch) then
       call psimax_writeResults()
       call psimax_destroy()
    end if

    if (mEpartSED .or. mCoC .or. mSed) then
       deallocate(asgn)
    endif
    if (mEpartSED)then
       call epart_write(iul,1)
       open(42,file=trim(basename)//'.see')
       call epart_write(42,2)
       close(42)
    endif

    if (mSED) then
       call sed_write()
       if (allocated(x_max)) deallocate(x_max,y_max,z_max)
    endif

    logmode = sOldLogmode

    if (mAccumulate) then
      ! remember the actual sample size
      mAccumulateSampleSize = getSampleSize(sample)

      ! replace the current sample with the "history" sample
      call destroySample(sample)
      sample = history

      call destroySample(history)

      n = getSampleSizeAllNodes(sample)
      if (MASTER .and. logmode >= 2) then
        write(iul,'(/a,i12)'), " sample accumulation: new total sample size is ", n
      endif
    endif

    if (mAutoCorrelation) then
      deallocate(Estep,AC,ACRingBuffer,ACresult)
      call autocorrStat%destroy()
      call autocorrEStat%destroy()
    endif

    if (mReconf .and. mRcf>2) then
      call qmc_reconfDestroy()
    endif

#ifdef WTIMER
    if (MASTER .and. logmode >= 2) then
      wtimer2 = omp_get_wtime()
      write(iul,'(/a)') " qmc wall clock times (master, in sec and relative to eloc):"
      twtime = wtimer2-wtimer1
      write(iul,'(a,f20.3)')       " total:    ",twtime
      write(iul,'(a,f20.3,f10.3)') " eloc:     ",wtimer(WTIMER_ELOC),wtimer(WTIMER_ELOC)/twtime
      ewtime = wtimer(WTIMER_ELOC)
      write(iul,'(a,f20.3,f10.3)') " phi:      ",wtimer(WTIMER_PHI),wtimer(WTIMER_PHI)/ewtime
      write(iul,'(a,f20.3,f10.3)') " jastrow:  ",wtimer(WTIMER_JAS),wtimer(WTIMER_JAS)/ewtime
      write(iul,'(a,f20.3,f10.3)') " pseudo:   ",wtimer(WTIMER_PP),wtimer(WTIMER_PP)/ewtime
      write(iul,'(a,f20.3,f10.3)') " AO:       ",wtimer(WTIMER_AO),wtimer(WTIMER_AO)/ewtime
      write(iul,'(a,f20.3,f10.3)') " MO:       ",wtimer(WTIMER_MO),wtimer(WTIMER_MO)/ewtime
      write(iul,'(a,f20.3,f10.3)') " AOMO:     ",wtimer(WTIMER_AOMO),wtimer(WTIMER_AOMO)/ewtime
      write(iul,'(a,f20.3,f10.3)') " mdet:     ",wtimer(WTIMER_MDET),wtimer(WTIMER_MDET)/ewtime
    endif
#endif

  contains

    subroutine internal_accumulateInit()
      ! if this is an accumulation run, the sample that was passed in is actual
      ! the history of previously seen samples. if this is the first accumulation
      ! run after creating a sample, store this size as the number of actual
      ! walkers used
      history = sample

      if (getSampleSize(sample) == currentSampleSize()) then
        mAccumulateSampleSize = getSampleSize(sample)
      endif

      if (MASTER .and. logmode >= 3) then
        write(iul,'(a,i8)') " accumulate run: sample copied to history. Initial (local) size: ", getSampleSize(history)
      endif

      ! the actual walkers used in this run are now only the last
      ! |mAccumulateSampleSize| samples in the history
      if (mAccumulateSampleSize < getSampleSize(sample)) then
         call reduceToLast(sample, mAccumulateSampleSize)
         if (MASTER .and. logmode >= 3) then
            write(iul,'(a,i8)') " accumulate: reduced actual sample to ", mAccumulateSampleSize
         endif
      endif

      ! pre-allocate the new history to avoid expensive reallocation in the
      ! qmc loop
      histSize = getSampleSize(history) + &
        (mSteps / mStepStride) * mAccumulateSampleSize
      call reallocateSample(history, histSize)
    end subroutine internal_accumulateInit

    subroutine internal_doStuff()

       if (mMaxSearch) then
          call psimax_opt(rwbp)
       endif

       if (mWalkerStatistics) then
          do i=1,bs
             rwp => rwbp(i)
             call rwAddToStat(rwp)
          enddo
       endif

       if (mAccumulate) then
          do i=1,bs
            call appendWalker(history, rwbp(i))
          enddo
       endif

       ! Assignment and  SED/EPart Statistics
       if (mCOC .or. mSED .or. mEpartSED) call assign_get(rwbp,asgn)
       if (mSED) call sed_add(rwbp,asgn)
       if (mEpartSED) then
          do i=1,bs
             rwp => rwbp(i)
             call epart_add(rwp,asgn(:,i),1)
          enddo
       endif
       if (mCOC) call coc_add(rwbp,asgn)

       ! collecting trajectory data
       if (mTrajectory .and. first) call display(rwbp(1),2,43)

    end subroutine internal_doStuff

    subroutine autocorrAdd()
       if (idx /= getSampleSize(sample)) call abortp("autocorrAdd: incorrect idx value")
       if (tStep > mAutoCorrMax) then
         do i=1,idx
            AC(0,i) = Estep(i)*Estep(i)
            do k=1,mAutoCorrMax
               j = mod(tStep-k-1,mAutoCorrMax) + 1
               AC(k,i) = ACRingBuffer(j,i)*Estep(i)
            end do
            j = mod(tStep-1,mAutoCorrMax) + 1
            ACRingBuffer(j,i) = Estep(i)
            call autocorrStat%add(AC(:,i))
            call autocorrEStat%add(Estep(i))
         end do
       else ! fill ring buffer
         j = mod(tStep-1,mAutoCorrMax) + 1
         do i=1,idx
            ACRingBuffer(j,i) = Estep(i)
         end do
       end if

    end subroutine autocorrAdd

  end subroutine qmc_run



  subroutine qmc_branch(sample)
  !---------------------------!

  ! different branching schemes for current walker sample

    real*8, parameter :: KILL_THRESH=1.d-12
    real*8, parameter :: SAMPLE_SIZE_FACTOR=0.9d0

    type(RWSample), intent(inout) :: sample

    integer  ps,sampleSize,nsplit,i,maxSize
    real*8   wNew,xi,w,wJoin
    type(RandomWalker), pointer :: rwp,rwpJoin


    ! Kill walkers with weight < KILL_THRESH
    rwp => getFirst(sample)
    do
       if (wgt(rwp) < KILL_THRESH) then
          call deleteCurrentWalker(sample)
          if (.not.isValid(sample)) exit     ! true after deleting last element
       else
          if (.not.isNext(sample)) exit
          rwp => getNext(sample)
       endif
    enddo

    if (mJoin) then  ! join two walkers with weight < mLThresh (Umrigar, 1993)

       rwpJoin => null()                     ! index of 1st walker
       rwp => getFirst(sample)                ! oder "=>" ???
       do
          w = wgt(rwp)
          if (w < mLThresh) then
             if (.not.associated(rwpJoin)) then          ! 1st walker for joining
                wJoin   = w
                rwpJoin  => rwp
                if (.not.isNext(sample)) exit
                rwp => getNext(sample)
             else                            ! 2nd found
                xi = myran()
                wNew = w + wJoin
                if (xi < w/wNew) then
                   call setWeight(rwp,wNew)
                   rwpJoin = rwp
                else
                   call setWeight(rwpJoin,wNew)
                endif
                call deleteCurrentWalker(sample)
                if (.not.isValid(sample)) exit
                rwpJoin => null()
             endif
          else
             if (.not.isNext(sample)) exit
             rwp => getNext(sample)
          endif
       enddo

    else if (mKill) then

       rwp => getFirst(sample)
       do
          w = wgt(rwp)
          if (w < mLThresh) then
             xi = myran()
             wNew = 2.d0*mLThresh
             if (xi < w/wNew) then
                ! walker survives
                call setWeight(rwp,wNew)
                if (.not.isNext(sample)) exit
                rwp => getNext(sample)
             else
                call deleteCurrentWalker(sample)
                if (.not.isValid(sample)) exit
             endif
          else
             if (.not.isNext(sample)) exit
             rwp => getNext(sample)
          endif
       enddo

    endif

    if (mSplit) then

       sampleSize = getSampleSize(sample)
       maxSize    = SAMPLE_SIZE_FACTOR*getMaxSampleSize(sample)
       rwp => getFirst(sample)
       ps = 0
       do
          w = wgt(rwp)
          if (w > mUThresh) then
             nsplit = int(w)
             ! do not split if maxSize reached
             if (logmode >= 3 .and. nsplit>3) write(iul,'(A,I6)') 'BRANCH: nsplit=',nsplit
             if (getSampleSize(sample)+nsplit-1 > maxSize) then
                write(999,*) "WARNING:qmc_branch: maxSize hit",getMyTaskId(),getSampleSize(sample)+nsplit-1,maxSize
                exit
             end if
             wNew = w/nsplit
             call setWeight(rwp,wNew)
             do i=1,nsplit-1
                call appendWalker(sample,rwp)
             enddo
          endif
          ps = ps +1
          if (ps == sampleSize) exit
          rwp => getNext(sample)
       enddo

    endif

  end subroutine qmc_branch



  subroutine qmc_rcf(sample)
  !-------------------------!
  ! reconfiguration of sample

  type(RWSample), intent(inout) :: sample
  type(RandomWalker), pointer :: rwp,rwpJoin

  select case (mRcf)
    case (1)
      call qmc_reconf1(sample)
    case (2)
      call qmc_reconf2(sample)
    case (3:4)
      call qmc_reconfNew(sample)
    case (5:6)
      ! kill persistent walkers by setting weight to zero
      rwp => getFirst(sample)
      do
        if (persist(rwp) > mPerMax) call setWeight(rwp,0.d0)
        if (.not.isNext(sample)) exit
        rwp => getNext(sample)
      enddo
      call qmc_reconfNew(sample)
    case default
      call abortp("qmc_rcf: illegal value for mRcf")
  end select

  end subroutine qmc_rcf


  subroutine pop_control(sample,sampleWeightAllNodes,ERef,st)
  !---------------------------------------------------------!

  ! modify reference energy for population control
  ! Popcontrl global = 1 ; local = 2

  type(RWSample), intent(inout) :: sample
  real*8,intent(out)            :: ERef
  integer*8,intent(in)          :: st
  real*8                        :: Samplesize
  real*8,intent(out)            :: sampleWeightAllNodes
  real*8                        :: sampleWeight


    if(mPopctrl == POP_GLOBAL) then
       Sampleweight = getSampleWeightAllNodes(sample)
       Samplesize   = mTargetSampleSizeAllNodes
       sampleWeightAllNodes = getSampleWeightAllNodes(sample)
    elseif(mPopctrl == POP_LOCAL) then
       Sampleweight = getSampleWeight(sample)
       Samplesize   = mTargetSampleSize
       sampleWeightAllNodes = getSampleWeightAllNodes(sample)
    endif
    ERef = sERef - mWFac*log(Sampleweight/Samplesize) !
    if (logmode >= 3) then
       write(iul,'(A,2F15.4,I6,2F15.4,I15)') 'POP:',Sampleweight,getSampleWeight(sample),   &
             getSampleSize(sample),ERef,sERef,st
       call flush(iul)
    endif
    if (sampleWeight > 2*Samplesize) then
       write(iul + mytid,'(A,G13.3)') ' sample overflow: all weights=',Sampleweight
       write(iul + mytid,'(A,2G18.6)') ' sample overflow: ERef,sERef=',ERef,sERef
       call displaySample(sample,iul + mytid)
       call myMPIBarrier()
       call abortp('QMC: walker overflow')
    endif

  end subroutine pop_control


  subroutine kill_persist(sample)
  !-----------------------------!
    ! kill all persistent walkers from sample
    type(RWSample), intent(inout) :: sample

    type(RandomWalker), pointer :: rwp

    rwp => getFirst(sample)
    do
       if (persist(rwp) > mPerMax) then
          call deleteCurrentWalker(sample)
          if (.not.isValid(sample)) exit
       else
          if (.not.isNext(sample)) exit
          rwp => getNext(sample)
       endif
    enddo

  end subroutine kill_persist



  subroutine qmc_adaptTau(accStat,tauFound)
  !---------------------------------------!

    ! adapting time step tau to the target acceptance ratio (AR)
    ! as accumulated in accStat
    ! note: adaptation is based on the model AR = exp(-a*tau)
    ! or tau = - log(AR)/a
    integer, parameter             :: minAdaptSteps = 500
    type(simpleStat),intent(inout) :: accStat      ! < acceptance ratio > of some steps for adapting tau
    logical, intent(out)           :: tauFound
    real*8 accRatio,oldTau,newTau

    tauFound = .false.

    if (dataCountAllNodes(accStat) < minAdaptSteps) return

    accRatio = meanAllNodes(accStat)

    ! adapt tau even if convergence reached
    oldTau = propagator_timeStep()
    newTau = oldTau * log(mTargetAR) / log(accRatio)
    call propagator_setTimeStep(newTau)

    if (abs(accRatio-mTargetAR) < mEpsAR) then
       mAdaptTau = .false.
       tauFound = .true.
    else
       call reset(accStat)
    endif

  end subroutine qmc_adaptTau



  subroutine initWalkerStat(lines,nl)
  !---------------------------------!

  ! read and initialize $walker_stat section

  integer                     :: nl
  character(len=120)          :: lines(nl)
  integer                     :: iflag
  integer                     :: blockStride
  integer                     :: stepStride
  character(len=1)            :: statType
  logical                     :: finda

  if (nproc > 1) call abortp('walkerstat only in serial runs')
  statType = 'c'
  if (finda(lines,nl,'spherical')) statType = 's'

  mWalkerStatistics = .true.
  mStatType = statType

  open(42,file=trim(baseName)//'.reg')

  end subroutine initWalkerStat



  subroutine initTrajectory(lines,nl)
  !---------------------------------!

  ! read and initialize $gen section

  use global
  integer                     :: nl
  character(len=120)          :: lines(nl)
  integer                     :: iflag
  integer                     :: blockStart,stepStride
  logical                     :: trajectory

  if (nproc > 1) call abortp('trajectory only in serial runs')
  open(43,file=trim(baseName)//'.trj')

  mTrajectory = .true.

  call trajectory_head(43)

  end subroutine initTrajectory

end MODULE qmc
