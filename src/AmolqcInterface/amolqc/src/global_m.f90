

MODULE global

! some global constants and data types for amolqc program
!
! $Id: global_m.f90,v 1.1.1.1 2007/04/25 13:42:20 luechow Exp $
!
! $Log: global_m.f90,v $
! Revision 1.1.1.1  2007/04/25 13:42:20  luechow
! QMC program amolqc. rewritten in 2006. AL
!


  use utilsmodule

! global constants

!!!  integer,parameter :: psmax=10000      ! max. no. of Random Walkers
!!!  integer,parameter :: nthrmax=64       ! max. no. of threads (parallel)
!!!  integer,parameter :: pwmax=40         ! max. additional weights/Elocs
!!!  integer,parameter :: evmax=5        ! max. saves for eigenvalue calc.
!!!  integer,parameter :: mpiArrayMax=1000 ! size of MPI-Buffer

  type ElecConfig
    real*8, allocatable :: x(:)
    real*8, allocatable :: y(:)
    real*8, allocatable :: z(:)
  end type ElecConfig

  type Coord
   real*8, pointer    :: x(:) => null()
   real*8, pointer    :: y(:) => null()
   real*8, pointer    :: z(:) => null()
  end type Coord

  integer,parameter :: nmax=202            ! max. # of dimensions(=elecs)
  integer,parameter :: ndmax=nmax/2+5      ! max. size of determinants
  integer,parameter :: MaxMPINodes = 2048  ! max. # of MPI-nodes
  integer,parameter :: MAXSUBLINES=50      ! max. # of lines in subroutines in .in script
  integer,parameter :: MAXLEN=120          ! max. line length of .in script lines
  integer,parameter :: MAXSUBS=5           ! max. of subroutine in .in script
  real*8,parameter  :: Pi= 3.1415926535897932d0
  real*8,parameter  :: bohr2angs = 0.529177249d0
  real*8,parameter  :: debye=2.541765d0    ! au -> Debye
  integer, parameter :: MAXSIZE=100        ! array size for saving E and var
  integer, parameter :: WTIMER_ELOC=1, WTIMER_PHI=2, WTIMER_AO=3, WTIMER_MO=4, WTIMER_AOMO=5, &
                        WTIMER_MDET=6, WTIMER_PP=7, WTIMER_JAS=8


! global variables
  integer      :: ne=0         ! number of electrons
  integer      :: nproc=1      ! number of parallel processors (1: serial)
  integer      :: mytid=0      ! task id (0 master or serial, 1..nproc-1)
  integer      :: logmode=2    ! output verbosity (0=none,1,2,..6)
  integer      :: iul=6        ! defaults to stdout (machine dependent)
  integer      :: iull=100     ! logging unit for node specific output
  integer      :: loopIter = 1 ! current loop index for loop in input file
  integer      :: optIter = 0  ! total number of cycles in all optimizations
  integer      :: savedData = 0! # of data in save storage
  integer      :: subLen(MAXSUBS) ! actual line length of subroutines
  real*8       :: E_trial=0    ! (initial) trial energy
  real*8       :: timesum1,timesum2,timesum3,timesum4,timesum5,timesum6
  real*8       :: currEnergy,currVar,currStdDev
  real*8       :: saveResults(3,MAXSIZE)
  real*8       :: wtimer(8)=0
  character*30 :: baseName=""  ! base name of .in file
  character*60 :: version=""   ! version number of program
  character*20 :: envvar="AMOLQC" ! environment variable pointing to installation
  character*180:: amolqcpath="toset"   ! path returned from envvar
  character(len=MAXLEN)       :: subLines(MAXSUBLINES,MAXSUBS)=''
  character(len=120)          :: subNames(MAXSUBS)=''
  logical      :: MASTER=.true.! .true. if MPI Master or serial
  logical      :: PARALLEL_RUN=.false. ! .t. if parallel run (nproc > 1)
  logical      :: useLAlib = .true.           ! use external BLAS LAPACK
  logical      :: maxDomainSampling = .false.

  private :: envvar,amolqcpath,loopIter,currEnergy,currVar,currStdDev,MAXSIZE,saveResults,savedData,optIter

contains

  integer pure function getNprocs()
     getNprocs = nproc
  end function

  subroutine setNprocs(n)
     integer, intent(in) :: n
     nproc = n
  end subroutine

  integer pure function getMyTaskId()
     getMyTaskId = mytid
  end function

  subroutine getMyTaskIdChar(mytaskid)
     character(len=*), intent(inout) :: mytaskid
     character(len=8) fs
     integer l
     l = len(mytaskid)
     fs = "(i"//char(48+l)//"."//char(48+l)//")"
     write(mytaskid,fs) mytid
  end subroutine

  subroutine setMyTaskId(n)
     integer, intent(in) :: n
     mytid = n
  end subroutine

  integer pure function getNElec()
     getNElec = ne
  end function

  subroutine setNElec(n)
     integer, intent(in) :: n
     ne = n
  end subroutine

  subroutine setEnvironmentVariableName(name)
     character(len=*), intent(in) :: name
     envvar = name
  end subroutine setEnvironmentVariableName

  subroutine getEnvironmentVariableName(name)
     character(len=*), intent(inout) :: name
     name = envvar
  end subroutine getEnvironmentVariableName

  subroutine setAmolqcPath(path)
     character(len=*), intent(in) :: path
     if (len(trim(path)) > len(amolqcpath)) then
        call abortp("global: path length exceeds amolqcpath def")
     endif
     amolqcpath = path
  end subroutine setAmolqcPath

  subroutine getAmolqcPath(path)
     character(len=*), intent(inout) :: path
     if (amolqcpath=="toset") then
        call abortp("global: amolqc path not set")
     endif
     if (len(trim(amolqcpath)) > len(path)) then
        call abortp("global: path length too small for amolqcpath def")
     endif
     path = amolqcpath
  end subroutine getAmolqcPath

  subroutine setCurrentLoopIdx(idx)
    integer, intent(in) :: idx
    loopIter = idx
  end subroutine setCurrentLoopIdx

  pure integer function getCurrentLoopIdx()
    getCurrentLoopIdx = loopIter
  end function getCurrentLoopIdx

  subroutine setOptIter(idx)
    integer, intent(in) :: idx
    optIter = idx
  end subroutine setOptIter

  subroutine getPlusOptIter(idx)
    integer, intent(out) :: idx
    optIter = optIter+1
    idx=optIter
  end subroutine getPlusOptIter

  pure integer function getOptIter()
    getOptIter = optIter
  end function getOptIter

  subroutine setCurrentResult(E,sigma,V)  ! do not allow setting individual values
    real*8, intent(in) :: E,sigma,V
    currEnergy = E
    currStdDev = sigma
    currVar = V
  end subroutine setCurrentResult

  pure real*8 function getCurrentEnergy()
    getCurrentEnergy = currEnergy
  end function getCurrentEnergy

  pure real*8 function getCurrentVariance()
    getCurrentVariance = currVar
  end function getCurrentVariance

  pure real*8 function getCurrentStdDev()
    getCurrentStdDev = currStdDev
  end function getCurrentStdDev

  subroutine initSaveResults()
    saveResults = 0
    savedData = 0
  end subroutine initSaveResults

  subroutine global_saveResult(lines,nl)
    character(len=120), intent(in) :: lines(:)
    integer, intent(in)           :: nl
    integer idx,n,iflag
    call getinta(lines,nl,'idx=$idx+',n,iflag)
    if (iflag == 0) then
       idx = getCurrentLoopIdx() + n
    else if (finda(lines,nl,'idx=$idx')) then
       idx = getCurrentLoopIdx()
    else
       idx = savedData + 1
       call getinta(lines,nl,'idx=',idx,iflag)
    end if
    if (idx > MAXSIZE) call abortp("addSaveResult: exceeding result storage size")
    saveResults(1,idx) = getCurrentEnergy()
    saveResults(2,idx) = getCurrentStdDev()
    saveResults(3,idx) = getCurrentVariance()
    savedData = max(idx,savedData)
  end subroutine global_saveResult

  subroutine global_printSavedResults(lines,nl)
    character(len=120), intent(in) :: lines(:)
    integer, intent(in)           :: nl
    integer i,i1,i2,iflag
    if (MASTER) then
       i1 = 1
       call getinta(lines,nl,'start=',i1,iflag)
       i1 = min(i1,savedData); i1 = max(i1,1)
       i2 = savedData
       call getinta(lines,nl,'end=',i2,iflag)
       i2 = min(i2,savedData); i2 = max(i2,1)

       write(iul,'(/a/)') 'table of results:'
       write(iul,*) '  i          E            stddev           var'
       write(iul,*) '-------------------------------------------------------'
       do i=i1,i2
          write(iul,'(i5,2f15.5,f15.3)') i,saveResults(1,i),saveResults(2,i),saveResults(3,i)
       end do
       write(iul,'(a/)') '-------------------------------------------------------'
     end if
  end subroutine global_printSavedResults

  logical function global_exitIf(lines,nl)
  !--------------------------------------!
    character(len=120), intent(in) :: lines(:)
    integer, intent(in)           :: nl
    real*8 emax,emin,vmax,vmin
    integer iflag

    global_exitIf = .false.

    call getdbla(lines,nl,'energy>',emax,iflag)
    if (iflag==0) then
       global_exitIf = global_exitIf .or. getCurrentEnergy() > emax
    end if

    call getdbla(lines,nl,'energy<',emin,iflag)
    if (iflag==0) then
       global_exitIf = global_exitIf .or. getCurrentEnergy() < emin
    end if

    call getdbla(lines,nl,'variance>',vmax,iflag)
    if (iflag==0) then
       global_exitIf = global_exitIf .or. getCurrentVariance() > vmax
    end if

    call getdbla(lines,nl,'variance<',vmin,iflag)
    if (iflag==0) then
       global_exitIf = global_exitIf .or. getCurrentVariance() < vmin
    end if

  end function global_exitIf


end module global

integer pure function getiul()
  use global
  implicit none
  getiul = iul
end function getiul


