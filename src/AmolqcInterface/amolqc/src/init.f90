!
! init.f90: contains initAmolqc finalizeAmolqc and parseInputFile
!
! $Id: init.f90,v 1.2 2008/02/05 22:21:35 luechow Exp $
!
! $Log: init.f90,v $
! Revision 1.2  2008/02/05 22:21:35  luechow
! minor bugs removed after compilation with sun compiler.
! parselib.f adapted to fortran standard
!
! Revision 1.1.1.1  2007/04/25 13:42:20  luechow
! QMC program amolqc. rewritten in 2006. AL
!
!

MODULE InitModule

  use global
  use MachModule, only: myGetArgLocal,myGetDate,myGetHost,mygetEnv
  use utilsmodule, only: init_ran
  use Utils, only: readFileParallel
  use versionModule, only: getVersion
  implicit none

CONTAINS

!-----------------------!
subroutine initAmolqc()
!-----------------------!

  !character(len=*), intent(in) :: vs
  character(len=40) :: name
  character(len=180) :: path
  integer ierr,nn

  call myMPIInitialize(ierr)
  call getVersion(version)
  call getEnvironmentVariableName(name)
  call myGetEnv(name,path)
  call setAmolqcPath(path)
  call initSaveResults()

end subroutine initAmolqc


!-------------------------!
subroutine finalizeAmolqc()
!-------------------------!


  integer ierr
  character(len=26) :: date

  if (MASTER) then
     call myGetDate(date)
     write(iul,'(/2A//A//)') 'Amolqc run finished on ',date,'Bye!'
     close(iul)
  endif

  ! Properly finish MPI
  call myMPIfinalize(ierr)
end subroutine finalizeAmolqc



!--------------------------------!
subroutine readInputFile(lines,nl)
!--------------------------------!


  character(len=*), intent(inout) :: lines(:)
  integer, intent(out)            :: nl
  character(len=40)               :: date,host
  character(len=3)                :: str
  character(len=180)              :: path
  integer counter
  logical fileExists
  character(len=1) cntstr

  integer id,io,iflag,ierr,nn

  nl = 0

  call myGetArgLocal(1,baseName,ierr)
  if (MASTER) then      ! only MASTER reads input file
      if (ierr /= 0) call abortp(" amolqc requires a command line argument (base name of .in file)")
      iul=9
      iull=100
      inquire(file=trim(baseName)//'.out',exist=fileExists)
      if (fileExists) then
         ! do not overwrite but find new out file name
         do counter=1,9
            write(cntstr,'(I1)') counter
            inquire(file=trim(baseName)//'.out-'//cntstr,exist=fileExists)
            if (.not.fileExists) exit
         end do
         open(iul,file=trim(baseName)//'.out-'//cntstr,status='new',iostat=io)
         call assert(io==0,'(readInputFile): opening output file '//trim(baseName)//'.out-'//cntstr//' failed')
      else
         open(iul,file=trim(baseName)//'.out',status='new',iostat=io)
         call assert(io==0,'(readInputFile): opening output file '//trim(baseName)//'.out failed')
      end if

      inquire(file=trim(baseName)//'.in',exist=fileExists)
      if (.not.fileExists) then
         call abortp('input file '//trim(baseName)//'.in not found')
      end if

      call myGetDate(date)
      call myGetHost(host)

      write(iul,'(/)')
      write(iul,*) '             __  __    ____    _         ____     _____   '
      write(iul,*) '     /\     |  \/  |  / __ \  | |       / __ \   / ____|  '
      write(iul,*) '    /  \    | \  / | | |  | | | |      | |  | | | |       '
      write(iul,*) '   / /\ \   | |\/| | | |  | | | |      | |  | | | |       '
      write(iul,*) '  / ____ \  | |  | | | |__| | | |____  | |__| | | |____   '
      write(iul,*) ' /_/    \_\ |_|  |_|  \____/  |______|  \___\_\  \_____|  '

      write(iul,'(//A/)') ' Atoms and Molecules with Quantum Monte Carlo -- electron structure code'
      write(iul,'(A)') ' initial version:'
      write(iul,'(A/)') '  Arne Luechow, Penn State University, 2/1996'
      write(iul,'(A)') ' main author:'
      write(iul,'(A/)') '  Arne Luechow, RWTH Aachen University, 52056 Aachen, Germany'
      write(iul,'(A)') ' with contributions from:'
      write(iul,'(A)') '  James B. Anderson'
      write(iul,'(A)') '  Sebastian Manten, Christian Diedrich, Annika Bande, Tony Scott,'
      write(iul,'(A)') '  Rene Petz, Raphael Berner, Alexander Sturm, Kaveh Haghighi Mood'
      write(iul,'(//2A/)') ' version ',version
      write(iul,'(/5A,I4,A)') ' run started on ',trim(host),' at ',trim(date), &
          ' on ',nproc,' processor(s)'
      call getAmolqcPath(path)
      write(iul,'(2A)') ' using path: ',trim(path)

  else
      iull = 100+getMyTaskId()
      iul = iull
  endif

  call readFileParallel(mytid,trim(baseName)//'.in',lines,nl)

end subroutine readInputFile



!--------------------------!
subroutine initGen(lines,nl)
!--------------------------!

  ! read and initialize $gen section

  integer                     :: nl
  character(len=120)          :: lines(nl)
  integer                     :: seed,iflag

  real*8                      :: dummy

  ! get seed and initialize random number generator

  call getinta(lines,nl,'seed=',seed,iflag)
  if (iflag /= 0 .or. seed <= 0) call abortp('$gen: seed > 0 required')
  seed = abs(seed) + mytid
  dummy = init_ran(seed)

  call getinta(lines,nl,'verbose=',logmode,iflag)
  if (.not. MASTER .and. logmode<3) logmode=0
  if (logmode>=2) then
     write(iul,'(/A/)') ' =======>      $gen - initializing RNG and setting general parameters       <======='
     write(iul,'(A,I7,A,I2)') ' seed =',seed,'     verbose level =',logmode
  endif

end subroutine initGen

end module InitModule

