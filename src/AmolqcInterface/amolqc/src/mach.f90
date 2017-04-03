!
!
!
! collection of machine-dependent utility routines
!
! adapt these to your compiler
!
! $Id: mach.f90,v 1.1.1.1 2007/04/25 13:42:20 luechow Exp $
!
! $Log: mach.f90,v $
! Revision 1.1.1.1  2007/04/25 13:42:20  luechow
! QMC program amolqc. rewritten in 2006. AL
!
!
module MachModule

contains

!---------------------------!
subroutine myGetArg(n,s,ierr)
!---------------------------!

! master reads arg broadcasts to other procs 
  use global
  implicit none
  integer, intent(in)              :: n    ! get n-th command line argument
  character(len=*), intent(out)    :: s    ! return argument
  integer, intent(out)             :: ierr ! error code: 1=no argument
  integer iargc


  ierr = 0
  if (mytid==0) then
     if (iargc() < n) then
        ierr = 1
     else
        call getarg(n,s)
     endif
  endif
  call myMPIBcastInteger(ierr,1)
  if (ierr == 0) call myMPIBcastString(s,len(s))

end subroutine myGetArg
  
!--------------------------------!
subroutine myGetArgLocal(n,s,ierr)
!--------------------------------!

  use global
  implicit none
  integer, intent(in)              :: n    ! get n-th command line argument
  character(len=*), intent(out)    :: s    ! return argument
  integer, intent(out)             :: ierr ! error code: 1=no argument
  integer iargc

  ierr = 0
  if (iargc() < n) then
     ierr = 1
  else
     call getarg(n,s)
  endif

end subroutine myGetArgLocal
  

!--------------------------
subroutine myGetHost(hostn)
!--------------------------
  ! returns the hostname in 'hostn'
  implicit none
  character(len=*) :: hostn

  call getenv('HOSTNAME',hostn)
  !call hostnm(hostn)    ! SGI, Alpha
  !ihost = hostnm_(hostn)   ! AIX
end subroutine myGetHost


!-----------------------------
subroutine myGetEnv(name,value)
!-----------------------------
  ! returns content of environment variable name
  use global
  implicit none
  character(len=*) name, value

  if (mytid==0) call getenv(name,value)
  call myMPIBcastString(value,len(value))

end subroutine myGetEnv

!--------------------------
subroutine myGetDate(date)
!--------------------------
  !returns current time and date in 'date'
  implicit none

  character(len=*) date
  call fdate(date)
  !call fdate_(date)     ! AIX
end subroutine myGetDate

end module MachModule
