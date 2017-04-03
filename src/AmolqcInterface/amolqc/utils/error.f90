
!------  error handling routines

! free subroutines (to avoid necessity to include modules in all simple programs)


  !-----------------------------------
  subroutine assert(logicalExp,string)
  !-----------------------------------
    
    logical, intent(in)                    :: logicalExp
    character(len=*), optional, intent(in) :: string
    integer ierr
    character(len=72) :: errstring

    if (.not. logicalExp) then
       if (present(string)) then
          errstring = 'Assertion violated: '//string
       else
          errstring = 'Assertion violated!'
       endif
       call abortp(errstring)
    end if

  end subroutine assert

  !----------------------------------------------!
  subroutine assertEqualAbsolute(t1,t2,tol,string)
  !----------------------------------------------!

  real*8, intent(in)                     :: t1
  real*8, intent(in)                     :: t2
  real*8, intent(in)                     :: tol
  character(len=*), optional, intent(in) :: string
  character(len=72) :: errstring

  if (abs(t1-t2) > tol) then
     if (present(string)) then
        errstring = 'Assertion violated: '//string
     else
        errstring = 'Assertion violated!'
     endif
     call abortp(errstring)
  end if

  end subroutine assertEqualAbsolute

  !----------------------------------------------!
  subroutine assertEqualRelative(t1,t2,tol,string)
  !----------------------------------------------!

  real*8, intent(in)                     :: t1
  real*8, intent(in)                     :: t2
  real*8, intent(in)                     :: tol
  character(len=*), optional, intent(in) :: string
  character(len=72) :: errstring

  if (abs((t1-t2)/t1) > tol) then
     if (present(string)) then
        errstring = 'Assertion violated: '//string
     else
        errstring = 'Assertion violated!'
     endif
     call abortp(errstring)
  end if

  end subroutine assertEqualRelative


  !-----------------!
  subroutine error(s)
  !-----------------!
    character(len=*), intent(in) :: s
    integer ierr
    
    call abortp(s)
    
  end subroutine error


