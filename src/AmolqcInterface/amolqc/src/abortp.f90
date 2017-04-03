!
! keep abortp outside module for compatability
!
  !-------------------
  subroutine abortp(s)
  !-------------------

    ! general abnormal thread termination

    implicit none
    
    character(*), intent(in) :: s
    integer ierr,iul,mytid,getMytid,getiul

    mytid = getMytid()
    iul = getiul()
    write(iul,*) ' *** ABNORMAL TERMINATION *** '
    write(iul,*) ' from rank ',mytid
    write(iul,*) s
    
    call myMPIabort(ierr)
    call myMPIfinalize(ierr)

    STOP
  end subroutine abortp



