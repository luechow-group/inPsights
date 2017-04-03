!
! keep abortp outside module for compatability
!
  !-------------------
  subroutine abortp(s)
  !-------------------

#ifdef PARALLEL
    !!!use mpi
#endif

    ! general abnormal thread termination

    implicit none
#ifdef PARALLEL
    include 'mpif.h'
#endif
    
    character(*), intent(in) :: s
    integer iexit,ierr
    write(*,*) ' *** ABNORMAL TERMINATION *** '
    write(*,*) s
#ifdef PARALLEL
    iexit = 1
    call MPI_ABORT(MPI_COMM_WORLD,iexit,ierr)
#endif
    STOP
  end subroutine abortp



