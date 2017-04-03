program test

  use statistics
  implicit none
  include 'mpif.h'
  integer, parameter :: MASTER=0
  integer ierr,mytid,nproc,i
  real*8  meanAll
  type(simpleStat)  :: stat

  ! MPI initialization
  call mpi_init(ierr)
  call mpi_comm_rank(MPI_COMM_WORLD, mytid, ierr)
  call mpi_comm_size(MPI_COMM_WORLD, nproc, ierr)


  call reset(stat)

  if (mytid==MASTER) print*,' running on ',nproc,' nodes'

  do i=1,100
     call addData(stat,mytid+1.0d-3*dble(i))
     meanAll = meanAllNodes(stat)
     if (mytid==MASTER) print*,meanAll
  enddo

  write(*,*) mytid,mean(stat),dataCount(stat)

end program test
