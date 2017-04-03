
! collection of routines that encapsulate MPI

! provide 2 implementations: one with MPI call, the other with dummy routines
! enable MPI version by linking

! for speed let compiler inline these routines

! Arne Luechow, 2006

integer function getMytid()
  use global
  implicit none

  getMytid = mytid
end function getMytid


subroutine myMPIInitialize(ierr)
  use global
  implicit none
  include 'mpif.h'
  integer ierr,nn

  ! MPI initialization
  call mpi_init(ierr)
  call mpi_comm_rank(MPI_COMM_WORLD, mytid, ierr)
  call mpi_comm_size(MPI_COMM_WORLD, nproc, ierr)
  call assert(nproc <= MaxMPINodes,'mpiMPIInitialize: recompile with increased MaxMPINodes')
  if (mytid > 0) MASTER = .false.
  if (nproc > 0) PARALLEL_RUN = .true.
end subroutine myMPIInitialize


subroutine myMPIFinalize(ierr)
  use global
  implicit none
  include 'mpif.h'
  integer ierr

  ! Properly finish MPI
  call mpi_finalize(ierr)
end subroutine myMPIFinalize


subroutine myMPIBarrier(ierr)
  use global
  implicit none
  include 'mpif.h'
  integer ierr

  call mpi_barrier(MPI_COMM_WORLD,ierr)
end subroutine myMPIBarrier


subroutine myMPIBcastString(string,n)
  use global
  implicit none
  include 'mpif.h'
  integer n,ierr
  character*(*) string

  call mpi_bcast(string,n,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)

end subroutine myMPIBcastString


subroutine myMPIBcastInteger(array,n)
  use global
  implicit none
  include 'mpif.h'
  integer n,ierr
  integer array(n)

  call mpi_bcast(array,n,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

end subroutine myMPIBcastInteger


subroutine myMPIBcastDouble(array,n)
  use global
  implicit none
  include 'mpif.h'
  integer n,ierr
  double precision array(n)

  call mpi_bcast(array,n,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

end subroutine myMPIBcastDouble


subroutine myMPIReduceSumDouble(sendbuf,recvbuf,n)
  ! summing all arrays into recvbuf on master (id=0)
  use global
  implicit none
  include 'mpif.h'
  integer n,ierr
  real*8 sendbuf(n),recvbuf(n)

  call mpi_reduce(sendbuf,recvbuf,n,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)

end subroutine myMPIReduceSumDouble


subroutine myMPIReduceSumInteger(sendbuf,recvbuf,n)
  ! summing all arrays into recvbuf on master (id=0)
  use global
  implicit none
  include 'mpif.h'
  integer n,ierr
  integer sendbuf(n),recvbuf(n)

  call mpi_reduce(sendbuf,recvbuf,n,MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,ierr)

end subroutine myMPIReduceSumInteger


subroutine myMPIAllReduceSumDouble(sendbuf,recvbuf,n)
  use global
  implicit none
  include 'mpif.h'
  integer n,ierr
  real*8 sendbuf(n),recvbuf(n)

  call mpi_allreduce(sendbuf,recvbuf,n,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)

end subroutine myMPIAllReduceSumDouble


subroutine myMPIAllReduceSumInteger(sendbuf,recvbuf,n)
  use global
  implicit none
  include 'mpif.h'
  integer n,ierr
  integer sendbuf(n),recvbuf(n)

  call mpi_allreduce(sendbuf,recvbuf,n,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)

end subroutine myMPIAllReduceSumInteger


subroutine myMPIGatherInteger(sendbuf,n,recvbuf,ierr)
  use global
  implicit none
  include 'mpif.h'

  integer n,ierr
  integer sendbuf(n),recvbuf(MaxMPINodes*n)

  call MPI_GATHER(sendbuf,n,MPI_INTEGER,recvbuf,n,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

end subroutine myMPIGatherInteger


subroutine myMPIScatterInteger(sendbuf,n,recvbuf,ierr)
  use global
  implicit none
  include 'mpif.h'

  integer n,ierr
  integer sendbuf(MaxMPINodes*n),recvbuf(n)

  call MPI_SCATTER(sendbuf,n,MPI_INTEGER,recvbuf,n,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

end subroutine myMPIScatterInteger


subroutine myMPIGatherDouble(sendbuf,n,recvbuf,ierr)
  use global
  implicit none
  include 'mpif.h'

  integer n,ierr
  real*8 sendbuf(n),recvbuf(MaxMPINodes*n)

  call MPI_GATHER(sendbuf,n,MPI_DOUBLE_PRECISION,recvbuf,n,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

end subroutine myMPIGatherDouble

! simple wrapper for MPI_GATHERV, which gathers a different amount of data from
! each process. this wrapper assumes that the data should be received in order,
! i.e. the "displacement" (array offset) in the receive buffer for a thread is
! the sum of all element counts of the previous threads.
subroutine myMPIGatherDoubleV(sendbuf,n,recvbuf,recvcnt,ierr)
  use global
  implicit none
  include 'mpif.h'

  integer n,ierr
  real*8 sendbuf(n),recvbuf(MaxMPINodes*n)
  integer displacements(nproc), recvcnt(nproc), displ, i

  displ = 0
  do i = 1, nproc
    displacements(i) = displ
    displ = displ + recvcnt(i)
  enddo

  call MPI_GATHERV(sendbuf,n,MPI_DOUBLE_PRECISION,recvbuf,recvcnt,displacements,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
end subroutine myMPIGatherDoubleV

subroutine myMPIAllGatherDouble(sendbuf,n,recvbuf,ierr)
  use global
  implicit none
  include 'mpif.h'

  integer n,ierr
  real*8 sendbuf(n),recvbuf(MaxMPINodes*n)

  call MPI_ALLGATHER(sendbuf,n,MPI_DOUBLE_PRECISION,recvbuf,n,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)

end subroutine myMPIAllGatherDouble


subroutine myMPIScatterDouble(sendbuf,n,recvbuf,ierr)
  use global
  implicit none
  include 'mpif.h'

  integer n,ierr
  real*8 sendbuf(MaxMPINodes*n),recvbuf(n)

  call MPI_SCATTER(sendbuf,n,MPI_DOUBLE_PRECISION,recvbuf,n,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

end subroutine myMPIScatterDouble


subroutine myMPISendDouble(mpiV,n,id,tag)
  ! send vector to node 'id'
  implicit none
  include 'mpif.h'
  integer                :: n        ! size of array
  real*8                 :: mpiV(n)
  integer                :: id       ! id of receiving node
  integer                :: tag      ! message tag
  integer ierr

  !write(*,*) 'mpis:',n,id,tag

  call MPI_SEND(mpiV,n,MPI_DOUBLE_PRECISION,id,tag,MPI_COMM_WORLD,ierr)
  !write(*,*) 'mpiss:',n,id,tag,ierr
end subroutine myMPISendDouble

subroutine myMPIReceiveDouble(mpiV,n,id,tag)
  ! receive vector from node 'id'
  implicit none
  include 'mpif.h'
  integer                :: n        ! size of array
  real*8                 :: mpiV(n)
  integer                :: id       ! id of sending node
  integer                :: tag      ! message tag
  integer ierr, status(MPI_STATUS_SIZE)

  !write(*,*) 'mpir:',n,id,tag

  call MPI_RECV(mpiV,n,MPI_DOUBLE_PRECISION,id,tag,MPI_COMM_WORLD,status,ierr)
  !write(*,*) 'mpirr:',n,id,tag,ierr
end subroutine myMPIReceiveDouble


subroutine myMPISendInteger(mpiV,n,id,tag)
  ! send vector to node 'id'
  implicit none
  include 'mpif.h'
  integer                :: n        ! size of array
  integer                :: mpiV(n)
  integer                :: id       ! id of receiving node
  integer                :: tag      ! message tag
  integer ierr

  call MPI_SEND(mpiV,n,MPI_INTEGER,id,tag,MPI_COMM_WORLD,ierr)
end subroutine myMPISendInteger

subroutine myMPIReceiveInteger(mpiV,n,id,tag)
  ! receive vector from node 'id'
  implicit none
  include 'mpif.h'
  integer                :: n        ! size of array
  integer                :: mpiV(n)
  integer                :: id       ! id of sending node
  integer                :: tag      ! message tag
  integer ierr, status(MPI_STATUS_SIZE)

  call MPI_RECV(mpiV,n,MPI_INTEGER,id,tag,MPI_COMM_WORLD,status,ierr)
end subroutine myMPIReceiveInteger


subroutine myMPISendString(string,id)
  implicit none
  include 'mpif.h'
  character(len=*)       :: string ! string to send (blocking)
  integer, intent(inout) :: id     ! mpi task id to send to
  integer ierr

  call MPI_SEND(string,len(string),MPI_CHARACTER,id,len(string), MPI_COMM_WORLD,ierr)
end subroutine myMPISendString

subroutine myMPIReceiveString(string)
  implicit none
  include 'mpif.h'
  character(len=*)       :: string ! string to receive (blocking) from master
  integer ierr, status(MPI_STATUS_SIZE)

  call MPI_SEND(string,len(string),MPI_CHARACTER,0,MPI_COMM_WORLD,status,ierr)
end subroutine myMPIReceiveString


subroutine myMPIabort(ierr)
 implicit none
 include 'mpif.h'
 integer,intent(out) :: ierr

 call MPI_ABORT(MPI_COMM_WORLD,ierr)

end subroutine myMPIabort

real*8 function myMPIWallTime()
   implicit none
   include 'mpif.h'

   myMPIWallTime = MPI_Wtime()

end function myMPIWallTime



