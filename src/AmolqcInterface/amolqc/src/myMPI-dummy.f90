
! collection of routines that encapsulate MPI

! this is the dummy version that provides empty interfaces only

! Arne Luechow, 2006

integer function getMytid()
  implicit none

  getMytid = 0
end function getMytid


subroutine myMPIInitialize(ierr)
  use global
  implicit none
  integer ierr

  ierr = 0
  mytid = 0
  nproc = 1
end subroutine myMPIInitialize


subroutine myMPIFinalize(ierr)
  implicit none
  integer ierr

  ierr = 0
end subroutine myMPIFinalize


subroutine myMPIBarrier(ierr)
  implicit none
  integer ierr

  ierr = 0
end subroutine myMPIBarrier


subroutine myMPIBcastString(string,n)
  use global
  implicit none
  integer n
  character*(*) string

end subroutine myMPIBcastString


subroutine myMPIBcastInteger(array,n)
  use global
  implicit none
  integer n,ierr
  integer array(n)

end subroutine myMPIBcastInteger


subroutine myMPIBcastDouble(array,n)
  use global
  implicit none
  integer n,ierr
  double precision array(n)

end subroutine myMPIBcastDouble


subroutine myMPIReduceSumDouble(sendbuf,recvbuf,n)
  use global
  implicit none
  integer n
  real*8 sendbuf(n),recvbuf(n)

  recvbuf = sendbuf

end subroutine myMPIReduceSumDouble


subroutine myMPIReduceSumInteger(sendbuf,recvbuf,n)
  use global
  implicit none
  integer n
  integer sendbuf(n),recvbuf(n)

  recvbuf = sendbuf

end subroutine myMPIReduceSumInteger

subroutine myMPIAllReduceSumDouble(sendbuf,recvbuf,n)
  use global
  implicit none
  integer n
  real*8 sendbuf(n),recvbuf(n)

  recvbuf = sendbuf

end subroutine myMPIAllReduceSumDouble


subroutine myMPIAllReduceSumInteger(sendbuf,recvbuf,n)
  use global
  implicit none
  integer n
  integer sendbuf(n),recvbuf(n)

  recvbuf = sendbuf

end subroutine myMPIAllReduceSumInteger


subroutine myMPIGatherInteger(sendbuf,n,recvbuf,ierr)
  implicit none

  integer n,ierr
  integer sendbuf(n),recvbuf(n)

  recvbuf = sendbuf
  ierr = 0
end subroutine myMPIGatherInteger


subroutine myMPIScatterInteger(sendbuf,n,recvbuf,ierr)
  implicit none

  integer n,ierr
  integer sendbuf(n),recvbuf(n)

  recvbuf = sendbuf
  ierr = 0
end subroutine myMPIScatterInteger


subroutine myMPIGatherDouble(sendbuf,n,recvbuf,ierr)
  implicit none

  integer n,ierr
  real*8 sendbuf(n),recvbuf(n)

  recvbuf = sendbuf
  ierr = 0
end subroutine myMPIGatherDouble

subroutine myMPIGatherDoubleV(sendbuf,n,recvbuf,recvcnt,ierr)
  implicit none

  integer n,ierr
  real*8 sendbuf(n),recvbuf(n)
  integer recvcnt(1)

  recvbuf = sendbuf
  ierr = 0
end subroutine myMPIGatherDoubleV

subroutine myMPIAllGatherDouble(sendbuf,n,recvbuf,ierr)
  implicit none

  integer n,ierr
  real*8 sendbuf(n),recvbuf(n)

  recvbuf = sendbuf
  ierr = 0
end subroutine myMPIAllGatherDouble

subroutine myMPIScatterDouble(sendbuf,n,recvbuf,ierr)
  implicit none

  integer n,ierr
  real*8 sendbuf(n),recvbuf(n)

  recvbuf = sendbuf
  ierr = 0
end subroutine myMPIScatterDouble


subroutine myMPISendDouble(mpiV,n,id,tag)
  ! send vector to node 'id'
  implicit none
  integer                :: n        ! size of array
  real*8                 :: mpiV(n)
  integer                :: id       ! id of receiving node
  integer                :: tag      ! message tag
  integer ierr

end subroutine myMPISendDouble

subroutine myMPIReceiveDouble(mpiV,n,id,tag)
  ! receive vector from node 'id'
  implicit none
  integer                :: n        ! size of array
  real*8                 :: mpiV(n)
  integer                :: id       ! id of receiving node
  integer                :: tag      ! message tag
  integer ierr, status

end subroutine myMPIReceiveDouble


subroutine myMPISendInteger(mpiV,n,id,tag)
  ! send vector to node 'id'
  implicit none
  integer                :: n        ! size of array
  integer                :: mpiV(n)
  integer                :: id       ! id of receiving node
  integer                :: tag      ! message tag
  integer ierr

end subroutine myMPISendInteger

subroutine myMPIReceiveInteger(mpiV,n,id,tag)
  ! receive vector from node 'id'
  implicit none
  integer                :: n        ! size of array
  integer                :: mpiV(n)
  integer                :: id       ! id of receiving node
  integer                :: tag      ! message tag
  integer ierr, status

end subroutine myMPIReceiveInteger


subroutine myMPISendString(string,id)
  implicit none
  character(len=*)       :: string ! string to send (blocking)
  integer, intent(inout) :: id     ! mpi task id to send to
  integer ierr

  ierr = 0
end subroutine myMPISendString

subroutine myMPIReceiveString(string)
  implicit none
  character(len=*)       :: string ! string to receive (blocking) from master
  integer ierr, status

  string = ''
end subroutine myMPIReceiveString

subroutine myMPIabort(ierr)
 implicit none
 integer,intent(out) :: ierr

 ierr = 1

end subroutine myMPIabort



real*8 function myMPIWallTime()
   implicit none
   integer*8 count,count_rate

   call system_clock(count,count_rate)
   myMPIWallTime = dble(count)/dble(count_rate)

end function myMPIWallTime

