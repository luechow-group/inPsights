! amolqcInterfaceModule provides a C interface to amolqc
! to be used as part of the amolqc lib.
! create lib with
! make src/libamolqc.a or make src/libamolqc_mpi.a (for parallel version, requires MPI)

! using the amolqc interface:
! call the appropriate functions from C
! compile C part
! link with libraries src/libamolqc.a utils/libutils.a -Lutils -lutils -L/usr/local/lib -llapack -lblas
! (assuming gcc and lapack/blas is installed in /usr/local/bin|lib)
! link at the end with fortran runtime lib when using C compiler for linkage


module amolqcInterfaceModule
   use iso_c_binding, only: c_double, c_int, c_ptr, c_f_pointer
   use global
   use InitModule
   use elocaldata
   use elocal, only: wf_init, eloc
   use econfigsModule
   use initialPositionsModule

   implicit none

   ! enum is not really used in this example
   enum, bind(c)
      enumerator :: GAUSSIAN, DENSITY, LMO
   end enum

contains

   subroutine amolqc_init() bind(c)
      integer dummy, seed
      character(len=80) :: lines(1)
      character(len=40) :: fname

      call initAmolqc()

      ! init rng and logmode (replacing init.f90:initGen)
      seed = 101 + mytid
      dummy = init_ran(seed)
      ! try supressing output
      logmode = 0
   end subroutine amolqc_init

   subroutine amolqc_set_wf(nelecs, natoms) bind(c)
      integer(c_int) :: nelecs   ! # of elecs
      integer(c_int) :: natoms   ! # of atoms
      integer nl
      character(len=80) :: lines(1)
      character(len=40) :: fname
      nl = 1
      fname = "'t.wf'"
      lines(1) = "$wf(read, file="//trim(fname)//")"
      print*, lines(1)
      call wf_init(lines,nl)
      !call setNumberOfElectrons(ne)    ! for randomwalker initialization
      !call setNumberOfCenters(ncenter)
      nelecs = ne
      natoms = ncenter
   end subroutine amolqc_set_wf

   subroutine amolqc_initial_positions(mode, n, x) bind(c)
      integer(c_int), value :: mode   ! enum above
      integer(c_int), value :: n   
      real(c_double)  :: x(3*n)
      real*8  xx(n), yy(n), zz(n)
      integer i

      if (n /= ne) call abortp("amolqc_initial_positions: n /= ne")

      do i = 1, n
         xx(i) = x(3*i-2)
         yy(i) = x(3*i-1)
         zz(i) = x(3*i)
      enddo

      call createRandomElectronPositions(mode,xx,yy,zz)

      do i = 1, n
         x(3*i-2) = xx(i) 
         x(3*i-1) = yy(i) 
         x(3*i)   = zz(i) 
      enddo
   end subroutine amolqc_initial_positions

   subroutine amolqc_eloc(x, n, phi, u, drift, elocal) bind(c)
      integer(c_int), value :: n
      real(c_double)        :: x(3*n)
      real(c_double)        :: phi
      real(c_double)        :: u
      real(c_double)        :: elocal
      real(c_double)        :: drift(3*n)
      real*8  xx(n), yy(n), zz(n)
      integer i
      type(eConfigArray) :: ec

      if (n /= ne) call abortp("amolqc_eloc: n /= ne")

      do i = 1, n
         xx(i) = x(3*i-2)
         yy(i) = x(3*i-1)
         zz(i) = x(3*i)
      enddo

      call eConfigArray_new(ec,n,1)
      call eConfigArray_set(ec,1,xx,yy,zz)
      call eloc(0,ec,'none')

      phi = elPhi(1)
      u = elU(1)
      elocal = elEloc(1)

      do i = 1, n
         drift(3*i-2) = elxDrift(i, 1)
         drift(3*i-1) = elyDrift(i, 1)
         drift(3*i)   = elzDrift(i, 1)
      enddo

      call eConfigArray_destroy(ec)
   end subroutine amolqc_eloc

end module amolqcInterfaceModule
