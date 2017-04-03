module utilsmodule
  use rannum
  use mrg
  use mt19937

! module to make explicit all interfaces of subroutines and functions
! in the utils library (that are not included in modules)
! use utilsmodules has to be included in routines that use routines from
! utils

implicit none

  interface
    pure subroutine abortp(s)
      character(*), intent(in) :: s
    end subroutine abortp
  end interface

! error.f90
  interface
    pure subroutine assert(logicalExp,string)
      logical, intent(in)                    :: logicalExp
      character(len=*), optional, intent(in) :: string
    end subroutine assert
  end interface

  interface
    subroutine assertEqualAbsolute(t1,t2,tol,string)
      real*8, intent(in)                     :: t1
      real*8, intent(in)                     :: t2
      real*8, intent(in)                     :: tol
      character(len=*), optional, intent(in) :: string
    end subroutine assertEqualAbsolute
  end interface

  interface
    subroutine assertEqualRelative(t1,t2,tol,string)
      real*8, intent(in)                     :: t1
      real*8, intent(in)                     :: t2
      real*8, intent(in)                     :: tol
      character(len=*), optional, intent(in) :: string
    end subroutine assertEqualRelative
  end interface

  interface
    subroutine error(s)
      character(len=*), intent(in) :: s
    end subroutine error
  end interface

! nelmead.f
  interface
    subroutine nelmead(x0,n0,maxiter,delta,eps,fmin,func1,iul,logmode)
      integer n0                    ! tatsaechliche Dimension
      real*8 x0(n0)                 ! zu optimierender Vektor
      integer maxiter               ! max Anzahl an Iterationen
      real*8 delta                  ! anfaengliche Auslenkung
      real*8 eps                    ! Zielgenauigkeit
      real*8 fmin                   ! minimierter Funktionswert
      interface                     ! external function
        real*8 function func1(x,n)
          integer n
          real*8 x(n)
        end function func1
      end interface
      integer iul                   ! open file unit for log messages
      integer logmode               ! verboseness level
    end subroutine nelmead
  end interface

! sdiag.f
  interface
    subroutine sdiag(A,MZA,MSA,D,MD,V,MZV,MSV,NB,IERR)
      integer MZA,MSA,MD,MZV,MSV,NB,IERR
      real*8 A(MZA,MSA),D(MD),V(MZV,MSV)
    end subroutine sdiag
  end interface

! nmtr
  interface
    subroutine dnmtr(n,x,f,g,H,ldH,frtol,fatol,fmin,task, &
                     delta,diag,scale,isave,dsave,xc,gc,s,gs,wa)
      character*60 task
      logical scale
      integer n, ldH
      integer isave(4)
      real*8 f, frtol, fatol, fmin, delta
      real*8 x(n), g(n), H(ldH,n), diag(n)
      real*8 dsave(7)
      real*8 xc(n), gc(n), s(n), gs(n), wa(n)
    end subroutine
  end interface

! parselib.f
  interface
    subroutine getblk(iu,itoken,ftoken,ldim,lines,nl)
      integer, intent(in)             :: iu
      character(len=*), intent(in)    :: itoken,ftoken
      integer, intent(in)             :: ldim
      character(len=*), intent(inout) :: lines(ldim)
      integer, intent(out)            :: nl
    end subroutine getblk
  end interface

  interface
    subroutine getNextBlock(allLines,nla,idx,itoken,ftoken,ctoken,ldim,lines,nl)
      integer, intent(in)          :: nla
      character(len=*), intent(in) :: allLines(nla)
      character(len=*), intent(in) :: itoken,ftoken,ctoken
      integer, intent(inout)       :: idx
      integer, intent(in)          :: ldim
      character(len=*), intent(inout) :: lines(ldim)
      integer, intent(out)         :: nl
    end subroutine getNextBlock
  end interface

  interface
    subroutine getdblf(iu,target,value,iflag)
      integer iu,iflag
      real*8 value
      character target*(*)
    end subroutine getdblf
  end interface

  interface
    subroutine getdbla(lines,nl,target,value,iflag)
      integer nl,iflag
      real*8 value
      character lines(nl)*80,target*(*)
    end subroutine getdbla
  end interface

  interface
    subroutine getintf(iu,target,value,iflag)
      integer iu,iflag
      integer value
      character target*(*)
    end subroutine getintf
  end interface

  interface
    subroutine getinta(lines,nl,target,value,iflag)
      integer nl,iflag
      integer value
      character lines(nl)*80,target*(*)
    end subroutine getinta
  end interface

  interface
    subroutine getint8f(iu,target,value,iflag)
      integer iu,iflag
      integer*8 value
      character target*(*)
    end subroutine getint8f
  end interface

  interface
    subroutine getint8a(lines,nl,target,value,iflag)
      integer nl,iflag
      integer*8 value
      character lines(nl)*80,target*(*)
    end subroutine getint8a
  end interface

  interface
    subroutine getstrf(iu,target,value,iflag)
      integer iu,iflag
      character*(*) value
      character target*(*)
    end subroutine getstrf
  end interface

  interface
    subroutine getstra(lines,nl,target,value,iflag)
      integer nl,iflag
      character*(*) value
      character lines(nl)*80,target*(*)
    end subroutine getstra
  end interface

  interface
    subroutine getlogf(iu,target,value,iflag)
      integer iu,iflag
      logical value
      character target*(*)
    end subroutine getlogf
  end interface

  interface
    subroutine getloga(lines,nl,target,value,iflag)
      integer nl,iflag
      logical value
      character lines(nl)*80,target*(*)
    end subroutine getloga
  end interface

  interface
    logical function findf(iu,target)
      integer iu
      character target*(*)
    end function findf
  end interface

  interface
    logical function finda(lines,nl,target)
      integer nl
      character lines(nl)*80,target*(*)
    end function finda
  end interface

  interface
    integer function ifinda(lines,nl,target)
      integer nl
      character lines(nl)*80,target*(*)
    end function ifinda
  end interface

! ran_m.f
  interface myran
#ifdef MRG
    module procedure mrg_ran
#elif MT
    module procedure mt_ran
#else
    module procedure ran2_ran
#endif
  end interface myran

  interface mygran
#ifdef MRG
    module procedure mrg_gran
#elif MT
    module procedure mt_gran
#else
    module procedure ran2_gran
#endif
  end interface mygran

  interface init_ran
#ifdef MRG
    module procedure init_mrgran
#elif MT
    module procedure init_mtran
#else
    module procedure init_ran2
#endif
  end interface init_ran

end module
