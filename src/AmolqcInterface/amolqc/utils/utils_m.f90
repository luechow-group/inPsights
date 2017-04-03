!==========!
module Utils
!==========!

  ! collection of simple general purpose routines

  ! $Id: utils_m.f90,v 1.1.1.1 2007/04/25 13:42:19 luechow Exp $

  ! $Log: utils_m.f90,v $
  ! Revision 1.1.1.1  2007/04/25 13:42:19  luechow
  ! QMC program amolqc. rewritten in 2006. AL
  !

#ifdef PARALLEL
  !!!use mpi
#endif

  implicit none
#ifdef PARALLEL
  include 'mpif.h'
#endif

contains


!------  string handling ----------------------------------


  !-------------------------------------!
  subroutine tokenize(line,token,nTokens)
  !-------------------------------------!
    ! splits line into "tokens" separated by blanks or comma (not white space!)
    integer, parameter                :: lSize=120
    character(len=*), intent(in)      :: line
    character(len=*), intent(inout)   :: token(:)
    integer, intent(out)              :: nTokens
    character(len=lSize) :: lineTmp
    integer              :: idx,idx1,idx2

    if (len_trim(line) > lSize) call error("tokenize: line too long")
    lineTmp = trim(line)
    nTokens = 0
    do
       lineTmp = adjustl(lineTmp)
       idx1 = index(lineTmp," ")
       idx2 = index(lineTmp,",")
       if (idx1==0) then
          idx = idx2
       else if (idx2 > 0) then
          idx = min(idx1,idx2)
       else
          idx = idx1
       endif
       if (idx==0 .or. len_trim(lineTmp)==0) exit
       nTokens = nTokens + 1
       if (nTokens > size(token)) call error("tokenize: too many tokens")
       token(nTokens) = lineTmp(1:idx-1)
       lineTmp = lineTmp(idx+1:)
    enddo
    lineTmp = adjustl(lineTmp)
    if (len_trim(lineTmp)>0) then
       nTokens = nTokens + 1
       if (nTokens > size(token)) call error("tokenize: too many tokens")
       token(nTokens) = lineTmp
    endif
  end subroutine tokenize


  !---------------------------------------------------------!
  subroutine split(string,splitString,leftString,rightString)
  !---------------------------------------------------------!
    ! split string at splitString (e.g. "=")
    character(len=*), intent(in)    :: string
    character(len=*), intent(in)    :: splitString
    character(len=*), intent(inout) :: leftString
    character(len=*), intent(inout) :: rightString
    integer idx

    idx = index(string,splitString)
    if (idx == 0) then
       leftString = string
       rightString = ""
    else
       leftString = string(1:idx-1)
       rightString = string(idx+len_trim(splitString):)
    endif
  end subroutine split

  !---------------------------------------------------!
  subroutine replace(string,needle,replacement,global_)
  !---------------------------------------------------!
    ! replace one (or all, if |global| = .true.) occurence of |needle| in
    ! |string| with |replacement|
    ! note: if the declared length of |string| is not large enough to contain
    !       the |string| with all occurences of |needle| replaced, it will be
    !       cut off at the end
    character(len=*), intent(inout) :: string
    character(len=*), intent(in)    :: needle
    character(len=*), intent(in)    :: replacement
    logical, optional               :: global_
    logical :: global = .false.
    integer :: i, slen
    if(present(global_)) global = global_

    slen = len(needle)
    i = 1
    do while(i < len(string) - slen + 1)
      if (string(i:i+slen-1) == needle) then
        string = string(1:i-1) // trim(replacement) // string(i+slen:)
        i = i + slen - 1

        if(.not. global) exit
      endif
      i = i + 1
    enddo
  end subroutine replace

  !---------------------------------!
  function strToUpper(in) result(res)
  !---------------------------------!
    ! converts all lowercase ASCII characters in |in| to their uppercase
    ! equivalents
    character(len=*), intent(in) :: in
    character(len=len(in)) :: res

    integer :: i, chr

    res(:) = in(:)
    do i = 1, len(in)
      chr = iachar(in(i:i))
      if (chr >= iachar("a") .and. chr <= iachar("z")) then
        res(i:i) = achar(chr - 32)
      endif
    enddo
  end function strToUpper


  !----------------------------------------------!
  function getToken(string,initString,finalString)
  !----------------------------------------------!
    ! get the string between the initString and finalString
    character(len=120)              :: getToken
    character(len=*)              :: string
    character(len=*)              :: initString
    character(len=*)              :: finalString
    integer idx1,idx2,i1,i2

    idx1 = index(string,initString)
    idx2 = index(string,finalString)

    if (idx2 > idx1) then
       i1  = idx1 + len_trim(initString)
       i2  = idx2 - 1
       getToken = string(i1:i2)
    else
       getToken = ""
    endif
  end function getToken


  !---------------------------------!
  integer function numberType(string)
  !---------------------------------!
    ! 0:"no number", 1:"integer", 2:"real or double"
    character(len=*), intent(in) :: string
    integer idx,n,io
    real*8 a

    idx = index(string,".")
    if (idx > 0) then
       ! possibly float
       read(string,*,iostat=io) a
       if (io==0) then
          numberType = 2
       else
          numberType = 0
       endif
    else
       ! possibly integer
       read(string,*,iostat=io) n
       if (io==0) then
          numberType = 1
       else
          numberType = 0
       endif
    endif
  end function numberType

  !---------------------------------------------------------
  subroutine expandMacro(allLines,nla,macroLines,mLines,idx)
  !---------------------------------------------------------

      ! insert macro in macroLines into string array allLines at position idx
      ! allLines(idx-1) holds the macro cmd
      integer, parameter           :: MAXTOKEN = 20
      ! allow exchange of values only when 1st character is schar
      character, parameter         :: xchar = 'x'        !

      character(len=*), intent(inout) :: allLines(:)     ! lines of in file containing macros
      integer, intent(inout)       :: nla                ! real # of lines
      character(len=*), intent(inout) :: macroLines(:)   ! string array holding macro
      integer, intent(in)          :: mLines             ! real # of macrolines
      integer, intent(in)       :: idx                   ! insert macro at this line

      character(len=120) :: strtmp,key,value
      integer i,k,t,io,nToken
      character(len=120) :: tokens(MAXTOKEN)

      ! replace all matching entries from lines(nl) i.e. the macro
      ! argument, in newLines

      strtmp = getToken(allLines(idx-1),'(',')')
      call tokenize(strtmp,tokens,nToken)

      do t=1,nToken
         call split(tokens(t),"=",key,value)
         do i=1,mLines
            if (macroLines(i)(1:1)/=xchar) cycle
            k = index(macroLines(i),trim(key)//"=")
            if (k>0) then
               call replaceEntry(macroLines(i),k,tokens(t))
            endif
         enddo
      enddo

      ! insert new lines
      call assert(nla+mLines<=size(allLines),"expandMacro: overflow in in-file lines")
      do i=nla,idx,-1
         allLines(i+mLines) = allLines(i)
      end do
      nla = nla + mLines
      do i=1,mLines
         if (macroLines(i)(1:1)==xchar) then
            allLines(idx+i-1) = macroLines(i)(2:)
         else
            allLines(idx+i-1) = macroLines(i)
         end if
      end do

  end subroutine expandMacro


  !---------------------------------------------!
  subroutine readFileLocal(fName,allLines,nLines)
  !---------------------------------------------!

   character(len=*), intent(in)    :: fName         ! file name of file to be read
   character(len=*), intent(inout) :: allLines(:)   ! string array to hold full file
   integer, intent(out)            :: nLines        ! output: number of lines read from file
   integer, parameter              :: iu=21
   integer                         :: io,ierr,nn
   character(len=72)               :: errMsg


   open(iu,file=fName,status='old',iostat=io)
   if (io /= 0) then
      errMsg = trim(fName) // ' could not be opened'
      call error(errMsg)
   end if
   nLines = 1
   do
      read(iu,'(A)',iostat=io) allLines(nLines)
      if (io /= 0) then
         nLines = nLines-1
         exit
      endif
      if (nLines==size(allLines)) then
         errMsg = trim(fName) // ' is longer than the lines array'
         call error(errMsg)
      end if
      nLines = nLines + 1
   end do
   close(iu)
  end subroutine readFileLocal


  !-------------------------------------------------------!
  subroutine readFileParallel(myrank,fName,allLines,nLines)
  !-------------------------------------------------------!

   integer, intent(in)             :: myrank        ! mpi rank
   character(len=*), intent(in)    :: fName         ! file name of file to be read
   character(len=*), intent(inout) :: allLines(:)   ! string array to hold full file
   integer, intent(out)            :: nLines        ! output: number of lines read from file
   integer, parameter              :: iu=21
   integer, parameter              :: MASTER=0
   integer                         :: io,ierr,nn
   character(len=72)               :: errMsg


   if (myrank == MASTER) then
      open(iu,file=fName,status='old',iostat=io)
      if (io /= 0) then
         errMsg = trim(fName) // ' could not be opened'
         call error(errMsg)
      end if
      nLines = 1
      do
         read(iu,'(A)',iostat=io) allLines(nLines)
         if (io /= 0) then
            nLines = nLines-1
            exit
         endif
         if (nLines==size(allLines)) then
            errMsg = trim(fName) // ' is longer than the lines array'
            call error(errMsg)
         end if
         nLines = nLines + 1
      end do
      close(iu)
   end if

#ifdef PARALLEL
   call mpi_bcast(nLines,1,MPI_INTEGER,MASTER,MPI_COMM_WORLD,ierr)
   call mpi_bcast(allLines,len(allLines(1))*nLines,MPI_CHARACTER,MASTER,MPI_COMM_WORLD,ierr)
#endif

  end subroutine readFileParallel



  !--------------------------------!
  real*8 function numDerivative(f,h)
  !--------------------------------!
  real*8, intent(in)  :: f(:)
  real*8, intent(in)  :: h
  select case (size(f))
  case (3)
     numDerivative =  (f(3) - f(1)) / (2.d0*h)
  case (5)
     numDerivative =  ( 8.d0*(f(4)-f(2)) - (f(5)-f(1)) ) / (12.d0*h)
  case (7)
     numDerivative =  ( 45.d0*(f(5)-f(3)) - 9.d0*(f(6)-f(2)) + (f(7)-f(1)) ) &
                   /  (60.d0*h)
  case default
     call abortp('numDeriv: wrong rule')
  end select
  end function numDerivative


  !---------------------------------!
  real*8 function num2Derivative(f,h)
  !---------------------------------!
  real*8, intent(in)  :: f(:)
  real*8, intent(in)  :: h
  select case (size(f))
  case (3)
     num2Derivative =  (f(1) - 2d0*f(2) + f(3)) / h**2
  case (5)
     num2Derivative =  ( -30d0*f(3) + 16d0*(f(2)+f(4)) - (f(1)+f(5)) ) / (12d0*h**2)
  case (7)
     num2Derivative =  ( -490d0*f(4) + 270d0*(f(3)+f(5)) - 27d0*(f(2)+f(6)) + 2d0*(f(1)+f(7)) ) &
                    /  (180d0*h**2)
  case default
     call abortp('num2Derivative: wrong rule')
  end select
  end function num2Derivative

  !--------------------------------------!
  function getTimeString(sec) result (str)
  !--------------------------------------!
   real*8, intent(in) :: sec
   character(len=17)  :: str
   integer            :: days, hours, minutes, secs, tsecs
   real*8             :: seconds
   seconds = sec
   days = int(seconds/86400)
   seconds = seconds - days*86400
   hours = int(seconds/3600)
   seconds = seconds - hours*3600
   minutes = int(seconds/60)
   seconds = seconds - minutes*60
   secs = int(seconds)
   seconds = seconds - secs
   tsecs = int(1000*seconds)
   write(str,'(i3,a2,i2.2,a1,i2.2,a1,i2.2,a1,i3.3)') days,'d ',hours,':',minutes,':',secs,'.',tsecs
  end function getTimeString

  !--------------------------------------!
  integer function findInList(list,n)
  !--------------------------------------!
     integer, intent(in) :: list(:)
     integer, intent(in) :: n
     integer i
     findInList = 0
     do i=1,size(list)
        if (list(i)==n) then
           findInList = i
           exit
        end if
     end do
  end function findInList


   logical function int_isInList(vec,value)
   !-------------------------------------!
      integer, intent(in) :: vec(:)
      integer, intent(in) :: value
      integer j
      int_isInList = .false.
      do j=1,size(vec)
         if (vec(j)==value) then
            int_isInList = .true.
            exit
         end if 
      end do
   end function int_isInList


   recursive function int_isInSortedList(vec,value) result(inList)
   !------------------------------------------------------------!
      ! binary search, expecting increasing order
      logical :: inList
      integer, intent(in) :: vec(:)
      integer, intent(in) :: value
      integer mid
      inList = .false.
      if (size(vec)==0) return
      mid = (size(vec)+1) / 2
      if (value == vec(mid)) then
         inList = .true.
      else if (value < vec(mid)) then
         inList = int_isInSortedList(vec(1:mid-1),value)
      else
         inList = int_isInSortedList(vec(mid+1:size(vec)),value)
      end if
   end function int_isInSortedList


   recursive subroutine double_findInSortedList(lb,ub,vec,value,left,right)
   !----------------------------------------------------------------------!
      ! binary search, expecting sorted elements in vec
      ! value is in [vec(left),vec(right)]
      ! output: smallest interval: right = left + 1
      integer, intent(in) :: lb, ub    ! bounds of vec!
      real*8, intent(in)  :: vec(lb:ub)
      real*8, intent(in)  :: value
      integer, intent(inout) :: left, right ! on entry/exit: search interval 
      integer mid
      if (right - left == 1) return
      mid = left + (right - left) / 2
      if (value <= vec(mid)) then
         right = mid
      else 
         left = mid
      end if
      call double_findInSortedList(lb,ub,vec,value,left,right)
   end subroutine double_findInSortedList


  function intToStr(int) result (str)
  !--------------------------------------!
     integer,intent(in) :: int
     character(len=10) :: str
     write(str,"(I0)") int
  end function intToStr

end module Utils
