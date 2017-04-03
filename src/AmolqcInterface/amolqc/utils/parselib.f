c parselib.f: 
c A collection of subroutines for parsing input files:
c
c $Id: parselib.f,v 1.2 2008/02/05 22:21:35 luechow Exp $
c
c $Log: parselib.f,v $
c Revision 1.2  2008/02/05 22:21:35  luechow
c minor bugs removed after compilation with sun compiler.
c parselib.f adapted to fortran standard
c
c Revision 1.1.1.1  2007/04/25 13:42:19  luechow
c QMC program amolqc. rewritten in 2006. AL
c
c Revision 1.1.1.1  2003/12/28 15:08:07  diedrch
c Initial Code 281203
c
c Revision 2.0  1999/08/18 16:13:44  luechow
c initial f90
c
c Revision 1.1  1998/03/05 13:14:58  luechow
c erste Version in Modulform
c
c
c
c
c     -------------------------------------------------
      subroutine getblk(iu,itoken,ftoken,ldim,lines,nl)
c     -------------------------------------------------

c getblk: get string array lines(1:nl) from open file with unit iu 
c _between_ the line with the initial token itoken and line with the 
c final token ftoken.
c Returns nl=0 if either the tokens have not been found or 
c no lines are in between.
c 

      implicit none

      integer, intent(in)             :: iu
      character(len=*), intent(in)    :: itoken,ftoken
      integer, intent(in)             :: ldim
      character(len=*), intent(inout) :: lines(ldim)
      integer, intent(out)            :: nl
      integer i,k,io
      character*120 line
      
      nl = 0
      rewind(iu)      
      do
        read(iu,'(A)',iostat=io) line
        k = index(line,itoken)
        if (k .gt. 0 .or. io .ne. 0) goto 101
      enddo
 101  continue
      if (io .eq. 0) then
         do
            read(iu,'(A)',iostat=io) line
            k = index(line,ftoken)
            if (k .gt. 0 .or. io .ne. 0) goto 201
            if (nl == ldim)
     &         call error("parselib:getblk: ldim too small")
            nl = nl + 1
            lines(nl) = line
         enddo
 201     continue
      endif
         
      return 
      end      

c====================================================

c     --------------------------------------------------------------
      subroutine getNextBlock(allLines,nla,idx,itoken,ftoken,ctoken,
     &                        ldim,lines,nl)
c     --------------------------------------------------------------

c getNextBlock: get string array lines(1:nl) from AllLines(1:nla)
c _starting_ from the line with the initial token itoken and up to the
c line with the final token ftoken.
c lines starting in column 1 with the comment token ctoken are ignored
c start search from idx. change idx to line of final token ftoken +1
c Returns nl=0 if either the tokens have not been found or
c number of lines in block 'lines'
c

      implicit none

      integer, intent(in)          :: nla
      character(len=*), intent(in) :: allLines(nla)
      character(len=*), intent(in) :: itoken,ftoken,ctoken
      integer, intent(inout)       :: idx
      integer, intent(in)          :: ldim
      character(len=*), intent(inout) :: lines(ldim)
      integer, intent(out)         :: nl

      character line*120
      integer i,k,n,idx0

      nl = 0
      idx0 =idx
      do i=idx,nla
         if (allLines(i)(1:1)==ctoken(1:1)) goto 100
         k = index(allLines(i),itoken)
         if ( k .gt. 0) goto 101
 100     continue
      enddo
 101  continue

      if (k .gt. 0) then
         n = i
         nl = 1
         do i=n,nla
            k = index(allLines(i),ftoken)
            if (nl > ldim) call error("getNextBlock: wrong dimension")
            lines(nl) = allLines(i)
            idx = i+1
            if (k .gt. 0) goto 201
            nl = nl + 1
         enddo
         nl = 0
         idx = idx0
 201     continue
      endif

      return
      end

c====================================================

c     ----------------------------------------
      subroutine getdblf(iu,target,value,iflag)
c     ----------------------------------------

c getdblf: get double precision value after string target in open file
c with unit iu
c iflag: 0 if target found, 1 if not
c

      implicit none

      integer i,iu,k,kf,io,iflag
      real*8 value
      character line*120,target*(*)
      
      rewind(iu)      
      do i=1,1000
        read(iu,'(A)',iostat=io) line
        k = index(line,target)
        if ( k .gt. 0 .or. io .ne. 0) goto 101
      enddo
 101  continue

      if (k .eq. 0) then
         iflag = 1                      ! target not found
      else
         iflag = 0                      ! target found
         k = k + len(target)
         kf = index(line(k:),')')
         if (kf == 0) then
            kf = len(line)
         else
            kf = k + kf - 2
         endif
         read(line(k:kf),*) value         ! internal file used for conversion
      endif

      return
      end      
                 
c============================================

c     ----------------------------------------------
      subroutine getdbla(lines,nl,target,value,iflag)
c     ----------------------------------------------

c getdbla: get double precision value after string target in string array
c lines(1:nl)
c iflag: 0 if target found, 1 if not
c
      implicit none
      
      integer i,k,kf,nl,iflag
      real*8 value
      character lines(nl)*120,target*(*)
      
      do i=1,nl
         k = index(lines(i),target)
         if ( k .gt. 0) goto 101
      enddo
 101  continue
      
      if (k .eq. 0) then
         iflag = 1                          ! target not found
      else
         iflag = 0                          ! target found
         k = k + len(target)
         kf = index(lines(i)(k:),')')
         if (kf == 0) then
            kf = len(lines(i))
         else
            kf = k + kf - 2
         endif
         read(lines(i)(k:kf),*) value         ! internal file used for conversion
      endif

      return
      end      

c====================================================

c     -----------------------------------------
      subroutine getintf(iu,target,value,iflag)
c     -----------------------------------------

c getintf: get integer value after string target in open file
c with unit iu
c iflag: 0 if target found, 1 if not
c
      implicit none

      integer i,iu,k,kf,io,value,iflag
      character line*120,target*(*)
      
      rewind(iu)      
      do i=1,1000
        read(iu,'(A)',iostat=io) line
        k = index(line,target)
        if ( k .gt. 0 .or. io .ne. 0) goto 101
      enddo
 101  continue

      if (k .eq. 0) then
         iflag = 1                      ! target not found
      else
         iflag = 0                      ! target found
         k = k + len(target)
         kf = index(line(k:),')')
         if (kf == 0) then
            kf = len(line)
         else
            kf = k + kf - 2
         endif
         read(line(k:kf),*) value         ! internal file used for conversion
      endif

      return
      end      
                 
c====================================================

c     ------------------------------------------
      subroutine getint8f(iu,target,value,iflag)
c     ------------------------------------------

c getintf: get integer value after string target in open file
c with unit iu
c iflag: 0 if target found, 1 if not
c
      implicit none

      integer*8 value
      integer i,iu,k,kf,io,iflag,n
      character line*120,target*(*)
      character str*80
      
      rewind(iu)      
      do i=1,1000
        read(iu,'(A)',iostat=io) line
        k = index(line,target)
        if ( k .gt. 0 .or. io .ne. 0) goto 101
      enddo
 101  continue

      if (k .eq. 0) then
         iflag = 1                      ! target not found
      else
         iflag = 0                      ! target found
         k = k + len(target)
         kf = index(line(k:),')')
         if (kf == 0) then
            kf = len(line)
         else
            kf = k + kf - 2
         endif
         read(line(k:kf),*,iostat=io) str
         if (io /= 0) then
            read(line(k:),*,iostat=io) str
            if (io /= 0) call abortp("cannot parse string")
         endif
         n = len(trim(str))
         if (str(n:n)=='k') then
            read(str(:n-1),*) value
            value = value * 1000
         else if (str(n:n)=='M') then
            read(str(:n-1),*) value
            value = value * 1000000
         else
            read(str,*) value
         endif            
      endif

      return
      end      
                 
c============================================

c     -----------------------------------------------
      subroutine getinta(lines,nl,target,value,iflag)
c     -----------------------------------------------

c getinta: get integer value after string target in string array
c lines(1:nl)
c iflag: 0 if target found, 1 if not
c
      implicit none

      integer i,k,kf,nl,value,iflag
      character lines(nl)*120,target*(*)
      
      do i=1,nl
         k = index(lines(i),target)
         if ( k .gt. 0) goto 101
      enddo
 101  continue

      if (k .eq. 0) then
         iflag = 1                        ! target not found
      else
         iflag = 0                        ! target found
         k = k + len(target)
         kf = index(lines(i)(k:),')')
         if (kf == 0) then
            kf = len(lines(i))
         else
            kf = k + kf - 2
         endif
         read(lines(i)(k:kf),*) value       ! internal file used for conversion
      endif
      
      return
      end      

c============================================

c     ------------------------------------------------
      subroutine getint8a(lines,nl,target,value,iflag)
c     ------------------------------------------------

c getint8a: get integer value after string target in string array
c lines(1:nl)
c iflag: 0 if target found, 1 if not
c
      implicit none

      integer*8 value
      integer i,k,kf,n,nl,iflag,io
      character lines(nl)*120,target*(*)
      character str*80
      
      do i=1,nl
         k = index(lines(i),target)
         if ( k .gt. 0) goto 101
      enddo
 101  continue

      if (k .eq. 0) then
         iflag = 1                        ! target not found
      else
         iflag = 0                        ! target found
         k = k + len(target)
         kf = index(lines(i)(k:),')')
         if (kf == 0) then
            kf = len(lines(i))
         else
            kf = k + kf - 2
         endif
         read(lines(i)(k:kf),*,iostat=io) str
         if (io /= 0) then
            read(lines(i)(k:),*,iostat=io) str
            if (io /= 0) call abortp("cannot parse string")
         endif
         n = len(trim(str))
         if (str(n:n)=='k') then
            read(str(:n-1),*) value
            value = value * 1000
         else if (str(n:n)=='M') then
            read(str(:n-1),*) value
            value = value * 1000000
         else
            read(str,*) value
         endif            
      endif
      
      return
      end      

c====================================================

c     ---------------------------------------
      subroutine getstrf(iu,target,str,iflag)
c     ---------------------------------------

c getstrf: get string str (with apostroph!) after string target in open file
c with unit iu
c iflag: 0 if target found, 1 if not
c
      implicit none

      integer i,iu,k,kf,io,iflag
      character line*120,target*(*),str*(*)
      
      rewind(iu)      
      do i=1,1000
        read(iu,'(A)',iostat=io) line
        k = index(line,target)
        if ( k .gt. 0 .or. io .ne. 0) goto 101
      enddo
 101  continue

      if (k .eq. 0) then
         iflag = 1                      ! target not found
      else
         iflag = 0                      ! target found
         k = k + len(target)
         read(line(k:),*,iostat=io) str         ! internal file used for conversion
         if (io /= 0) then
            kf = index(line(k:),')')
            if (kf == 0) then
               kf = len(line)
            else
               kf = k + kf - 2
            endif
            read(line(k:kf),*) str         ! internal file used for conversion
         endif
      endif

      return
      end      
                 
c============================================

c     ---------------------------------------------
      subroutine getstra(lines,nl,target,str,iflag)
c     ---------------------------------------------

c getstra: get string str (with apostrophs) after string target in string array
c lines(1:nl)
c iflag: 0 if target found, 1 if not
c
      implicit none

      integer i,k,kf,nl,iflag,io
      character lines(nl)*120,target*(*),str*(*)
      
      do i=1,nl
         k = index(lines(i),target)
         if ( k .gt. 0) goto 101
      enddo
 101  continue

      if (k .eq. 0) then
         iflag = 1                        ! target not found
      else
         iflag = 0                        ! target found
         k = k + len(target)
         kf = index(lines(i)(k:),')')
         if (kf == 0) then
            kf = len(lines(i))
         else
            kf = k + kf - 2
         endif
         read(lines(i)(k:kf),*,iostat=io) str
         if (io /= 0) then
            read(lines(i)(k:),*,iostat=io) str
            if (io /= 0) call abortp("cannot parse string")
         endif
      endif
      
      return
      end      

c====================================================

c     -----------------------------------------
      subroutine getlogf(iu,target,value,iflag)
c     -----------------------------------------

c getlogf: get logical value after string target in open file
c with unit iu
c iflag: 0 if target found, 1 if not
c
      implicit none

      integer i,iu,k,kf,io,iflag
      character line*120,target*(*)
      logical value
      
      rewind(iu)      
      do i=1,1000
        read(iu,'(A)',iostat=io) line
        k = index(line,target)
        if ( k .gt. 0 .or. io .ne. 0) goto 101
      enddo
 101  continue

      if (k .eq. 0) then
         iflag = 1                      ! target not found
      else
         iflag = 0                      ! target found
         k = k + len(target)
         kf = index(line,')')
         if (kf == 0) then
            kf = len(line)
         else
            kf = k + kf - 2
         endif
         read(line(k:kf),*) value         ! internal file used for conversion
      endif

      return
      end      
                 
c============================================

c     -----------------------------------------------
      subroutine getloga(lines,nl,target,value,iflag)
c     -----------------------------------------------

c getloga: get logical value after string target in string array
c lines(1:nl)
c iflag: 0 if target found, 1 if not
c
      implicit none

      integer i,k,kf,nl,iflag
      character lines(nl)*120,target*(*)
      logical value
      
      do i=1,nl
         k = index(lines(i),target)
         if ( k .gt. 0) goto 101
      enddo
 101  continue

      if (k .eq. 0) then
         iflag = 1                        ! target not found
      else
         iflag = 0                        ! target found
         k = k + len(target)
         kf = index(lines(i)(k:),')')
         if (kf == 0) then
            kf = len(lines(i))
         else
            kf = k + kf - 2
         endif
         read(lines(i)(k:kf),*) value       ! internal file used for conversion
      endif
      
      return
      end      


c====================================================

c     ---------------------------------
      logical function findf(iu,target)
c     ---------------------------------

c findf: find string target in open file
c with unit iu
c
      implicit none

      integer i,iu,k,io
      character line*120,target*(*)
      
      rewind(iu)      
c      do i=1,1000
c        read(iu,'(A)',iostat=io) line
c        k = index(line,target)
c        if ( k .gt. 0 .or. io .ne. 0) goto 101
c      enddo
c 101  continue

ccc modified 24.02.03 by CD

      io=0
      do while (io.eq.0)
        read(iu,'(A)',iostat=io) line
        k = index(line,target)
        if ( k .gt. 0 ) io=1
      enddo

ccc end CD

      if (k .eq. 0) then
         findf = .false.                ! target not found
      else
         findf = .true.                 ! target found
      endif

      return
      end      
                 
c============================================

c     ---------------------------------------
      logical function finda(lines,nl,target)
c     ---------------------------------------

c finda: find string target in string array
c lines(1:nl)
c
      implicit none

      integer i,k,nl
      character lines(nl)*120,target*(*)
      
      do i=1,nl
         k = index(lines(i),target)
         if ( k .gt. 0) goto 101
      enddo
 101  continue

      if (k .eq. 0) then
         finda = .false.                ! target not found
      else
         finda = .true.                 ! target found
      endif
      
      return
      end      

c============================================

c     ----------------------------------------
      integer function ifinda(lines,nl,target)
c     ----------------------------------------

c ifinda: find string target in string array
c lines(1:nl). return line number if found,
c else zero
c
      implicit none

      integer i,k,nl
      character lines(nl)*120,target*(*)
      
      do i=1,nl
         k = index(lines(i),target)
         if ( k .gt. 0) goto 101
      enddo
 101  continue

      if (k .eq. 0) then
         ifinda = 0                ! target not found
      else
         ifinda = i                 ! target found
      endif
      
      return
      end


c============================================

c     -----------------------------------------------
      subroutine replaceEntry(line,entryIdx,newEntry)
c     -----------------------------------------------

      implicit none
      integer nl,entryIdx
      character line*(*),newEntry*(*)
      character str*40
      integer nnew,i,io,nold,offset

      nnew = len_trim(newEntry)
      read(line(entryIdx:),*,iostat=io) str
      nold = len_trim(str)
      if (str(nold:nold)==')') nold = nold-1

      offset = nnew - nold
      call shiftTail(line,entryIdx,offset)
      line(entryIdx:entryIdx+nnew-1) = newEntry(1:nnew)
      return
      end



c============================================

c     -------------------------------------
      subroutine shiftTail(line,idx,offset)
c     -------------------------------------

      ! shift tail of line starting from idx by 'offset' characters
      ! offset > 0: insert 'offset' spaces at 'idx' position
      ! offset < 0: remove 'offset' characters starting from 'idx'

      implicit none
      character line*120
      integer idx,offset
      integer newlen,oldlen,i,ii

      oldlen = len_trim(line)
      newlen = oldlen + offset

      if (newlen > 120) then
         newlen = 120
         oldlen = newlen - offset
      endif

      if (offset > 0) then ! shift to right
         do i=oldlen,idx,-1
            ii = i+offset
            line(ii:ii) = line(i:i)
         enddo
         do i=idx,idx+offset-1
            line(i:i) = ' '
         enddo
      else if (offset < 0) then ! shift to left
         do i=idx,newlen
            ii = i-offset
            line(i:i) = line(ii:ii)
         enddo
         do i=newlen+1,oldlen
            line(i:i) = ' '
         enddo
      endif

      return
      end
