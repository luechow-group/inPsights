c 
c myLinAlg.f: contains miscellaneous linear algebra subroutines
c
c depends on BLAS and on UMFPACK

c $Id: myLinAlg.f,v 1.1.1.1 2007/04/25 13:42:20 luechow Exp $

c $Log: myLinAlg.f,v $
c Revision 1.1.1.1  2007/04/25 13:42:20  luechow
c QMC program amolqc. rewritten in 2006. AL
c


c     ---------------------------------------
      subroutine invcupd(col,i,n,nmax,a1,det)
c     ---------------------------------------

c invupd updates the inverse matrix a1 and the determinant det
c for a matrix with a new column i col.
c
c Both BLAS and NON-BLAS versions in one routine
c
c input: 
c     col(nmax) : new column in matrix
c             i : column index
c             n : actual dimension of matrix a1 and vector col
c          nmax : definition of matrix and vector
c  a1(nmax,nmax): old inverse matrix
c           det : old determinant 

c output:
c  a1(nmax,nmax): updated inverse matrix
c           det : updated determinant
cTS
      use global, only: useLAlib
      implicit none
      integer n,nmax,i,j,k,l
      real*8 alpha,f(nmax)
      real*8 col(nmax),a1(nmax,nmax),det,r,tmp,one,zero
c
      if (useLAlib) then
        one = 1.0d0
        zero = 0.0d0
        call DGEMV('N',n,n,one,a1,nmax,col,1,zero,f,1) 
        r =f(i)
        det = r*det
        alpha=one/r

c       // update row i of inverse
c       a1(1,i) is address of first element of i-th column
c       In FORTRAN, matrices are stored column-wise as arrays
        call DSCAL(n,alpha,a1(i,1),nmax)

c       // update all other columns
        do j=1,n
           if (j.eq.i) goto 202
           alpha = -f(j)
           call DAXPY(n,alpha,a1(i,1),nmax,a1(j,1),nmax)
 202       continue
        enddo
      else
        r = 0d0
        do k=1,n
           r = r + col(k)*a1(i,k)
        enddo
        det = r*det
c       // update row i 
        do k=1,n
           a1(i,k) = a1(i,k)/r
        enddo

c       // update all other rows
        do j=1,n
           if (j.eq.i) goto 102
           tmp = 0d0
           do l=1,n
              tmp = tmp + col(l)*a1(j,l)
           enddo
           do k=1,n
              a1(j,k) = a1(j,k) - a1(i,k)*tmp
           enddo
 102       continue
        enddo
      endif

      end

c     ---------------------------------------
      subroutine invdetcalc(col,n,nmax,a1,det)
c     ---------------------------------------
c calculates the determinant for an updated inverse matrix (following the
c matrix determinant lemma). note that no matrices are modified by this,
c only the new determinant is calculated

      use global, only: useLAlib
      implicit none
      integer n,nmax,i
      real*8 col(nmax),a1(nmax),det
      real*8 ddot
c
      if (useLAlib) then
        det = DDOT(n,col(1:n),1,a1(1:n),1) * det
      else
        det = dot_product(col(1:n), a1(1:n)) * det
      endif

      end

c================================================

c     ---------------------------------------
      subroutine invrupd(row,i,n,nmax,a1,det)
c     ---------------------------------------

c invupd updates the inverse matrix a1 and the determinant det
c for a matrix with a new row i row(.).
c
c Both BLAS and NON-BLAS versions in one routine
c
c input: 
c     row(nmax) : new row in matrix
c             i : row index
c             n : actual dimension of matrix a1 and vector col
c          nmax : definition of matrix and vector
c  a1(nmax,nmax): old inverse matrix
c           det : old determinant 

c output:
c  a1(nmax,nmax): updated inverse matrix
c           det : updated determinant
c
cTS
      use global, only: useLAlib
      implicit none
      integer n,nmax,i,j,k,l
      real*8 alpha,f(nmax),row(nmax),a1(nmax,nmax),det,r,tmp,zero,one
c
      if (useLAlib) then
        zero = 0.d0
        one = 1.d0
        call DGEMV('T',n,n,one,a1,nmax,row,1,zero,f,1) 
        r =f(i)
        det = r*det
        alpha=one/r

c       // update column i of inverse
c       a1(1,i) is address of first element of i-th column
c       In FORTRAN, matrices are stored column-wise as arrays
        call DSCAL(n,alpha,a1(1,i),1)

c       // update all other columns
        do j=1,n
          if (j.eq.i) goto 202
          alpha = -f(j)
          call DAXPY(n,alpha,a1(1,i),1,a1(1,j),1)
 202      continue
        enddo

      else
        r = 0d0
        do k=1,n
           r = r + row(k)*a1(k,i)
        enddo
        det = r*det

c       // update column i of inverse
        do k=1,n
           a1(k,i) = a1(k,i)/r
        enddo

c       // update all other columns
        do j=1,n
           if (j.eq.i) goto 302
           tmp = 0d0
           do l=1,n
            tmp = tmp + row(l)*a1(l,j)
           enddo
           do k=1,n
              a1(k,j) = a1(k,j) - a1(k,i)*tmp
           enddo
 302       continue
        enddo
      endif

      end


c================================================

c     --------------------------------
      subroutine inv1(a,n,nmax,a1,det)
c     --------------------------------

c calculates the inverse matrix of a by updating the unit matrix

c input: 
c     col(nmax) : new column in matrix
c             i : column index
c             n : actual dimension of matrix a1 and vector col
c          nmax : definition of matrix and vector
c  a1(nmax,nmax): old inverse matrix
c           det : old determinant 

c output:
c  a1(nmax,nmax): updated inverse matrix
c           det : updated determinant


      implicit none
      integer n,nmax,i,j,k,l
      real*8 a(nmax,nmax),a1(nmax,nmax),det,r,tmp

c     // unit matrix is self-inverse with det=1
      do j=1,n
         do i=1,n
            a1(i,j) = 0d0
         enddo
         a1(j,j) = 1d0
      enddo
      det = 1d0

c     // loop over columns i for updating unit matrix
      do i=1,n

         r = 0d0
         do k=1,n
            r = r + a(k,i)*a1(i,k)
         enddo
         det = r*det

c        // update row i 
         do k=1,n
            a1(i,k) = a1(i,k)/r
         enddo

c        // update all other rows
         do j=1,n
            if (j.eq.i) goto 102
            tmp = 0d0
            do l=1,n
               tmp = tmp + a(l,i)*a1(j,l)
            enddo
            do k=1,n
               a1(j,k) = a1(j,k) - a1(i,k)*tmp
            enddo
 102        continue
         enddo

      enddo

      end

c================================================


c     -------------------------------
      subroutine inv(a,n,nmax,a1,det)
c     -------------------------------
            
c Inversion and calculation of determinant of matrix A

c input:
c a    : n x n matrix defined as a(nmax,nmax)

c output:
c a    : LU decomposed matrix A
c a1   : A**-1 inverse matrix of A
c det  : Det(A)

      use global, only: timesum6

      integer nmax,n,i,j
      integer indx(nmax)
      real*8 a(nmax,nmax),a1(nmax,nmax),det
      !real*8 vor,nach,cputime  

      call ludcmp(a,n,nmax,indx,det)
      do i=1,n
        det = det*a(i,i)
      enddo
      !vor = cputime()      
      do i=1,n
        do j=1,n
          a1(i,j) = 0d0
        enddo
        a1(i,i) = 1d0
      enddo
      do i=1,n
        call lubksb(a,n,nmax,indx,a1(1,i))
      enddo
      !nach = cputime()
      !timesum6 = timesum6 + (nach-vor)  

      return
      end




c ================================================

c     -------------------------------
      subroutine ludcmp(a,n,np,indx,d)
c     -------------------------------

c LU decomposition from Numerical Recipes (W.H.Press, 1989,p.35ff)

c input:
c a   : n x n matrix with defined as a(np,np) 

c output:
c a   : decomposed matrix
c indx: vector that records row permutation from pivoting
c d   : +/- 1 dep. on whether even/odd number of row changes

      integer n,np,indx(np)
      real*8 d,a(np,np),tiny
      parameter (tiny=1.0d-50)
      integer i,imax,j,k
      real*8 aamax,dum,sum,vv(np)

      d=1d0
      do 12 i=1,n
        aamax=0d0
        do 11 j=1,n
          if (abs(a(i,j)).gt.aamax) aamax=abs(a(i,j))
11      continue
        if (aamax.eq.0d0) call abortp('singular matrix in ludcmp')
        vv(i)=1d0/aamax
12    continue
      do 19 j=1,n
        do 14 i=1,j-1
          sum=a(i,j)
          do 13 k=1,i-1
            sum=sum-a(i,k)*a(k,j)
13        continue
          a(i,j)=sum
14      continue
        aamax=0d0
        do 16 i=j,n
          sum=a(i,j)
          do 15 k=1,j-1
            sum=sum-a(i,k)*a(k,j)
15        continue
          a(i,j)=sum
          dum=vv(i)*abs(sum)
          if (dum.ge.aamax) then
            imax=i
            aamax=dum
          endif
16      continue
        if (j.ne.imax)then
          do 17 k=1,n
            dum=a(imax,k)
            a(imax,k)=a(j,k)
            a(j,k)=dum
17        continue
          d=-d
          vv(imax)=vv(j)
        endif
        indx(j)=imax
        if(a(j,j).eq.0d0)a(j,j)=tiny
        if(j.ne.n)then
          dum=1d0/a(j,j)
          do 18 i=j+1,n
            a(i,j)=a(i,j)*dum
18        continue
        endif
19    continue

      return
      end
c  (c) copr. 1986-92 numerical recipes software 0#).


c ================================================

c     -------------------------------
      subroutine lubksb(a,n,np,indx,b)
c     -------------------------------

c backsubstitution from Numerical Recipes (W.H.Press, 1989,p.35ff)

c input:
c a   : n x n matrix with defined as a(np,np) (LU decomposed) 
c indx: vector that records row permutation from pivoting
c b   : rhs from A x = b

c output:
c a,n,indx: not modified
c b   : contains solution of A x = b

      integer n,np,indx(n)
      real*8 a(np,np),b(n)
      integer i,ii,j,ll
      real*8 sum
      ii=0
      do 12 i=1,n
        ll=indx(i)
        sum=b(ll)
        b(ll)=b(i)
        if (ii.ne.0)then
          do 11 j=ii,i-1
            sum=sum-a(i,j)*b(j)
11        continue
        else if (sum.ne.0d0) then
          ii=i
        endif
        b(i)=sum
12    continue
      do 14 i=n,1,-1
        sum=b(i)
        do 13 j=i+1,n
          sum=sum-a(i,j)*b(j)
13      continue
        b(i)=sum/a(i,i)
14    continue
      return
      end
c  (c) copr. 1986-92 numerical recipes software 0#).

      
