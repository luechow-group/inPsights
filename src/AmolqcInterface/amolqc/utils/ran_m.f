


      MODULE rannum

      ! based on ran2 from Numerical Recipes
      implicit none
      private
      public :: init_ran2, ran2_ran, ran2_gran, ran2_vran, ran2_vgran

      INTEGER idum
      INTEGER IM1,IM2,IMM1,IA1,IA2,IQ1,IQ2,IR1,IR2,NTAB,NDIV
      REAL*8 AM,EPS,RNMX
      PARAMETER (IM1=2147483563,IM2=2147483399,AM=1.d0/IM1,IMM1=IM1-1,
     *IA1=40014,IA2=40692,IQ1=53668,IQ2=52774,IR1=12211,IR2=3791,
     *NTAB=32,NDIV=1+IMM1/NTAB,EPS=1.2d-14,RNMX=1.d0-EPS)
      INTEGER idum2
      INTEGER iv(NTAB)
      INTEGER iy
      INTEGER iset
      INTEGER j,k
      REAL*8  gset

      CONTAINS

c     ------------------------------
      REAL*8 function init_ran2(seed)
c     ------------------------------

      INTEGER, intent(in) :: seed

      idum2 = 123456789
      iv(NTAB) = 0
      iy = 0
      iset = 0
      if (seed <= 0) then
         print*," seed must be positive"
         stop
      endif
      idum = -seed
      idum=max(-idum,1)
      idum2=idum
      do 11 j=NTAB+8,1,-1
         k=idum/IQ1
         idum=IA1*(idum-k*IQ1)-k*IR1
         if (idum.lt.0) idum=idum+IM1
         if (j.le.NTAB) iv(j)=idum
 11   continue
      iy=iv(1)

      init_ran2 = ran2_ran()

      end function init_ran2

c     ---------------------
      REAL*8 FUNCTION ran2_ran()
c     ---------------------

c slightly modified version of ran2 random number generator of Num. Recipes

      k=idum/IQ1
      idum=IA1*(idum-k*IQ1)-k*IR1
      if (idum.lt.0) idum=idum+IM1
      k=idum2/IQ2
      idum2=IA2*(idum2-k*IQ2)-k*IR2
      if (idum2.lt.0) idum2=idum2+IM2
      j=1+iy/NDIV
      iy=iv(j)-idum2
      iv(j)=idum
      if(iy.lt.1)iy=iy+IMM1
      ran2_ran=min(AM*iy,RNMX)
c      print*,'DBG:',ran2_ran
      end function ran2_ran
C  (C) Copr. 1986-92 Numerical Recipes Software 0#).


c     ----------------------
      REAL*8 FUNCTION ran2_gran()
c     ----------------------

c slightly modified version of gasdev (Gaussian deviates) from Num. Recipes

      REAL*8 fac,rsq,v1,v2

      if (iset.eq.0) then
1       v1=2.*ran2_ran()-1.
        v2=2.*ran2_ran()-1.
        rsq=v1**2+v2**2
        if(rsq.ge.1..or.rsq.eq.0.)goto 1
        fac=sqrt(-2.*log(rsq)/rsq)
        gset=v1*fac
        ran2_gran=v2*fac
        iset=1
      else
        ran2_gran=gset
        iset=0
      endif
      end function ran2_gran
C  (C) Copr. 1986-92 Numerical Recipes Software 0#).


c     --------------------
      subroutine ran2_vran(n,v)
c     --------------------

      integer, intent(in)   :: n
      real*8, intent(inout) :: v(:)
      integer i

      do i=1,n
         v(i) = ran2_ran()
      enddo

      end subroutine ran2_vran

c     ---------------------
      subroutine ran2_vgran(n,v)
c     ---------------------

      integer, intent(in)   :: n
      real*8, intent(inout) :: v(:)
      integer i

      do i=1,n
         v(i) = ran2_gran()
      enddo

      end subroutine ran2_vgran

      END MODULE rannum

