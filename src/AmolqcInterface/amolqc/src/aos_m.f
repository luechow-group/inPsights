
c F90 module for handling the atomic basis functions (AO's)

c use allocatable arrays for uao etc.
c factored out basis input and norm factors
c allocation, input and output moved to aosdata
c here only calculations based on aosdata
c AL, 2011

c combined with SM (26.6.00) version for linear scaling TS 2005/12/25

c $Id: aos_m.f,v 1.4 2008/02/07 16:58:33 luechow Exp $

c $Log: aos_m.f,v $
c Revision 1.4  2008/02/07 16:58:33  luechow
c correction
c
c Revision 1.3  2008/02/07 16:51:25  luechow
c doppelte Uebergabe (use und parameter) entfernt
c
c Revision 1.2  2007/05/31 08:28:11  luechow
c bug fix in aos_m.f and rwStatistics_m.f90. removed unused files
c
c Revision 1.1.1.1  2007/04/25 13:42:20  luechow
c QMC program amolqc. rewritten in 2006. AL
c
c Revision 1.12  2005/09/29 12:33:18  annika
c automatic basis set input from basis set library (aoinputex), automatic "so counting" (getso)
c
c Revision 1.11  2005/05/17 10:34:41  luechow
c reorganization of cubic spline and cusp correction.
c Both are in own modul and no longer part or aos
c Slight change in aospline to ignore unused so(bf) entries
c
c Revision 1.10  2005/03/11 13:33:11  luechow
c adding 'linscal' routines using sparse matrix algorithms with umfpack lib
c
c Revision 1.9  2004/10/11 11:51:28  annika
c Formatierungsfehler beseitigt
c
c Revision 1.8  2004/09/13 17:29:45  diedrch
c beautified output
c
c Revision 1.7  2004/08/12 17:37:10  diedrch
c Added routine ao1splcalc --> support fot ECP localisation with splines
c TS adapted for LINSCAL
c
c Revision 1.6  2004/06/25 14:22:04  diedrch
c reinserted the subroutines aocalc mocalc and aosplcalc from revision 1-0.
c Different names for the MO coefficient arrays prevent conflicts with aomo_calc
c
c Revision 1.5  2004/04/02 09:38:24  diedrch
c removed aocalc and ao1calc (now part of aomo_calc in aomo_m.f)
c
c Revision 1.4  2004/03/31 15:03:11  diedrch
c fixed segmentation fault for calculations without cutoff
c
c Revision 1.3  2004/03/30 15:01:08  diedrch
c allocate AO - cutoff array dynamically in order to keep it as small
c as possible.
c
c Revision 1.2  2004/03/25 17:34:37  diedrch
c - Some optimisation of the AO routines (precalculation of some products
c which are the same for each GTO in a contraction)
c - added support for AO - cutoff (aocalc, ao1calc)
c
c Revision 1.1.1.1  2003/12/28 15:08:07  diedrch
c Initial Code 281203
c
c Revision 2.0  1999/08/18 16:13:44  luechow
c initial f90
c

c Modifications:
c     09.09.1999 SM
c     Gaussianausgaben als AO-Eingabe moeglich
c     Basisfunktionen aus einem GTO werden nicht mehr ueber splines
c     berechnet sondern analytisch, das ist effektiver weil hier keine
c     Schleife ueber die Kontraktion benoetigt wird
c     Nur noch eine Spline-Tabelle fuer identische Basisfunktionen bei
c     gleichen Kernen, Zugriff erfolgt ueber Zuordnungstabelle so(basmax)


      MODULE aos

      use wfdata
      use aosdata
      use cubicspline
      use eConfigsModule
      implicit none

      CONTAINS


c     -----------------------------
      subroutine aocalc(ie,ec,rrai)
c     -----------------------------

c aocalc calculates all atomic orbitals for position vectors x,y,z
c (in configurations ec) including all required derivatives.
c This is the version that calculates several electron configurations
c Currently sequential calculation of electron configurations.

c on entry: requires rai distance matrix properly dimensioned and set
c           for all input position vector x,y,z.
c           for ie>0 position vector must be modified only compared to
c           last call to aocalc only at electron ie
c on exit:  updates AO data structure in aos_d (i.e. AOs and derivatives)

c Version 2.1 (28.4.98)    references to walker eliminated
c Version 2.0 (20.3.98)    module version
c Version 1.1 (6/25/1996)  deleted 5D code; one orb coeff for
c                          all p,d,f orbs; simultaneous evaluation
c                          for p,d,f orbs. contracted GTO's added.
c Version 1.0 (1/26/1996)

c input parameters:
      integer, intent(in)               :: ie              ! if >0 only AO's for electron ie recalculated
      type(eConfigArray), intent(inout) :: ec              ! electron configurations (intent is in)
      real*8, intent(in)                :: rrai(:,:,:)     ! r_ai electron-nucleus distances
c constants:
      real*8 sqr3,sqr5,sqr7
      parameter (sqr3=1.73205080756887729d0,sqr5=2.236067977499789696d0,
     &           sqr7=2.645751311064591d0)
c variables
      integer bf,a,i,i1,i2,j,nn,al,ii,ic,m,n,k,l,w
      integer fxxx,fyyy,fzzz,fxyy,fxxy,fxxz,fxzz,fyzz,fyyz,fxyz
      real*8 rai(ncenter,ne)
      real*8, pointer :: x(:),y(:),z(:)
      real*8 xx,yy,zz,rr,r2,alp,nrm,u,ux,dx,dy,dz,tmp,uux,
     .       dx2,dy2,dz2,dxyz,tmp1,tmp2
      logical gaussFOrder       ! .t.: Gaussian order for f function
                                ! .f.: Gamess==Turbomole order used

c bf refers to the degenerate set of cartesian
c basis function (S:1,P:3,D:6,F:10) as input, which may be of type STO
c or contracted GTO.
c al refers to the individual basis function, as used in LCAO-MO's.
c (composed in subroutine mdetwf)
c i refers to the current electron.

c-----Calculation of the AO's and their derivatives

      mAOElecConfigs = eConfigArray_size(ec)
      do w=1,mAOElecConfigs

      call eConfigArray_getPtr(ec,w,x,y,z)
      rai = rrai(:,:,w)

      ! check input data for NaN or Inf
      call assert(all(rai<huge(1.d0)),"aocalc: illegal rai values")
      call assert(all(abs(x(1:ne))<huge(1.d0)) .and.
     .  all(abs(y(1:ne))<huge(1.d0)) .and. all(abs(z(1:ne))<huge(1.d0)),
     .  "aocalc: illegal x,y,z coords in ec")


      if (evfmt=='gau' .or. evfmt=='mol' ) then
         gaussFOrder = .true.
      else
         gaussFOrder = .false.
      endif
      if (gaussFOrder) then
c                    // order: f_xxx, f_yyy, f_zzz, f_xyy, f_xxy, f_xxz,
c                    //   f_xzz, f_yzz, f_yyz, f_xyz  (like Gaussian)
         fxxx=0; fyyy=1; fzzz=2
         fxyy=3; fxxy=4; fxxz=5
         fxzz=6; fyzz=7; fyyz=8; fxyz=9
      else
c                    // order: f_xxx, f_yyy, f_zzz, f_xxy, f_xxz, f_yyx,
c                    //   f_yyz, f_zzx, f_zzy, f_xyz  (like GAMESS)
         fxxx=0; fyyy=1; fzzz=2
         fxxy=3; fxxz=4; fxyy=5
         fyyz=6; fxzz=7; fyzz=8; fxyz=9
      endif

      if (ie .eq. 0) then                     ! AO's for all electrons
         i1 = 1
         i2 = ne
      else
         i1 = ie                              ! only AO for electron ie
         i2 = ie
      endif

      do i=i1,i2                              ! loop over electrons

         xx = x(i)
         yy = y(i)
         zz = z(i)

         al = 1

         ! loop over basis functions:
         do n = 1,nbasf
           bf = n
           a = bc(bf)                         ! center of AO
           rr = rai(a,i)                      ! r_ai
           nn = bn(bf)                        ! n quantum no. of AO

           if (typ(bf) .eq. 'STO') then      ! STO-type basis function

            alp = bzet(bf)                    ! orbital exponent
            nrm = norm(bf)                    ! normalization constant

            if (nn .eq. 1) then               ! 1s orbital
               u  = nrm*exp(-alp*rr)
               ux = -alp*u/rr
               uao(al,i,w) = u
               uxao(al,i,w) = ux*(xx-atoms(a)%cx)
               uyao(al,i,w) = ux*(yy-atoms(a)%cy)
               uzao(al,i,w) = ux*(zz-atoms(a)%cz)
               u2ao(al,i,w) = alp*(alp-2d0/rr)*u
               al = al +1

            else if (nn .eq. 2) then
               if (bl(bf) .eq. 'S') then           ! 2s orbital
c                 // 2s orbital
                  u  = nrm*exp(-alp*rr)
                  ux = (1d0-alp*rr)*u/rr
                  uao(al,i,w) = rr*u
                  uxao(al,i,w) = ux*(xx-atoms(a)%cx)
                  uyao(al,i,w) = ux*(yy-atoms(a)%cy)
                  uzao(al,i,w) = ux*(zz-atoms(a)%cz)
                  u2ao(al,i,w) = alp*(alp*rr-2d0)*u + 2d0*ux
                  al = al +1
               else if (bl(bf) .eq. 'P') then      ! 2p orbital
                  u  = nrm*exp(-alp*rr)
                  dx = xx-atoms(a)%cx
                  dy = yy-atoms(a)%cy
                  dz = zz-atoms(a)%cz
                  ux = -alp*u/rr
                  uao(al,i,w)    = u*dx
                  uao(al+1,i,w)  = u*dy
                  uao(al+2,i,w)  = u*dz
                  uxao(al,i,w)   = ux*dx*dx + u
                  uxao(al+1,i,w) = ux*dx*dy
                  uxao(al+2,i,w) = ux*dx*dz
                  uyao(al,i,w)   = uxao(al+1,i,w)
                  uyao(al+1,i,w) = ux*dy*dy + u
                  uyao(al+2,i,w) = ux*dy*dz
                  uzao(al,i,w)   = uxao(al+2,i,w)
                  uzao(al+1,i,w) = uyao(al+2,i,w)
                  uzao(al+2,i,w) = ux*dz*dz + u
                  tmp          = alp*(alp-4d0/rr)*u
                  u2ao(al,i,w)   = tmp*dx
                  u2ao(al+1,i,w) = tmp*dy
                  u2ao(al+2,i,w) = tmp*dz
                  al = al +3
               else
                  call abortp(' getaos: n=2 and l .gt. 2')
               endif

            else if (nn .eq. 3) then
               if (bl(bf) .eq. 'S') then
c                 // 3s orbital
                  u  = nrm*exp(-alp*rr)
                  ux = (2d0-alp*rr)*u
                  uao(al,i,w) = rr*rr*u
                  uxao(al,i,w) = ux*(xx-atoms(a)%cx)
                  uyao(al,i,w) = ux*(yy-atoms(a)%cy)
                  uzao(al,i,w) = ux*(zz-atoms(a)%cz)
                  tmp = (2d0 + alp*rr*(-4d0 + alp*rr))*u
                  u2ao(al,i,w) = tmp + 2d0*ux
                  al = al +1
               else if (bl(bf) .eq. 'P') then
c                 //3p orbital
                  u = nrm*exp(-alp*rr)
                  dx = xx-atoms(a)%cx
                  dx2 = dx*dx
                  dy = yy-atoms(a)%cy
                  dy2= dy*dy
                  dz = zz-atoms(a)%cz
                  dz2 = dz*dz
                  uao(al,i,w)    = u*rr*dx
                  uao(al+1,i,w)  = u*rr*dy
                  uao(al+2,i,w)  = u*rr*dz
                  uxao(al,i,w)   = -alp*u*dx2 + u*dx2/rr + u*rr
                  uxao(al+1,i,w) = u*dx*dy/rr - alp*u*dx*dy
                  uxao(al+2,i,w) = u*dx*dz/rr - alp*u*dx*dz
                  uyao(al,i,w)   = uxao(al+1,i,w)
                  uyao(al+1,i,w) = -alp*u*dy2 + u*dy2/rr + u*rr
                  uyao(al+2,i,w) = u*dy*dz/rr - alp*u*dy*dz
                  uzao(al,i,w)   = uxao(al+2,i,w)
                  uzao(al+1,i,w) = uyao(al+2,i,w)
                  uzao(al+2,i,w) = -alp*u*dz2 + u*dz2/rr + u*rr
                  tmp1 = alp*alp*u/rr
                  tmp2 = tmp1*dx2 + tmp1*dy2 + tmp1*dz2 - 6*alp*u + 4*u/rr
                  u2ao(al,i,w)   = tmp2*dx
                  u2ao(al+1,i,w) = tmp2*dy
                  u2ao(al+2,i,w) = tmp2*dz
                  al = al +3
               else if (bl(bf) .eq. 'D') then         ! 3d orbital (6D)
c                 // Norm is different for d_xx and d_xy !
                  u   = nrm*exp(-alp*rr)
                  dx  = xx - atoms(a)%cx
                  dx2 = dx*dx
                  dy  = yy - atoms(a)%cy
                  dy2 = dy*dy
                  dz  = zz - atoms(a)%cz
                  dz2 = dz*dz
                  ux  = -alp*u/rr
                  tmp = alp*(alp-6d0/rr)

                  uao(al,i,w)    = dx2*u
                  uao(al+1,i,w)  = dy2*u
                  uao(al+2,i,w)  = dz2*u
                  uxao(al,i,w)   = dx*(2d0*u + ux*dx2)
                  uxao(al+1,i,w) = ux*dy2*dx
                  uxao(al+2,i,w) = ux*dz2*dx
                  uyao(al,i,w)   = ux*dx2*dy
                  uyao(al+1,i,w) = dy*(2d0*u + ux*dy2)
                  uyao(al+2,i,w) = ux*dz2*dy
                  uzao(al,i,w)   = ux*dx2*dz
                  uzao(al+1,i,w) = ux*dy2*dz
                  uzao(al+2,i,w) = dz*(2d0*u + ux*dz2)
                  u2ao(al,i,w)   = u*(2d0+tmp*dx2)
                  u2ao(al+1,i,w) = u*(2d0+tmp*dy2)
                  u2ao(al+2,i,w) = u*(2d0+tmp*dz2)

                  u = sqr3*u                   ! correction of norm for last 3
                  ux = sqr3*ux

                  uao(al+3,i,w)  = u*dx*dy
                  uao(al+4,i,w)  = u*dx*dz
                  uao(al+5,i,w)  = u*dy*dz

                  uxao(al+3,i,w) = dy*(u + ux*dx2)
                  uxao(al+4,i,w) = dz*(u + ux*dx2)
                  uxao(al+5,i,w) = ux*dx*dy*dz
                  uyao(al+3,i,w) = dx*(u + ux*dy2)
                  uyao(al+4,i,w) = uxao(al+5,i,w)
                  uyao(al+5,i,w) = dz*(u + ux*dy2)
                  uzao(al+3,i,w) = uxao(al+5,i,w)
                  uzao(al+4,i,w) = dx*(u + ux*dz2)
                  uzao(al+5,i,w) = dy*(u + ux*dz2)
                  u2ao(al+3,i,w) = u*dx*dy*tmp
                  u2ao(al+4,i,w) = u*dx*dz*tmp
                  u2ao(al+5,i,w) = u*dy*dz*tmp
                  al = al +6
               endif     ! bl

            else if (nn == 4) then
               if (bl(bf) == 'F') then ! 4f orbital (10F)
                  u   = nrm*exp(-alp*rr)
                  dx  = xx - atoms(a)%cx
                  dx2 = dx*dx
                  dy  = yy - atoms(a)%cy
                  dy2 = dy*dy
                  dz  = zz - atoms(a)%cz
                  dz2 = dz*dz
                  dxyz = dx*dy*dz
                  ux  = -alp*u/rr

c                 // f_xxx, f_yyy, f_zzz
                  uao(al+fxxx,i,w)  = dx2*dx*u
                  uao(al+fyyy,i,w)  = dy2*dy*u
                  uao(al+fzzz,i,w)  = dz2*dz*u
                  uxao(al+fxxx,i,w) = (3d0*u + ux*dx2)*dx2
                  uxao(al+fyyy,i,w) = dy2*dy*ux*dx
                  uxao(al+fzzz,i,w) = dz2*dz*ux*dx
                  uyao(al+fxxx,i,w) = dx2*dx*ux*dy
                  uyao(al+fyyy,i,w) = (3d0*u + ux*dy2)*dy2
                  uyao(al+fzzz,i,w) = dz2*dz*ux*dy
                  uzao(al+fxxx,i,w) = dx2*dx*ux*dz
                  uzao(al+fyyy,i,w) = dy2*dy*ux*dz
                  uzao(al+fzzz,i,w) = (3d0*u + ux*dz2)*dz2

                  tmp = 8*ux + alp*alp*u

                  u2ao(al+fxxx,i,w) = dx*(dx2*tmp + 6*u)
                  u2ao(al+fyyy,i,w) = dy*(dy2*tmp + 6*u)
                  u2ao(al+fzzz,i,w) = dz*(dz2*tmp + 6*u)

c                 // f_xxy, f_xxz, f_yyx, f_yyz, f_zzx, f_zzy
                  u = sqr5*u                   ! correction of norm
                  ux = sqr5*ux
                  uux = alp*alp*u + 8*ux

                  tmp = dx2*dy
                  uao(al+fxxy,i,w)  = tmp*u
                  uxao(al+fxxy,i,w) = ux*tmp*dx + 2*dx*dy*u
                  uyao(al+fxxy,i,w) = ux*tmp*dy + dx2*u
                  uzao(al+fxxy,i,w) = ux*tmp*dz
                  u2ao(al+fxxy,i,w) = tmp*uux + 2*dy*u

                  tmp = dx2*dz
                  uao(al+fxxz,i,w)  = tmp*u
                  uxao(al+fxxz,i,w) = ux*tmp*dx + 2*dx*dz*u
                  uyao(al+fxxz,i,w) = ux*tmp*dy
                  uzao(al+fxxz,i,w) = ux*tmp*dz + dx2*u
                  u2ao(al+fxxz,i,w) = tmp*uux + 2*dz*u

                  tmp = dy2*dx
                  uao(al+fxyy,i,w)  = tmp*u
                  uxao(al+fxyy,i,w) = ux*tmp*dx + dy2*u
                  uyao(al+fxyy,i,w) = ux*tmp*dy + 2*dy*dx*u
                  uzao(al+fxyy,i,w) = ux*tmp*dz
                  u2ao(al+fxyy,i,w) = tmp*uux + 2*dx*u

                  tmp = dy2*dz
                  uao(al+fyyz,i,w)  = tmp*u
                  uxao(al+fyyz,i,w) = ux*tmp*dx
                  uyao(al+fyyz,i,w) = ux*tmp*dy + 2*dy*dz*u
                  uzao(al+fyyz,i,w) = ux*tmp*dz + dy2*u
                  u2ao(al+fyyz,i,w) = tmp*uux + 2*dz*u

                  tmp = dz2*dx
                  uao(al+fxzz,i,w)  = tmp*u
                  uxao(al+fxzz,i,w) = ux*tmp*dx + dz2*u
                  uyao(al+fxzz,i,w) = ux*tmp*dy
                  uzao(al+fxzz,i,w) = ux*tmp*dz + 2*dx*dz*u
                  u2ao(al+fxzz,i,w) = tmp*uux + 2*dx*u

                  tmp = dz2*dy
                  uao(al+fyzz,i,w)  = tmp*u
                  uxao(al+fyzz,i,w) = ux*tmp*dx
                  uyao(al+fyzz,i,w) = ux*tmp*dy + dz2*u
                  uzao(al+fyzz,i,w) = ux*tmp*dz + 2*dy*dz*u
                  u2ao(al+fyzz,i,w) = tmp*uux + 2*dy*u

c                 // f_xyz
                  u = sqr3*u                  ! correction of norm
                  ux = sqr3*ux

                  uao(al+fxyz,i,w)  = dxyz*u
                  uxao(al+fxyz,i,w) = dxyz*(ux*dx + u/dx)
                  uyao(al+fxyz,i,w) = dxyz*(ux*dy + u/dy)
                  uzao(al+fxyz,i,w) = dxyz*(ux*dz + u/dz)
                  u2ao(al+fxyz,i,w) = dxyz*(alp*alp*u+8*ux)
                  al = al + 10
               else
                  call abortp("aocalc: currently only 4f implemented")
               endif     ! bl

            endif        ! nn


c          // Contracted GTO's as basis function (AO)
           else

            r2 = rr*rr

c           // only primitive cartesian gaussians: 1s,2p,3d,4f,5g
c           // i.e. no r factor. Thus nn is not used here.

            if (bl(bf) .eq. 'S') then                 ! 1s GTO
               uao(al,i,w) = 0d0
               uxao(al,i,w) = 0d0
               uyao(al,i,w) = 0d0
               uzao(al,i,w) = 0d0
               u2ao(al,i,w) = 0d0
               do ic=1,ngto(bf)                     ! loop over contraction
                  alp = cntrctn(1,ic,bf)
                  u = cntrctn(2,ic,bf) * exp(-alp*r2)
                  ux = -2d0*alp*u
                  uao(al,i,w) = uao(al,i,w) + u
                  uxao(al,i,w) = uxao(al,i,w) + ux*(xx-atoms(a)%cx)
                  uyao(al,i,w) = uyao(al,i,w) + ux*(yy-atoms(a)%cy)
                  uzao(al,i,w) = uzao(al,i,w) + ux*(zz-atoms(a)%cz)
                  u2ao(al,i,w) = u2ao(al,i,w) + ux*(3d0-2d0*alp*r2)
               enddo
               al = al +1

            else if (bl(bf) .eq. 'P') then             ! 2p GTO's
c              // do all 3 P simultaneously (same exponent is required)
c              // order p_x,p_y,p_z
               do ii=0,2
                  uao(al+ii,i,w) = 0d0
                  uxao(al+ii,i,w) = 0d0
                  uyao(al+ii,i,w) = 0d0
                  uzao(al+ii,i,w) = 0d0
                  u2ao(al+ii,i,w) = 0d0
               enddo
               do ic=1,ngto(bf)                      ! loop over contraction
                  alp = cntrctn(1,ic,bf)
                  u = cntrctn(2,ic,bf) * exp(-alp*r2)
                  dx = xx-atoms(a)%cx
                  dy = yy-atoms(a)%cy
                  dz = zz-atoms(a)%cz
                  ux = -2d0*alp*u
                  uao(al,i,w) = uao(al,i,w) + dx*u
                  uao(al+1,i,w) = uao(al+1,i,w) + dy*u
                  uao(al+2,i,w) = uao(al+2,i,w) + dz*u
                  uxao(al,i,w) = uxao(al,i,w) +  u + ux*dx*dx
                  uxao(al+1,i,w) = uxao(al+1,i,w) + ux*dx*dy
                  uxao(al+2,i,w) = uxao(al+2,i,w) + ux*dx*dz
                  uyao(al,i,w) = uyao(al,i,w) + ux*dx*dy
                  uyao(al+1,i,w) = uyao(al+1,i,w) + u + ux*dy*dy
                  uyao(al+2,i,w) = uyao(al+2,i,w) + ux*dy*dz
                  uzao(al,i,w) = uzao(al,i,w) + ux*dx*dz
                  uzao(al+1,i,w) = uzao(al+1,i,w) + ux*dy*dz
                  uzao(al+2,i,w) = uzao(al+2,i,w) + u + ux*dz*dz
                  tmp = (5d0-2d0*alp*r2)*ux
                  u2ao(al,i,w) = u2ao(al,i,w) + tmp*dx
                  u2ao(al+1,i,w) = u2ao(al+1,i,w) + tmp*dy
                  u2ao(al+2,i,w) = u2ao(al+2,i,w) + tmp*dz
               enddo
	           al = al +3

            else if (bl(bf) .eq. 'D') then ! 3d GTO
c              // do all 6 D simultaneously (same exponent is required)
c              // order: d_xx, d_yy, d_zz, d_xy, d_xz, d_yz  (like GAMESS)
               do ii=0,5
                  uao(al+ii,i,w) = 0d0
                  uxao(al+ii,i,w) = 0d0
                  uyao(al+ii,i,w) = 0d0
                  uzao(al+ii,i,w) = 0d0
                  u2ao(al+ii,i,w) = 0d0
               enddo
               do ic=1,ngto(bf)                      ! loop over contraction
                  alp = cntrctn(1,ic,bf)
                  u = cntrctn(2,ic,bf) * exp(-alp*r2)
                  dx = xx-atoms(a)%cx
                  dx2 = dx*dx
                  dy = yy-atoms(a)%cy
                  dy2 = dy*dy
                  dz = zz-atoms(a)%cz
                  dz2 = dz*dz
                  ux = -2d0*alp*u

                  uao(al,i,w)    = uao(al,i,w)    + dx2*u
                  uao(al+1,i,w)  = uao(al+1,i,w)  + dy2*u
                  uao(al+2,i,w)  = uao(al+2,i,w)  + dz2*u

                  uxao(al,i,w)   = uxao(al,i,w)   + (2d0*u + ux*dx2)*dx
                  uxao(al+1,i,w) = uxao(al+1,i,w) + dy2*ux*dx
                  uxao(al+2,i,w) = uxao(al+2,i,w) + dz2*ux*dx
                  uyao(al,i,w)   = uyao(al,i,w)   + dx2*ux*dy
                  uyao(al+1,i,w) = uyao(al+1,i,w) + (2d0*u + ux*dy2)*dy
                  uyao(al+2,i,w) = uyao(al+2,i,w) + dz2*ux*dy
                  uzao(al,i,w)   = uzao(al,i,w)   + dx2*ux*dz
                  uzao(al+1,i,w) = uzao(al+1,i,w) + dy2*ux*dz
                  uzao(al+2,i,w) = uzao(al+2,i,w) + (2d0*u + ux*dz2)*dz
                  tmp          = (7d0 - 2d0*alp*r2)*ux
                  u2ao(al,i,w)   = u2ao(al,i,w)   + 2d0*u + dx2*tmp
                  u2ao(al+1,i,w) = u2ao(al+1,i,w) + 2d0*u + dy2*tmp
                  u2ao(al+2,i,w) = u2ao(al+2,i,w) + 2d0*u + dz2*tmp

                  u = sqr3*u                   ! correction of norm for last 3
                  ux = sqr3*ux

                  uao(al+3,i,w)  = uao(al+3,i,w)  + dx*dy*u
                  uao(al+4,i,w)  = uao(al+4,i,w)  + dx*dz*u
                  uao(al+5,i,w)  = uao(al+5,i,w)  + dy*dz*u

                  tmp = ux*dx*dy*dz
                  uxao(al+3,i,w) = uxao(al+3,i,w) + (u + ux*dx2)*dy
                  uxao(al+4,i,w) = uxao(al+4,i,w) + (u + ux*dx2)*dz
                  uxao(al+5,i,w) = uxao(al+5,i,w) + tmp
                  uyao(al+3,i,w) = uyao(al+3,i,w) + (u + ux*dy2)*dx
                  uyao(al+4,i,w) = uyao(al+4,i,w) + tmp
                  uyao(al+5,i,w) = uyao(al+5,i,w) + (u + ux*dy2)*dz
                  uzao(al+3,i,w) = uzao(al+3,i,w) + tmp
                  uzao(al+4,i,w) = uzao(al+4,i,w) + (u + ux*dz2)*dx
                  uzao(al+5,i,w) = uzao(al+5,i,w) + (u + ux*dz2)*dy

                  tmp = (7d0 - 2d0*alp*r2)*ux
                  u2ao(al+3,i,w) = u2ao(al+3,i,w) + tmp*dx*dy
                  u2ao(al+4,i,w) = u2ao(al+4,i,w) + tmp*dx*dz
                  u2ao(al+5,i,w) = u2ao(al+5,i,w) + tmp*dy*dz
               enddo
               al = al +6

            else if (bl(bf)=='F') then     ! 4f GTO
c              // do all 10 F simultaneously (same exponent is required)
c              // order: f_xxx, f_yyy, f_zzz, fd_xxy, f_xxz, f_yyx,
c              //   f_yyz, f_zzx, f_zzy, f_xyz  (like GAMESS)
               do ii=0,9
                  uao(al+ii,i,w) = 0d0
                  uxao(al+ii,i,w) = 0d0
                  uyao(al+ii,i,w) = 0d0
                  uzao(al+ii,i,w) = 0d0
                  u2ao(al+ii,i,w) = 0d0
               enddo
               do ic=1,ngto(bf)                      ! loop over contraction
                  alp = cntrctn(1,ic,bf)
                  u = cntrctn(2,ic,bf) * exp(-alp*r2)
                  dx = xx-atoms(a)%cx
                  dx2 = dx*dx
                  dy = yy-atoms(a)%cy
                  dy2 = dy*dy
                  dz = zz-atoms(a)%cz
                  dz2 = dz*dz
                  dxyz = dx*dy*dz
                  ux = -2d0*alp*u

c                 // f_xxx, f_yyy, f_zzz
                  uao(al+fxxx,i,w)    = uao(al+fxxx,i,w)    + dx2*dx*u
                  uao(al+fyyy,i,w)  = uao(al+fyyy,i,w)  + dy2*dy*u
                  uao(al+fzzz,i,w)  = uao(al+fzzz,i,w)  + dz2*dz*u

                  uxao(al+fxxx,i,w)   = uxao(al+fxxx,i,w)   + (3d0*u + ux*dx2)*dx2
                  uxao(al+fyyy,i,w) = uxao(al+fyyy,i,w) + dy2*dy*ux*dx
                  uxao(al+fzzz,i,w) = uxao(al+fzzz,i,w) + dz2*dz*ux*dx
                  uyao(al+fxxx,i,w)   = uyao(al+fxxx,i,w)   + dx2*dx*ux*dy
                  uyao(al+fyyy,i,w) = uyao(al+fyyy,i,w) + (3d0*u + ux*dy2)*dy2
                  uyao(al+fzzz,i,w) = uyao(al+fzzz,i,w) + dz2*dz*ux*dy
                  uzao(al+fxxx,i,w)   = uzao(al+fxxx,i,w)   + dx2*dx*ux*dz
                  uzao(al+fyyy,i,w) = uzao(al+fyyy,i,w) + dy2*dy*ux*dz
                  uzao(al+fzzz,i,w) = uzao(al+fzzz,i,w) + (3d0*u + ux*dz2)*dz2
                  tmp          = (9d0 - 2d0*alp*r2)*ux
                  u2ao(al+fxxx,i,w)   = u2ao(al+fxxx,i,w)   + (6d0*u + dx2*tmp)*dx
                  u2ao(al+fyyy,i,w) = u2ao(al+fyyy,i,w) + (6d0*u + dy2*tmp)*dy
                  u2ao(al+fzzz,i,w) = u2ao(al+fzzz,i,w) + (6d0*u + dz2*tmp)*dz

c                 // f_xxy, f_xxz, f_yyx, f_yyz, f_zzx, f_zzy
                  u = sqr5*u                   ! correction of norm
                  ux = sqr5*ux

                  uao(al+fxxy,i,w)  = uao(al+fxxy,i,w)  + dx2*dy*u
                  uao(al+fxxz,i,w)  = uao(al+fxxz,i,w)  + dx2*dz*u
                  uao(al+fxyy,i,w)  = uao(al+fxyy,i,w)  + dy2*dx*u
                  uao(al+fyyz,i,w)  = uao(al+fyyz,i,w)  + dy2*dz*u
                  uao(al+fxzz,i,w)  = uao(al+fxzz,i,w)  + dz2*dx*u
                  uao(al+fyzz,i,w)  = uao(al+fyzz,i,w)  + dz2*dy*u

                  tmp = ux*dxyz
                  uxao(al+fxxy,i,w) = uxao(al+fxxy,i,w) + (2d0*u + ux*dx2)*dx*dy
                  uxao(al+fxxz,i,w) = uxao(al+fxxz,i,w) + (2d0*u + ux*dx2)*dx*dz
                  uxao(al+fxyy,i,w) = uxao(al+fxyy,i,w) + (u + ux*dx2)*dy2
                  uxao(al+fyyz,i,w) = uxao(al+fyyz,i,w) + tmp*dy
                  uxao(al+fxzz,i,w) = uxao(al+fxzz,i,w) + (u + ux*dx2)*dz2
                  uxao(al+fyzz,i,w) = uxao(al+fyzz,i,w) + tmp*dz
                  uyao(al+fxxy,i,w) = uyao(al+fxxy,i,w) + (u + ux*dy2)*dx2
                  uyao(al+fxxz,i,w) = uyao(al+fxxz,i,w) + tmp*dx
                  uyao(al+fxyy,i,w) = uyao(al+fxyy,i,w) + (2d0*u + ux*dy2)*dx*dy
                  uyao(al+fyyz,i,w) = uyao(al+fyyz,i,w) + (2d0*u + ux*dy2)*dy*dz
                  uyao(al+fxzz,i,w) = uyao(al+fxzz,i,w) + tmp*dz
                  uyao(al+fyzz,i,w) = uyao(al+fyzz,i,w) + (u + ux*dy2)*dz2
                  uzao(al+fxxy,i,w) = uzao(al+fxxy,i,w) + tmp*dx
                  uzao(al+fxxz,i,w) = uzao(al+fxxz,i,w) + (u + ux*dz2)*dx2
                  uzao(al+fxyy,i,w) = uzao(al+fxyy,i,w) + tmp*dy
                  uzao(al+fyyz,i,w) = uzao(al+fyyz,i,w) + (u + ux*dz2)*dy2
                  uzao(al+fxzz,i,w) = uzao(al+fxzz,i,w) + (2d0*u + ux*dz2)*dx*dz
                  uzao(al+fyzz,i,w) = uzao(al+fyzz,i,w) + (2d0*u + ux*dz2)*dy*dz

                  tmp = (9d0 - 2d0*alp*r2)*ux
                  u2ao(al+fxxy,i,w) = u2ao(al+fxxy,i,w) + (2d0*u + dx2*tmp)*dy
                  u2ao(al+fxxz,i,w) = u2ao(al+fxxz,i,w) + (2d0*u + dx2*tmp)*dz
                  u2ao(al+fxyy,i,w) = u2ao(al+fxyy,i,w) + (2d0*u + dy2*tmp)*dx
                  u2ao(al+fyyz,i,w) = u2ao(al+fyyz,i,w) + (2d0*u + dy2*tmp)*dz
                  u2ao(al+fxzz,i,w) = u2ao(al+fxzz,i,w) + (2d0*u + dz2*tmp)*dx
                  u2ao(al+fyzz,i,w) = u2ao(al+fyzz,i,w) + (2d0*u + dz2*tmp)*dy

c                 // f_xyz
                  u = sqr3*u                  ! correction of norm
                  ux = sqr3*ux

                  uao(al+fxyz,i,w)  = uao(al+fxyz,i,w)  + dxyz*u

                  uxao(al+fxyz,i,w) = uxao(al+fxyz,i,w) + (u + ux*dx2)*dy*dz
                  uyao(al+fxyz,i,w) = uyao(al+fxyz,i,w) + (u + ux*dy2)*dx*dz
                  uzao(al+fxyz,i,w) = uzao(al+fxyz,i,w) + (u + ux*dz2)*dx*dy

                  tmp = (9d0 - 2d0*alp*r2)*ux
                  u2ao(al+fxyz,i,w) = u2ao(al+fxyz,i,w) + dxyz*tmp
               enddo
               al = al + 10

            else if (bl(bf)=='G') then     ! 5g GTO
c              // do all 15 cartesian G simultaneously (same exponent is required)

               uao(al:al+14,i,w) = 0d0
               uxao(al:al+14,i,w) = 0d0
               uyao(al:al+14,i,w) = 0d0
               uzao(al:al+14,i,w) = 0d0
               u2ao(al:al+14,i,w) = 0d0

               if (gaussFOrder) then
                  call internal_GaussianOrderGFunctions()
               else
                  call internal_GamessOrderGFunctions()
               endif

               al = al + 15

            else
               call abortp('(getaos): wrong GTO')
            endif  ! bl
           endif  ! STO/GTO
         enddo  ! bf-loop over basis functions
      enddo  ! i-loop over electrons

      enddo  ! w-loop over electron configurations


      CONTAINS

         subroutine internal_GamessOrderGFunctions()

            do ic=1,ngto(bf)                      ! loop over contraction
               alp = cntrctn(1,ic,bf)
               u = cntrctn(2,ic,bf) * exp(-alp*r2)
               dx = xx-atoms(a)%cx
               dx2 = dx*dx
               dy = yy-atoms(a)%cy
               dy2 = dy*dy
               dz = zz-atoms(a)%cz
               dz2 = dz*dz
               dxyz = dx*dy*dz
               ux = -2d0*alp*u

c              // g_xxxx, g_yyyy, g_zzzz
               uao(al,i,w)    = uao(al,i,w)    + dx2*dx2*u
               uao(al+1,i,w)  = uao(al+1,i,w)  + dy2*dy2*u
               uao(al+2,i,w)  = uao(al+2,i,w)  + dz2*dz2*u

               uxao(al,i,w)   = uxao(al,i,w)   + (4d0*u + ux*dx2)*dx2*dx
               uxao(al+1,i,w) = uxao(al+1,i,w) + dy2*dy2*ux*dx
               uxao(al+2,i,w) = uxao(al+2,i,w) + dz2*dz2*ux*dx
               uyao(al,i,w)   = uyao(al,i,w)   + dx2*dx2*ux*dy
               uyao(al+1,i,w) = uyao(al+1,i,w) + (4d0*u + ux*dy2)*dy2*dy
               uyao(al+2,i,w) = uyao(al+2,i,w) + dz2*dz2*ux*dy
               uzao(al,i,w)   = uzao(al,i,w)   + dx2*dx2*ux*dz
               uzao(al+1,i,w) = uzao(al+1,i,w) + dy2*dy2*ux*dz
               uzao(al+2,i,w) = uzao(al+2,i,w) + (4d0*u + ux*dz2)*dz2*dz
               tmp          = (11d0 - 2d0*alp*r2)*ux
               u2ao(al,i,w)   = u2ao(al,i,w)   + (12d0*u + dx2*tmp)*dx2
               u2ao(al+1,i,w) = u2ao(al+1,i,w) + (12d0*u + dy2*tmp)*dy2
               u2ao(al+2,i,w) = u2ao(al+2,i,w) + (12d0*u + dz2*tmp)*dz2

c              // g_xxxy, g_xxxz, g_yyyx, g_yyyz, g_zzzx, g_zzzy
               u = sqr7*u                   ! correction of norm
               ux = sqr7*ux

               uao(al+3,i,w)  = uao(al+3,i,w)  + dx2*dx*dy*u
               uao(al+4,i,w)  = uao(al+4,i,w)  + dx2*dx*dz*u
               uao(al+5,i,w)  = uao(al+5,i,w)  + dy2*dy*dx*u
               uao(al+6,i,w)  = uao(al+6,i,w)  + dy2*dy*dz*u
               uao(al+7,i,w)  = uao(al+7,i,w)  + dz2*dz*dx*u
               uao(al+8,i,w)  = uao(al+8,i,w)  + dz2*dz*dy*u

               tmp = ux*dxyz
               uxao(al+3,i,w) = uxao(al+3,i,w) + (3d0*u + ux*dx2)*dx2*dy
               uxao(al+4,i,w) = uxao(al+4,i,w) + (3d0*u + ux*dx2)*dx2*dz
               uxao(al+5,i,w) = uxao(al+5,i,w) + (u + ux*dx2)*dy2*dy
               uxao(al+6,i,w) = uxao(al+6,i,w) + tmp*dy2
               uxao(al+7,i,w) = uxao(al+7,i,w) + (u + ux*dx2)*dz2*dz
               uxao(al+8,i,w) = uxao(al+8,i,w) + tmp*dz2
               uyao(al+3,i,w) = uyao(al+3,i,w) + (u + ux*dy2)*dx2*dx
               uyao(al+4,i,w) = uyao(al+4,i,w) + tmp*dx2
               uyao(al+5,i,w) = uyao(al+5,i,w) + (3d0*u + ux*dy2)*dy2*dx
               uyao(al+6,i,w) = uyao(al+6,i,w) + (3d0*u + ux*dy2)*dy2*dz
               uyao(al+7,i,w) = uyao(al+7,i,w) + tmp*dz2
               uyao(al+8,i,w) = uyao(al+8,i,w) + (u + ux*dy2)*dz2*dz
               uzao(al+3,i,w) = uzao(al+3,i,w) + tmp*dx2
               uzao(al+4,i,w) = uzao(al+4,i,w) + (u + ux*dz2)*dx2*dx
               uzao(al+5,i,w) = uzao(al+5,i,w) + tmp*dy2
               uzao(al+6,i,w) = uzao(al+6,i,w) + (u + ux*dz2)*dy2*dy
               uzao(al+7,i,w) = uzao(al+7,i,w) + (3d0*u + ux*dz2)*dz2*dx
               uzao(al+8,i,w) = uzao(al+8,i,w) + (3d0*u + ux*dz2)*dz2*dy

               tmp = (11d0 - 2d0*alp*r2)*ux
               u2ao(al+3,i,w) = u2ao(al+3,i,w) + (6d0*u + dx2*tmp)*dx*dy
               u2ao(al+4,i,w) = u2ao(al+4,i,w) + (6d0*u + dx2*tmp)*dx*dz
               u2ao(al+5,i,w) = u2ao(al+5,i,w) + (6d0*u + dy2*tmp)*dy*dx
               u2ao(al+6,i,w) = u2ao(al+6,i,w) + (6d0*u + dy2*tmp)*dy*dz
               u2ao(al+7,i,w) = u2ao(al+7,i,w) + (6d0*u + dz2*tmp)*dz*dx
               u2ao(al+8,i,w) = u2ao(al+8,i,w) + (6d0*u + dz2*tmp)*dz*dy

c              // g_xxyy, g_xxzz, g_yyzz
               u = sqr5 / sqr3 * u          ! correction of norm
               ux = sqr5 / sqr3 * ux

               uao(al+9,i,w)   = uao(al+9,i,w)  + dx2*dy2*u
               uao(al+10,i,w)  = uao(al+10,i,w)  + dx2*dz2*u
               uao(al+11,i,w)  = uao(al+11,i,w)  + dy2*dz2*u

               tmp = ux*dxyz
               uxao(al+9,i,w)  = uxao(al+9,i,w)  + (2d0*u + ux*dx2)*dx*dy2
               uxao(al+10,i,w) = uxao(al+10,i,w) + (2d0*u + ux*dx2)*dx*dz2
               uxao(al+11,i,w) = uxao(al+11,i,w) + tmp*dy*dz
               uyao(al+9,i,w)  = uyao(al+9,i,w)  + (2d0*u + ux*dy2)*dy*dx2
               uyao(al+10,i,w) = uyao(al+10,i,w) + tmp*dx*dz
               uyao(al+11,i,w) = uyao(al+11,i,w) + (2d0*u + ux*dy2)*dy*dz2
               uzao(al+9,i,w)  = uzao(al+9,i,w)  + tmp*dx*dy
               uzao(al+10,i,w) = uzao(al+10,i,w) + (2d0*u + ux*dz2)*dz*dx2
               uzao(al+11,i,w) = uzao(al+11,i,w) + (2d0*u + ux*dz2)*dz*dy2

               tmp = (11d0 - 2d0*alp*r2)*ux
               u2ao(al+9,i,w)  = u2ao(al+9,i,w)  + 2d0*u*(dx2+dy2) + dx2*dy2*tmp
               u2ao(al+10,i,w) = u2ao(al+10,i,w) + 2d0*u*(dx2+dz2) + dx2*dz2*tmp
               u2ao(al+11,i,w) = u2ao(al+11,i,w) + 2d0*u*(dy2+dz2) + dy2*dz2*tmp


c              // g_xxyz, g_yyxz, g_zzxy
               u = sqr3*u                  ! correction of norm
               ux = sqr3*ux

               uao(al+12,i,w)  = uao(al+12,i,w)  + dx*dxyz*u
               uao(al+13,i,w)  = uao(al+13,i,w)  + dy*dxyz*u
               uao(al+14,i,w)  = uao(al+14,i,w)  + dz*dxyz*u

               uxao(al+12,i,w) = uxao(al+12,i,w) + (2d0*u + ux*dx2)*dxyz
               uxao(al+13,i,w) = uxao(al+13,i,w) + (u + ux*dx2)*dy2*dz
               uxao(al+14,i,w) = uxao(al+14,i,w) + (u + ux*dx2)*dz2*dy
               uyao(al+12,i,w) = uyao(al+12,i,w) + (u + ux*dy2)*dx2*dz
               uyao(al+13,i,w) = uyao(al+13,i,w) + (2d0*u + ux*dy2)*dxyz
               uyao(al+14,i,w) = uyao(al+14,i,w) + (u + ux*dy2)*dz2*dx
               uzao(al+12,i,w) = uzao(al+12,i,w) + (u + ux*dz2)*dx2*dy
               uzao(al+13,i,w) = uzao(al+13,i,w) + (u + ux*dz2)*dy2*dx
               uzao(al+14,i,w) = uzao(al+14,i,w) + (2d0*u + ux*dz2)*dxyz

               tmp = (11d0 - 2d0*alp*r2)*ux
               u2ao(al+12,i,w) = u2ao(al+12,i,w) + (2d0*u + tmp*dx2)*dy*dz
               u2ao(al+13,i,w) = u2ao(al+13,i,w) + (2d0*u + tmp*dy2)*dx*dz
               u2ao(al+14,i,w) = u2ao(al+14,i,w) + (2d0*u + tmp*dz2)*dx*dy
            enddo

         end subroutine internal_GamessOrderGFunctions


         subroutine internal_GaussianOrderGFunctions()

            do ic=1,ngto(bf)                      ! loop over contraction
               alp = cntrctn(1,ic,bf)
               u = cntrctn(2,ic,bf) * exp(-alp*r2)
               dx = xx-atoms(a)%cx
               dx2 = dx*dx
               dy = yy-atoms(a)%cy
               dy2 = dy*dy
               dz = zz-atoms(a)%cz
               dz2 = dz*dz
               dxyz = dx*dy*dz
               ux = -2d0*alp*u

c              // g_xxxx, g_yyyy, g_zzzz
               uao(al+14,i,w)    = uao(al+14,i,w)    + dx2*dx2*u
               uao(al+4,i,w)  = uao(al+4,i,w)  + dy2*dy2*u
               uao(al,i,w)  = uao(al,i,w)  + dz2*dz2*u

               uxao(al+14,i,w)   = uxao(al+14,i,w)   + (4d0*u + ux*dx2)*dx2*dx
               uxao(al+4,i,w) = uxao(al+4,i,w) + dy2*dy2*ux*dx
               uxao(al,i,w) = uxao(al,i,w) + dz2*dz2*ux*dx
               uyao(al+14,i,w)   = uyao(al+14,i,w)   + dx2*dx2*ux*dy
               uyao(al+4,i,w) = uyao(al+4,i,w) + (4d0*u + ux*dy2)*dy2*dy
               uyao(al,i,w) = uyao(al,i,w) + dz2*dz2*ux*dy
               uzao(al+14,i,w)   = uzao(al+14,i,w)   + dx2*dx2*ux*dz
               uzao(al+4,i,w) = uzao(al+4,i,w) + dy2*dy2*ux*dz
               uzao(al,i,w) = uzao(al,i,w) + (4d0*u + ux*dz2)*dz2*dz
               tmp          = (11d0 - 2d0*alp*r2)*ux
               u2ao(al+14,i,w)   = u2ao(al+14,i,w)   + (12d0*u + dx2*tmp)*dx2
               u2ao(al+4,i,w) = u2ao(al+4,i,w) + (12d0*u + dy2*tmp)*dy2
               u2ao(al,i,w) = u2ao(al,i,w) + (12d0*u + dz2*tmp)*dz2

c              // g_xxxy, g_xxxz, g_yyyx, g_yyyz, g_zzzx, g_zzzy
               u = sqr7*u                   ! correction of norm
               ux = sqr7*ux

               uao(al+13,i,w)  = uao(al+13,i,w)  + dx2*dx*dy*u
               uao(al+12,i,w)  = uao(al+12,i,w)  + dx2*dx*dz*u
               uao(al+8,i,w)  = uao(al+8,i,w)  + dy2*dy*dx*u
               uao(al+3,i,w)  = uao(al+3,i,w)  + dy2*dy*dz*u
               uao(al+5,i,w)  = uao(al+5,i,w)  + dz2*dz*dx*u
               uao(al+1,i,w)  = uao(al+1,i,w)  + dz2*dz*dy*u

               tmp = ux*dxyz
               uxao(al+13,i,w) = uxao(al+13,i,w) + (3d0*u + ux*dx2)*dx2*dy
               uxao(al+12,i,w) = uxao(al+12,i,w) + (3d0*u + ux*dx2)*dx2*dz
               uxao(al+8,i,w) = uxao(al+8,i,w) + (u + ux*dx2)*dy2*dy
               uxao(al+3,i,w) = uxao(al+3,i,w) + tmp*dy2
               uxao(al+5,i,w) = uxao(al+5,i,w) + (u + ux*dx2)*dz2*dz
               uxao(al+1,i,w) = uxao(al+1,i,w) + tmp*dz2
               uyao(al+13,i,w) = uyao(al+13,i,w) + (u + ux*dy2)*dx2*dx
               uyao(al+12,i,w) = uyao(al+12,i,w) + tmp*dx2
               uyao(al+8,i,w) = uyao(al+8,i,w) + (3d0*u + ux*dy2)*dy2*dx
               uyao(al+3,i,w) = uyao(al+3,i,w) + (3d0*u + ux*dy2)*dy2*dz
               uyao(al+5,i,w) = uyao(al+5,i,w) + tmp*dz2
               uyao(al+1,i,w) = uyao(al+1,i,w) + (u + ux*dy2)*dz2*dz
               uzao(al+13,i,w) = uzao(al+13,i,w) + tmp*dx2
               uzao(al+12,i,w) = uzao(al+12,i,w) + (u + ux*dz2)*dx2*dx
               uzao(al+8,i,w) = uzao(al+8,i,w) + tmp*dy2
               uzao(al+3,i,w) = uzao(al+3,i,w) + (u + ux*dz2)*dy2*dy
               uzao(al+5,i,w) = uzao(al+5,i,w) + (3d0*u + ux*dz2)*dz2*dx
               uzao(al+1,i,w) = uzao(al+1,i,w) + (3d0*u + ux*dz2)*dz2*dy

               tmp = (11d0 - 2d0*alp*r2)*ux
               u2ao(al+13,i,w) = u2ao(al+13,i,w) + (6d0*u + dx2*tmp)*dx*dy
               u2ao(al+12,i,w) = u2ao(al+12,i,w) + (6d0*u + dx2*tmp)*dx*dz
               u2ao(al+8,i,w) = u2ao(al+8,i,w) + (6d0*u + dy2*tmp)*dy*dx
               u2ao(al+3,i,w) = u2ao(al+3,i,w) + (6d0*u + dy2*tmp)*dy*dz
               u2ao(al+5,i,w) = u2ao(al+5,i,w) + (6d0*u + dz2*tmp)*dz*dx
               u2ao(al+1,i,w) = u2ao(al+1,i,w) + (6d0*u + dz2*tmp)*dz*dy

c              // g_xxyy, g_xxzz, g_yyzz
               u = sqr5 / sqr3 * u          ! correction of norm
               ux = sqr5 / sqr3 * ux

               uao(al+11,i,w)   = uao(al+11,i,w)  + dx2*dy2*u
               uao(al+9,i,w)  = uao(al+9,i,w)  + dx2*dz2*u
               uao(al+2,i,w)  = uao(al+2,i,w)  + dy2*dz2*u

               tmp = ux*dxyz
               uxao(al+11,i,w)  = uxao(al+11,i,w)  + (2d0*u + ux*dx2)*dx*dy2
               uxao(al+9,i,w) = uxao(al+9,i,w) + (2d0*u + ux*dx2)*dx*dz2
               uxao(al+2,i,w) = uxao(al+2,i,w) + tmp*dy*dz
               uyao(al+11,i,w)  = uyao(al+11,i,w)  + (2d0*u + ux*dy2)*dy*dx2
               uyao(al+9,i,w) = uyao(al+9,i,w) + tmp*dx*dz
               uyao(al+2,i,w) = uyao(al+2,i,w) + (2d0*u + ux*dy2)*dy*dz2
               uzao(al+11,i,w)  = uzao(al+11,i,w)  + tmp*dx*dy
               uzao(al+9,i,w) = uzao(al+9,i,w) + (2d0*u + ux*dz2)*dz*dx2
               uzao(al+2,i,w) = uzao(al+2,i,w) + (2d0*u + ux*dz2)*dz*dy2

               tmp = (11d0 - 2d0*alp*r2)*ux
               u2ao(al+11,i,w)  = u2ao(al+11,i,w)  + 2d0*u*(dx2+dy2) + dx2*dy2*tmp
               u2ao(al+9,i,w) = u2ao(al+9,i,w) + 2d0*u*(dx2+dz2) + dx2*dz2*tmp
               u2ao(al+2,i,w) = u2ao(al+2,i,w) + 2d0*u*(dy2+dz2) + dy2*dz2*tmp


c              // g_xxyz, g_yyxz, g_zzxy
               u = sqr3*u                  ! correction of norm
               ux = sqr3*ux

               uao(al+10,i,w)  = uao(al+10,i,w)  + dx*dxyz*u
               uao(al+7,i,w)  = uao(al+7,i,w)  + dy*dxyz*u
               uao(al+6,i,w)  = uao(al+6,i,w)  + dz*dxyz*u

               uxao(al+10,i,w) = uxao(al+10,i,w) + (2d0*u + ux*dx2)*dxyz
               uxao(al+7,i,w) = uxao(al+7,i,w) + (u + ux*dx2)*dy2*dz
               uxao(al+6,i,w) = uxao(al+6,i,w) + (u + ux*dx2)*dz2*dy
               uyao(al+10,i,w) = uyao(al+10,i,w) + (u + ux*dy2)*dx2*dz
               uyao(al+7,i,w) = uyao(al+7,i,w) + (2d0*u + ux*dy2)*dxyz
               uyao(al+6,i,w) = uyao(al+6,i,w) + (u + ux*dy2)*dz2*dx
               uzao(al+10,i,w) = uzao(al+10,i,w) + (u + ux*dz2)*dx2*dy
               uzao(al+7,i,w) = uzao(al+7,i,w) + (u + ux*dz2)*dy2*dx
               uzao(al+6,i,w) = uzao(al+6,i,w) + (2d0*u + ux*dz2)*dxyz

               tmp = (11d0 - 2d0*alp*r2)*ux
               u2ao(al+10,i,w) = u2ao(al+10,i,w) + (2d0*u + tmp*dx2)*dy*dz
               u2ao(al+7,i,w) = u2ao(al+7,i,w) + (2d0*u + tmp*dy2)*dx*dz
               u2ao(al+6,i,w) = u2ao(al+6,i,w) + (2d0*u + tmp*dz2)*dx*dy
            enddo
         end subroutine internal_GaussianOrderGFunctions


      end subroutine aocalc

c ==================================================
cBH

c     --------------------------------
      subroutine aosplcalc(ie,ec,rrai)
c     --------------------------------

c aosplcalc calculates all atomic orbitals for position vector x,y,z
c including all required derivatives, using splines.

c on entry: requires rai distance matrix properly dimensioned and set
c           for input position vector x,y,z.
c           for ie>0 position vector must be modified only compared to
c           last call to aocalc only at electron ie
c on exit:  updates AO data structure in aos.h (i.e. AOs and derivatives)

c Version 1.1 (28.4.98):  reference to walker eliminated.
c Version 1.0 (22.11.97): AO calculation using splined radial parts
c                         modified Version 2.0 of aocalc (1.1 of getaos)
c                         splines are used only for contracted GTO's
      integer, intent(in)             :: ie                 ! if >0 only AO's for electron ie recalculated
      type(eConfigArray), intent(inout)  :: ec                 ! electron configurations
      real*8, intent(in)              :: rrai(:,:,:)     ! r_ai electron-nucleus distances
c constants:
      real*8 sqr3,sqr5,sqr7
      parameter (sqr3=1.73205080756887729d0,sqr5=2.236067977499789696d0,
     &           sqr7=2.645751311064591d0)
c variables
      integer bf,a,i,j,i1,i2,nn,al,ispl,m,n,k,l,w,aa,ii
      real*8 rai(ncenter,ne)
      real*8, pointer :: x(:),y(:),z(:)
      real*8 xx,yy,zz,rr,r2,alp,nrm,u,ux,uxx,u2,dx,dy,dz,tmp,
     .       dx2,dy2,dz2,dxyz,df
      logical gaussFOrder       ! .t.: Gaussian order for f function
                                ! .f.: Gamess==Turbomole order used

      integer :: cnums, cnump, cnumd, cnumf


c bf refers to the degenerate set of cartesian
c basis function (S:1,P:3,D:6,F:10) as input, which may be of typ STO
c or contracted GTO.
c al refers to the individual basis function, as used in LCAO-MO's.
c (composed in subroutine mdetwf)
c i refers to the current electron.
c-----Calculation of the AO's and their derivatives

      mAOElecConfigs = eConfigArray_size(ec)
      do w=1,mAOElecConfigs

      call eConfigArray_getPtr(ec,w,x,y,z)
      rai = rrai(:,:,w)

      ! check input data for NaN or Inf
      call assert(all(rai<huge(1.d0)),"aosplcalc: illegal rai values")
      call assert(all(abs(x(1:ne))<huge(1.d0)) .and.
     .  all(abs(y(1:ne))<huge(1.d0)) .and. all(abs(z(1:ne))<huge(1.d0)),
     .  "aosplcalc: illegal x,y,z coords in ec")

      if (evfmt=='gau' .or. evfmt=='mol' ) then
         gaussFOrder = .true.
      else
         gaussFOrder = .false.
      endif

      if (ie .eq. 0) then                     ! AO's for all electrons
         i1 = 1
         i2 = ne
      else
         i1 = ie                              ! only AO for electron ie
         i2 = ie
      endif

      do i=i1,i2                              ! loop over electrons

         xx = x(i)
         yy = y(i)
         zz = z(i)

         al = 1

         cnums = 0
         cnump = 0
         cnumd = 0
         cnumf = 0


         ! loop over basis functions:
         do n = 1,nbasf
           bf = n
           a = bc(bf)                         ! center of AO
           rr = rai(a,i)                      ! r_ai
           nn = bn(bf)                        ! n quantum no. of AO

           if (typ(bf) == 'STO') then      ! STO-type basis function

            alp = bzet(bf)                    ! orbital exponent
            nrm = norm(bf)                    ! normalization constant

            if (nn == 1) then                    ! 1s orbital
               u  = nrm*exp(-alp*rr)
               ux = -alp*u/rr
               uao(al,i,w) = u
               uxao(al,i,w) = ux*(xx-atoms(a)%cx)
               uyao(al,i,w) = ux*(yy-atoms(a)%cy)
               uzao(al,i,w) = ux*(zz-atoms(a)%cz)
               u2ao(al,i,w) = alp*(alp-2d0/rr)*u
               al = al+1
            else if (nn == 2) then
               if (bl(bf) == 'S') then           ! 2s orbital
c                 // 2s orbital
                  u  = nrm*exp(-alp*rr)
                  ux = (1d0-alp*rr)*u/rr
                  uao(al,i,w) = rr*u
                  uxao(al,i,w) = ux*(xx-atoms(a)%cx)
                  uyao(al,i,w) = ux*(yy-atoms(a)%cy)
                  uzao(al,i,w) = ux*(zz-atoms(a)%cz)
                  u2ao(al,i,w) = alp*(alp*rr-2d0)*u + 2d0*ux
                  al = al+1
               else if (bl(bf) == 'P') then      ! 2p orbital
                  u  = nrm*exp(-alp*rr)
                  dx = xx-atoms(a)%cx
                  dy = yy-atoms(a)%cy
                  dz = zz-atoms(a)%cz
                  ux = -alp*u/rr
                  uao(al,i,w)    = u*dx
                  uao(al+1,i,w)  = u*dy
                  uao(al+2,i,w)  = u*dz
                  uxao(al,i,w)   = ux*dx*dx + u
                  uxao(al+1,i,w) = ux*dx*dy
                  uxao(al+2,i,w) = ux*dx*dz
                  uyao(al,i,w)   = uxao(al+1,i,w)
                  uyao(al+1,i,w) = ux*dy*dy + u
                  uyao(al+2,i,w) = ux*dy*dz
                  uzao(al,i,w)   = uxao(al+2,i,w)
                  uzao(al+1,i,w) = uyao(al+2,i,w)
                  uzao(al+2,i,w) = ux*dz*dz + u
                  tmp          = alp*(alp-4d0/rr)*u
                  u2ao(al,i,w)   = tmp*dx
                  u2ao(al+1,i,w) = tmp*dy
                  u2ao(al+2,i,w) = tmp*dz
                  al = al+3
               else
                  call abortp(' aossplcalc: n=2 and l > 2')
               endif

            else if (nn == 3) then
               if (bl(bf) == 'S') then
c                 // 3s orbital
                  u  = nrm*exp(-alp*rr)
                  ux = (2d0-alp*rr)*u
                  uao(al,i,w) = rr*rr*u
                  uxao(al,i,w) = ux*(xx-atoms(a)%cx)
                  uyao(al,i,w) = ux*(yy-atoms(a)%cy)
                  uzao(al,i,w) = ux*(zz-atoms(a)%cz)
                  tmp = (2d0 + alp*rr*(-4d0 + alp*rr))*u
                  u2ao(al,i,w) = tmp + 2d0*ux
                  al = al+1
               else if (bl(bf) == 'P') then
                  call abortp(' aosplcalc: 3p orbitals not implemented')
               else if (bl(bf) == 'D') then         ! 3d orbital (6D)
c                 // Norm is different for d_xx and d_xy !

                  u   = nrm*exp(-alp*rr)
                  dx  = xx - atoms(a)%cx
                  dx2 = dx*dx
                  dy  = yy - atoms(a)%cy
                  dy2 = dy*dy
                  dz  = zz - atoms(a)%cz
                  dz2 = dz*dz
                  ux  = -alp*u/rr
                  tmp = alp*(alp-6d0/rr)

                  uao(al,i,w)    = dx2*u
                  uao(al+1,i,w)  = dy2*u
                  uao(al+2,i,w)  = dz2*u
                  uxao(al,i,w)   = dx*(2d0*u + ux*dx2)
                  uxao(al+1,i,w) = ux*dy2*dx
                  uxao(al+2,i,w) = ux*dz2*dx
                  uyao(al,i,w)   = ux*dx2*dy
                  uyao(al+1,i,w) = dy*(2d0*u + ux*dy2)
                  uyao(al+2,i,w) = ux*dz2*dy
                  uzao(al,i,w)   = ux*dx2*dz
                  uzao(al+1,i,w) = ux*dy2*dz
                  uzao(al+2,i,w) = dz*(2d0*u + ux*dz2)
                  u2ao(al,i,w)   = u*(2d0+tmp*dx2)
                  u2ao(al+1,i,w) = u*(2d0+tmp*dy2)
                  u2ao(al+2,i,w) = u*(2d0+tmp*dz2)

                  u = sqr3*u                   ! correction of norm for last 3
                  ux = sqr3*ux

                  uao(al+3,i,w)  = u*dx*dy
                  uao(al+4,i,w)  = u*dx*dz
                  uao(al+5,i,w)  = u*dy*dz

                  uxao(al+3,i,w) = dy*(u + ux*dx2)
                  uxao(al+4,i,w) = dz*(u + ux*dx2)
                  uxao(al+5,i,w) = ux*dx*dy*dz
                  uyao(al+3,i,w) = dx*(u + ux*dy2)
                  uyao(al+4,i,w) = uxao(al+5,i,w)
                  uyao(al+5,i,w) = dz*(u + ux*dy2)
                  uzao(al+3,i,w) = uxao(al+5,i,w)
                  uzao(al+4,i,w) = dx*(u + ux*dz2)
                  uzao(al+5,i,w) = dy*(u + ux*dz2)

                  u2ao(al+3,i,w) = u*dx*dy*tmp
                  u2ao(al+4,i,w) = u*dx*dz*tmp
                  u2ao(al+5,i,w) = u*dy*dz*tmp
                  al = al+6
               endif     ! bl
            endif        ! nn


c          // Contracted GTO's as basis function (AO)
c          // Evaluate with splines
           else
           if (so(bf)==0) then !only 1 GTO in contraction, no splines used !

            r2 = rr*rr

c           // only primitive cartesian gaussians: 1s,2p,3d,4f
c           // i.e. no r factor. Thus nn is not used here.

            if (bl(bf) == 'S') then                 ! 1s GTO
               alp = cntrctn(1,1,bf)
               u   = cntrctn(2,1,bf) * exp(-alp*r2)
               ux = -2d0*alp*u
               uao(al,i,w)  = u
               uxao(al,i,w) = ux*(xx-atoms(a)%cx)
               uyao(al,i,w) = ux*(yy-atoms(a)%cy)
               uzao(al,i,w) = ux*(zz-atoms(a)%cz)
               u2ao(al,i,w) = ux*(3d0-2d0*alp*r2)
               al = al+1

            else if (bl(bf) == 'P') then             ! 2p GTO's
c              // do all 3 P simultaneously (same exponent is required)
c              // order p_x,p_y,p_z
               alp = cntrctn(1,1,bf)
               u = cntrctn(2,1,bf) * exp(-alp*r2)
               dx = xx-atoms(a)%cx
               dy = yy-atoms(a)%cy
               dz = zz-atoms(a)%cz
               ux = -2d0*alp*u
               uao(al,i,w)    = dx*u
               uao(al+1,i,w)  = dy*u
               uao(al+2,i,w)  = dz*u
               uxao(al,i,w)   = u + ux*dx*dx
               uxao(al+1,i,w) = ux*dx*dy
               uxao(al+2,i,w) = ux*dx*dz
               uyao(al,i,w)   = ux*dx*dy
               uyao(al+1,i,w) = u + ux*dy*dy
               uyao(al+2,i,w) = ux*dy*dz
               uzao(al,i,w)   = ux*dx*dz
               uzao(al+1,i,w) = ux*dy*dz
               uzao(al+2,i,w) = u + ux*dz*dz
               tmp = (5d0-2d0*alp*r2)*ux
               u2ao(al,i,w)   = tmp*dx
               u2ao(al+1,i,w) = tmp*dy
               u2ao(al+2,i,w) = tmp*dz
               al = al+3

            else if (bl(bf) == 'D') then         ! 3d GTO
c              // do all 6 D simultaneously (same exponent is required)
c              // order: d_xx, d_yy, d_zz, d_xy, d_xz, d_yz  (like GAMESS)
               alp = cntrctn(1,1,bf)
               u = cntrctn(2,1,bf) * exp(-alp*r2)
               dx = xx-atoms(a)%cx
               dx2 = dx*dx
               dy = yy-atoms(a)%cy
               dy2 = dy*dy
               dz = zz-atoms(a)%cz
               dz2 = dz*dz
               ux = -2d0*alp*u
               uao(al,i,w)    = dx2*u
               uao(al+1,i,w)  = dy2*u
               uao(al+2,i,w)  = dz2*u
               uxao(al,i,w)   = (2d0*u + ux*dx2)*dx
               uxao(al+1,i,w) = dy2*ux*dx
               uxao(al+2,i,w) = dz2*ux*dx
               uyao(al,i,w)   = dx2*ux*dy
               uyao(al+1,i,w) = (2d0*u + ux*dy2)*dy
               uyao(al+2,i,w) = dz2*ux*dy
               uzao(al,i,w)   = dx2*ux*dz
               uzao(al+1,i,w) = dy2*ux*dz
               uzao(al+2,i,w) = (2d0*u + ux*dz2)*dz
               tmp          = (7d0 - 2d0*alp*r2)*ux
               u2ao(al,i,w)   = 2d0*u + dx2*tmp
               u2ao(al+1,i,w) = 2d0*u + dy2*tmp
               u2ao(al+2,i,w) = 2d0*u + dz2*tmp
               u = sqr3*u                   ! correction of norm for last 3
               ux = sqr3*ux
               uao(al+3,i,w)  = dx*dy*u
               uao(al+4,i,w)  = dx*dz*u
               uao(al+5,i,w)  = dy*dz*u
               tmp = ux*dx*dy*dz
               uxao(al+3,i,w) = (u + ux*dx2)*dy
               uxao(al+4,i,w) = (u + ux*dx2)*dz
               uxao(al+5,i,w) = tmp
               uyao(al+3,i,w) = (u + ux*dy2)*dx
               uyao(al+4,i,w) = tmp
               uyao(al+5,i,w) = (u + ux*dy2)*dz
               uzao(al+3,i,w) = tmp
               uzao(al+4,i,w) = (u + ux*dz2)*dx
               uzao(al+5,i,w) = (u + ux*dz2)*dy
               tmp = (7d0 - 2d0*alp*r2)*ux
               u2ao(al+3,i,w) = tmp*dx*dy
               u2ao(al+4,i,w) = tmp*dx*dz
               u2ao(al+5,i,w) = tmp*dy*dz
               al = al+6

            else if (bl(bf)=='F'.and..not.gaussFOrder) then     ! 3f GTO
c              // do all 10 F simultaneously (same exponent is required)
c              // order: f_xxx, f_yyy, f_zzz, fd_xxy, f_xxz, f_yyx,
c              //   f_yyz, f_zzx, f_zzy, f_xyz  (like GAMESS)
               alp = cntrctn(1,1,bf)
               u = cntrctn(2,1,bf) * exp(-alp*r2)
               dx = xx-atoms(a)%cx
               dx2 = dx*dx
               dy = yy-atoms(a)%cy
               dy2 = dy*dy
               dz = zz-atoms(a)%cz
               dz2 = dz*dz
               dxyz = dx*dy*dz
               ux = -2d0*alp*u
c              // f_xxx, f_yyy, f_zzz
               uao(al,i,w)    = dx2*dx*u
               uao(al+1,i,w)  = dy2*dy*u
               uao(al+2,i,w)  = dz2*dz*u
               uxao(al,i,w)   = (3d0*u + ux*dx2)*dx2
               uxao(al+1,i,w) = dy2*dy*ux*dx
               uxao(al+2,i,w) = dz2*dz*ux*dx
               uyao(al,i,w)   = dx2*dx*ux*dy
               uyao(al+1,i,w) = (3d0*u + ux*dy2)*dy2
               uyao(al+2,i,w) = dz2*dz*ux*dy
               uzao(al,i,w)   = dx2*dx*ux*dz
               uzao(al+1,i,w) = dy2*dy*ux*dz
               uzao(al+2,i,w) = (3d0*u + ux*dz2)*dz2
               tmp          = (9d0 - 2d0*alp*r2)*ux
               u2ao(al,i,w)   = (6d0*u + dx2*tmp)*dx
               u2ao(al+1,i,w) = (6d0*u + dy2*tmp)*dy
               u2ao(al+2,i,w) = (6d0*u + dz2*tmp)*dz
c              // f_xxy, f_xxz, f_yyx, f_yyz, f_zzx, f_zzy
               u = sqr5*u                   ! correction of norm
               ux = sqr5*ux
               uao(al+3,i,w)  = dx2*dy*u
               uao(al+4,i,w)  = dx2*dz*u
               uao(al+5,i,w)  = dy2*dx*u
               uao(al+6,i,w)  = dy2*dz*u
               uao(al+7,i,w)  = dz2*dx*u
               uao(al+8,i,w)  = dz2*dy*u

               tmp = ux*dxyz
               uxao(al+3,i,w) = (2d0*u + ux*dx2)*dx*dy
               uxao(al+4,i,w) = (2d0*u + ux*dx2)*dx*dz
               uxao(al+5,i,w) = (u + ux*dx2)*dy2
               uxao(al+6,i,w) = tmp*dy
               uxao(al+7,i,w) = (u + ux*dx2)*dz2
               uxao(al+8,i,w) = tmp*dz
               uyao(al+3,i,w) = (u + ux*dy2)*dx2
               uyao(al+4,i,w) = tmp*dx
               uyao(al+5,i,w) = (2d0*u + ux*dy2)*dx*dy
               uyao(al+6,i,w) = (2d0*u + ux*dy2)*dy*dz
               uyao(al+7,i,w) = tmp*dz
               uyao(al+8,i,w) = (u + ux*dy2)*dz2
               uzao(al+3,i,w) = tmp*dx
               uzao(al+4,i,w) = (u + ux*dz2)*dx2
               uzao(al+5,i,w) = tmp*dy
               uzao(al+6,i,w) = (u + ux*dz2)*dy2
               uzao(al+7,i,w) = (2d0*u + ux*dz2)*dx*dz
               uzao(al+8,i,w) = (2d0*u + ux*dz2)*dy*dz

               tmp = (9d0 - 2d0*alp*r2)*ux
               u2ao(al+3,i,w) = (2d0*u + dx2*tmp)*dy
               u2ao(al+4,i,w) = (2d0*u + dx2*tmp)*dz
               u2ao(al+5,i,w) = (2d0*u + dy2*tmp)*dx
               u2ao(al+6,i,w) = (2d0*u + dy2*tmp)*dz
               u2ao(al+7,i,w) = (2d0*u + dz2*tmp)*dx
               u2ao(al+8,i,w) = (2d0*u + dz2*tmp)*dy
c              // f_xyz
               u = sqr3*u                  ! correction of norm
               ux = sqr3*ux
               uao(al+9,i,w)  = dxyz*u
               uxao(al+9,i,w) = (u + ux*dx2)*dy*dz
               uyao(al+9,i,w) = (u + ux*dy2)*dx*dz
               uzao(al+9,i,w) = (u + ux*dz2)*dx*dy
               tmp = (9d0 - 2d0*alp*r2)*ux
               u2ao(al+9,i,w) = dxyz*tmp
               al = al + 10

            else if (bl(bf)=='F'.and.gaussFOrder) then     ! 3f GTO
c              // do all 10 F simultaneously (same exponent is required)
c              // order: f_xxx, f_yyy, f_zzz, f_xyy, f_xxy, f_xxz,
c              //   f_xzz, f_yzz, f_yyz, f_xyz  (like Gaussian)
               alp = cntrctn(1,1,bf)
               u = cntrctn(2,1,bf) * exp(-alp*r2)
               dx = xx-atoms(a)%cx
               dx2 = dx*dx
               dy = yy-atoms(a)%cy
               dy2 = dy*dy
               dz = zz-atoms(a)%cz
               dz2 = dz*dz
               dxyz = dx*dy*dz
               ux = -2d0*alp*u
c              // f_xxx, f_yyy, f_zzz
               uao(al,i,w)    = dx2*dx*u
               uao(al+1,i,w)  = dy2*dy*u
               uao(al+2,i,w)  = dz2*dz*u
               uxao(al,i,w)   = (3d0*u + ux*dx2)*dx2
               uxao(al+1,i,w) = dy2*dy*ux*dx
               uxao(al+2,i,w) = dz2*dz*ux*dx
               uyao(al,i,w)   = dx2*dx*ux*dy
               uyao(al+1,i,w) = (3d0*u + ux*dy2)*dy2
               uyao(al+2,i,w) = dz2*dz*ux*dy
               uzao(al,i,w)   = dx2*dx*ux*dz
               uzao(al+1,i,w) = dy2*dy*ux*dz
               uzao(al+2,i,w) = (3d0*u + ux*dz2)*dz2
               tmp          = (9d0 - 2d0*alp*r2)*ux
               u2ao(al,i,w)   = (6d0*u + dx2*tmp)*dx
               u2ao(al+1,i,w) = (6d0*u + dy2*tmp)*dy
               u2ao(al+2,i,w) = (6d0*u + dz2*tmp)*dz
c              // f_xxy, f_xxz, f_yyx, f_yyz, f_zzx, f_zzy
               u = sqr5*u                   ! correction of norm
               ux = sqr5*ux

               uao(al+4,i,w)  = dx2*dy*u
               uao(al+5,i,w)  = dx2*dz*u
               uao(al+3,i,w)  = dy2*dx*u
               uao(al+8,i,w)  = dy2*dz*u
               uao(al+6,i,w)  = dz2*dx*u
               uao(al+7,i,w)  = dz2*dy*u

               tmp = ux*dxyz
               uxao(al+4,i,w) = (2d0*u + ux*dx2)*dx*dy
               uxao(al+5,i,w) = (2d0*u + ux*dx2)*dx*dz
               uxao(al+3,i,w) = (u + ux*dx2)*dy2
               uxao(al+8,i,w) = tmp*dy
               uxao(al+6,i,w) = (u + ux*dx2)*dz2
               uxao(al+7,i,w) = tmp*dz
               uyao(al+4,i,w) = (u + ux*dy2)*dx2
               uyao(al+5,i,w) = tmp*dx
               uyao(al+3,i,w) = (2d0*u + ux*dy2)*dx*dy
               uyao(al+8,i,w) = (2d0*u + ux*dy2)*dy*dz
               uyao(al+6,i,w) = tmp*dz
               uyao(al+7,i,w) = (u + ux*dy2)*dz2
               uzao(al+4,i,w) = tmp*dx
               uzao(al+5,i,w) = (u + ux*dz2)*dx2
               uzao(al+3,i,w) = tmp*dy
               uzao(al+8,i,w) = (u + ux*dz2)*dy2
               uzao(al+6,i,w) = (2d0*u + ux*dz2)*dx*dz
               uzao(al+7,i,w) = (2d0*u + ux*dz2)*dy*dz

               tmp = (9d0 - 2d0*alp*r2)*ux
               u2ao(al+4,i,w) = (2d0*u + dx2*tmp)*dy
               u2ao(al+5,i,w) = (2d0*u + dx2*tmp)*dz
               u2ao(al+3,i,w) = (2d0*u + dy2*tmp)*dx
               u2ao(al+8,i,w) = (2d0*u + dy2*tmp)*dz
               u2ao(al+6,i,w) = (2d0*u + dz2*tmp)*dx
               u2ao(al+7,i,w) = (2d0*u + dz2*tmp)*dy

c              // f_xyz
               u = sqr3*u                  ! correction of norm
               ux = sqr3*ux
               uao(al+9,i,w)  = dxyz*u
               uxao(al+9,i,w) = (u + ux*dx2)*dy*dz
               uyao(al+9,i,w) = (u + ux*dy2)*dx*dz
               uzao(al+9,i,w) = (u + ux*dz2)*dx*dy
               tmp = (9d0 - 2d0*alp*r2)*ux
               u2ao(al+9,i,w) = dxyz*tmp
               al = al + 10

            else if (bl(bf)=='G') then     ! 5g GTO
c              // do all 15 cartesian G simultaneously (same exponent is required)

               uao(al:al+14,i,w) = 0d0
               uxao(al:al+14,i,w) = 0d0
               uyao(al:al+14,i,w) = 0d0
               uzao(al:al+14,i,w) = 0d0
               u2ao(al:al+14,i,w) = 0d0

               if (gaussFOrder) then
                  call internal_GaussianOrderGFunctions()
               else
                  call internal_GamessOrderGFunctions()
               endif

               al = al + 15

            else
               call abortp('(aossplcalc): wrong GTO')
            endif  ! bl


           else	 !more then 1 GTO in contraction, splines used !
            r2 = rr*rr
            j  = (csplnpnt-1)*rr/(csalpha+rr)  + 1
            df = rr - csplx(j)

c           // only primitive cartesian gaussians: 1s,2p,3d,4f
c           // i.e. no r factor. Thus nn is not used here.

            if (bl(bf) == 'S') then                 ! 1s GTO

               ispl       = 3*so(bf)-2
               uao(al,i,w)  = cspla(ispl,j) + df*(csplb(ispl,j)
     .                    + df*(csplc(ispl,j) + df*cspld(ispl,j)))
               ispl       = ispl + 1
               ux         = cspla(ispl,j) + df*(csplb(ispl,j)
     .                    + df*(csplc(ispl,j) + df*cspld(ispl,j)))
               ispl       = ispl + 1
               u2         = cspla(ispl,j) + df*(csplb(ispl,j)
     .                    + df*(csplc(ispl,j) + df*cspld(ispl,j)))
               uxao(al,i,w) = ux*(xx-atoms(a)%cx)/rr
               uyao(al,i,w) = ux*(yy-atoms(a)%cy)/rr
               uzao(al,i,w) = ux*(zz-atoms(a)%cz)/rr
               u2ao(al,i,w) = u2 + 2*ux/rr

               al = al+1


            else if (bl(bf) == 'P') then             ! 2p GTO's

c              // do all 3 P simultaneously (same exponent is required)
c              // order p_x,p_y,p_z
               ispl = 3*so(bf)-2
               u    = cspla(ispl,j) + df*(csplb(ispl,j)
     .              + df*(csplc(ispl,j) + df*cspld(ispl,j)))
               ispl = ispl + 1
               ux   = cspla(ispl,j) + df*(csplb(ispl,j)
     .              + df*(csplc(ispl,j) + df*cspld(ispl,j)))
               ispl = ispl + 1
               uxx  = cspla(ispl,j) + df*(csplb(ispl,j)
     .              + df*(csplc(ispl,j) + df*cspld(ispl,j)))

               dx = xx-atoms(a)%cx
               dy = yy-atoms(a)%cy
               dz = zz-atoms(a)%cz
               uao(al,i,w) = dx*u
               uao(al+1,i,w) = dy*u
               uao(al+2,i,w) = dz*u
               uxao(al,i,w)   = u + ux*dx*dx
               uxao(al+1,i,w) = ux*dx*dy
               uxao(al+2,i,w) = ux*dx*dz
               uyao(al,i,w)   = ux*dx*dy
               uyao(al+1,i,w) = u + ux*dy*dy
               uyao(al+2,i,w) = ux*dy*dz
               uzao(al,i,w)   = ux*dx*dz
               uzao(al+1,i,w) = ux*dy*dz
               uzao(al+2,i,w) = u + ux*dz*dz
               u2ao(al,i,w)   = uxx*dx
               u2ao(al+1,i,w) = uxx*dy
               u2ao(al+2,i,w) = uxx*dz
               al = al+3

            else if (bl(bf) == 'D') then         ! 3d GTO

c              // do all 6 D simultaneously (same exponent is required)
c              // order: d_xx, d_yy, d_zz, d_xy, d_xz, d_yz  (like GAMESS)
               ispl = 3*so(bf)-2
               u    = cspla(ispl,j) + df*(csplb(ispl,j)
     .              + df*(csplc(ispl,j) + df*cspld(ispl,j)))
               ispl = ispl + 1
               ux   = cspla(ispl,j) + df*(csplb(ispl,j)
     .              + df*(csplc(ispl,j) + df*cspld(ispl,j)))
               ispl = ispl + 1
               uxx   = cspla(ispl,j) + df*(csplb(ispl,j)
     .              + df*(csplc(ispl,j) + df*cspld(ispl,j)))

               dx = xx-atoms(a)%cx
               dx2 = dx*dx
               dy = yy-atoms(a)%cy
               dy2 = dy*dy
               dz = zz-atoms(a)%cz
               dz2 = dz*dz
               uao(al,i,w)    = dx2*u
               uao(al+1,i,w)  = dy2*u
               uao(al+2,i,w)  = dz2*u
               uxao(al,i,w)   = (2d0*u + ux*dx2)*dx
               uxao(al+1,i,w) = dy2*ux*dx
               uxao(al+2,i,w) = dz2*ux*dx
               uyao(al,i,w)   = dx2*ux*dy
               uyao(al+1,i,w) = (2d0*u + ux*dy2)*dy
               uyao(al+2,i,w) = dz2*ux*dy
               uzao(al,i,w)   = dx2*ux*dz
               uzao(al+1,i,w) = dy2*ux*dz
               uzao(al+2,i,w) = (2d0*u + ux*dz2)*dz
               u2ao(al,i,w)   = 2d0*u + dx2*uxx
               u2ao(al+1,i,w) = 2d0*u + dy2*uxx
               u2ao(al+2,i,w) = 2d0*u + dz2*uxx

               u = sqr3*u                   ! correction of norm for last 3
               ux = sqr3*ux
               uxx = sqr3*uxx

               uao(al+3,i,w)  = dx*dy*u
               uao(al+4,i,w)  = dx*dz*u
               uao(al+5,i,w)  = dy*dz*u

               tmp = ux*dx*dy*dz
               uxao(al+3,i,w) = (u + ux*dx2)*dy
               uxao(al+4,i,w) = (u + ux*dx2)*dz
               uxao(al+5,i,w) = tmp
               uyao(al+3,i,w) = (u + ux*dy2)*dx
               uyao(al+4,i,w) = tmp
               uyao(al+5,i,w) = (u + ux*dy2)*dz
               uzao(al+3,i,w) = tmp
               uzao(al+4,i,w) = (u + ux*dz2)*dx
               uzao(al+5,i,w) = (u + ux*dz2)*dy

               u2ao(al+3,i,w) = uxx*dx*dy
               u2ao(al+4,i,w) = uxx*dx*dz
               u2ao(al+5,i,w) = uxx*dy*dz
               al = al+6

            else if (bl(bf)=='F'.and..not.gaussFOrder) then     ! 3f GTO
c              // do all 10 F simultaneously (same exponent is required)
c              // order: f_xxx, f_yyy, f_zzz, fd_xxy, f_xxz, f_yyx,
c              //   f_yyz, f_zzx, f_zzy, f_xyz  (like GAMESS)

               ispl = 3*so(bf)-2
               u    = cspla(ispl,j) + df*(csplb(ispl,j)
     .              + df*(csplc(ispl,j) + df*cspld(ispl,j)))
               ispl = ispl + 1
               ux   = cspla(ispl,j) + df*(csplb(ispl,j)
     .              + df*(csplc(ispl,j) + df*cspld(ispl,j)))
               ispl = ispl + 1
               uxx   = cspla(ispl,j) + df*(csplb(ispl,j)
     .              + df*(csplc(ispl,j) + df*cspld(ispl,j)))

               dx = xx-atoms(a)%cx
               dx2 = dx*dx
               dy = yy-atoms(a)%cy
               dy2 = dy*dy
               dz = zz-atoms(a)%cz
               dz2 = dz*dz
               dxyz = dx*dy*dz

c              // f_xxx, f_yyy, f_zzz
               uao(al,i,w)    = dx2*dx*u
               uao(al+1,i,w)  = dy2*dy*u
               uao(al+2,i,w)  = dz2*dz*u

               uxao(al,i,w)   = (3d0*u + ux*dx2)*dx2
               uxao(al+1,i,w) = dy2*dy*ux*dx
               uxao(al+2,i,w) = dz2*dz*ux*dx
               uyao(al,i,w)   = dx2*dx*ux*dy
               uyao(al+1,i,w) = (3d0*u + ux*dy2)*dy2
               uyao(al+2,i,w) = dz2*dz*ux*dy
               uzao(al,i,w)   = dx2*dx*ux*dz
               uzao(al+1,i,w) = dy2*dy*ux*dz
               uzao(al+2,i,w) = (3d0*u + ux*dz2)*dz2
               u2ao(al,i,w)   = (6d0*u + dx2*uxx)*dx
               u2ao(al+1,i,w) = (6d0*u + dy2*uxx)*dy
               u2ao(al+2,i,w) = (6d0*u + dz2*uxx)*dz

c              // f_xxy, f_xxz, f_yyx, f_yyz, f_zzx, f_zzy
               u = sqr5*u                   ! correction of norm
               ux = sqr5*ux
               uxx = sqr5*uxx

               uao(al+3,i,w)  = dx2*dy*u
               uao(al+4,i,w)  = dx2*dz*u
               uao(al+5,i,w)  = dy2*dx*u
               uao(al+6,i,w)  = dy2*dz*u
               uao(al+7,i,w)  = dz2*dx*u
               uao(al+8,i,w)  = dz2*dy*u
               tmp = ux*dxyz
               uxao(al+3,i,w) = (2d0*u + ux*dx2)*dx*dy
               uxao(al+4,i,w) = (2d0*u + ux*dx2)*dx*dz
               uxao(al+5,i,w) = (u + ux*dx2)*dy2
               uxao(al+6,i,w) = tmp*dy
               uxao(al+7,i,w) = (u + ux*dx2)*dz2
               uxao(al+8,i,w) = tmp*dz
               uyao(al+3,i,w) = (u + ux*dy2)*dx2
               uyao(al+4,i,w) = tmp*dx
               uyao(al+5,i,w) = (2d0*u + ux*dy2)*dx*dy
               uyao(al+6,i,w) = (2d0*u + ux*dy2)*dy*dz
               uyao(al+7,i,w) = tmp*dz
               uyao(al+8,i,w) = (u + ux*dy2)*dz2
               uzao(al+3,i,w) = tmp*dx
               uzao(al+4,i,w) = (u + ux*dz2)*dx2
               uzao(al+5,i,w) = tmp*dy
               uzao(al+6,i,w) = (u + ux*dz2)*dy2
               uzao(al+7,i,w) = (2d0*u + ux*dz2)*dx*dz
               uzao(al+8,i,w) = (2d0*u + ux*dz2)*dy*dz

               u2ao(al+3,i,w) = (2d0*u + dx2*uxx)*dy
               u2ao(al+4,i,w) = (2d0*u + dx2*uxx)*dz
               u2ao(al+5,i,w) = (2d0*u + dy2*uxx)*dx
               u2ao(al+6,i,w) = (2d0*u + dy2*uxx)*dz
               u2ao(al+7,i,w) = (2d0*u + dz2*uxx)*dx
               u2ao(al+8,i,w) = (2d0*u + dz2*uxx)*dy

c              // f_xyz
               u = sqr3*u                  ! correction of norm
               ux = sqr3*ux
               uxx = sqr3*uxx

               uao(al+9,i,w)  = dxyz*u

               uxao(al+9,i,w) = (u + ux*dx2)*dy*dz
               uyao(al+9,i,w) = (u + ux*dy2)*dx*dz
               uzao(al+9,i,w) = (u + ux*dz2)*dx*dy
               u2ao(al+9,i,w) = dxyz*uxx
               al = al +10

            else if (bl(bf)=='F'.and.gaussFOrder) then     ! f GTO
c              // do all 10 F simultaneously (same exponent is required)
c              // order: f_xxx, f_yyy, f_zzz, f_xyy, f_xxy, f_xxz,
c              //   f_xzz, f_yzz, f_yyz, f_xyz  (like Gaussian)

               ispl = 3*so(bf)-2
               u    = cspla(ispl,j) + df*(csplb(ispl,j)
     .              + df*(csplc(ispl,j) + df*cspld(ispl,j)))
               ispl = ispl + 1
               ux   = cspla(ispl,j) + df*(csplb(ispl,j)
     .              + df*(csplc(ispl,j) + df*cspld(ispl,j)))
               ispl = ispl + 1
               uxx   = cspla(ispl,j) + df*(csplb(ispl,j)
     .              + df*(csplc(ispl,j) + df*cspld(ispl,j)))

               dx = xx-atoms(a)%cx
               dx2 = dx*dx
               dy = yy-atoms(a)%cy
               dy2 = dy*dy
               dz = zz-atoms(a)%cz
               dz2 = dz*dz
               dxyz = dx*dy*dz

c              // f_xxx, f_yyy, f_zzz
               uao(al,i,w)    = dx2*dx*u
               uao(al+1,i,w)  = dy2*dy*u
               uao(al+2,i,w)  = dz2*dz*u

               uxao(al,i,w)   = (3d0*u + ux*dx2)*dx2
               uxao(al+1,i,w) = dy2*dy*ux*dx
               uxao(al+2,i,w) = dz2*dz*ux*dx
               uyao(al,i,w)   = dx2*dx*ux*dy
               uyao(al+1,i,w) = (3d0*u + ux*dy2)*dy2
               uyao(al+2,i,w) = dz2*dz*ux*dy
               uzao(al,i,w)   = dx2*dx*ux*dz
               uzao(al+1,i,w) = dy2*dy*ux*dz
               uzao(al+2,i,w) = (3d0*u + ux*dz2)*dz2
               u2ao(al,i,w)   = (6d0*u + dx2*uxx)*dx
               u2ao(al+1,i,w) = (6d0*u + dy2*uxx)*dy
               u2ao(al+2,i,w) = (6d0*u + dz2*uxx)*dz

c              // f_xxy, f_xxz, f_yyx, f_yyz, f_zzx, f_zzy
               u = sqr5*u                   ! correction of norm
               ux = sqr5*ux
               uxx = sqr5*uxx

               uao(al+4,i,w)  = dx2*dy*u
               uao(al+5,i,w)  = dx2*dz*u
               uao(al+3,i,w)  = dy2*dx*u
               uao(al+8,i,w)  = dy2*dz*u
               uao(al+6,i,w)  = dz2*dx*u
               uao(al+7,i,w)  = dz2*dy*u
               tmp = ux*dxyz
               uxao(al+4,i,w) = (2d0*u + ux*dx2)*dx*dy
               uxao(al+5,i,w) = (2d0*u + ux*dx2)*dx*dz
               uxao(al+3,i,w) = (u + ux*dx2)*dy2
               uxao(al+8,i,w) = tmp*dy
               uxao(al+6,i,w) = (u + ux*dx2)*dz2
               uxao(al+7,i,w) = tmp*dz
               uyao(al+4,i,w) = (u + ux*dy2)*dx2
               uyao(al+5,i,w) = tmp*dx
               uyao(al+3,i,w) = (2d0*u + ux*dy2)*dx*dy
               uyao(al+8,i,w) = (2d0*u + ux*dy2)*dy*dz
               uyao(al+6,i,w) = tmp*dz
               uyao(al+7,i,w) = (u + ux*dy2)*dz2
               uzao(al+4,i,w) = tmp*dx
               uzao(al+5,i,w) = (u + ux*dz2)*dx2
               uzao(al+3,i,w) = tmp*dy
               uzao(al+8,i,w) = (u + ux*dz2)*dy2
               uzao(al+6,i,w) = (2d0*u + ux*dz2)*dx*dz
               uzao(al+7,i,w) = (2d0*u + ux*dz2)*dy*dz

               u2ao(al+4,i,w) = (2d0*u + dx2*uxx)*dy
               u2ao(al+5,i,w) = (2d0*u + dx2*uxx)*dz
               u2ao(al+3,i,w) = (2d0*u + dy2*uxx)*dx
               u2ao(al+8,i,w) = (2d0*u + dy2*uxx)*dz
               u2ao(al+6,i,w) = (2d0*u + dz2*uxx)*dx
               u2ao(al+7,i,w) = (2d0*u + dz2*uxx)*dy

c              // f_xyz
               u = sqr3*u                  ! correction of norm
               ux = sqr3*ux
               uxx = sqr3*uxx

               uao(al+9,i,w)  = dxyz*u

               uxao(al+9,i,w) = (u + ux*dx2)*dy*dz
               uyao(al+9,i,w) = (u + ux*dy2)*dx*dz
               uzao(al+9,i,w) = (u + ux*dz2)*dx*dy
               u2ao(al+9,i,w) = dxyz*uxx
               al = al +10

            else
               call abortp('(aosplcalc): wrong GTO')
            endif       ! bl
            endif	! simple GTO, no splines
           endif        ! STO/GTO
         enddo          ! bf-loop over basis functions
      enddo             ! i-loop over electrons

      enddo ! w-loop over elec configs


      CONTAINS

         subroutine internal_GamessOrderGFunctions()

            alp = cntrctn(1,1,bf)
            u = cntrctn(2,1,bf) * exp(-alp*r2)
            dx = xx-atoms(a)%cx
            dx2 = dx*dx
            dy = yy-atoms(a)%cy
            dy2 = dy*dy
            dz = zz-atoms(a)%cz
            dz2 = dz*dz
            dxyz = dx*dy*dz
            ux = -2d0*alp*u

c           // g_xxxx, g_yyyy, g_zzzz
            uao(al,i,w)    = uao(al,i,w)    + dx2*dx2*u
            uao(al+1,i,w)  = uao(al+1,i,w)  + dy2*dy2*u
            uao(al+2,i,w)  = uao(al+2,i,w)  + dz2*dz2*u

            uxao(al,i,w)   = uxao(al,i,w)   + (4d0*u + ux*dx2)*dx2*dx
            uxao(al+1,i,w) = uxao(al+1,i,w) + dy2*dy2*ux*dx
            uxao(al+2,i,w) = uxao(al+2,i,w) + dz2*dz2*ux*dx
            uyao(al,i,w)   = uyao(al,i,w)   + dx2*dx2*ux*dy
            uyao(al+1,i,w) = uyao(al+1,i,w) + (4d0*u + ux*dy2)*dy2*dy
            uyao(al+2,i,w) = uyao(al+2,i,w) + dz2*dz2*ux*dy
            uzao(al,i,w)   = uzao(al,i,w)   + dx2*dx2*ux*dz
            uzao(al+1,i,w) = uzao(al+1,i,w) + dy2*dy2*ux*dz
            uzao(al+2,i,w) = uzao(al+2,i,w) + (4d0*u + ux*dz2)*dz2*dz
            tmp          = (11d0 - 2d0*alp*r2)*ux
            u2ao(al,i,w)   = u2ao(al,i,w)   + (12d0*u + dx2*tmp)*dx2
            u2ao(al+1,i,w) = u2ao(al+1,i,w) + (12d0*u + dy2*tmp)*dy2
            u2ao(al+2,i,w) = u2ao(al+2,i,w) + (12d0*u + dz2*tmp)*dz2

c           // g_xxxy, g_xxxz, g_yyyx, g_yyyz, g_zzzx, g_zzzy
            u = sqr7*u                   ! correction of norm
            ux = sqr7*ux

            uao(al+3,i,w)  = uao(al+3,i,w)  + dx2*dx*dy*u
            uao(al+4,i,w)  = uao(al+4,i,w)  + dx2*dx*dz*u
            uao(al+5,i,w)  = uao(al+5,i,w)  + dy2*dy*dx*u
            uao(al+6,i,w)  = uao(al+6,i,w)  + dy2*dy*dz*u
            uao(al+7,i,w)  = uao(al+7,i,w)  + dz2*dz*dx*u
            uao(al+8,i,w)  = uao(al+8,i,w)  + dz2*dz*dy*u

            tmp = ux*dxyz
            uxao(al+3,i,w) = uxao(al+3,i,w) + (3d0*u + ux*dx2)*dx2*dy
            uxao(al+4,i,w) = uxao(al+4,i,w) + (3d0*u + ux*dx2)*dx2*dz
            uxao(al+5,i,w) = uxao(al+5,i,w) + (u + ux*dx2)*dy2*dy
            uxao(al+6,i,w) = uxao(al+6,i,w) + tmp*dy2
            uxao(al+7,i,w) = uxao(al+7,i,w) + (u + ux*dx2)*dz2*dz
            uxao(al+8,i,w) = uxao(al+8,i,w) + tmp*dz2
            uyao(al+3,i,w) = uyao(al+3,i,w) + (u + ux*dy2)*dx2*dx
            uyao(al+4,i,w) = uyao(al+4,i,w) + tmp*dx2
            uyao(al+5,i,w) = uyao(al+5,i,w) + (3d0*u + ux*dy2)*dy2*dx
            uyao(al+6,i,w) = uyao(al+6,i,w) + (3d0*u + ux*dy2)*dy2*dz
            uyao(al+7,i,w) = uyao(al+7,i,w) + tmp*dz2
            uyao(al+8,i,w) = uyao(al+8,i,w) + (u + ux*dy2)*dz2*dz
            uzao(al+3,i,w) = uzao(al+3,i,w) + tmp*dx2
            uzao(al+4,i,w) = uzao(al+4,i,w) + (u + ux*dz2)*dx2*dx
            uzao(al+5,i,w) = uzao(al+5,i,w) + tmp*dy2
            uzao(al+6,i,w) = uzao(al+6,i,w) + (u + ux*dz2)*dy2*dy
            uzao(al+7,i,w) = uzao(al+7,i,w) + (3d0*u + ux*dz2)*dz2*dx
            uzao(al+8,i,w) = uzao(al+8,i,w) + (3d0*u + ux*dz2)*dz2*dy

            tmp = (11d0 - 2d0*alp*r2)*ux
            u2ao(al+3,i,w) = u2ao(al+3,i,w) + (6d0*u + dx2*tmp)*dx*dy
            u2ao(al+4,i,w) = u2ao(al+4,i,w) + (6d0*u + dx2*tmp)*dx*dz
            u2ao(al+5,i,w) = u2ao(al+5,i,w) + (6d0*u + dy2*tmp)*dy*dx
            u2ao(al+6,i,w) = u2ao(al+6,i,w) + (6d0*u + dy2*tmp)*dy*dz
            u2ao(al+7,i,w) = u2ao(al+7,i,w) + (6d0*u + dz2*tmp)*dz*dx
            u2ao(al+8,i,w) = u2ao(al+8,i,w) + (6d0*u + dz2*tmp)*dz*dy

c           // g_xxyy, g_xxzz, g_yyzz
            u = sqr5 / sqr3 * u          ! correction of norm
            ux = sqr5 / sqr3 * ux

            uao(al+9,i,w)   = uao(al+9,i,w)  + dx2*dy2*u
            uao(al+10,i,w)  = uao(al+10,i,w)  + dx2*dz2*u
            uao(al+11,i,w)  = uao(al+11,i,w)  + dy2*dz2*u

            tmp = ux*dxyz
            uxao(al+9,i,w)  = uxao(al+9,i,w)  + (2d0*u + ux*dx2)*dx*dy2
            uxao(al+10,i,w) = uxao(al+10,i,w) + (2d0*u + ux*dx2)*dx*dz2
            uxao(al+11,i,w) = uxao(al+11,i,w) + tmp*dy*dz
            uyao(al+9,i,w)  = uyao(al+9,i,w)  + (2d0*u + ux*dy2)*dy*dx2
            uyao(al+10,i,w) = uyao(al+10,i,w) + tmp*dx*dz
            uyao(al+11,i,w) = uyao(al+11,i,w) + (2d0*u + ux*dy2)*dy*dz2
            uzao(al+9,i,w)  = uzao(al+9,i,w)  + tmp*dx*dy
            uzao(al+10,i,w) = uzao(al+10,i,w) + (2d0*u + ux*dz2)*dz*dx2
            uzao(al+11,i,w) = uzao(al+11,i,w) + (2d0*u + ux*dz2)*dz*dy2

            tmp = (11d0 - 2d0*alp*r2)*ux
            u2ao(al+9,i,w)  = u2ao(al+9,i,w)  + 2d0*u*(dx2+dy2) + dx2*dy2*tmp
            u2ao(al+10,i,w) = u2ao(al+10,i,w) + 2d0*u*(dx2+dz2) + dx2*dz2*tmp
            u2ao(al+11,i,w) = u2ao(al+11,i,w) + 2d0*u*(dy2+dz2) + dy2*dz2*tmp


c           // g_xxyz, g_yyxz, g_zzxy
            u = sqr3*u                  ! correction of norm
            ux = sqr3*ux

            uao(al+12,i,w)  = uao(al+12,i,w)  + dx*dxyz*u
            uao(al+13,i,w)  = uao(al+13,i,w)  + dy*dxyz*u
            uao(al+14,i,w)  = uao(al+14,i,w)  + dz*dxyz*u

            uxao(al+12,i,w) = uxao(al+12,i,w) + (2d0*u + ux*dx2)*dxyz
            uxao(al+13,i,w) = uxao(al+13,i,w) + (u + ux*dx2)*dy2*dz
            uxao(al+14,i,w) = uxao(al+14,i,w) + (u + ux*dx2)*dz2*dy
            uyao(al+12,i,w) = uyao(al+12,i,w) + (u + ux*dy2)*dx2*dz
            uyao(al+13,i,w) = uyao(al+13,i,w) + (2d0*u + ux*dy2)*dxyz
            uyao(al+14,i,w) = uyao(al+14,i,w) + (u + ux*dy2)*dz2*dx
            uzao(al+12,i,w) = uzao(al+12,i,w) + (u + ux*dz2)*dx2*dy
            uzao(al+13,i,w) = uzao(al+13,i,w) + (u + ux*dz2)*dy2*dx
            uzao(al+14,i,w) = uzao(al+14,i,w) + (2d0*u + ux*dz2)*dxyz

            tmp = (11d0 - 2d0*alp*r2)*ux
            u2ao(al+12,i,w) = u2ao(al+12,i,w) + (2d0*u + tmp*dx2)*dy*dz
            u2ao(al+13,i,w) = u2ao(al+13,i,w) + (2d0*u + tmp*dy2)*dx*dz
            u2ao(al+14,i,w) = u2ao(al+14,i,w) + (2d0*u + tmp*dz2)*dx*dy

         end subroutine internal_GamessOrderGFunctions


         subroutine internal_GaussianOrderGFunctions()

            alp = cntrctn(1,1,bf)
            u = cntrctn(2,1,bf) * exp(-alp*r2)
            dx = xx-atoms(a)%cx
            dx2 = dx*dx
            dy = yy-atoms(a)%cy
            dy2 = dy*dy
            dz = zz-atoms(a)%cz
            dz2 = dz*dz
            dxyz = dx*dy*dz
            ux = -2d0*alp*u

c           // g_xxxx, g_yyyy, g_zzzz
            uao(al+14,i,w)    = uao(al+14,i,w)    + dx2*dx2*u
            uao(al+4,i,w)  = uao(al+4,i,w)  + dy2*dy2*u
            uao(al,i,w)  = uao(al,i,w)  + dz2*dz2*u

            uxao(al+14,i,w)   = uxao(al+14,i,w)   + (4d0*u + ux*dx2)*dx2*dx
            uxao(al+4,i,w) = uxao(al+4,i,w) + dy2*dy2*ux*dx
            uxao(al,i,w) = uxao(al,i,w) + dz2*dz2*ux*dx
            uyao(al+14,i,w)   = uyao(al+14,i,w)   + dx2*dx2*ux*dy
            uyao(al+4,i,w) = uyao(al+4,i,w) + (4d0*u + ux*dy2)*dy2*dy
            uyao(al,i,w) = uyao(al,i,w) + dz2*dz2*ux*dy
            uzao(al+14,i,w)   = uzao(al+14,i,w)   + dx2*dx2*ux*dz
            uzao(al+4,i,w) = uzao(al+4,i,w) + dy2*dy2*ux*dz
            uzao(al,i,w) = uzao(al,i,w) + (4d0*u + ux*dz2)*dz2*dz
            tmp          = (11d0 - 2d0*alp*r2)*ux
            u2ao(al+14,i,w)   = u2ao(al+14,i,w)   + (12d0*u + dx2*tmp)*dx2
            u2ao(al+4,i,w) = u2ao(al+4,i,w) + (12d0*u + dy2*tmp)*dy2
            u2ao(al,i,w) = u2ao(al,i,w) + (12d0*u + dz2*tmp)*dz2

c           // g_xxxy, g_xxxz, g_yyyx, g_yyyz, g_zzzx, g_zzzy
            u = sqr7*u                   ! correction of norm
            ux = sqr7*ux

            uao(al+13,i,w)  = uao(al+13,i,w)  + dx2*dx*dy*u
            uao(al+12,i,w)  = uao(al+12,i,w)  + dx2*dx*dz*u
            uao(al+8,i,w)  = uao(al+8,i,w)  + dy2*dy*dx*u
            uao(al+3,i,w)  = uao(al+3,i,w)  + dy2*dy*dz*u
            uao(al+5,i,w)  = uao(al+5,i,w)  + dz2*dz*dx*u
            uao(al+1,i,w)  = uao(al+1,i,w)  + dz2*dz*dy*u

            tmp = ux*dxyz
            uxao(al+13,i,w) = uxao(al+13,i,w) + (3d0*u + ux*dx2)*dx2*dy
            uxao(al+12,i,w) = uxao(al+12,i,w) + (3d0*u + ux*dx2)*dx2*dz
            uxao(al+8,i,w) = uxao(al+8,i,w) + (u + ux*dx2)*dy2*dy
            uxao(al+3,i,w) = uxao(al+3,i,w) + tmp*dy2
            uxao(al+5,i,w) = uxao(al+5,i,w) + (u + ux*dx2)*dz2*dz
            uxao(al+1,i,w) = uxao(al+1,i,w) + tmp*dz2
            uyao(al+13,i,w) = uyao(al+13,i,w) + (u + ux*dy2)*dx2*dx
            uyao(al+12,i,w) = uyao(al+12,i,w) + tmp*dx2
            uyao(al+8,i,w) = uyao(al+8,i,w) + (3d0*u + ux*dy2)*dy2*dx
            uyao(al+3,i,w) = uyao(al+3,i,w) + (3d0*u + ux*dy2)*dy2*dz
            uyao(al+5,i,w) = uyao(al+5,i,w) + tmp*dz2
            uyao(al+1,i,w) = uyao(al+1,i,w) + (u + ux*dy2)*dz2*dz
            uzao(al+13,i,w) = uzao(al+13,i,w) + tmp*dx2
            uzao(al+12,i,w) = uzao(al+12,i,w) + (u + ux*dz2)*dx2*dx
            uzao(al+8,i,w) = uzao(al+8,i,w) + tmp*dy2
            uzao(al+3,i,w) = uzao(al+3,i,w) + (u + ux*dz2)*dy2*dy
            uzao(al+5,i,w) = uzao(al+5,i,w) + (3d0*u + ux*dz2)*dz2*dx
            uzao(al+1,i,w) = uzao(al+1,i,w) + (3d0*u + ux*dz2)*dz2*dy

            tmp = (11d0 - 2d0*alp*r2)*ux
            u2ao(al+13,i,w) = u2ao(al+13,i,w) + (6d0*u + dx2*tmp)*dx*dy
            u2ao(al+12,i,w) = u2ao(al+12,i,w) + (6d0*u + dx2*tmp)*dx*dz
            u2ao(al+8,i,w) = u2ao(al+8,i,w) + (6d0*u + dy2*tmp)*dy*dx
            u2ao(al+3,i,w) = u2ao(al+3,i,w) + (6d0*u + dy2*tmp)*dy*dz
            u2ao(al+5,i,w) = u2ao(al+5,i,w) + (6d0*u + dz2*tmp)*dz*dx
            u2ao(al+1,i,w) = u2ao(al+1,i,w) + (6d0*u + dz2*tmp)*dz*dy

c           // g_xxyy, g_xxzz, g_yyzz
            u = sqr5 / sqr3 * u          ! correction of norm
            ux = sqr5 / sqr3 * ux

            uao(al+11,i,w)   = uao(al+11,i,w)  + dx2*dy2*u
            uao(al+9,i,w)  = uao(al+9,i,w)  + dx2*dz2*u
            uao(al+2,i,w)  = uao(al+2,i,w)  + dy2*dz2*u

            tmp = ux*dxyz
            uxao(al+11,i,w)  = uxao(al+11,i,w)  + (2d0*u + ux*dx2)*dx*dy2
            uxao(al+9,i,w) = uxao(al+9,i,w) + (2d0*u + ux*dx2)*dx*dz2
            uxao(al+2,i,w) = uxao(al+2,i,w) + tmp*dy*dz
            uyao(al+11,i,w)  = uyao(al+11,i,w)  + (2d0*u + ux*dy2)*dy*dx2
            uyao(al+9,i,w) = uyao(al+9,i,w) + tmp*dx*dz
            uyao(al+2,i,w) = uyao(al+2,i,w) + (2d0*u + ux*dy2)*dy*dz2
            uzao(al+11,i,w)  = uzao(al+11,i,w)  + tmp*dx*dy
            uzao(al+9,i,w) = uzao(al+9,i,w) + (2d0*u + ux*dz2)*dz*dx2
            uzao(al+2,i,w) = uzao(al+2,i,w) + (2d0*u + ux*dz2)*dz*dy2

            tmp = (11d0 - 2d0*alp*r2)*ux
            u2ao(al+11,i,w)  = u2ao(al+11,i,w)  + 2d0*u*(dx2+dy2) + dx2*dy2*tmp
            u2ao(al+9,i,w) = u2ao(al+9,i,w) + 2d0*u*(dx2+dz2) + dx2*dz2*tmp
            u2ao(al+2,i,w) = u2ao(al+2,i,w) + 2d0*u*(dy2+dz2) + dy2*dz2*tmp


c           // g_xxyz, g_yyxz, g_zzxy
            u = sqr3*u                  ! correction of norm
            ux = sqr3*ux

            uao(al+10,i,w)  = uao(al+10,i,w)  + dx*dxyz*u
            uao(al+7,i,w)  = uao(al+7,i,w)  + dy*dxyz*u
            uao(al+6,i,w)  = uao(al+6,i,w)  + dz*dxyz*u

            uxao(al+10,i,w) = uxao(al+10,i,w) + (2d0*u + ux*dx2)*dxyz
            uxao(al+7,i,w) = uxao(al+7,i,w) + (u + ux*dx2)*dy2*dz
            uxao(al+6,i,w) = uxao(al+6,i,w) + (u + ux*dx2)*dz2*dy
            uyao(al+10,i,w) = uyao(al+10,i,w) + (u + ux*dy2)*dx2*dz
            uyao(al+7,i,w) = uyao(al+7,i,w) + (2d0*u + ux*dy2)*dxyz
            uyao(al+6,i,w) = uyao(al+6,i,w) + (u + ux*dy2)*dz2*dx
            uzao(al+10,i,w) = uzao(al+10,i,w) + (u + ux*dz2)*dx2*dy
            uzao(al+7,i,w) = uzao(al+7,i,w) + (u + ux*dz2)*dy2*dx
            uzao(al+6,i,w) = uzao(al+6,i,w) + (2d0*u + ux*dz2)*dxyz

            tmp = (11d0 - 2d0*alp*r2)*ux
            u2ao(al+10,i,w) = u2ao(al+10,i,w) + (2d0*u + tmp*dx2)*dy*dz
            u2ao(al+7,i,w) = u2ao(al+7,i,w) + (2d0*u + tmp*dy2)*dx*dz
            u2ao(al+6,i,w) = u2ao(al+6,i,w) + (2d0*u + tmp*dz2)*dx*dy

         end subroutine internal_GaussianOrderGFunctions


      end subroutine aosplcalc


c================================================
cBH

c     -------------------------------
      subroutine ao1calc(ie,x,y,z,rai)
c     -------------------------------

c aocalc calculates all atomic orbitals for position vector x,y,z
c without derivatives. For ie>0 update only

c on entry: requires rai distance matrix properly dimensioned and set
c           for input position vector x,y,z.
c           for ie>0 position vector must be modified only compared to
c           last call only at electron ie
c on exit:  updates AO data structure in aos.h (i.e. AOs)

c Version 1.3 (28.4.98)    references to walker eliminated
c Version 1.1 (6/26/1996) simultaneous p,d calculation
c Version 1.0 (5/29/1996)

c input parameters:
      integer ie                 ! if >0 only AO's for electron ie recalculated
      real*8   x(:),          ! x,y,z coordinates of position vector
     .         y(:),
     .         z(:),
     .         rai(:,:)    ! r_ai electron-nucleus distances
c output (to local data structure (aos.h) in commons block):
c     uao  : AO-array uao(n,i,w) for nth AO at electron i
cEH
c constants:
      real*8 sqr3,sqr5,sqr7
      parameter (sqr3=1.73205080756887729d0,sqr5=2.236067977499789696d0,
     &           sqr7=2.645751311064591d0)
c variables
c     integer al,bf,a,i,ii,i1,i2,nn,ic
      integer al,bf,a,i,ii,i1,i2,nn,ic,m,n,k,l
      real*8 xx,yy,zz,rr,alp,nrm,u,dx,dy,dz,r2,dx2,dy2,dz2,dxyz
cTS
      integer nend
      logical gaussFOrder       ! .t.: Gaussian order for f function
                                ! .f.: Gamess==Turbomole order used


      mAOElecConfigs = 1

c-----Calculation of the AO's

      if (evfmt=='gau' .or. evfmt=='mol' ) then
         gaussFOrder = .true.
      else
         gaussFOrder = .false.
      endif

      if (ie .eq. 0) then                     ! AO's for all electrons
         i1 = 1
         i2 = ne
      else
         i1 = ie                              ! only AO for electron ie
         i2 = ie
      endif

      do i=i1,i2                              ! loop over electrons

         xx = x(i)
         yy = y(i)
         zz = z(i)

         al = 1

         do n = 1,nbasf                   ! loop over basis functions
           bf = n
           a = bc(bf)                         ! center of AO
           rr = rai(a,i)                      ! r_ai
           nn = bn(bf)                        ! n quantum no. of AO

           if (typ(bf) == 'STO') then      ! STO-type basis function

            alp = bzet(bf)                    ! orbital exponent
            nrm = norm(bf)                    ! normalization constant

            if (nn == 1) then               ! 1s orbital
               u  = nrm*exp(-alp*rr)
               uao(al,i,1) = u
               al = al+1
            else if (nn == 2) then
               if (bl(bf) == 'S') then           ! 2s orbital
c                 // 2s orbital
                  u  = nrm*exp(-alp*rr)
                  uao(al,i,1) = rr*u
                  al = al+1
               else if (bl(bf) == 'P') then      ! 2p orbital
                  u  = nrm*exp(-alp*rr)
                  dx = xx-atoms(a)%cx
                  dy = yy-atoms(a)%cy
                  dz = zz-atoms(a)%cz
                  uao(al,i,1)    = u*dx
                  uao(al+1,i,1)  = u*dy
                  uao(al+2,i,1)  = u*dz
                  al = al+3

               else
                  call abortp(' ao1calc: n=2 and l > 2')
               endif

            else if (nn == 3) then
               if (bl(bf) == 'S') then
c                 // 3s orbital
                  u  = nrm*exp(-alp*rr)
                  uao(al,i,1) = rr*rr*u
                  al = al+1
               else if (bl(bf) == 'P') then
                  call abortp('ao1calc: 3p orbitals not implemented')

               else if (bl(bf) == 'D') then         ! 3d orbital (6D)
c                 // Norm is different for d_xx and d_xy !
                  u   = nrm*exp(-alp*rr)
                  dx  = xx - atoms(a)%cx
                  dx2 = dx*dx
                  dy  = yy - atoms(a)%cy
                  dy2 = dy*dy
                  dz  = zz - atoms(a)%cz
                  dz2 = dz*dz

                  uao(al,i,1)    = dx2*u
                  uao(al+1,i,1)  = dy2*u
                  uao(al+2,i,1)  = dz2*u

                  u = sqr3*u                   ! correction of norm for last 3
                  uao(al+3,i,1)  = u*dx*dy
                  uao(al+4,i,1)  = u*dx*dz
                  uao(al+5,i,1)  = u*dy*dz
                  al = al+6
               endif     ! bl
            endif        ! nn


c          // Contracted GTO's as basis function (AO)
           else

            r2 = rr*rr

c           // only primitive cartesian gaussians: 1s,2p,3d,4f
c           // i.e. no r factor. Thus nn is not used here.

            if (bl(bf) == 'S') then                 ! 1s GTO
               uao(al,i,1) = 0d0
               do ic=1,ngto(bf)                     ! loop over contraction
                  alp = cntrctn(1,ic,bf)
                  u = cntrctn(2,ic,bf) * exp(-alp*r2)
                  uao(al,i,1) = uao(al,i,1) + u
               enddo
               al = al+1

            else if (bl(bf) == 'P') then             ! 2p GTO's
c              // do all 3 P simultaneously (same exponent is required)
c              // order p_x,p_y,p_z
               do ii=0,2
                  uao(al+ii,i,1) = 0d0
               enddo
               do ic=1,ngto(bf)                      ! loop over contraction
                  alp = cntrctn(1,ic,bf)
                  u = cntrctn(2,ic,bf) * exp(-alp*r2)
                  dx = xx-atoms(a)%cx
                  dy = yy-atoms(a)%cy
                  dz = zz-atoms(a)%cz
                  uao(al,i,1) = uao(al,i,1) + dx*u
                  uao(al+1,i,1) = uao(al+1,i,1) + dy*u
                  uao(al+2,i,1) = uao(al+2,i,1) + dz*u
               enddo
               al = al+3

            else if (bl(bf) == 'D') then         ! 3d GTO
c              // do all 6 D simultaneously (same exponent is required)
c              // order: d_xx, d_yy, d_zz, d_xy, d_xz, d_yz  (like GAMESS)
               do ii=0,5
                  uao(al+ii,i,1) = 0d0
               enddo
               do ic=1,ngto(bf)                      ! loop over contraction
                  alp = cntrctn(1,ic,bf)
                  u = cntrctn(2,ic,bf) * exp(-alp*r2)
                  dx = xx-atoms(a)%cx
                  dx2 = dx*dx
                  dy = yy-atoms(a)%cy
                  dy2 = dy*dy
                  dz = zz-atoms(a)%cz
                  dz2 = dz*dz

                  uao(al,i,1)    = uao(al,i,1)    + dx2*u
                  uao(al+1,i,1)  = uao(al+1,i,1)  + dy2*u
                  uao(al+2,i,1)  = uao(al+2,i,1)  + dz2*u

                  u = sqr3*u                   ! correction of norm for last 3

                  uao(al+3,i,1)  = uao(al+3,i,1)  + dx*dy*u
                  uao(al+4,i,1)  = uao(al+4,i,1)  + dx*dz*u
                  uao(al+5,i,1)  = uao(al+5,i,1)  + dy*dz*u
               enddo
               al = al+6

            else if (bl(bf)=='F'.and..not.gaussFOrder) then     ! f GTO
c              // do all 10 F simultaneously (same exponent is required)
c              // order: f_xxx, f_yyy, f_zzz, fd_xxy, f_xxz, f_yyx,
c              //   f_yyz, f_zzx, f_zzy, f_xyz  (like GAMESS)
               do ii=0,9
                  uao(al+ii,i,1) = 0d0
               enddo
               do ic=1,ngto(bf)                      ! loop over contraction
                  alp = cntrctn(1,ic,bf)
                  u = cntrctn(2,ic,bf) * exp(-alp*r2)
                  dx = xx-atoms(a)%cx
                  dx2 = dx*dx
                  dy = yy-atoms(a)%cy
                  dy2 = dy*dy
                  dz = zz-atoms(a)%cz
                  dz2 = dz*dz
                  dxyz = dx*dy*dz

c                 // f_xxx, f_yyy, f_zzz
                  uao(al,i,1)    = uao(al,i,1)    + dx2*dx*u
                  uao(al+1,i,1)  = uao(al+1,i,1)  + dy2*dy*u
                  uao(al+2,i,1)  = uao(al+2,i,1)  + dz2*dz*u

c                 // f_xxy, f_xxz, f_yyx, f_yyz, f_zzx, f_zzy
                  u = sqr5*u                   ! correction of norm

                  uao(al+3,i,1)  = uao(al+3,i,1)  + dx2*dy*u
                  uao(al+4,i,1)  = uao(al+4,i,1)  + dx2*dz*u
                  uao(al+5,i,1)  = uao(al+5,i,1)  + dy2*dx*u
                  uao(al+6,i,1)  = uao(al+6,i,1)  + dy2*dz*u
                  uao(al+7,i,1)  = uao(al+7,i,1)  + dz2*dx*u
                  uao(al+8,i,1)  = uao(al+8,i,1)  + dz2*dy*u
c                 // f_xyz
                  u = sqr3*u                  ! correction of norm
                  uao(al+9,i,1)  = uao(al+9,i,1)  + dxyz*u
               enddo
               al = al+10

            else if (bl(bf)=='F'.and.gaussFOrder) then     ! 3f GTO
c              // do all 10 F simultaneously (same exponent is required)
c              // order: f_xxx, f_yyy, f_zzz, f_xyy, f_xxy, f_xxz,
c              //   f_xzz, f_yzz, f_yyz, f_xyz  (like Gaussian)
               do ii=0,9
                  uao(al+ii,i,1) = 0d0
               enddo
               do ic=1,ngto(bf)                      ! loop over contraction
                  alp = cntrctn(1,ic,bf)
                  u = cntrctn(2,ic,bf) * exp(-alp*r2)
                  dx = xx-atoms(a)%cx
                  dx2 = dx*dx
                  dy = yy-atoms(a)%cy
                  dy2 = dy*dy
                  dz = zz-atoms(a)%cz
                  dz2 = dz*dz
                  dxyz = dx*dy*dz

c                 // f_xxx, f_yyy, f_zzz
                  uao(al,i,1)    = uao(al,i,1)    + dx2*dx*u
                  uao(al+1,i,1)  = uao(al+1,i,1)  + dy2*dy*u
                  uao(al+2,i,1)  = uao(al+2,i,1)  + dz2*dz*u

c                 // f_xxy, f_xxz, f_yyx, f_yyz, f_zzx, f_zzy
                  u = sqr5*u                   ! correction of norm

                  uao(al+3,i,1)  = uao(al+3,i,1)  + dy2*dx*u
                  uao(al+4,i,1)  = uao(al+4,i,1)  + dx2*dy*u
                  uao(al+5,i,1)  = uao(al+5,i,1)  + dx2*dz*u
                  uao(al+6,i,1)  = uao(al+6,i,1)  + dz2*dx*u
                  uao(al+7,i,1)  = uao(al+7,i,1)  + dz2*dy*u
                  uao(al+8,i,1)  = uao(al+8,i,1)  + dy2*dz*u
c                 // f_xyz
                  u = sqr3*u                  ! correction of norm
                  uao(al+9,i,1)  = uao(al+9,i,1)  + dxyz*u
               enddo
               al = al + 10

            else if (bl(bf)=='G') then     ! 5g GTO
c              // do all 15 cartesian G simultaneously (same exponent is required)

               uao(al:al+14,i,1) = 0d0

               if (gaussFOrder) then
                  call internal_GaussianOrderGFunctions()
               else
                  call internal_GamessOrderGFunctions()
               endif

               al = al + 15


            else
               call abortp('(ao1calc): wrong GTO')
            endif  ! bl
           endif  ! STO/GTO
         enddo  ! bf-loop over basis functions
      enddo  ! i-loop over electrons


      CONTAINS

         subroutine internal_GamessOrderGFunctions()

            do ic=1,ngto(bf)                      ! loop over contraction
               alp = cntrctn(1,ic,bf)
               u = cntrctn(2,ic,bf) * exp(-alp*r2)
               dx = xx-atoms(a)%cx
               dx2 = dx*dx
               dy = yy-atoms(a)%cy
               dy2 = dy*dy
               dz = zz-atoms(a)%cz
               dz2 = dz*dz
               dxyz = dx*dy*dz

c              // g_xxxx, g_yyyy, g_zzzz
               uao(al,i,1)    = uao(al,i,1)    + dx2*dx2*u
               uao(al+1,i,1)  = uao(al+1,i,1)  + dy2*dy2*u
               uao(al+2,i,1)  = uao(al+2,i,1)  + dz2*dz2*u

c              // g_xxxy, g_xxxz, g_yyyx, g_yyyz, g_zzzx, g_zzzy
               u = sqr7*u                   ! correction of norm

               uao(al+3,i,1)  = uao(al+3,i,1)  + dx2*dx*dy*u
               uao(al+4,i,1)  = uao(al+4,i,1)  + dx2*dx*dz*u
               uao(al+5,i,1)  = uao(al+5,i,1)  + dy2*dy*dx*u
               uao(al+6,i,1)  = uao(al+6,i,1)  + dy2*dy*dz*u
               uao(al+7,i,1)  = uao(al+7,i,1)  + dz2*dz*dx*u
               uao(al+8,i,1)  = uao(al+8,i,1)  + dz2*dz*dy*u

c              // g_xxyy, g_xxzz, g_yyzz
               u = sqr5 / sqr3 * u          ! correction of norm

               uao(al+9,i,1)   = uao(al+9,i,1)  + dx2*dy2*u
               uao(al+10,i,1)  = uao(al+10,i,1)  + dx2*dz2*u
               uao(al+11,i,1)  = uao(al+11,i,1)  + dy2*dz2*u

c              // g_xxyz, g_yyxz, g_zzxy
               u = sqr3*u                  ! correction of norm

               uao(al+12,i,1)  = uao(al+12,i,1)  + dx*dxyz*u
               uao(al+13,i,1)  = uao(al+13,i,1)  + dy*dxyz*u
               uao(al+14,i,1)  = uao(al+14,i,1)  + dz*dxyz*u
            enddo

         end subroutine internal_GamessOrderGFunctions


         subroutine internal_GaussianOrderGFunctions()
               do ic=1,ngto(bf)                      ! loop over contraction
               alp = cntrctn(1,ic,bf)
               u = cntrctn(2,ic,bf) * exp(-alp*r2)
               dx = xx-atoms(a)%cx
               dx2 = dx*dx
               dy = yy-atoms(a)%cy
               dy2 = dy*dy
               dz = zz-atoms(a)%cz
               dz2 = dz*dz
               dxyz = dx*dy*dz

c              // g_xxxx, g_yyyy, g_zzzz
               uao(al+14,i,1)    = uao(al+14,i,1)    + dx2*dx2*u
               uao(al+4,i,1)  = uao(al+4,i,1)  + dy2*dy2*u
               uao(al,i,1)  = uao(al,i,1)  + dz2*dz2*u

c              // g_xxxy, g_xxxz, g_yyyx, g_yyyz, g_zzzx, g_zzzy
               u = sqr7*u                   ! correction of norm

               uao(al+13,i,1)  = uao(al+13,i,1)  + dx2*dx*dy*u
               uao(al+12,i,1)  = uao(al+12,i,1)  + dx2*dx*dz*u
               uao(al+8,i,1)  = uao(al+8,i,1)  + dy2*dy*dx*u
               uao(al+3,i,1)  = uao(al+3,i,1)  + dy2*dy*dz*u
               uao(al+6,i,1)  = uao(al+6,i,1)  + dz2*dz*dx*u
               uao(al+1,i,1)  = uao(al+1,i,1)  + dz2*dz*dy*u

c              // g_xxyy, g_xxzz, g_yyzz
               u = sqr5 / sqr3 * u          ! correction of norm

               uao(al+11,i,1)   = uao(al+11,i,1)  + dx2*dy2*u
               uao(al+9,i,1)  = uao(al+9,i,1)  + dx2*dz2*u
               uao(al+2,i,1)  = uao(al+2,i,1)  + dy2*dz2*u

c              // g_xxyz, g_yyxz, g_zzxy
               u = sqr3*u                  ! correction of norm

               uao(al+10,i,1)  = uao(al+10,i,1)  + dx*dxyz*u
               uao(al+7,i,1)  = uao(al+7,i,1)  + dy*dxyz*u
               uao(al+6,i,1)  = uao(al+6,i,1)  + dz*dxyz*u
            enddo
            call abortp("gaussian order g functions are not yet implemented")
         end subroutine internal_GaussianOrderGFunctions

      end subroutine ao1calc

c ==================================================
cBH

c     ----------------------------------
      subroutine ao1splcalc(ie,x,y,z,rai)
c     ----------------------------------

c ao1splcalc calculates all atomic orbitals for position vector x,y,z


c on entry: requires rai distance matrix properly dimensioned and set
c           for input position vector x,y,z.
c           for ie>0 position vector must be modified only compared to
c           last call to aocalc only at electron ie

c input parameters:
      integer ie                 ! if >0 only AO's for electron ie recalculated
      real*8   x(:),          ! x,y,z coordinates of position vector
     .         y(:),
     .         z(:),
     .        rai(:,:)     ! r_ai electron-nucleus distances

cEH
c constants:
      real*8 sqr3,sqr5,sqr7
      parameter (sqr3=1.73205080756887729d0,sqr5=2.236067977499789696d0,
     &           sqr7=2.645751311064591d0)
c variables
c     integer bf,a,i,j,i1,i2,nn,al,ispl
      integer al,bf,a,i,j,ii,i1,i2,nn,ic,m,n,k,l,ispl

      real*8 xx,yy,zz,rr,r2,alp,nrm,u,ux,uxx,u2,dx,dy,dz,tmp,
     .       dx2,dy2,dz2,dxyz,df
cTS
      integer nend
      logical gaussFOrder       ! .t.: Gaussian order for f function
                                ! .f.: Gamess==Turbomole order used


c bf refers to the degenerate set of cartesian
c basis function (S:1,P:3,D:6,F:10) as input, which may be of type STO
c or contracted GTO.
c al refers to the individual basis function, as used in LCAO-MO's.
c (composed in subroutine mdetwf)
c i refers to the current electron.
c-----Calculation of the AO's and their derivatives

      mAOElecConfigs = 1

      if (evfmt=='gau' .or. evfmt=='mol' ) then
         gaussFOrder = .true.
      else
         gaussFOrder = .false.
      endif

      if (ie .eq. 0) then                     ! AO's for all electrons
         i1 = 1
         i2 = ne
      else
         i1 = ie                              ! only AO for electron ie
         i2 = ie
      endif

      do i=i1,i2                              ! loop over electrons

         xx = x(i)
         yy = y(i)
         zz = z(i)

         al = 1

         ! loop over basis functions:
         do n = 1,nbasf
           bf = n
           a = bc(bf)                         ! center of AO
           rr = rai(a,i)                      ! r_ai
           nn = bn(bf)                        ! n quantum no. of AO

           if (typ(bf) == 'STO') then      ! STO-type basis function
             call abortp('STOs make no sense in ao1splcalc')

c          // Contracted GTO's as basis function (AO)
c          // Evaluate with splines
           else
           if (so(bf)==0) then !only 1 GTO in contraction, no splines used !

            r2 = rr*rr

c           // only primitive cartesian gaussians: 1s,2p,3d,4f
c           // i.e. no r factor. Thus nn is not used here.

            if (bl(bf) == 'S') then                 ! 1s GTO
               alp = cntrctn(1,1,bf)
               u   = cntrctn(2,1,bf) * exp(-alp*r2)
               uao(al,i,1)  = u
               al = al + 1
            else if (bl(bf) == 'P') then             ! 2p GTO's
c              // do all 3 P simultaneously (same exponent is required)
c              // order p_x,p_y,p_z
               alp = cntrctn(1,1,bf)
               u = cntrctn(2,1,bf) * exp(-alp*r2)
               dx = xx-atoms(a)%cx
               dy = yy-atoms(a)%cy
               dz = zz-atoms(a)%cz
               uao(al,i,1)    = dx*u
               uao(al+1,i,1)  = dy*u
               uao(al+2,i,1)  = dz*u

               al = al + 3

            else if (bl(bf) == 'D') then         ! 3d GTO
c              // do all 6 D simultaneously (same exponent is required)
c              // order: d_xx, d_yy, d_zz, d_xy, d_xz, d_yz  (like GAMESS)
               alp = cntrctn(1,1,bf)
               u = cntrctn(2,1,bf) * exp(-alp*r2)
               dx = xx-atoms(a)%cx
               dx2 = dx*dx
               dy = yy-atoms(a)%cy
               dy2 = dy*dy
               dz = zz-atoms(a)%cz
               dz2 = dz*dz
               uao(al,i,1)    = dx2*u
               uao(al+1,i,1)  = dy2*u
               uao(al+2,i,1)  = dz2*u

               u = sqr3*u                   ! correction of norm for last 3

               uao(al+3,i,1)  = dx*dy*u
               uao(al+4,i,1)  = dx*dz*u
               uao(al+5,i,1)  = dy*dz*u

               al = al + 6

            else if (bl(bf)=='F'.and..not.gaussFOrder) then     ! 3f GTO
c              // do all 10 F simultaneously (same exponent is required)
c              // order: f_xxx, f_yyy, f_zzz, fd_xxy, f_xxz, f_yyx,
c              //   f_yyz, f_zzx, f_zzy, f_xyz  (like GAMESS)
               alp = cntrctn(1,1,bf)
               u = cntrctn(2,1,bf) * exp(-alp*r2)
               dx = xx-atoms(a)%cx
               dx2 = dx*dx
               dy = yy-atoms(a)%cy
               dy2 = dy*dy
               dz = zz-atoms(a)%cz
               dz2 = dz*dz
               dxyz = dx*dy*dz

c              // f_xxx, f_yyy, f_zzz
               uao(al,i,1)    = dx2*dx*u
               uao(al+1,i,1)  = dy2*dy*u
               uao(al+2,i,1)  = dz2*dz*u

c              // f_xxy, f_xxz, f_yyx, f_yyz, f_zzx, f_zzy
               u = sqr5*u                   ! correction of norm
               uao(al+3,i,1)  = dx2*dy*u
               uao(al+4,i,1)  = dx2*dz*u
               uao(al+5,i,1)  = dy2*dx*u
               uao(al+6,i,1)  = dy2*dz*u
               uao(al+7,i,1)  = dz2*dx*u
               uao(al+8,i,1)  = dz2*dy*u

c              // f_xyz
               u = sqr3*u                  ! correction of norm
               uao(al+9,i,1)  = dxyz*u

               al = al + 10

            else if (bl(bf)=='F'.and.gaussFOrder) then         ! 3f GTO
c              // do all 10 F simultaneously (same exponent is required)
c              // order: f_xxx, f_yyy, f_zzz, f_xyy, f_xxy, f_xxz,
c              //   f_xzz, f_yzz, f_yyz, f_xyz  (like Gaussian)
               alp = cntrctn(1,1,bf)
               u = cntrctn(2,1,bf) * exp(-alp*r2)
               dx = xx-atoms(a)%cx
               dx2 = dx*dx
               dy = yy-atoms(a)%cy
               dy2 = dy*dy
               dz = zz-atoms(a)%cz
               dz2 = dz*dz
               dxyz = dx*dy*dz

c              // f_xxx, f_yyy, f_zzz
               uao(al,i,1)    = dx2*dx*u
               uao(al+1,i,1)  = dy2*dy*u
               uao(al+2,i,1)  = dz2*dz*u

c              // f_xyy, f_xxy, f_xxz, f_xzz, f_yzz, f_yyz
               u = sqr5*u                   ! correction of norm
               uao(al+3,i,1)  = dy2*dx*u
               uao(al+4,i,1)  = dx2*dy*u
               uao(al+5,i,1)  = dx2*dz*u
               uao(al+6,i,1)  = dz2*dx*u
               uao(al+7,i,1)  = dz2*dy*u
               uao(al+8,i,1)  = dy2*dz*u

c              // f_xyz
               u = sqr3*u                  ! correction of norm
               uao(al+9,i,1)  = dxyz*u

               al = al + 10

            else if (bl(bf)=='G') then     ! 5g GTO
c              // do all 15 cartesian G simultaneously (same exponent is required)
               uao(al:al+14,i,1) = 0d0

               if (gaussFOrder) then
                  call internal_GaussianOrderGFunctions()
               else
                  call internal_GamessOrderGFunctions()
               endif

               al = al + 15

            else
               call abortp('(getaos): wrong GTO')
            endif  ! bl


           else	 !more then 1 GTO in contraction, splines used !
            r2 = rr*rr
            j  = (csplnpnt-1)*rr/(csalpha+rr)  + 1
            df = rr - csplx(j)

c           // only primitive cartesian gaussians: 1s,2p,3d,4f
c           // i.e. no r factor. Thus nn is not used here.

            if (bl(bf) == 'S') then                 ! 1s GTO
               ispl       = 3*so(bf)-2
               uao(al,i,1)  = cspla(ispl,j) + df*(csplb(ispl,j)
     .                    + df*(csplc(ispl,j) + df*cspld(ispl,j)))

               al = al + 1

            else if (bl(bf) == 'P') then             ! 2p GTO's
c              // do all 3 P simultaneously (same exponent is required)
c              // order p_x,p_y,p_z
               ispl = 3*so(bf)-2
               u    = cspla(ispl,j) + df*(csplb(ispl,j)
     .              + df*(csplc(ispl,j) + df*cspld(ispl,j)))

               dx = xx-atoms(a)%cx
               dy = yy-atoms(a)%cy
               dz = zz-atoms(a)%cz
               uao(al,i,1) = dx*u
               uao(al+1,i,1) = dy*u
               uao(al+2,i,1) = dz*u

               al = al + 3

            else if (bl(bf) == 'D') then         ! 3d GTO
c              // do all 6 D simultaneously (same exponent is required)
c              // order: d_xx, d_yy, d_zz, d_xy, d_xz, d_yz  (like GAMESS)
               ispl = 3*so(bf)-2
               u    = cspla(ispl,j) + df*(csplb(ispl,j)
     .              + df*(csplc(ispl,j) + df*cspld(ispl,j)))

               dx = xx-atoms(a)%cx
               dx2 = dx*dx
               dy = yy-atoms(a)%cy
               dy2 = dy*dy
               dz = zz-atoms(a)%cz
               dz2 = dz*dz
               uao(al,i,1)    = dx2*u
               uao(al+1,i,1)  = dy2*u
               uao(al+2,i,1)  = dz2*u

               u = sqr3*u                   ! correction of norm for last 3

               uao(al+3,i,1)  = dx*dy*u
               uao(al+4,i,1)  = dx*dz*u
               uao(al+5,i,1)  = dy*dz*u

               al = al + 6

            else if (bl(bf)=='F'.and..not.gaussFOrder) then  ! 3f GTO
c              // do all 10 F simultaneously (same exponent is required)
c              // order: f_xxx, f_yyy, f_zzz, fd_xxy, f_xxz, f_yyx,
c              //   f_yyz, f_zzx, f_zzy, f_xyz  (like GAMESS)

               ispl = 3*so(bf)-2
               u    = cspla(ispl,j) + df*(csplb(ispl,j)
     .              + df*(csplc(ispl,j) + df*cspld(ispl,j)))

               dx = xx-atoms(a)%cx
               dx2 = dx*dx
               dy = yy-atoms(a)%cy
               dy2 = dy*dy
               dz = zz-atoms(a)%cz
               dz2 = dz*dz
               dxyz = dx*dy*dz

c              // f_xxx, f_yyy, f_zzz
               uao(al,i,1)    = dx2*dx*u
               uao(al+1,i,1)  = dy2*dy*u
               uao(al+2,i,1)  = dz2*dz*u

c              // f_xxy, f_xxz, f_yyx, f_yyz, f_zzx, f_zzy
               u = sqr5*u                   ! correction of norm

               uao(al+3,i,1)  = dx2*dy*u
               uao(al+4,i,1)  = dx2*dz*u
               uao(al+5,i,1)  = dy2*dx*u
               uao(al+6,i,1)  = dy2*dz*u
               uao(al+7,i,1)  = dz2*dx*u
               uao(al+8,i,1)  = dz2*dy*u

c              // f_xyz
               u = sqr3*u                  ! correction of norm

               uao(al+9,i,1)  = dxyz*u

               al = al + 10

            else if (bl(bf)=='F'.and.gaussFOrder) then  ! 3f GTO
c              // do all 10 F simultaneously (same exponent is required)
c              // order: f_xxx, f_yyy, f_zzz, f_xyy, f_xxy, f_xxz,
c              //   f_xzz, f_yzz, f_yyz, f_xyz  (like Gaussian)

               ispl = 3*so(bf)-2
               u    = cspla(ispl,j) + df*(csplb(ispl,j)
     .              + df*(csplc(ispl,j) + df*cspld(ispl,j)))

               dx = xx-atoms(a)%cx
               dx2 = dx*dx
               dy = yy-atoms(a)%cy
               dy2 = dy*dy
               dz = zz-atoms(a)%cz
               dz2 = dz*dz
               dxyz = dx*dy*dz

c              // f_xxx, f_yyy, f_zzz
               uao(al,i,1)    = dx2*dx*u
               uao(al+1,i,1)  = dy2*dy*u
               uao(al+2,i,1)  = dz2*dz*u

c              // f_xyy, f_xxy, f_xxz, f_xzz, f_yzz, f_yyz
               u = sqr5*u                   ! correction of norm

               uao(al+3,i,1)  = dy2*dx*u
               uao(al+4,i,1)  = dx2*dy*u
               uao(al+5,i,1)  = dx2*dz*u
               uao(al+6,i,1)  = dz2*dx*u
               uao(al+7,i,1)  = dz2*dy*u
               uao(al+8,i,1)  = dy2*dz*u

c              // f_xyz
               u = sqr3*u                  ! correction of norm

               uao(al+9,i,1)  = dxyz*u

               al = al + 10

            else
               call abortp('(getaos): wrong GTO')
            endif ! bl
            endif	! simple GTO, no splines
           endif  ! STO/GTO
         enddo  ! bf-loop over basis functions
      enddo  ! i-loop over electrons

      CONTAINS

         subroutine internal_GamessOrderGFunctions()

            alp = cntrctn(1,1,bf)
            u = cntrctn(2,1,bf) * exp(-alp*r2)
            dx = xx-atoms(a)%cx
            dx2 = dx*dx
            dy = yy-atoms(a)%cy
            dy2 = dy*dy
            dz = zz-atoms(a)%cz
            dz2 = dz*dz
            dxyz = dx*dy*dz

c           // g_xxxx, g_yyyy, g_zzzz
            uao(al,i,1)    = uao(al,i,1)    + dx2*dx2*u
            uao(al+1,i,1)  = uao(al+1,i,1)  + dy2*dy2*u
            uao(al+2,i,1)  = uao(al+2,i,1)  + dz2*dz2*u

c           // g_xxxy, g_xxxz, g_yyyx, g_yyyz, g_zzzx, g_zzzy
            u = sqr7*u                   ! correction of norm

            uao(al+3,i,1)  = uao(al+3,i,1)  + dx2*dx*dy*u
            uao(al+4,i,1)  = uao(al+4,i,1)  + dx2*dx*dz*u
            uao(al+5,i,1)  = uao(al+5,i,1)  + dy2*dy*dx*u
            uao(al+6,i,1)  = uao(al+6,i,1)  + dy2*dy*dz*u
            uao(al+7,i,1)  = uao(al+7,i,1)  + dz2*dz*dx*u
            uao(al+8,i,1)  = uao(al+8,i,1)  + dz2*dz*dy*u

c           // g_xxyy, g_xxzz, g_yyzz
            u = sqr5 / sqr3 * u          ! correction of norm

            uao(al+9,i,1)   = uao(al+9,i,1)  + dx2*dy2*u
            uao(al+10,i,1)  = uao(al+10,i,1)  + dx2*dz2*u
            uao(al+11,i,1)  = uao(al+11,i,1)  + dy2*dz2*u

c           // g_xxyz, g_yyxz, g_zzxy
            u = sqr3*u                  ! correction of norm

            uao(al+12,i,1)  = uao(al+12,i,1)  + dx*dxyz*u
            uao(al+13,i,1)  = uao(al+13,i,1)  + dy*dxyz*u
            uao(al+14,i,1)  = uao(al+14,i,1)  + dz*dxyz*u

         end subroutine internal_GamessOrderGFunctions


         subroutine internal_GaussianOrderGFunctions()

            alp = cntrctn(1,1,bf)
            u = cntrctn(2,1,bf) * exp(-alp*r2)
            dx = xx-atoms(a)%cx
            dx2 = dx*dx
            dy = yy-atoms(a)%cy
            dy2 = dy*dy
            dz = zz-atoms(a)%cz
            dz2 = dz*dz
            dxyz = dx*dy*dz

c           // g_xxxx, g_yyyy, g_zzzz
            uao(al+14,i,1)    = uao(al+14,i,1)    + dx2*dx2*u
            uao(al+4,i,1)  = uao(al+4,i,1)  + dy2*dy2*u
            uao(al,i,1)  = uao(al,i,1)  + dz2*dz2*u

c           // g_xxxy, g_xxxz, g_yyyx, g_yyyz, g_zzzx, g_zzzy
            u = sqr7*u                   ! correction of norm

            uao(al+13,i,1)  = uao(al+13,i,1)  + dx2*dx*dy*u
            uao(al+12,i,1)  = uao(al+12,i,1)  + dx2*dx*dz*u
            uao(al+8,i,1)  = uao(al+8,i,1)  + dy2*dy*dx*u
            uao(al+3,i,1)  = uao(al+3,i,1)  + dy2*dy*dz*u
            uao(al+5,i,1)  = uao(al+5,i,1)  + dz2*dz*dx*u
            uao(al+1,i,1)  = uao(al+1,i,1)  + dz2*dz*dy*u

c           // g_xxyy, g_xxzz, g_yyzz
            u = sqr5 / sqr3 * u          ! correction of norm

            uao(al+11,i,1)   = uao(al+11,i,1)  + dx2*dy2*u
            uao(al+9,i,1)  = uao(al+9,i,1)  + dx2*dz2*u
            uao(al+2,i,1)  = uao(al+2,i,1)  + dy2*dz2*u

c           // g_xxyz, g_yyxz, g_zzxy
            u = sqr3*u                  ! correction of norm

            uao(al+10,i,1)  = uao(al+10,i,1)  + dx*dxyz*u
            uao(al+7,i,1)  = uao(al+7,i,1)  + dy*dxyz*u
            uao(al+6,i,1)  = uao(al+6,i,1)  + dz*dxyz*u

         end subroutine internal_GaussianOrderGFunctions


      end subroutine ao1splcalc



      END MODULE aos
