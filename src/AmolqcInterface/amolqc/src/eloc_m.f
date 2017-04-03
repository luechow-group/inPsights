c f90 module eloc -- local energy calculation

c $Id: eloc_m.f,v 1.2 2008/02/05 22:21:35 luechow Exp $

c $Log: eloc_m.f,v $
c Revision 1.2  2008/02/05 22:21:35  luechow
c minor bugs removed after compilation with sun compiler.
c parselib.f adapted to fortran standard
c
c Revision 1.1.1.1  2007/04/25 13:42:20  luechow
c QMC program amolqc. rewritten in 2006. AL
c
c Revision 1.13  2005/09/29 12:48:02  annika
c use_ecp as logical variable instead of ecptype
c
c Revision 1.12  2005/03/11 13:33:11  luechow
c adding 'linscal' routines using sparse matrix algorithms with umfpack lib
c
c Revision 1.11  2005/02/16 12:44:55  luechow
c adaptation to new Jastrow interface
c
c Revision 1.10  2005/02/08 14:51:12  diedrch
c added subroutine cutyn; used to get last action of the 'elcutoff' routine
c
c Revision 1.9  2004/09/13 17:29:02  diedrch
c Added support for the routines in aomocut_m.f
c
c Revision 1.8  2004/08/20 14:27:07  diedrch
c adapted eloc_m.f and mos_m.f for new version of aomo_m.f
c
c Revision 1.7  2004/08/12 17:38:11  diedrch
c psit and psit_new now call ao1splcalc instead of ao1calc if splines are used
c
c Revision 1.6  2004/06/25 14:23:09  diedrch
c Now allows both methods of AO/MO evaluation (aomo_calc as well as aocalc/mocalc)
c
c Revision 1.5  2004/04/02 09:37:37  diedrch
c now calling aomo_calc and aomo1_calc instead of aocalc, mocalc, ao1calc
c and mo1calc
c
c Revision 1.4  2004/03/25 17:37:34  diedrch
c modified 'jasp' call because 'jasp' now has seperate initialisation option
c
c Revision 1.3  2004/01/06 18:11:30  diedrch
c Avoid persistend walker in subroutine sample
c
c Revision 1.2  2003/12/31 14:36:43  diedrch
c loop structure for ecp localization rearranged. ppint and
c ecpcalc combined to one subroutine ecpcalc. Data objects
c which change during ECP localisation are now stored and restored
c afterwords. This avoids the additional psit call in ppinti
c
c Revision 1.1.1.1  2003/12/28 15:08:07  diedrch
c Initial Code 281203
c
c Revision 2.3  2002/10/14 18:44:00  luechow
c final stable version from Sebastian. Oct 02
c
c Revision 2.0  1999/08/18 16:13:44  luechow
c initial f90

c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

c---------------------------------------------------------------

      MODULE elocal

      use omp_lib
      use wfdata
      use wfModule, only: readwriteWF
      use wfParamsModule
      use jastrowparamdata
      use elocaldata
      use eConfigsModule
      use jastrow, only: jasCalcWDerivs, JasCalc
      use aos, only: aocalc, aosplcalc
      use mos, only: mocalc
      use mdet, only: mdetcalc
      use mdetparam
      use moparam
      use ecpModule
      use ecpio, only: initecpparams
      use aomo
      use RdataUpdateModule

      implicit none


      type(WFType) :: mWF

      private
      public :: wf_init, wf_ecp_init, wf_getNCoreElecs, eloc, vpotcalc, elCutOff, elocUpdateDistAndPot


      CONTAINS

      ! first new routines for wf objects, currently only one, hidden from QMC
      ! stored here


      subroutine wf_init(lines,nl)
         integer, intent(in)          :: nl
         character(len=*), intent(in) :: lines(nl)

         call readwriteWF(lines,nl,mWF)
      end subroutine wf_init


      subroutine wf_ecp_init(lines,nl)
         integer, intent(in)          :: nl
         character(len=*), intent(in) :: lines(nl)
         call initecpparams(lines,nl,mWF%ecp)
      end subroutine wf_ecp_init


      function wf_getNCoreElecs(a) result(res)
         integer :: res
         integer, optional, intent(in) :: a
         if (present(a)) then
            res = mWF%ecp%NCoreElecs(atom=a)
         else
            res = mWF%ecp%NCoreElecs()
         endif
      end function wf_getNCoreElecs


c     -------------------------------------------------------------------------
      subroutine eloc(ie, ec, optType, wfpDef, wfpDT, twoLevelStep, doCalc,
     .                tmove, tau, isMoved)
c     -------------------------------------------------------------------------

c TODO:
c free format!
c remove optType (is part of wfpDef)
c check two-level propagator (requires saving intermediate results)
c remove EConfigArray! replace by EConfig




c ELOC calculates psi, grad(psi), kin. energy ekin,
c and pot. energy Vpot, for given electron configurations (x,y,z) in ec
c one-electron update for ie>0. Results are stored in module
c data structure elocdata (eloc_d.f90) for fast and flexible access.

c Modifications/Versions:
c This version allows efficient calculation of several electron
c configurations simultaneously. ec contains several configurations
c (max is mMaxElecConfigs), AL, 2011
c Version 1.5 (28.4.98) "module" version. No reference to walker, Just
c                       evaluate E_local and other stuff (stored in
c                       local data structure in eloc.h) for some position
c                       vector.
c Version 1.4 (9.3.98)  include spline evaluation of radial part of AO's
C Version 1.3 (sept 97) module version: use aos,mos,jastrow,mdet modules
c Version 1.2 (5/31/96) include ECP's
c Version 1.1 (5/20/96) allow single-electron moves; calc local energy
c                       for the orbital part only (COMMON eltest)
c Version 1.0 (1/29/96)
      integer, intent(in)                              :: ie       ! 0: update all electron; i: update only electron i
      type(eConfigArray), intent(inout)                :: ec       ! contains electron configurations: T move on exit
      character(len=*), intent(in)                     :: optType  ! transferred to call to JasWDerivs (REMOVE! REPLACED BY WFPDEF)
      type(WFParamDef), intent(in), optional           :: wfpDef   ! definition of wf parameters; triggers calculation of
                                                                   !                                    param derivatives
      type(WFParamDerivTerms), intent(inout), optional :: wfpDT    ! contains on output parameter derivatives (of Psi and Elocal)
      integer, intent(in), optional                    :: twoLevelStep     ! two Level propagator step
      logical, intent(in), optional                    :: doCalc(eConfigArray_size(ec)) ! for twoLevel
      integer, intent(in), optional                    :: tmove    ! triggers T move calculation
      real*8, intent(in), optional                     :: tau      ! time step for T moves
      logical, intent(inout), optional                 :: isMoved(:) ! T move accepted? for each EConfig

      integer i,ii,a,itmp,nd,idx,iull,n,nElecConfigs,ierr
      real*8 ekin,tmp,vpot
      real*8 wtimer_i(8),wtimer_f(8),wtimerPhi(4)
      ! automatic arrays for actual electron number!
      real*8 x(ne),y(ne),z(ne)
      real*8 rai(ncenter,ne),rij(ne,ne)
      real*8 rrai(ncenter,ne,eConfigArray_size(ec))    ! rai for all elec configs
      real*8 rrij(ne,ne,eConfigArray_size(ec))         ! rij for all elec configs
      real*8 u,ugrad(3*ne),ulapl,ulapli(ne)
      real*8 uu(eConfigArray_size(ec))
      real*8 uugrad(3*ne,eConfigArray_size(ec))
      real*8 uulapl(eConfigArray_size(ec))
      real*8 uulapli(ne,eConfigArray_size(ec))
      real*8 phi(eConfigArray_size(ec)),
     .       fgrad(3*ne,eConfigArray_size(ec)),
     .       flapl(eConfigArray_size(ec)),
     .       flapli(ne,eConfigArray_size(ec))
      real*8 ecpLocal, ecpNonlocal
      real*8, allocatable :: ecpNLocalk(:)
      type(TMoveDataType) :: TMoveData


#ifdef WTIMER
      if (MASTER) wtimer_i(WTIMER_ELOC) = omp_get_wtime()
#endif

      nElecConfigs = eConfigArray_size(ec)

      call assert(nElecConfigs==1,"eloc: many EConfigs disabled")

      iull = mytid + 80  ! for write output for each node

      rrai = 0; rrij = 0
      rai = 0; rij = 0


      if (present(twoLevelStep)) then
         call elocUpdateDistAndPot(ec,ie,nElecConfigs,rai,rij,rrai,rrij)
         if (twoLevelStep == 1) then

#ifdef WTIMER
            if (MASTER) wtimer_i(WTIMER_PHI) = omp_get_wtime()
#endif

            call resetEloc()
            call elocCalcPhi(ie,nElecconfigs,rrai,ec,wtimerPhi)

#ifdef WTIMER
            if (MASTER) wtimer_f(WTIMER_PHI) = omp_get_wtime()
#endif

         else

#ifdef WTIMER
            if (MASTER) wtimer_i(WTIMER_JAS) = omp_get_wtime()
#endif

            call elocCalcJas(ie,ec,optType,nElecConfigs,uu,uugrad,uulapl,uulapli,rrai,rrij,twoLevelStep,doCalc)

#ifdef WTIMER
            if (MASTER) wtimer_f(WTIMER_JAS) = omp_get_wtime()
#endif
            if (twoLevelStep == 2) return

         endif

      else

         call resetEloc()
         call elocUpdateDistAndPot(ec,ie,nElecConfigs,rai,rij,rrai,rrij)

#ifdef WTIMER
         if (MASTER) wtimer_i(WTIMER_PHI) = omp_get_wtime()
#endif

         call elocCalcPhi(ie,nElecconfigs,rrai,ec,wtimerPhi)

#ifdef WTIMER
         if (MASTER) wtimer_f(WTIMER_PHI) = omp_get_wtime()
         if (MASTER) wtimer_i(WTIMER_JAS) = omp_get_wtime()
#endif

         call elocCalcJas(ie,ec,optType,nElecConfigs,uu,uugrad,uulapl,uulapli,rrai,rrij)

#ifdef WTIMER
         if (MASTER) wtimer_f(WTIMER_JAS) = omp_get_wtime()
#endif
      endif

      if (mWF%ecp%isInitialised()) then

#ifdef WTIMER
         if (MASTER) wtimer_i(WTIMER_PP) = omp_get_wtime()
#endif
         if (present(wfpDef) .and. present(wfpDT)) then
            if (wfpDT%ELiCalc .and. (wfpDT%noCalc.eqv. .false.) )then
               call elocCalcECP(ie, ec, nElecConfigs, rrai, rrij, ecpLocal,
     .                       ecpNonlocal, ecpNLocalk=ecpNLocalk, wfpDef=wfpDef)
            else
               call elocCalcECP(ie, ec, nElecConfigs, rrai, rrij, ecpLocal,
     .                       ecpNonlocal)

            end if
         else if (present(tmove)) then
            TMoveData%tmove = tmove
            TMoveData%tau = tau
            call elocCalcECP(ie, ec, nElecConfigs, rrai, rrij, ecpLocal,
     .                       ecpNonlocal, TMoveData=TMoveData)
            isMoved(1) = TMoveData%isMoved
         else
            call elocCalcECP(ie,ec,nElecConfigs,rrai,rrij,ecpLocal,ecpNonlocal)
         endif

         elVPot(1) = elVPot(1) + ecpLocal + ecpNonlocal
         elECPPot = ecpLocal
         elECPPotNl = ecpNonlocal

#ifdef WTIMER
         if (MASTER) wtimer_f(WTIMER_PP) = omp_get_wtime()
#endif
      endif

      n = 1   !!! only one EConfig !!!

      if (present(doCalc)) then
         if (.not. doCalc(1)) return
      endif

c     // Grad(Psi_G) / psi, Laplacian and kinetic energy
      phi(n) = elPhi(n)
      fgrad(:,n) = elFgrad(:,n)
      flapl(n) = elFlapl(n)
      flapli(:,n) = elFlapli(:,n)
      ii = 1
      do i=1,ne
         elxDrift(i,n) = fgrad(ii,n)/phi(n) + uugrad(ii,n)
         ii = ii+1
         elyDrift(i,n) = fgrad(ii,n)/phi(n) + uugrad(ii,n)
         ii = ii+1
         elzDrift(i,n) = fgrad(ii,n)/phi(n) + uugrad(ii,n)
         ii = ii+1
      enddo
c     // kinetic energy for Psi_G
      ii = 1
      do i=1,ne
         elEkini(i,n) = -0.5d0*( flapli(i,n)/phi(n) + uulapli(i,n)
     .        + 2d0*fgrad(ii,n)/phi(n) * uugrad(ii,n)
     .        + 2d0*fgrad(ii+1,n)/phi(n) * uugrad(ii+1,n)
     .        + 2d0*fgrad(ii+2,n)/phi(n) * uugrad(ii+2,n)
     .        + uugrad(ii,n) * uugrad(ii,n)
     .        + uugrad(ii+1,n) * uugrad(ii+1,n)
     .        + uugrad(ii+2,n) * uugrad(ii+2,n) )
         ii = ii+3
      enddo

      if(do_epart) elEkin(1:ne,n) = elEkini(1:ne,n)

      ekin = 0d0                ! local kinetic energy with Jastrow
      do i=1,ne
         ekin = ekin + elEkini(i,n)
      enddo

c      // local energies
      elEloc(n) = ekin + elVPot(n) + vpot0

      if (logmode >= 5) then
         write(iul,'(A,5G20.10)') "eloc:",elEloc(n),ekin,elVPot(n)+vpot0,elECPPot,elECPPotNl
         write(iul,'(A,2G20.10)') "phi,U:", elPhi(n),uu
      endif

      if (mElocCutOff) then
         call elCutOff(n)
         if (logmode >= 5) write(iul,'(A,4G15.6)') "eloc (after cutoff):",elEloc(n),ekin,elVPot(n)+vpot0,elECPPot
      endif

      if (present(wfpDef) .and. present(wfpDT)) then
           wfpDT%eloc = elEloc(1)
           wfpDT%phi  = phi(1)
           wfpDT%U    = uu(1)
         ! calculate contributions to wf parameter derivative terms
         if(wfpDT%noCalc .eqv. .false.) call internal_calcDerivContribs()
      endif



#ifdef WTIMER
      if (MASTER) then
         wtimer_f(WTIMER_ELOC) = omp_get_wtime()
         wtimer(WTIMER_ELOC) = wtimer(WTIMER_ELOC) + (wtimer_f(WTIMER_ELOC)-wtimer_i(WTIMER_ELOC))
         wtimer(WTIMER_PHI) = wtimer(WTIMER_PHI) + (wtimer_f(WTIMER_PHI)-wtimer_i(WTIMER_PHI))
         wtimer(WTIMER_PP) = wtimer(WTIMER_PP) + (wtimer_f(WTIMER_PP)-wtimer_i(WTIMER_PP))
         wtimer(WTIMER_JAS) = wtimer(WTIMER_JAS) + (wtimer_f(WTIMER_JAS)-wtimer_i(WTIMER_JAS))
         wtimer(WTIMER_AO:WTIMER_MDET) = wtimer(WTIMER_AO:WTIMER_MDET) + wtimerPhi
      endif
#endif


      contains

         subroutine internal_calcDerivContribs()

         integer np, npJ, npMO, npCI, k, l, i
         real*8 tmp

         call wfparams_getNJastrowParams(wfpDef,npJ)
         call wfparams_getNDetParams(wfpDef,npCI,npMO)
         np   = npJ + npCI + npMO

         if (.not.allocated(wfpDT%Eli)) then
            allocate(wfpDT%fi(np), wfpDT%Eli(np), wfpDT%fij(np,np))
         else
            if (size(wfpDT%Eli) /= np) call abortp("eloc:internal_calcDerivContribs: illegal sizes")
         endif

C          wfpDT%eloc = elEloc(1)
C          wfpDT%phi  = phi(1)
C          wfpDT%U    = uu(1)

         ! parameter order: Jastrow, MO, CI

         ! Jastrow block
         if (npJ > 0) then

            if (wfpDT%fiCalc) then
               do k = 1, npJ
                  wfpDT%fi(k) = uk(k)
               end do
            end if

            if (wfpDT%ELiCalc) then
               do k = 1, npJ
                  tmp = 0.d0
                  do i = 1, ne
                     tmp = tmp + elxDrift(i,1)*ukgrad(3*i-2,k) + elyDrift(i,1)*ukgrad(3*i-1,k)
     &                         + elzDrift(i,1)*ukgrad(3*i,k)
                  end do
                  wfpDT%ELi(k) = -0.5d0*uklapl(k) - tmp
               end do
               if (mWF%ecp%isInitialised()) then
                  wfpDT%ELi(1:npJ) = wfpDT%ELi(1:npJ)
     &                             + ecpNLocalk(1:npJ) - ecpNonlocal*uk(1:npJ)
               end if
            end if

            if (wfpDT%fijCalc) then
               do k = 1, npJ
                  do l = 1, k
                     wfpDT%fij(k,l) = uk(l) * uk(k)
                  end do
               end do
            end if

         endif

         ! MO block
         if (npMO > 0) then
            call moparam_calcderivs(wfpDT%ELiCalc)

            if (wfpDT%fiCalc) then
               do i = 1, npMO
                  wfpDT%fi(npJ+i) = mok(i)/elPhi(1)
               end do
            end if

            if (wfpDT%ELiCalc) then
               do i = 1, npMO
                  tmp = moklapl(i) + 2*dot_product(mokgrad(:,i),elUGrad(1:3*ne,1))
     &                 + mok(i)*dot_product(elUGrad(1:3*ne,1),elUGrad(1:3*ne,1)) + mok(i)*elULapl(1)
                  tmp = tmp * (-0.5d0) / elPhi(1)
                  wfpDT%Eli(npJ+i) = tmp + (elVpot(1)+vpot0 - elEloc(1))*mok(i)/elPhi(1)
               end do
               if (mWF%ecp%isInitialised()) then
                  wfpDT%ELi(npJ+1:npJ+npMO) = wfpDT%ELi(npJ+1:npJ+npMO) +
     &              ecpNLocalk(npJ+1:npJ+npMO) - ecpNonlocal*mok(1:npMO)/elPhi(1)
               end if
            end if

            if (wfpDT%fijCalc) then
               wfpDT%fij(npJ+1:npJ+npMO,npJ+1:npJ+npMO) = 0.d0
            end if

         endif

         ! CI block
         if (npCI > 0) then
            call mdetcalcderivs()

            if (wfpDT%fiCalc) then
               do i = 1, npCI
                  wfpDT%fi(npMO+npJ+i) = fk(i)/elPhi(1)
               end do
            end if

            if (wfpDT%ELiCalc) then
               do i = 1, npCI
                  tmp = fklapl(i) + 2*dot_product(fkgrad(:,i),elUGrad(1:3*ne,1))
     &                 + fk(i)*dot_product(elUGrad(1:3*ne,1),elUGrad(1:3*ne,1)) + fk(i)*elULapl(1)
                  tmp = tmp * (-0.5d0) / elPhi(1)
                  wfpDT%ELi(npMO+npJ+i) = tmp + (elVpot(1)+vpot0 - elEloc(1))*fk(i)/elPhi(1)
               end do
               if (mWF%ecp%isInitialised()) then
                  wfpDT%ELi(npJ+npMO+1:np) = wfpDT%ELi(npJ+npMO+1:np) +
     &               ecpNLocalk(npJ+npMO+1:np) - ecpNonlocal*fk(1:npCI)/elPhi(1)
               end if
            end if

            if (wfpDT%fijCalc) then
               wfpDT%fij(npMO+npJ+1:np,npMO+npJ+1:np) = 0.d0
            end if
         end if


         ! mixed blocks
         if (wfpDT%fijCalc) then
            ! jas-MO
            do k = 1, npJ
               do i = 1, npMO
                  wfpDT%fij(npJ+i,k)  =  mok(i)/elPhi(1) * uk(k)
               end do
            end do
            ! jas-CI
            do k = 1, npJ
               do i = 1, npCI
                  wfpDT%fij(npJ+npMO+i,k)  =  fk(i)/elPhi(1) * uk(k)
               end do
            end do
            ! MO-CI
            do k = 1, npMO
               do i = 1, npCI
                  wfpDT%fij(npJ+npMO+i,npJ+k)  =  fk(i)/elPhi(1) * mok(k)
               end do
            end do
         end if

         end subroutine internal_calcDerivContribs
      end subroutine eloc


c===========================================================

c     --------------------------!
      subroutine resetEloc()
c     --------------------------!
        elPhi = 0
        elFgrad = 0
        elFlapl = 0
        elU = 0
        elUgrad = 0
        elUlapl = 0
        elxDrift = 0; elyDrift = 0; elzDrift = 0
        elECPPot = 0
        elECPPotNl=0
        elVpot = 0
        elEkini = 0
        elEloc = 0
        elSingularity = .false.
      end subroutine resetEloc

c===========================================================

c     -----------------------------------------!
      subroutine vpotcalc(ie,x,y,z,rai,rij,vpot,n)
c     -----------------------------------------!

c vpotcalc calculates the potential energy for given position vector
c NOTE: rij is calculated ONLY FOR i<j
c ie == 0: new calculation of vpot,rai and rij
c ie == n: update only for n-th electron
c on entry: for one-electron update the vectors x,y,z must be changed compared
c           to the last call to this routine only at electron ie.
c           rai and rij contain the old distances.
c           vpot is then required have 'old' value of vpot

      integer, intent(in),optional :: n
      integer, intent(in)   :: ie    ! ie==0: recalculate rai,rij,vpot,
                                     !        else update electron ie
      real*8, intent(in)    ::   x(:),y(:),z(:)
      real*8, intent(inout) ::   rai(:,:),rij(:,:)
      real*8, intent(inout) ::   vpot  ! input: for ie>0 old potential
                                       ! output: vpot for x,y,z
      integer a,i,j
      real*8 vpoti(nmax)             ! individual potential energy
      save vpoti

c     // Potential energy vpot and distances
      if (do_epart) then
      if (ie == 0) then
c        // Recalculate all distances and potential energy
         vpot = 0.d0
         do i=1,ne
            vpoti(i) = 0.d0
            do a=1,ncenter
               rai(a,i) = sqrt((x(i)-atoms(a)%cx)**2
     .              + (y(i)-atoms(a)%cy)**2 + (z(i)-atoms(a)%cz)**2)
               vpoti(i) = vpoti(i) - atoms(a)%za/rai(a,i)
               elVne(a,i,n) = -atoms(a)%za/rai(a,i)
            enddo
            vpot = vpot + vpoti(i)
            do j=i+1,ne
               rij(i,j) = sqrt((x(i)-x(j))**2 +
     .              (y(i)-y(j))**2 + (z(i)-z(j))**2)
               vpot = vpot + 1.d0/rij(i,j)
               elVee(i,j,n) = 1.d0/rij(i,j)
            enddo
         enddo
      else
c        // Update distances and potential energy for electron ie
         i = ie
         vpot = vpot - vpoti(i)
         vpoti(i) = 0.d0
         do a=1,ncenter
            rai(a,i) = sqrt((x(i)-atoms(a)%cx)**2
     .           + (y(i)-atoms(a)%cy)**2 + (z(i)-atoms(a)%cz)**2)
            vpoti(i) = vpoti(i) - atoms(a)%za/rai(a,i)
            elVne(a,i,n) = -atoms(a)%za/rai(a,i)
         enddo
         vpot = vpot + vpoti(i)
         do j=1,i-1
            vpot = vpot - 1d0/rij(j,i)
            rij(j,i) = sqrt((x(i)-x(j))**2 +
     .           (y(i)-y(j))**2 + (z(i)-z(j))**2)
            vpot = vpot + 1d0/rij(j,i)
            elVee(j,i,n) = 1d0/rij(j,i) !is this correct?
         enddo
         do j=i+1,ne
            vpot = vpot - 1d0/rij(i,j)
            rij(i,j) = sqrt((x(i)-x(j))**2 +
     .           (y(i)-y(j))**2 + (z(i)-z(j))**2)
            vpot = vpot + 1d0/rij(i,j)
            elVee(i,j,n) = 1d0/rij(i,j) !is this correct?
         enddo
      endif !! (ie==0)
      else  !! (.not.do_epart)
      if (ie == 0) then
c        // Recalculate all distances and potential energy
         vpot = 0.d0
         do i=1,ne
            vpoti(i) = 0.d0
            do a=1,ncenter
               rai(a,i) = sqrt((x(i)-atoms(a)%cx)**2
     .              + (y(i)-atoms(a)%cy)**2 + (z(i)-atoms(a)%cz)**2)
               vpoti(i) = vpoti(i) - atoms(a)%za/rai(a,i)
            enddo
            vpot = vpot + vpoti(i)
            do j=i+1,ne
               rij(i,j) = sqrt((x(i)-x(j))**2 +
     .              (y(i)-y(j))**2 + (z(i)-z(j))**2)
               vpot = vpot + 1.d0/rij(i,j)
            enddo
         enddo
      else
c        // Update distances and potential energy for electron ie
         i = ie
         vpot = vpot - vpoti(i)
         vpoti(i) = 0.d0
         do a=1,ncenter
            rai(a,i) = sqrt((x(i)-atoms(a)%cx)**2
     .           + (y(i)-atoms(a)%cy)**2 + (z(i)-atoms(a)%cz)**2)
            vpoti(i) = vpoti(i) - atoms(a)%za/rai(a,i)
         enddo
         vpot = vpot + vpoti(i)
         do j=1,i-1
            vpot = vpot - 1d0/rij(j,i)
            rij(j,i) = sqrt((x(i)-x(j))**2 +
     .           (y(i)-y(j))**2 + (z(i)-z(j))**2)
            vpot = vpot + 1d0/rij(j,i)
         enddo
         do j=i+1,ne
            vpot = vpot - 1d0/rij(i,j)
            rij(i,j) = sqrt((x(i)-x(j))**2 +
     .           (y(i)-y(j))**2 + (z(i)-z(j))**2)
            vpot = vpot + 1d0/rij(i,j)
         enddo
      endif !! (ie==0)
      endif !! (do_epart)


      end subroutine vpotcalc



c     ------------------------
      subroutine elCutOff(n)
c     ------------------------

c  Cutoff of Rothstein/Vrbik (see Umrigar 1993 paper)
c  cutoff for both local energy and drift. Cutoff vanishes
c  for tau->0. Reference for local energy cutoff is 'evar'.
c  Careful: 'evar' must have an appropriate value
c  Alternatively: use best current estimate for <E>
c  Now: CUTOFF for E_local(Psi_T) only, not for Psi_G!
c  but: for Drift(Psi_G)

      integer, intent(in) :: n    ! electron configuration
      integer i
      real*8 frac
      logical driftCut

      frac = mCutOffFactor/sqrt(mTauCutOff)

      if (abs(elEloc(n)-E_trial) > 2d0*frac) then
         if (elEloc(n) >= E_trial) then
            elEloc = E_trial + 2d0*frac
         else
            elEloc = E_trial - 2d0*frac
         endif
         mElocCut = mElocCut + 1
      endif
      mElocCutCount = mElocCutCount + 1

      frac = mCutOffFactor/mTauCutOff
      driftCut = .false.
      do i=1,ne
         if (abs(elxDrift(i,n)) > frac) then
            driftCut = .true.
            if (elxDrift(i,n) >= 0d0) then
               elxDrift(i,n) = frac
            else
               elxDrift(i,n) = -frac
            endif
         endif
         if (abs(elyDrift(i,n)) > frac) then
            driftCut = .true.
            if (elyDrift(i,n) >= 0d0) then
               elyDrift(i,n) = frac
            else
               elyDrift(i,n) = -frac
            endif
         endif
         if (abs(elzDrift(i,n)) > frac) then
            driftCut = .true.
            if (elzDrift(i,n) >= 0d0) then
               elzDrift(i,n) = frac
            else
               elzDrift(i,n) = -frac
            endif
         endif
      enddo
      if (driftCut) then
         mDriftCut = mDriftCut + 1
      end if
      mDriftCutCount = mDriftCutCount + 1

      end subroutine elCutOff


c     ------------------------
      subroutine elocUpdateDistAndPot(ec,ie,nElecConfigs,rai,rij,rrai,rrij,doCalc)
c     ------------------------
      integer, intent(in)                :: ie       ! 0: update all electron; i: update only electron i
      type(eConfigArray), intent(inout)  :: ec  ! contains electron configurations
      real*8 x(ne),y(ne),z(ne),vpot
      integer :: n
      integer,intent(in) :: nElecConfigs
      real*8,intent(inout) :: rai(ncenter,ne),rij(ne,ne)
      real*8,intent(inout) :: rrai(ncenter,ne,eConfigArray_size(ec))    ! rai for all elec configs
      real*8,intent(inout) :: rrij(ne,ne,eConfigArray_size(ec))         ! rij for all elec configs
      logical,intent(in),optional :: doCalc(eConfigArray_size(ec))

      do n=1,nElecConfigs
         if(present(doCalc)) then
          if(.not. doCalc(n)) cycle
         endif
         call eConfigArray_get(ec,n,x,y,z)
         ! check input data for NaN or Inf
         call assert(all(abs(x)<huge(1.d0)) .and.
     .     all(abs(y)<huge(1.d0)) .and. all(abs(z)<huge(1.d0)),
     .     "eloc: illegal x,y,z coords in ec")
c        // calculate/update particle distances and potential energy
         call vpotcalc(0,x,y,z,rai,rij,vpot,n)
         call assert(all(rai<huge(1.d0)),"eloc: illegal rai values")
         call assert(all(rij<huge(1.d0)),"eloc: illegal rij values")
         elVpot(n) = vpot
         rrai(:,:,n) = rai
         rrij(:,:,n) = rij
       end do
       return
       end subroutine elocUpdateDistAndPot


c     ----------------------------------------------------------
       subroutine elocCalcPhi(ie,nElecconfigs,rrai,ec,wtimerPhi)
c     ----------------------------------------------------------
      integer, intent(in)                :: ie       ! 0: update all electron; i: update only electron i
      integer,intent(in) :: nElecConfigs
      real*8,intent(in) :: rrai(:,:,:)    ! rai for all elec configs
      real*8,intent(out):: wtimerPhi(4)
      type(eConfigArray), intent(inout)  :: ec  ! contains electron configurations
      real*8 x(ne),y(ne),z(ne)
      integer :: n,ierr,i
      real*8 phi(eConfigArray_size(ec)),
     .       fgrad(3*ne,eConfigArray_size(ec)),
     .       flapl(eConfigArray_size(ec)),
     .       flapli(ne,eConfigArray_size(ec))
      real*8 wtimer1,wtimer2,wtimer3,wtimer4

      wtimerPhi = 0.d0

      if (aomocomb) then ! use combined AO/MO calculation

         ECLOOP: do n=1,nElecConfigs

         call eConfigArray_get(ec,n,x,y,z)
         !rai = rrai(:,:,n)
#ifdef WTIMER
         if (MASTER) wtimer1 = omp_get_wtime()
#endif
         if (spline) then
            if (cutmo) then
               call aomocutspl_calc(ie,x,y,z,rrai(:,:,n))
            else
               call aomospl_calc(ie,x,y,z,rrai(:,:,n))
            endif
         else
            if (cutmo) then
               call aomocut_calc(ie,x,y,z,rrai(:,:,n))
            else
               call aomo_calc(ie,x,y,z,rrai(:,:,n))
            endif
         endif
#ifdef WTIMER
         if (MASTER) wtimer2 = omp_get_wtime()
#endif
         phi(n) = 0
         fgrad(:,n) = 0
         flapli(:,n) = 0
         flapl(n) = 0
         call mdetcalc(ie,1,phi(n:n),fgrad(:,n:n),flapli(:,n:n),
     .                 flapl(n:n),ierr)
#ifdef WTIMER
         if (MASTER) wtimer3 = omp_get_wtime()
#endif
         if (ierr > 0) then
            if (mStopAtSingularity) then
               write(iull,*) 'mdetcalc failed for:'
               do i=1,ne
                  write(iull,'(i5,3f14.7)') i,x(i),y(i),z(i)
               end do
               call abortp("eloc: singularity in determinant")
            else
               elSingularity = .true.
            end if
         end if
         elPhi(n) = phi(n)
         elFgrad(:,n) = fgrad(:,n)
         elFlapl(n) = flapl(n)
         elFlapli(:,n) = flapli(:,n)

#ifdef WTIMER
         if (MASTER) then
            wtimerPhi(3) = wtimerPhi(3) + (wtimer2 - wtimer1)
            wtimerPhi(4) = wtimerPhi(4) + (wtimer3 - wtimer2)
         endif
#endif

         enddo ECLOOP

      else

#ifdef WTIMER
        if (MASTER) wtimer1 = omp_get_wtime()
#endif
        if (spline) then
          call aosplcalc(ie,ec,rrai)
        else
          call aocalc(ie,ec,rrai)
        endif
c       // calculate MO's for current AO's
#ifdef WTIMER
        if (MASTER) wtimer2 = omp_get_wtime()
#endif
        call mocalc(ie)
#ifdef WTIMER
        if (MASTER) wtimer3 = omp_get_wtime()
#endif

        phi = 0
        fgrad = 0
        flapli = 0
        flapl = 0
        call mdetcalc(ie,nElecconfigs,phi,fgrad,flapli,flapl,ierr)
#ifdef WTIMER
        if (MASTER) wtimer4 = omp_get_wtime()
#endif
        if (ierr > 0) then
            if (mStopAtSingularity) then
               write(iull,*) 'mdetcalc failed for:'
               do i=1,ne
                  write(iull,'(i5,3f14.7)') i,x(i),y(i),z(i)
               end do
               call abortp("eloc: singularity in determinant")
            else
               elSingularity = .true.
            end if
        end if

        do n = 1, nElecConfigs
          elPhi(n) = phi(n)
          elFgrad(:,n) = fgrad(:,n)
          elFlapl(n) = flapl(n)
          elFlapli(:,n) = flapli(:,n)
        enddo

#ifdef WTIMER
        if (MASTER) then
           wtimerPhi(1) = (wtimer2 - wtimer1)
           wtimerPhi(2) = (wtimer3 - wtimer2)
           wtimerPhi(4) = (wtimer4 - wtimer3)
        endif
#endif

      endif
      end subroutine elocCalcPhi


c     ----------------------------------------------------------------------------------------------------------
      subroutine elocCalcJas(ie,ec,optType,nElecConfigs,uu,uugrad,uulapl,uulapli,rrai,rrij,twoLevelStep,calcJas)
c     ----------------------------------------------------------------------------------------------------------
      integer, intent(in)                :: ie,nElecConfigs       ! 0: update all electron; i: update only electron i
      type(eConfigArray), intent(inout)  :: ec  ! contains electron configurations
      character(len=*), intent(in)       :: optType ! transferred to call to JasWDerivs

      integer :: n
      ! automatic arrays for actual electron number!
      real*8             :: x(ne),y(ne),z(ne)
      real*8,intent(in)  :: rrai(ncenter,ne,eConfigArray_size(ec))
      real*8,intent(in)  :: rrij(ne,ne,eConfigArray_size(ec))
      real*8             :: u,ugrad(3*ne),ulapl,ulapli(ne)
      real*8,intent(out) :: uu(eConfigArray_size(ec))
      real*8,intent(out) :: uugrad(3*ne,eConfigArray_size(ec))
      real*8,intent(out) :: uulapl(eConfigArray_size(ec))
      real*8,intent(out) :: uulapli(ne,eConfigArray_size(ec))
      integer,intent(in),optional         :: twoLevelStep
      logical,intent(in),optional         :: calcJas(eConfigArray_size(ec))
      integer :: step
      logical :: calc

      step = 3
      if(present(twoLevelStep)) step = twoLevelStep

      do n = 1, nElecConfigs
      if (step == 2 .and. jastype /= 'none') then
        ! calculate only jastrow value without derivatives
        u = 0
        ugrad = 0
        ulapl = 0
        ulapli = 0
        call jasCalc(.false., ie, rrai(:,:,n), rrij(:,:,n), u)
      else if (step == 3 .and. jastype /= 'none') then
        ! calculate jastrow
        call eConfigArray_get(ec,n,x,y,z)
        u = 0
        ugrad = 0
        ulapl = 0
        ulapli = 0
        calc = .true.
        if(present(calcJas)) calc = calcJas(n)

        if(calc) then
          call jasCalcWDerivs(ie,x,y,z,rrai(:,:,n),rrij(:,:,n),optType,u,ugrad,ulapl,ulapli,n)
        endif
      else
        ! assign default values
        u = 0
        ugrad = 0
        ulapl = 0
        ulapli = 0
      endif
      uu(n) = u
      uugrad(:,n) = ugrad
      uulapl(n) = ulapl
      uulapli(:,n) = ulapli
      elU(n) = u
      elUGrad(:,n) = ugrad(:)
      elULapl(n) = ulapl
      enddo ! nElecConfigs

      end subroutine elocCalcJas

c     -----------------------------------------------------------------------------------------------------
      subroutine elocCalcECP(ie, ec, nElecConfigs, rrai, rrij, ecpLocal,
     .                       ecpNonlocal, ecpNLocalk, wfpDef, TMoveData)
c     -----------------------------------------------------------------------------------------------------
         integer, intent(in)                :: ie            ! ie==0: all electron, ie>0: update ie only
         integer, intent(in)                :: nElecConfigs
         type(eConfigArray), intent(inout)  :: ec            ! contains electron configurations
         real*8, intent(in)                 :: rrai(:,:,:)   ! nuc-elec distances for all EConfigs
         real*8, intent(in)                 :: rrij(:,:,:)   ! elec-elec distances for all EConfigs
         real*8, intent(out)                :: ecpLocal      ! local ECP contribution
         real*8, intent(out)                :: ecpNonlocal   ! localized nonlocal ECP contribution
         real*8, intent(out), allocatable, optional      :: ecpNLocalk(:) ! parameter derivative of nonlocal contrib
         type(WFParamDef), intent(in), optional :: wfpDef    ! wfParamDefinition with info about Jastrow terms
         type(TMoveDataType), intent(inout), optional :: TMoveData ! triggers T move calculation

         real*8             :: x(size(rrij,1)),y(size(rrij,1)),z(size(rrij,1))
         integer            :: n,np,npJ,npJ1,npJ2,npJnl,npCI,npMO, i
         type(RdataUpdate)  :: Rdu

         !
         ! check if  arg 'ie' makes any sense. One-electron moves? Really implemented?
         !

         call assert(nElecConfigs==1,"elocCalcECP: fast ECP Jastrow currently only for walker_block==1")
         ! careful: ec used for aos access for anisotropic Jastrow
         ! also: initialization of Jastrow terms for ECP updates relying on walker_block==1
         !
         ecpLocal = 0.d0
         ecpNonlocal = 0.d0
         call eConfigArray_get(ec,1,x,y,z)

         !!!write(iul,'(a)') 'DBG:elocCalcECP:x,y,z:'
         !!!do i=1,ne
         !!!   write(iul,'(i3,3g20.10)') i,x(i),y(i),z(i)
         !!!enddo

         if (present(ecpNLocalk) .and. present(wfpDef)) then
            ecpNLocalk = 0.d0
            call wfparams_getNJastrowParams(wfpDef,npJ,npJ1,npJ2,npJnl)
            call wfparams_getNDetParams(wfpDef,npCI,npMO)
            np = npJ + npCI + npMO
            call Rdu%initWParamDerivs(x,y,z,npJ1,npJ2,npCI,npMO,rrai(:,:,1),rrij(:,:,1))  ! check if new (i.e. allocation) here OK
            Rdu%phi0 = elPhi(1)
            if (.not.allocated(ecpNLocalk)) then
               allocate(ecpNLocalk(np))
            else
               if (size(ecpNLocalk) /=np) call abortp("elocCalcECP: inconsistent sizes")
            end if
            call mWF%ecp%calculate(Rdu, ecpLocal, ecpNonlocal, ecpNLocalk)
         else if (present(TMoveData)) then
            ! possible T move
            call Rdu%init(x,y,z,rrai(:,:,1),rrij(:,:,1))                   !!! check if new (i.e. allocation) here OK
            Rdu%phi0 = elPhi(1)
            call mWF%ecp%calculate(Rdu, ecpLocal, ecpNonlocal, TMoveData=TMoveData)
            if (TMoveData%isMoved) then
               call eConfigArray_set(ec, 1, Rdu%x, Rdu%y, Rdu%z)
            end if
         else
            call Rdu%init(x,y,z,rrai(:,:,1),rrij(:,:,1))                   !!! check if new (i.e. allocation) here OK
            Rdu%phi0 = elPhi(1)
            call mWF%ecp%calculate(Rdu, ecpLocal, ecpNonlocal)
         end if

         call Rdu%delete()


      end subroutine elocCalcECP

      END MODULE elocal
