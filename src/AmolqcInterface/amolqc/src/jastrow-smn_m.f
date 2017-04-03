c
c module jastrowSM for calculating the Jastrow part
c of the wavefunction using the Schmidt/Moskowitz form
c
c $Id: jastrow-smn_m.f,v 1.1.1.1 2007/04/25 13:42:20 luechow Exp $
c
c $Log: jastrow-smn_m.f,v $
c Revision 1.1.1.1  2007/04/25 13:42:20  luechow
c QMC program amolqc. rewritten in 2006. AL
c
c Revision 1.8  2005/02/16 12:42:49  luechow
c new generic Jastrow module allowing to select specific Jastrow implementations
c using new parameter 'jastype'
c
c
c SM 04.11.1999
c
c using scaled distance  [a*r/(1-a*r)]
c                                                   power
c "term types" :              type  opt  1-nc L-nc  n k l
c                       cusp    1    0     0    0   0 0 0
c                       ee      1    0     0    0   2 0 0
c                       en      2    0     1    4   0 2 0
c                       een 1   3    0     1    4   0 1 1
c                       een 2   4    0     1    4   1 1 0

c (en=electron-nucleus, enn= electron-electron-nucleus,
c  1-nc=first nucleus, L-nc=last nucleus)
c
c note: nuclei that shall use the same en oder enn term
c have to be entered sequentially in the geometry part.
c I.e. The geometry listing MUST be ordered according to
c the atoms.
c
c opt: type of allowed parameter-range during the optimization

c JASSM calculates the Jastrow part and its derivatives for the actual
c psip ps. The Jastrow exponent used is identical to Schmidt/Moskowitz
c (JCP 93, 4172 (1990)).


c     ----------------
      MODULE jastrowSM
c     ----------------

      use atomModule
      use wfdata
      use jastrowparamdata
      use RdataUpdateModule
      use Utils, only: tokenize
      implicit none

      private
      public :: jasinput_sm, jasinput_sm_new, jasoutput_sm, jas_shortoutput_sm,
     .          jas_shortoutput_sm_new,
     .          jasoutput_sm_new, jasChangeType_sm, jassmall, jassmallOnlyUk,
     .          jassmp, getVectorLenSM,getVectorSM, putVectorSM, jas_diffeecusp_sm,
     .          jassmInit, jassmUpdate, jassmInitWithUk, jassmUpdateWithUk

      integer jmode             ! jmode only used in getvec putvec
                                ! jmode = 0 ,only bpg opt.
                                ! jmode = 1 ,+    apg

      integer kpmax             ! max power of k terms
      parameter (kpmax=4)

      integer lmax              ! real max power el-nucl
      integer mmax              ! real max power el-el
      integer nclast            ! last nucl. with en- and een-Terms
                                ! normally only H-atoms follows

      integer k0                ! no. of lin. Jastrow terms
      integer tcontrol(kmax,10)   ! term-control-unit,
                                ! first column:  term-type,
                                ! 2nd col.: optimization parameter
                                ! 3-4   column:  nucleus term belongs to
                                ! 5-10  column:  term-parameters
      real*8 :: apg(amax) = 1.d0   ! el-nucl,  Jas. parameter
      real*8 :: bpg = 1.d0         ! el-el. Jas. parameter
      real*8 :: cjas(kmax) = 0.d0  ! lin. Jastrow coeff.
      real*8  lincut1,lincut2	  ! max. Werte der lin. und nicht-lin.
     .  ,nonlcut,varmax,elmax	  ! Parameter bei Optimierung
      logical :: diffeecusp  ! Params for differnet CUSP condition for spin-like

c     required for OEM = one electron move
      real*8  ijuu(nmax,nmax)
      real*8  igradx(nmax,nmax),igrady(nmax,nmax),
     .        igradz(nmax,nmax),ilapla(nmax,nmax)
      real*8  oijuu(nmax,nmax)
      real*8  oigradx(nmax,nmax),oigrady(nmax,nmax),
     .        oigradz(nmax,nmax),oilapla(nmax,nmax)
      real*8  ju,uu1(nmax),ugrad1(3*nmax),ulapli1(nmax)
      real*8  oju,uu2(nmax),ugrad2(3*nmax),ulapli2(nmax)
      real*8  ouu1(nmax),ougrad1(3*nmax),oulapli1(nmax)
      real*8  ouu2(nmax),ougrad2(3*nmax),oulapli2(nmax)
      real*8  rijbar(kpmax,nmax,nmax),raibar(kpmax,amax,nmax)
      real*8  derivi(kpmax,nmax,nmax),fderiva(kpmax,amax,nmax)
      real*8  sderiva(kpmax,amax,nmax)
      real*8  orijbar(kpmax,nmax,nmax),oraibar(kpmax,amax,nmax)
      real*8  oderivi(kpmax,nmax,nmax),ofderiva(kpmax,amax,nmax)
      real*8  osderiva(kpmax,amax,nmax)

#ifdef NEWINIT
c     required for efficient ECP update initialization
      real*8, allocatable :: mFij(:,:), mGi(:)
      real*8              :: mJu = 0
#endif

      CONTAINS

c     --------------------------------
      subroutine jasinput_sm(lines,nl)
c     --------------------------------

c     jasinput reads Jastrow related input from 'lines'

      character(len=*), intent(in) :: lines(:)! lines array
      integer, intent(in)          :: nl      ! actual # of lines
      integer a,k,tpar,alstat,io

      read(lines(2),*) k0
      if (k0 /= 0) then
        read(lines(2),*,iostat=io) k0,jmode,lmax,mmax,nclast
        call assert(io==0,
     .     '(jasinput_sm) format: k0 jmode lmax mmax nclast')
        read(lines(3),*,iostat=io) lincut1,lincut2,nonlcut,varmax,elmax
        call assert(io==0,
     .     '(jasinput_sm) format: lincut1, lincut2, ncut, var, el')
        read(lines(5),*,iostat=io) bpg,(apg(a),a=1,atoms(nclast)%sa)
        call assert(io==0,
     .     '(jasinput_sm) format: b a_1 .. a_sa(nclast)')
        call assert(k0<=kmax,'(jasinput_sm): kmax too small')

        do k=1,k0
           read(lines(k+6),*) (tcontrol(k,tpar),tpar = 1,7),cjas(k)
        enddo
      else
        nclast = 0
        bpg = 1d0
      endif

#ifdef NEWINIT
      allocate(mFij(ne,ne), mGi(ne))
#endif


      end subroutine jasinput_sm


c=====================================================================

c     ------------------------------------
      subroutine jasinput_sm_new(lines,nl)
c     ------------------------------------

c     jasinput reads Jastrow related input from 'lines'

      character(len=*), intent(in) :: lines(:)! lines array
      integer, intent(in)          :: nl      ! actual # of lines
      integer a,k,alstat,io,idx,kk
      integer idxStart,idxEnd,i,j,line,pwr,enStart,enEnd
      integer jasCenter
      character(len=10) :: word(20)
      integer ::nWords

      lincut1 = 20.0
      lincut2 = 20.0
      nonlcut = 5.0
      varmax  = 100.0
      elmax   = 1.0
      jmode   = 1
      nclast = atoms_getNCLast(atoms)

      if (nclast==ncenter) then
         jasCenter = nscenter
      else if (nclast < ncenter) then
         jasCenter = nscenter - 1       ! ignore hydrogens for jastrow
      else
         call abortp("jasinput_sm_new: illegal nclast value")
      end if

      if (jastype=='sm0') then
         k0     = 1
      else
         if (jastype=='sm1') then
            k0     = 2 + jasCenter
         else if (jastype=='sm11') then
            k0     = 2 + 2*jasCenter
         else if (jastype=='sm2') then
            k0     = 4 + 3*jasCenter
         else if (jastype=='sm21') then
            k0     = 4 + 4*jasCenter
         else if (jastype=='sm3') then
            k0     = 4 + 5*jasCenter
         else if (jastype=='sm31') then
            k0     = 4 + 6*jasCenter
         else if (jastype=='sm4') then
            k0     = 4 + 15*jasCenter
         else if (jastype=='sm41') then
            k0     = 4 + 20*jasCenter
         else
            call abortp("jastrow (SM): this jastrow type not yet implemented")
         endif
      endif

      idx=2
      diffeecusp=.FALSE.
      call tokenize(lines(idx),word,nWords)
      if (nWords >= 1) then
         if (word(1)=='F' .or. word(1)=='T') then
            read(lines(idx),"(L1)", iostat=io) diffeecusp
            call assert(io==0,
     .      '(jasinput_sm_new) format: <additional logical parameters>')
            idx=idx+1
         endif
      endif

      call assert(k0<=kmax,'(jasinput_sm_new): kmax too small')

      if (jastype=='sm0') then
         k0     = 1
         lmax   = 0
         mmax   = 1
         tcontrol(1,1:7) = 0
         bpg = 1.d0
         cjas(1) = 0.5d0     ! cusp condition
         ! e-e terms
         if (diffeecusp) then
            !type 7 used for spin like and spin unlike CUSP
            tcontrol(1,1) = 7
         else
            tcontrol(1,1) = 1
         end if
         tcontrol(1,5)   = 1
         line = 1
      else
         lmax   = 4
         mmax   = 4
         if (jastype=='sm1' .or. jastype=='sm11') then
            lmax = 2
            mmax = 2
         endif
         read(lines(idx),*,iostat=io) bpg,(apg(a),a=1,atoms(nclast)%sa)
         call assert(io==0,
     .      '(jasinput_sm_new) format: b a_1 .. a_sa(nclast)')
         kk=0; idx=idx+1
         do
            read(lines(idx),*,iostat=io) (cjas(k),k=kk+1,min(k0,kk+5))
            call assert(io==0,
     .         '(jasinput_sm_new) format: c_1 .. c_k0, 5 per row')
            kk=kk+5; idx=idx+1
         if (kk>=k0) exit
         end do
         tcontrol(1:k0,1:7) = 0
         ! e-e terms
         if (diffeecusp) then
         !type 7 used for spin like and spin unlike CUSP
            tcontrol(1,1) = 7
         else
            tcontrol(1,1) = 1
         end if
         tcontrol(1,5) = 1
         tcontrol(2,1) = 1
         tcontrol(2,5) = 2
         tcontrol(3,1) = 1
         tcontrol(3,5) = 3
         tcontrol(4,1) = 1
         tcontrol(4,5) = 4
         idxStart = 1
         line = 4
         if (jastype=='sm1'.or.jastype=='sm11') then
            line = 2
         endif
         do i=1,atoms(nclast)%sa
            idxEnd = ncenter
            do a=idxStart+1,ncenter
               if (atoms(a)%sa /= atoms(idxStart)%sa) then
                  idxEnd = a-1
                  exit
               endif
            enddo
            ! e-n terms for same atom type sa(idxStart)
            enEnd=4
            if (jastype=='sm1'.or.jastype=='sm11') enEnd=2
            if (jastype=='sm1'.or.jastype=='sm2'.or.jastype=='sm3'
     &          .or.jastype=='sm4') enStart=2
            if (jastype=='sm11'.or.jastype=='sm21'.or.jastype=='sm31'
     &         .or.jastype=='sm41')  enStart=1
            do pwr=enStart,enEnd
               line = line + 1
               tcontrol(line,1) = 2
               tcontrol(line,3) = idxStart
               tcontrol(line,4) = idxEnd
               tcontrol(line,6) = pwr
            enddo
            if (jastype(1:3)=='sm3') then
               line = line + 1
               tcontrol(line,1) = 3
               tcontrol(line,3) = idxStart
               tcontrol(line,4) = idxEnd
               tcontrol(line,6) = 2
               tcontrol(line,7) = 2
               line = line + 1
               tcontrol(line,1) = 4
               tcontrol(line,3) = idxStart
               tcontrol(line,4) = idxEnd
               tcontrol(line,5) = 2
               tcontrol(line,6) = 2
            else if (jastype(1:3)=='sm4') then
               do pwr=enStart,enEnd
                  line = line + 1
                  tcontrol(line,1) = 3
                  tcontrol(line,3) = idxStart
                  tcontrol(line,4) = idxEnd
                  tcontrol(line,6) = pwr
                  tcontrol(line,7) = pwr
               enddo
               do pwr=enStart,enEnd
                  line = line + 1
                  tcontrol(line,1) = 4
                  tcontrol(line,3) = idxStart
                  tcontrol(line,4) = idxEnd
                  tcontrol(line,5) = 2
                  tcontrol(line,6) = pwr
               enddo
               do pwr=enStart,enEnd
                  line = line + 1
                  tcontrol(line,1) = 4
                  tcontrol(line,3) = idxStart
                  tcontrol(line,4) = idxEnd
                  tcontrol(line,5) = 3
                  tcontrol(line,6) = pwr
               enddo
               do pwr=enStart,enEnd
                  line = line + 1
                  tcontrol(line,1) = 4
                  tcontrol(line,3) = idxStart
                  tcontrol(line,4) = idxEnd
                  tcontrol(line,5) = 4
                  tcontrol(line,6) = pwr
               enddo
            endif
            idxStart = idxEnd + 1
         enddo
         if (line /= k0)
     .       call abortp("jastrow_sm_new: error in jastrow input")
      endif

#ifdef NEWINIT
      allocate(mFij(ne,ne), mGi(ne))
#endif

      end subroutine jasinput_sm_new


c=====================================================================

c     -------------------------------
      subroutine jasChangeType_sm(jt)
c     -------------------------------

      character(len=*),intent(in) :: jt

      integer a,iu,k
      integer idxStart,idxEnd,i,j,line,pwr,enStart,enEnd
      character(len=9) :: oldJasType
      integer jasCenter

      oldJasType = jastype
      jastype = jt

      diffeecusp=.FALSE.

      lincut1 = 20.0
      lincut2 = 20.0
      nonlcut = 5.0
      varmax  = 100.0
      elmax   = 1.0
      jmode   = 1

      nclast = atoms_getNCLast(atoms)

      if (nclast==ncenter) then
         jasCenter = nscenter
      else if (nclast < ncenter) then
         jasCenter = nscenter - 1       ! ignore hydrogens for jastrow
      else
         call abortp("jasinput_sm_new: illegal nclast value")
      end if

      if (jastype=='sm0') then
         k0     = 1
         lmax   = 0
         mmax   = 1
         tcontrol(1,1:7) = 0
         bpg = 1.d0
         cjas(1) = 0.5d0     ! cusp condition
         ! e-e terms
         if (diffeecusp) then
            !type 7 used for spin like and spin unlike CUSP
            tcontrol(1,1) = 7
         else
            tcontrol(1,1) = 1
         end if
         tcontrol(1,5)   = 1
         line = 1
      else
         if (jastype=='sm1') then
            k0     = 2 + jasCenter
         else if (jastype=='sm11') then
            k0     = 2 + 2*jasCenter
         else if (jastype=='sm2') then
            k0     = 4 + 3*jasCenter
         else if (jastype=='sm21') then
            k0     = 4 + 4*jasCenter
         else if (jastype=='sm3') then
            k0     = 4 + 5*jasCenter
         else if (jastype=='sm31') then
            k0     = 4 + 6*jasCenter
         else if (jastype=='sm4') then
            k0     = 4 + 15*jasCenter
         else if (jastype=='sm41') then
            k0     = 4 + 20*jasCenter
         else
            call abortp("jastrow (SM): this jastrow type not yet implemented")
         endif

         lmax   = 4
         mmax   = 4
         if (jastype=='sm1' .or. jastype=='sm11') then
            lmax = 2
            mmax = 2
         endif
         tcontrol(1:k0,1:7) = 0
         ! e-e terms
         if (diffeecusp) then
            !type 7 used for spin like and spin unlike CUSP
            tcontrol(1,1) = 7
         else
            tcontrol(1,1) = 1
         end if
         tcontrol(1,5) = 1
         tcontrol(2,1) = 1
         tcontrol(2,5) = 2
         tcontrol(3,1) = 1
         tcontrol(3,5) = 3
         tcontrol(4,1) = 1
         tcontrol(4,5) = 4
         idxStart = 1
         line = 4
         if (jastype=='sm1' .or. jastype=='sm11') then
            line = 2
         endif
         do i=1,atoms(nclast)%sa
            idxEnd = ncenter
            do a=idxStart+1,ncenter
               if (atoms(a)%sa /= atoms(idxStart)%sa) then
                  idxEnd = a-1
                  exit
               endif
            enddo
            ! e-n terms for same atom type sa(idxStart)
            enEnd=4
            if (jastype=='sm1'.or.jastype=='sm11') enEnd=2
            if (jastype=='sm1'.or.jastype=='sm2'.or.jastype=='sm3'
     &          .or.jastype=='sm4') enStart=2
            if (jastype=='sm11'.or.jastype=='sm21'.or.jastype=='sm31'
     &         .or.jastype=='sm41')  enStart=1
            do pwr=enStart,enEnd
               line = line + 1
               tcontrol(line,1) = 2
               tcontrol(line,3) = idxStart
               tcontrol(line,4) = idxEnd
               tcontrol(line,6) = pwr
            enddo
            if (jastype(1:3)=='sm3') then
               line = line + 1
               tcontrol(line,1) = 3
               tcontrol(line,3) = idxStart
               tcontrol(line,4) = idxEnd
               tcontrol(line,6) = 2
               tcontrol(line,7) = 2
               line = line + 1
               tcontrol(line,1) = 4
               tcontrol(line,3) = idxStart
               tcontrol(line,4) = idxEnd
               tcontrol(line,5) = 2
               tcontrol(line,6) = 2
            else if (jastype(1:3)=='sm4') then
               do pwr=enStart,enEnd
                  line = line + 1
                  tcontrol(line,1) = 3
                  tcontrol(line,3) = idxStart
                  tcontrol(line,4) = idxEnd
                  tcontrol(line,6) = pwr
                  tcontrol(line,7) = pwr
               enddo
               do pwr=enStart,enEnd
                  line = line + 1
                  tcontrol(line,1) = 4
                  tcontrol(line,3) = idxStart
                  tcontrol(line,4) = idxEnd
                  tcontrol(line,5) = 2
                  tcontrol(line,6) = pwr
               enddo
               do pwr=enStart,enEnd
                  line = line + 1
                  tcontrol(line,1) = 4
                  tcontrol(line,3) = idxStart
                  tcontrol(line,4) = idxEnd
                  tcontrol(line,5) = 3
                  tcontrol(line,6) = pwr
               enddo
               do pwr=enStart,enEnd
                  line = line + 1
                  tcontrol(line,1) = 4
                  tcontrol(line,3) = idxStart
                  tcontrol(line,4) = idxEnd
                  tcontrol(line,5) = 4
                  tcontrol(line,6) = pwr
               enddo
            endif
            idxStart = idxEnd + 1
         enddo
         if (line /= k0) call abortp("jastrow_reset_sm: wrong input")
      endif

      if (oldJasType(1:2)=='no'.or.oldJasType=='sm0') then
         bpg = 1.d0
         apg = 1.d0
         cjas(1:k0) = 0.d0
         if (JPD_enforceECusp) cjas(1) = 0.5d0
      else if (oldJasType=='sm1') then
         if (jastype=='sm2') then
            do a=jasCenter-1,0,-1
               cjas(4+3*a+1) = cjas(2+a+1)
               cjas(4+3*a+2:4+3*a+3) = 0.d0
            end do
            cjas(3:4) = 0.d0
         else if (jastype=='sm3') then
            do a=jasCenter-1,0,-1
               cjas(4+5*a+1) = cjas(2+a+1)
               cjas(4+5*a+2:4+3*a+5) = 0.d0
            end do
            cjas(3:4) = 0.d0
         else
            call abortp("jasChangeType: this change not implemented")
         end if
      end if

      if (logmode >= 2) then
         write(iul,*) ' changing SM Jastrow from ',trim(oldJasType),
     &        ' to type ',trim(jastype)
         write(iul,*) ' with ',k0,' parameters:'
         write(iul,'(10G12.4)') bpg,apg(1:jasCenter)
         write(iul,'(10G12.4)') cjas(1:k0)
         if (logmode >=3) then
            write(iul,*) 'new tcontrol:'
            do k=1,k0
               write(iul,'(7I5)') tcontrol(k,1:7)
            enddo
         endif
      endif

      call flush(iul)

      end subroutine jasChangeType_sm


c=======================================================================


c     ---------------------------
      subroutine jasoutput_sm(iu)
c     ---------------------------

c     jasoutput writes Jastrow in input format to file unit 'iu'

      integer a,k,iu

      if (k0 /= 0) then
         write(iu,'(5i5)') k0,jmode,lmax,mmax,nclast
         write(iu,'(5f13.3)') lincut1,lincut2,nonlcut,varmax,elmax
         write(iu,'(5F14.7)') bpg,(apg(a),a=1,atoms(nclast)%sa)
         do k=1,k0
            write(iu,'(7I3,3X,F14.7)') tcontrol(k,1),tcontrol(k,2),
     .           tcontrol(k,3),tcontrol(k,4),
     .           tcontrol(k,5),tcontrol(k,6),
     .           tcontrol(k,7),cjas(k)
         enddo
      endif
      end subroutine jasoutput_sm

c=======================================================================


c     ---------------------------------
      subroutine jas_shortoutput_sm(iu)
c     ---------------------------------

      integer k,iu

      if (k0 /= 0) then
         write(iu,'(i5,a)') k0," generic Schmidt-Moskowitz terms"
      endif
      end subroutine jas_shortoutput_sm

c=======================================================================


c     -------------------------------
      subroutine jasoutput_sm_new(iu)
c     -------------------------------

      integer, intent(in) :: iu

c     jasoutput writes Jastrow terms to file unit 'iu'
c     suitable for wf-file

      integer a,k

      if (diffeecusp) then
         write(iu,'(L1)') diffeecusp
      endif
      if (k0.ne.0) then
         write(iu,'(5F10.5)') bpg,(apg(a),a=1,atoms(nclast)%sa)
         write(iu,'(5F10.5)') (cjas(k),k=1,k0)
      endif
      end subroutine jasoutput_sm_new

c=======================================================================


c     -------------------------------------
      subroutine jas_shortoutput_sm_new(iu)
c     -------------------------------------

      integer, intent(in) :: iu

      integer a,k
      if (k0 /= 0) then
         write(iu,'(i5,2a)') k0," generic Schmidt-Moskowitz terms of type ",trim(jastype)
      endif
      if (diffeecusp) then
         write(iu,'(a)')"with different cusp term for spin like and spin unlike e-"
      endif
      end subroutine jas_shortoutput_sm_new

c=======================================================================

c     ----------------------------------------------------------------
      subroutine jassmall(x,y,z,rai,rij,optType,ju,jud,julapl,julapli,uuk)
c     ----------------------------------------------------------------

c rijbar:  (rij/(1+brij))**m
c raibar:  (rai/(1+arai))**l
c derivi:  term used for calculation of gradient of (rij/(1+brij))**m upper triangle,
c          term used for calculation of laplacian of (rij/(1+brij))**m lower triangle
c fderiva: term used for calculation of gradient of (rai/(1+arai))**l ,
c sdiriva: term used for calculation of laplacian of (rai/(1+arai))**l


      real*8, intent(in)            ::  x(:),y(:),z(:)     ! cartesian coords
      real*8, intent(in)            ::  rai(:,:),rij(:,:)  ! distances r_ai, r_ij
      character(len=*)              ::  optType            ! parameter optimization
      real*8, intent(out)           ::  ju                 ! in exp(U)
      real*8, intent(out)           ::  jud(:)             ! \nabla U
      real*8, intent(out)           ::  julapl             ! laplacian U
      real*8, intent(out)           ::  julapli(:)         ! \nabla_i^2 U  ((x_i,y_i,z_i))
      real*8, optional, intent(out) ::  uuk(:)
      integer i,j,a,k
      integer s,v
      real*8  uu(kmax),ugrad(3*nmax,kmax),ulapl(kmax),ulapli(nmax,kmax)
      real*8  htmp2a,htmp2b,htmp3
      real*8  xij(nmax,nmax),yij(nmax,nmax),zij(nmax,nmax)
      real*8  xai(amax,nmax),yai(amax,nmax),zai(amax,nmax)
      real*8  xyzi,xyzj
      real*8  tmp,tmp1,tmp2,t

#ifdef NEWINIT
            mFij = 0.d0
            mGi  = 0.d0
#endif


c---------Enforce cusp condition for unlike electrons--------------
c         ((here to allow modification of bpg for optimization))
      if (JPD_enforceECusp) then
         cjas(1) = 0.5d0/bpg
      endif

c---------Calculation of Jastrow factor and its deriviatives-------
      do k = 1,k0
         uu(k) = 0
         ulapl(k) = 0
         do i = 1,3*ne
            ugrad(i,k)  = 0
         enddo
         do i = 1,ne
            ulapli(i,k) = 0
         enddo
      enddo

c---------Precalculation of Terms-------------------------------
      t    = ne-1

      do i = 1,ne
         do a = 1,nclast
            xai(a,i) = x(i)-atoms(a)%cx
            yai(a,i) = y(i)-atoms(a)%cy
            zai(a,i) = z(i)-atoms(a)%cz
            tmp1           = 1/(1+apg(atoms(a)%sa)*rai(a,i))
            tmp2           = tmp1*apg(atoms(a)%sa)*rai(a,i)
            raibar(1,a,i)  = tmp2
            fderiva(1,a,i) = tmp2*tmp1/(rai(a,i)*rai(a,i))
            sderiva(1,a,i) = fderiva(1,a,i)*tmp1*2d0
            do s = 2,lmax
              raibar(s,a,i)  = tmp2*raibar(s-1,a,i)
              fderiva(s,a,i) = s*raibar(s-1,a,i)*fderiva(1,a,i)
              sderiva(s,a,i) = fderiva(s,a,i)*(s+1)*tmp1
            enddo
         enddo
         do j = i+1,ne
            xij(j,i) = x(i)-x(j)
            yij(j,i) = y(i)-y(j)
            zij(j,i) = z(i)-z(j)
            tmp1          = 1/(1+bpg*rij(i,j))
            tmp2          = tmp1*bpg*rij(i,j)
            rijbar(1,j,i) = tmp2
            derivi(1,j,i) = tmp2*tmp1/(rij(i,j)*rij(i,j))
            derivi(1,i,j) = derivi(1,j,i)*tmp1*2d0
            do v = 2,mmax
              rijbar(v,j,i)  = tmp2*rijbar(v-1,j,i)
              derivi(v,j,i)  = v*rijbar(v-1,j,i)*derivi(1,j,i)
              derivi(v,i,j)  = derivi(v,j,i)*(v+1)*tmp1
            enddo
         enddo
      enddo

c---------main part-----------------------------------------------------

      do k = 1,k0

c------------Electron-Electron-Correlation Terms------------------------
      if (tcontrol(k,1)==1) then
      v = tcontrol(k,5)
      do i = 1,ne
         do j = i+1,ne
            uu(k)          =  uu(k)          + rijbar(v,j,i)
#ifdef NEWINIT
            mFij(j,i) = mFij(j,i) + cjas(k)*rijbar(v,j,i)
#endif
            ugrad(3*i-2,k) =  ugrad(3*i-2,k) + derivi(v,j,i)*xij(j,i)
            ugrad(3*i-1,k) =  ugrad(3*i-1,k) + derivi(v,j,i)*yij(j,i)
            ugrad(3*i,k)   =  ugrad(3*i,k)   + derivi(v,j,i)*zij(j,i)
            ugrad(3*j-2,k) =  ugrad(3*j-2,k) - derivi(v,j,i)*xij(j,i)
            ugrad(3*j-1,k) =  ugrad(3*j-1,k) - derivi(v,j,i)*yij(j,i)
            ugrad(3*j,k)   =  ugrad(3*j,k)   - derivi(v,j,i)*zij(j,i)
            ulapli(i,k)    =  ulapli(i,k)    + derivi(v,i,j)
            ulapli(j,k)    =  ulapli(j,k)    + derivi(v,i,j)
         enddo
      enddo
      endif


      ! satisfy cusp for same and different spin
      if (tcontrol(k,1)==7) then
      v = 1                         ! only linear r_ij term!
      do i = 1,nalpha
         do j = i+1,nalpha
            uu(k)          =  uu(k)          + 0.5d0*rijbar(v,j,i)
#ifdef NEWINIT
            mFij(j,i) = mFij(j,i) + cjas(k)*0.5d0*rijbar(v,j,i)
#endif
            ugrad(3*i-2,k) =  ugrad(3*i-2,k) + 0.5d0*derivi(v,j,i)*xij(j,i)
            ugrad(3*i-1,k) =  ugrad(3*i-1,k) + 0.5d0*derivi(v,j,i)*yij(j,i)
            ugrad(3*i,k)   =  ugrad(3*i,k)   + 0.5d0*derivi(v,j,i)*zij(j,i)
            ugrad(3*j-2,k) =  ugrad(3*j-2,k) - 0.5d0*derivi(v,j,i)*xij(j,i)
            ugrad(3*j-1,k) =  ugrad(3*j-1,k) - 0.5d0*derivi(v,j,i)*yij(j,i)
            ugrad(3*j,k)   =  ugrad(3*j,k)   - 0.5d0*derivi(v,j,i)*zij(j,i)
            ulapli(i,k)    =  ulapli(i,k)    + 0.5d0*derivi(v,i,j)
            ulapli(j,k)    =  ulapli(j,k)    + 0.5d0*derivi(v,i,j)
         enddo
         do j = nalpha+1,ne
            uu(k)          =  uu(k)          + rijbar(v,j,i)
#ifdef NEWINIT
            mFij(j,i) = mFij(j,i) + cjas(k)*rijbar(v,j,i)
#endif
            ugrad(3*i-2,k) =  ugrad(3*i-2,k) + derivi(v,j,i)*xij(j,i)
            ugrad(3*i-1,k) =  ugrad(3*i-1,k) + derivi(v,j,i)*yij(j,i)
            ugrad(3*i,k)   =  ugrad(3*i,k)   + derivi(v,j,i)*zij(j,i)
            ugrad(3*j-2,k) =  ugrad(3*j-2,k) - derivi(v,j,i)*xij(j,i)
            ugrad(3*j-1,k) =  ugrad(3*j-1,k) - derivi(v,j,i)*yij(j,i)
            ugrad(3*j,k)   =  ugrad(3*j,k)   - derivi(v,j,i)*zij(j,i)
            ulapli(i,k)    =  ulapli(i,k)    + derivi(v,i,j)
            ulapli(j,k)    =  ulapli(j,k)    + derivi(v,i,j)
         enddo
      enddo
      do i = nalpha+1,ne
         do j = i+1,ne
            uu(k)          =  uu(k)          + 0.5d0*rijbar(v,j,i)
#ifdef NEWINIT
            mFij(j,i) = mFij(j,i) + cjas(k)*0.5d0*rijbar(v,j,i)
#endif
            ugrad(3*i-2,k) =  ugrad(3*i-2,k) + 0.5d0*derivi(v,j,i)*xij(j,i)
            ugrad(3*i-1,k) =  ugrad(3*i-1,k) + 0.5d0*derivi(v,j,i)*yij(j,i)
            ugrad(3*i,k)   =  ugrad(3*i,k)   + 0.5d0*derivi(v,j,i)*zij(j,i)
            ugrad(3*j-2,k) =  ugrad(3*j-2,k) - 0.5d0*derivi(v,j,i)*xij(j,i)
            ugrad(3*j-1,k) =  ugrad(3*j-1,k) - 0.5d0*derivi(v,j,i)*yij(j,i)
            ugrad(3*j,k)   =  ugrad(3*j,k)   - 0.5d0*derivi(v,j,i)*zij(j,i)
            ulapli(i,k)    =  ulapli(i,k)    + 0.5d0*derivi(v,i,j)
            ulapli(j,k)    =  ulapli(j,k)    + 0.5d0*derivi(v,i,j)
         enddo
      enddo
      endif


c---------Electron-Nucleus-Correlation-Terms----------------------------------

      if (tcontrol(k,1)==2) then
      s   = tcontrol(k,6)
      do i = 1,ne
         do a = tcontrol(k,3),tcontrol(k,4)
            uu(k)          =  uu(k)          + t*raibar(s,a,i)
#ifdef NEWINIT
            mGi(i) = mGi(i) + cjas(k)*t*raibar(s,a,i)
#endif
            ugrad(3*i-2,k) =  ugrad(3*i-2,k) + t*fderiva(s,a,i)*xai(a,i)
            ugrad(3*i-1,k) =  ugrad(3*i-1,k) + t*fderiva(s,a,i)*yai(a,i)
            ugrad(3*i,k)   =  ugrad(3*i,k)   + t*fderiva(s,a,i)*zai(a,i)
            ulapli(i,k)    =  ulapli(i,k)    + t*sderiva(s,a,i)
     	 enddo
      enddo
      endif

c---------Electron-Electron-Nucleus-Correlation-Terms----------------------

      if (tcontrol(k,1)==3) then
      s  = tcontrol(k,6)
      do i = 1,ne
        do j = i+1,ne
          do a = tcontrol(k,3),tcontrol(k,4)

            htmp2a  = raibar(s,a,j)*fderiva(s,a,i)
            htmp2b  = raibar(s,a,i)*fderiva(s,a,j)

            uu(k)          =  uu(k)       + raibar(s,a,i)*raibar(s,a,j)
#ifdef NEWINIT
            mFij(j,i) = mFij(j,i) + cjas(k)*raibar(s,a,i)*raibar(s,a,j)
#endif
            ugrad(3*i-2,k) =  ugrad(3*i-2,k) + htmp2a*xai(a,i)
            ugrad(3*i-1,k) =  ugrad(3*i-1,k) + htmp2a*yai(a,i)
            ugrad(3*i,k)   =  ugrad(3*i,k)   + htmp2a*zai(a,i)
            ugrad(3*j-2,k) =  ugrad(3*j-2,k) + htmp2b*xai(a,j)
            ugrad(3*j-1,k) =  ugrad(3*j-1,k) + htmp2b*yai(a,j)
            ugrad(3*j,k)   =  ugrad(3*j,k)   + htmp2b*zai(a,j)
            ulapli(i,k)    =  ulapli(i,k) + raibar(s,a,j)*sderiva(s,a,i)
            ulapli(j,k)    =  ulapli(j,k) + raibar(s,a,i)*sderiva(s,a,j)
          enddo
        enddo
      enddo
      endif
c-----------------------------------------------------------------------

      if (tcontrol(k,1)==4) then
      v  = tcontrol(k,5)
      s  = tcontrol(k,6)
      do i = 1,ne
        do j = i+1,ne
          do a = tcontrol(k,3),tcontrol(k,4)

            xyzi =xij(j,i)*xai(a,i)+yij(j,i)*yai(a,i)+zij(j,i)*zai(a,i)
            xyzj =xij(j,i)*xai(a,j)+yij(j,i)*yai(a,j)+zij(j,i)*zai(a,j)

            tmp1    = raibar(s,a,i)+raibar(s,a,j)
            tmp2    = tmp1*derivi(v,i,j)

            htmp2a  = rijbar(v,j,i)*fderiva(s,a,i)
            htmp2b  = rijbar(v,j,i)*fderiva(s,a,j)
            htmp3   = tmp1*derivi(v,j,i)

            uu(k)         =uu(k)         +tmp1*rijbar(v,j,i)
#ifdef NEWINIT
            mFij(j,i) = mFij(j,i) + cjas(k)*tmp1*rijbar(v,j,i)
#endif
            ugrad(3*i-2,k)=ugrad(3*i-2,k)+htmp2a*xai(a,i)+htmp3*xij(j,i)
            ugrad(3*i-1,k)=ugrad(3*i-1,k)+htmp2a*yai(a,i)+htmp3*yij(j,i)
            ugrad(3*i,k)  =ugrad(3*i,k)  +htmp2a*zai(a,i)+htmp3*zij(j,i)
            ugrad(3*j-2,k)=ugrad(3*j-2,k)+htmp2b*xai(a,j)-htmp3*xij(j,i)
            ugrad(3*j-1,k)=ugrad(3*j-1,k)+htmp2b*yai(a,j)-htmp3*yij(j,i)
            ugrad(3*j,k)  =ugrad(3*j,k)  +htmp2b*zai(a,j)-htmp3*zij(j,i)
            ulapli(i,k)   =ulapli(i,k)   +tmp2
     .                    +2*derivi(v,j,i)*fderiva(s,a,i)*xyzi
     .                    +rijbar(v,j,i)*sderiva(s,a,i)
            ulapli(j,k)   =ulapli(j,k)   +tmp2
     .                    -2*derivi(v,j,i)*fderiva(s,a,j)*xyzj
     .                    +rijbar(v,j,i)*sderiva(s,a,j)
          enddo
        enddo
      enddo
      endif
c-----------------------------------------------------------------------

      enddo ! K-Schleife
c--------Calculation of U and its Deriviatives--------------------------

      ju = 0d0
      do i=1,3*ne
         jud(i) = 0d0
      enddo
      do i=1,ne
         julapli(i) = 0d0
      enddo

      do k=1,k0
        ju = ju + cjas(k)*uu(k)
      enddo

#ifdef NEWINIT
      mJu = ju
#endif

      do k=1,k0
        do i=1,3*ne
          jud(i) = jud(i) + cjas(k)*ugrad(i,k)
        enddo
        do i=1,ne
          julapli(i)  = julapli(i)  + cjas(k)*ulapli(i,k)
          ulapl(k)    = ulapl(k)    +         ulapli(i,k)
        enddo
      enddo

      julapl = sum(julapli(1:ne))

      if (optType == 'jastrow' .or. optType == 'jas+ci' .or. optType == 'jas+mo' .or. optType == 'jas+mo+ci') then
         call assert(allocated(uk),
     .               'jassmall: jastrowParamData not allocated')
         call assert(uk_params == k0-1,'jassmall: mismatch of # params')
         uk = uu(2:k0)
         ukgrad = ugrad(1:3*ne,2:k0)
         uklapl = ulapl(2:k0)
         uklapli = ulapli(1:ne,2:k0)
      endif
      if (present(uuk)) then
            do,i=1,k0-1
               uuk(i)=uu(1+i)
            enddo
       endif
      end subroutine jassmall


c=======================================================================

c     -------------------------------
      subroutine jassmInit(Rdu,iMode)
c     -------------------------------

c this version uses RdataUpdate to keep track of auxiliary data for one electron updates
c initialize auxiliary data
c no electron derivatives here!

c rijbar:  (rij/(1+b*rij))**m
c raibar:  (rai/(1+a*rai))**l

      type(RdataUpdate), intent(inout) :: Rdu             ! data structure for electron update calculations
      integer, intent(in)              :: iMode

      integer i,j,a,k
      integer s,v
      !real*8  uu(k0)
      real*8  tmp,tmp1,tmp2,t


      call Rdu%initENSize(lmax,nclast)

c---------Enforce cusp condition for unlike electrons--------------
c         ((here to allow modification of bpg for optimization))
      if (JPD_enforceECusp) then
         cjas(1) = 0.5d0/bpg
      endif

c---------Precalculation of Terms-------------------------------

      !uu = 0
      t    = ne-1

      if (iMode == 0) then
        do i = 1,ne
           do a = 1,nclast
              tmp1          = apg(atoms(a)%sa)*Rdu%rai(a,i)/(1.d0+apg(atoms(a)%sa)*Rdu%rai(a,i))
              raibar(1,a,i) = tmp1
              do s = 2,lmax
                 raibar(s,a,i)  = raibar(s-1,a,i)*tmp1
              enddo
           enddo
           do j = i+1,ne
              tmp2          = bpg*Rdu%rij(i,j)/(1.d0+bpg*Rdu%rij(i,j))
              rijbar(1,j,i) = tmp2
              do v = 2,mmax
                 rijbar(v,j,i)  = rijbar(v-1,j,i)*tmp2
              enddo
           enddo
        enddo
      endif

c---------main part-----------------------------------------------------

      if (iMode <= 1) then

      Rdu%Fij = 0.d0
      Rdu%Gi = 0.d0

      do k = 1,k0

         ! Electron-Electron-Correlation Terms------------------------
         if (tcontrol(k,1)==1) then
            v = tcontrol(k,5)
            do i = 1,ne
               do j = i+1,ne
                  !uu(k) = uu(k) + rijbar(v,j,i)
                  Rdu%Fij(j,i) = Rdu%Fij(j,i) + cjas(k)*rijbar(v,j,i)
               enddo
            enddo
         endif

         ! satisfy cusp for same and different spin
         if (tcontrol(k,1)==7) then
            v = 1                         ! only linear r_ij term!
            do i = 1,nalpha
               do j = i+1,nalpha
                  !uu(k) =  uu(k) + 0.5d0*rijbar(v,j,i)
                  Rdu%Fij(j,i) = Rdu%Fij(j,i) + cjas(k)*0.5d0*rijbar(v,j,i)
               enddo
               do j = nalpha+1,ne
                  !uu(k) =  uu(k) + rijbar(v,j,i)
                  Rdu%Fij(j,i) = Rdu%Fij(j,i) + cjas(k)*rijbar(v,j,i)
               enddo
            enddo
            do i = nalpha+1,ne
               do j = i+1,ne
                  !uu(k) =  uu(k) + 0.5d0*rijbar(v,j,i)
                  Rdu%Fij(j,i) = Rdu%Fij(j,i) + cjas(k)*0.5d0*rijbar(v,j,i)
               enddo
            enddo
         endif

         ! Electron-Nucleus-Correlation-Terms----------------------------------

         if (tcontrol(k,1)==2) then
            s = tcontrol(k,6)
            do i = 1,ne
               do a = tcontrol(k,3),tcontrol(k,4)
                  !uu(k) = uu(k) + t*raibar(s,a,i)
                  Rdu%Gi(i) = Rdu%Gi(i) + cjas(k)*t*raibar(s,a,i)
               enddo
            enddo
         endif

         ! Electron-Electron-Nucleus-Correlation-Terms----------------------

         if (tcontrol(k,1)==3) then
            s  = tcontrol(k,6)
            do i = 1,ne
              do j = i+1,ne
                do a = tcontrol(k,3),tcontrol(k,4)
                  !uu(k) =  uu(k) + raibar(s,a,i)*raibar(s,a,j)
                  Rdu%Fij(j,i) = Rdu%Fij(j,i) + cjas(k)*raibar(s,a,i)*raibar(s,a,j)
                enddo
              enddo
            enddo
         endif
   
         if (tcontrol(k,1)==4) then
            v  = tcontrol(k,5)
            s  = tcontrol(k,6)
            do i = 1,ne
              do j = i+1,ne
                do a = tcontrol(k,3),tcontrol(k,4)
                  tmp = raibar(s,a,i) + raibar(s,a,j)
                  !uu(k) =uu(k) + tmp*rijbar(v,j,i)
                  Rdu%Fij(j,i) = Rdu%Fij(j,i) + cjas(k)*tmp*rijbar(v,j,i)
                enddo
              enddo
            enddo
         endif
   
      enddo ! K-Schleife

      ! Calculation of U

      ju = 0d0
      do i=1,ne
        ju = ju + Rdu%Gi(i)
        do j=i+1,ne
          ju = ju + Rdu%Fij(j,i)
        enddo
      enddo
      !!!write(iul,'(a,g20.10)') 'DBG:jassmInit1',ju

      !!ju = 0d0
      !!do k=1,k0
      !!  ju = ju + cjas(k)*uu(k)
      !!enddo

      Rdu%U = ju
      Rdu%U0 = ju
      Rdu%ieJasold = 0

C       write(iul,'(a,g20.10)') 'DBG:jassmInit2',ju
C       write(iul,'(a)') 'DBG:uu1:F'
C       do i=1,ne
C         write(iul,'(i3)') i
C         write(iul,'(5g20.10)') (Rdu%Fij(j,i),j=1,ne)
C       enddo
C       write(iul,'(a)') 'DBG:uu2:G'
C       write(iul,'(5g20.10)') (Rdu%Gi(i),i=1,ne)

#ifdef NEWINIT
      else if (iMode == 2) then

         ! use data from previous jassmall call

         Rdu%Fij = mFij
         Rdu%Gi = mGi
         Rdu%U = mJu
         Rdu%U0 = mJu
         Rdu%ieJasold = 0
#endif

      else 

         call abortp('jassmInit: illegal mode')

      endif

C      write(iul,'(a,g20.10)') 'DBG:jassmInit:U',Rdu%U

      call Rdu%markJastrowValid()

      end subroutine jassmInit


c=======================================================================

c     -------------------------------
      subroutine jassmInitWithUk(Rdu)
c     -------------------------------

c this version uses RdataUpdate to keep track of auxiliary data for one electron updates
c initialize auxiliary data
c no electron derivatives here!

c rijbar:  (rij/(1+b*rij))**m
c raibar:  (rai/(1+a*rai))**l

      type(RdataUpdate), intent(inout) :: Rdu             ! data structure for electron update calculations

      integer i,j,a,k,k1,k2,km1
      integer s,v
      !real*8  uu(k0)
      real*8  tmp,tmp1,tmp2,t

      call assert(Rdu%npJ1+Rdu%npJ2==k0-1,"jassmInitWithUk: inconsistent jastrow parameter numbers")

      call Rdu%initENSize(lmax,nclast)

c---------Enforce cusp condition for unlike electrons--------------
c         ((here to allow modification of bpg for optimization))
      if (JPD_enforceECusp) then
         cjas(1) = 0.5d0/bpg
      endif

c---------Precalculation of Terms-------------------------------

      !uu = 0
      t    = ne-1

      do i = 1,ne
         do a = 1,nclast
            tmp1          = apg(atoms(a)%sa)*Rdu%rai(a,i)/(1.d0+apg(atoms(a)%sa)*Rdu%rai(a,i))
            raibar(1,a,i) = tmp1
            do s = 2,lmax
               raibar(s,a,i)  = raibar(s-1,a,i)*tmp1
            enddo
         enddo
         do j = i+1,ne
            tmp2          = bpg*Rdu%rij(i,j)/(1.d0+bpg*Rdu%rij(i,j))
            rijbar(1,j,i) = tmp2
            do v = 2,mmax
               rijbar(v,j,i)  = rijbar(v-1,j,i)*tmp2
            enddo
         enddo
      enddo

c---------main part-----------------------------------------------------

      Rdu%Fij = 0.d0
      Rdu%Gi = 0.d0
      Rdu%Fijk = 0.d0
      Rdu%Gki = 0.d0

      Rdu%Uk = 0.d0

      k1 = 0
      k2 = 0

      do k = 1,k0

         km1 = k - 1

         ! Electron-Electron-Correlation Terms------------------------
         if (tcontrol(k,1)==1) then
            v = tcontrol(k,5)
            if (v==1) then     ! cusp term, fixed parameter
               do i = 1,ne
                  do j = i+1,ne
                     Rdu%Fij(j,i)    = Rdu%Fij(j,i) + cjas(k)*rijbar(v,j,i)
                  enddo
               enddo
            else
               k2 = k2 + 1
               do i = 1,ne
                  do j = i+1,ne
                     Rdu%Fij(j,i)    = Rdu%Fij(j,i) + cjas(k)*rijbar(v,j,i)
                     Rdu%Fijk(j,i,k2) = Rdu%Fijk(j,i,k2) + rijbar(v,j,i)
                     Rdu%Uk(km1) = Rdu%Uk(km1) + rijbar(v,j,i)
                  enddo
               enddo
            endif
         endif

         ! satisfy cusp for same and different spin
         if (tcontrol(k,1)==7) then
            v = 1                         ! only linear r_ij term! fixed parameter!
            do i = 1,nalpha
               do j = i+1,nalpha
                  Rdu%Fij(j,i)    = Rdu%Fij(j,i) + cjas(k)*0.5d0*rijbar(v,j,i)
               enddo
               do j = nalpha+1,ne
                  Rdu%Fij(j,i) = Rdu%Fij(j,i) + cjas(k)*rijbar(v,j,i)
               enddo
            enddo
            do i = nalpha+1,ne
               do j = i+1,ne
                  Rdu%Fij(j,i) = Rdu%Fij(j,i) + cjas(k)*0.5d0*rijbar(v,j,i)
               enddo
            enddo
         endif

         ! Electron-Nucleus-Correlation-Terms----------------------------------

         if (tcontrol(k,1)==2) then
            s = tcontrol(k,6)
            k1 = k1 + 1
            do i = 1,ne
               do a = tcontrol(k,3),tcontrol(k,4)
                  tmp = t*raibar(s,a,i)
                  Rdu%Gi(i)    = Rdu%Gi(i) + cjas(k)*tmp
                  Rdu%Gki(k1,i) = Rdu%Gki(k1,i) + tmp
                  Rdu%Uk(km1) = Rdu%Uk(km1) + tmp
               enddo
            enddo
         endif

         ! Electron-Electron-Nucleus-Correlation-Terms----------------------

         if (tcontrol(k,1)==3) then
            s  = tcontrol(k,6)
            k2 = k2 + 1
            do i = 1,ne
              do j = i+1,ne
                do a = tcontrol(k,3),tcontrol(k,4)
                  tmp2 = raibar(s,a,i)*raibar(s,a,j)
                  Rdu%Fij(j,i)    = Rdu%Fij(j,i) + cjas(k)*tmp2
                  Rdu%Fijk(j,i,k2) = Rdu%Fijk(j,i,k2) + tmp2
                  Rdu%Uk(km1) = Rdu%Uk(km1) + tmp2
                enddo
              enddo
            enddo
         endif
   
         if (tcontrol(k,1)==4) then
            v  = tcontrol(k,5)
            s  = tcontrol(k,6)
            k2 = k2 + 1
            do i = 1,ne
              do j = i+1,ne
                do a = tcontrol(k,3),tcontrol(k,4)
                  tmp1 = rijbar(v,j,i)*(raibar(s,a,i) + raibar(s,a,j))
                  Rdu%Fij(j,i)    = Rdu%Fij(j,i) + cjas(k)*tmp1
                  Rdu%Fijk(j,i,k2) = Rdu%Fijk(j,i,k2) + tmp1
                  Rdu%Uk(km1) = Rdu%Uk(km1) + tmp1
                enddo
              enddo
            enddo
         endif
   
      enddo ! K-Schleife

      call assert(k1 == size(Rdu%Gki,1) .and. k2 == size(Rdu%Fijk,3) .and. k1+k2==k0-1,
     &   "jassmInitWithUk: k1,k2 illegal values")

      ! Calculation of U

      Rdu%U = 0d0
      do i=1,ne
        Rdu%U = Rdu%U + Rdu%Gi(i)
        do j=i+1,ne
          Rdu%U = Rdu%U + Rdu%Fij(j,i)
        enddo
      enddo
      !!!write(iul,'(a,g20.10)') 'DBG:jassmInit1',ju

      Rdu%U0 = Rdu%U
      Rdu%Uk0 = Rdu%Uk
      Rdu%ieJasold = 0

      !!!write(iul,'(a,g20.10)') 'DBG:jassmInitwpd:U,Uk',Rdu%U
      !!!write(iul,'(5g20.10)') Rdu%Uk 

C       write(iul,'(a)') 'Fijk:'
C       do k=1,k0-1
C         do i=1,ne
C           write(iul,'(2i4,4x,g20.10)') k,i,Rdu%Gki(k,i)
C           do j=i+1,ne
C             write(iul,'(3i4,g20.10)') k,i,j,Rdu%Fijk(j,i,k)
C           enddo
C         enddo
C       enddo



C       write(iul,'(a,g20.10)') 'DBG:jassmInit2',ju
C       write(iul,'(a)') 'DBG:uu1:F'
C       do i=1,ne
C         write(iul,'(i3)') i
C         write(iul,'(5g20.10)') (Rdu%Fij(j,i),j=1,ne)
C       enddo
C       write(iul,'(a)') 'DBG:uu2:G'
C       write(iul,'(5g20.10)') (Rdu%Gi(i),i=1,ne)

      call Rdu%markJastrowValid()

      end subroutine jassmInitWithUk


c=======================================================================

c     ------------------------------
      subroutine jassmUpdate(Rdu,ie)
c     ------------------------------

c this version uses RdataUpdate to keep track of auxiliary data for one electron updates
c one electron update of U: restore old electron (ieold) data when new electron (ie /= ieold)
c no electron derivatives here!

c rijbar:  (rij/(1+b*rij))**m
c raibar:  (rai/(1+a*rai))**l

      type(RdataUpdate), intent(inout) :: Rdu    ! data structure for electron update calculations
      integer, intent(in)              :: ie     ! update for electron ie

      integer i,j,a,k
      integer s,v
      !real*8  uu(k0)
      real*8  tmp,tmp1,tmp2,t,Fsum

      ! cusp enforcement (s. jassmall/Init) should not be necessary

      if (ie /= Rdu%ieJasold) then
         if (Rdu%ieJasold /= 0) then
           ! necessary to restore one-electron terms for previous electron ieold
           ! also restore U to original U
           Rdu%Fij(Rdu%ieJasold,1:Rdu%ieJasold-1) = Rdu%Fijold(1:Rdu%ieJasold-1) 
           Rdu%Fij(Rdu%ieJasold+1:ne,Rdu%ieJasold) = Rdu%Fijold(Rdu%ieJasold+1:ne)
           Rdu%Gi(Rdu%ieJasold) = Rdu%Giold
           raibar(1:Rdu%endim1,1:Rdu%endim2,Rdu%ieJasold) = Rdu%raibarOld(:,:)
           Rdu%U = Rdu%U0
         endif
         Rdu%Fijold(1:ie-1) = Rdu%Fij(ie,1:ie-1)
         Rdu%Fijold(ie+1:ne) = Rdu%Fij(ie+1:ne,ie)
         Rdu%Giold = Rdu%Gi(ie)
         Rdu%raibarOld(:,:) = raibar(1:Rdu%endim1,1:Rdu%endim2,ie)
         Rdu%ieJasold = ie
      ! else same electron as before, simply update ie
      endif

      !!!write(iul,'(ai3,g20.10)') 'DBG:jassmUpdate:ie,Uold:',ie,Rdu%U

      ! subtract current=previous values the three Jastrow terms
      Fsum = 0.d0
      do i=1,ie-1
         Fsum = Fsum + Rdu%Fij(ie,i)
      enddo
      do j=ie+1,ne
         Fsum = Fsum + Rdu%Fij(j,ie)
      enddo
      Rdu%U = Rdu%U - Rdu%Gi(ie) - Fsum

      !!!write(iul,'(a)') 'DBG:Fsumold,Giold,Uold-Giold-Fsumold'
      !!!write(iul,'(5g20.10)') Fsum,Rdu%Gi(ie),Rdu%U


c---------Precalculation of Terms-------------------------------

      !uu = 0
      t    = ne-1

      i = ie ! only ie !

      do a = 1,nclast
         tmp1          = apg(atoms(a)%sa)*Rdu%rai(a,i)/(1.d0+apg(atoms(a)%sa)*Rdu%rai(a,i))
         raibar(1,a,i) = tmp1
         do s = 2,lmax
            raibar(s,a,i)  = raibar(s-1,a,i)*tmp1
         enddo
      enddo
      do j = 1,i-1
         tmp2          = bpg*Rdu%rij(j,i)/(1.d0+bpg*Rdu%rij(j,i))
         rijbar(1,i,j) = tmp2
         do v = 2,mmax
            rijbar(v,i,j)  = rijbar(v-1,i,j)*tmp2
         enddo
      enddo
      do j = i+1,ne
         tmp2          = bpg*Rdu%rij(i,j)/(1.d0+bpg*Rdu%rij(i,j))
         rijbar(1,j,i) = tmp2
         do v = 2,mmax
            rijbar(v,j,i)  = rijbar(v-1,j,i)*tmp2
         enddo
      enddo

c---------main part-----------------------------------------------------

      ! recalculate the ie terms   
      Rdu%Fij(i,1:i-1) = 0.d0
      Rdu%Fij(i+1:ne,i) = 0.d0
      Rdu%Gi(i) = 0.d0

      do k = 1,k0

         ! Electron-Electron-Correlation Terms------------------------
         if (tcontrol(k,1)==1) then
            v = tcontrol(k,5)
            do j = 1,i-1
               Rdu%Fij(i,j) = Rdu%Fij(i,j) + cjas(k)*rijbar(v,i,j)
               !!!write(iul,'(a,2i3,g20.10)') 'DBG:---->',j,i,rijbar(v,i,j)               
            enddo
            do j = i+1,ne
               Rdu%Fij(j,i) = Rdu%Fij(j,i) + cjas(k)*rijbar(v,j,i)
            enddo
         endif

         ! satisfy cusp for same and different spin
         if (tcontrol(k,1)==7) then
               !!!write(iul,'(a)') 'DBG: OH OH OH!'               
            v = 1                         ! only linear r_ij term!
            if (i <= nalpha) then
               do j = 1,i-1
                  !uu(k) =  uu(k) + 0.5d0*rijbar(v,i,j)
                  Rdu%Fij(i,j) = Rdu%Fij(i,j) + cjas(k)*0.5d0*rijbar(v,i,j)
               enddo
               do j = i+1,nalpha
                  !uu(k) =  uu(k) + 0.5d0*rijbar(v,j,i)
                  Rdu%Fij(j,i) = Rdu%Fij(j,i) + cjas(k)*0.5d0*rijbar(v,j,i)
               enddo
               do j = nalpha+1,ne
                  !uu(k) =  uu(k) + rijbar(v,j,i)
                  Rdu%Fij(j,i) = Rdu%Fij(j,i) + cjas(k)*rijbar(v,j,i)
               enddo
            else
               do j = i,nalpha
                  !uu(k) =  uu(k) + 0.5d0*rijbar(v,j,i)
                  Rdu%Fij(i,j) = Rdu%Fij(i,j) + cjas(k)*rijbar(v,i,j)
               enddo
               do j = nalpha+1,i-1
                  !uu(k) =  uu(k) + 0.5d0*rijbar(v,j,i)
                  Rdu%Fij(i,j) = Rdu%Fij(i,j) + cjas(k)*0.5d0*rijbar(v,i,j)
               enddo
               do j = i+1,ne
                  !uu(k) =  uu(k) + 0.5d0*rijbar(v,j,i)
                  Rdu%Fij(j,i) = Rdu%Fij(j,i) + cjas(k)*0.5d0*rijbar(v,j,i)
               enddo
            endif
         endif

         ! Electron-Nucleus-Correlation-Terms----------------------------------

         if (tcontrol(k,1)==2) then
            s = tcontrol(k,6)
            do a = tcontrol(k,3),tcontrol(k,4)
               !uu(k) = uu(k) + t*raibar(s,a,i)
               Rdu%Gi(i) = Rdu%Gi(i) + cjas(k)*t*raibar(s,a,i)
            enddo
         endif

         ! Electron-Electron-Nucleus-Correlation-Terms----------------------

         if (tcontrol(k,1)==3) then
            s  = tcontrol(k,6)
            do j = 1,i-1
              do a = tcontrol(k,3),tcontrol(k,4)
                !uu(k) =  uu(k) + raibar(s,a,i)*raibar(s,a,j)
                Rdu%Fij(i,j) = Rdu%Fij(i,j) + cjas(k)*raibar(s,a,i)*raibar(s,a,j)
                !!!write(iul,'(a,3i3,2g20.10)') 'DBG1:---->',j,i,a,raibar(s,a,i),raibar(s,a,j)          
              enddo
            enddo
            do j = i+1,ne
              do a = tcontrol(k,3),tcontrol(k,4)
                !uu(k) =  uu(k) + raibar(s,a,i)*raibar(s,a,j)
                Rdu%Fij(j,i) = Rdu%Fij(j,i) + cjas(k)*raibar(s,a,i)*raibar(s,a,j)
              enddo
            enddo
         endif
   
         if (tcontrol(k,1)==4) then
            v  = tcontrol(k,5)
            s  = tcontrol(k,6)
            do j = 1,i-1
              do a = tcontrol(k,3),tcontrol(k,4)
                tmp1 = raibar(s,a,i) + raibar(s,a,j)
                !uu(k) =uu(k) + tmp*rijbar(v,j,i)
                Rdu%Fij(i,j) = Rdu%Fij(i,j) + cjas(k)*tmp1*rijbar(v,i,j)
                !!!write(iul,'(a,3i3,3g20.10)') 'DBG2:---->',j,i,a,raibar(s,a,i),raibar(s,a,j),rijbar(v,i,j)
              enddo
            enddo
            do j = i+1,ne
              do a = tcontrol(k,3),tcontrol(k,4)
                tmp1 = raibar(s,a,i) + raibar(s,a,j)
                !uu(k) =uu(k) + tmp*rijbar(v,j,i)
                Rdu%Fij(j,i) = Rdu%Fij(j,i) + cjas(k)*tmp1*rijbar(v,j,i)
              enddo
            enddo
         endif
   
      enddo ! K-Schleife


      ! update of U
      Fsum = 0.d0
      do i=1,ie-1
         Fsum = Fsum + Rdu%Fij(ie,i)
      enddo
      do j=ie+1,ne
         Fsum = Fsum + Rdu%Fij(j,ie)
      enddo
      Rdu%U = Rdu%U + Rdu%Gi(i) + Fsum


C       write(iul,'(a)') 'DBG:JassmUpdate:NewFij'
C       do i=1,ne
C         write(iul,'(i3)') i
C         write(iul,'(5g20.10)') (Rdu%Fij(j,i),j=1,ne)
C       enddo
C       write(iul,'(a)') 'DBG:JassmUpdate:New:Fsum,Ginew,Unew'
C       write(iul,'(5g20.10)') Fsum,Rdu%Gi(ie),Rdu%U


      call Rdu%markJastrowValid()

      end subroutine jassmUpdate


c=======================================================================

c     ------------------------------------
      subroutine jassmUpdateWithUk(Rdu,ie)
c     ------------------------------------

c this version uses RdataUpdate to keep track of auxiliary data for one electron updates
c one electron update of U: restore old electron (ieold) data when new electron (ie /= ieold)
c no electron derivatives here!
c but parameter deriv Uk updates

c rijbar:  (rij/(1+b*rij))**m
c raibar:  (rai/(1+a*rai))**l

      type(RdataUpdate), intent(inout) :: Rdu    ! data structure for electron update calculations
      integer, intent(in)              :: ie     ! update for electron ie

      integer i,j,a,k,k1,k2,km1
      integer s,v
      !real*8  uu(k0)
      real*8  tmp,tmp1,tmp2,t,Fsum
      real*8 :: Fksum(size(Rdu%Fijk,3))

      ! cusp enforcement (s. jassmall/Init) should not be necessary

      if (ie /= Rdu%ieJasold) then
         if (Rdu%ieJasold /= 0) then
           ! necessary to restore one-electron terms for previous electron ieold
           ! also restore U to original U
           Rdu%Fij(Rdu%ieJasold,1:Rdu%ieJasold-1) = Rdu%Fijold(1:Rdu%ieJasold-1) 
           Rdu%Fij(Rdu%ieJasold+1:ne,Rdu%ieJasold) = Rdu%Fijold(Rdu%ieJasold+1:ne)
           Rdu%Fijk(Rdu%ieJasold,1:Rdu%ieJasold-1,:) = Rdu%Fijkold(1:Rdu%ieJasold-1,:) 
           Rdu%Fijk(Rdu%ieJasold+1:ne,Rdu%ieJasold,:) = Rdu%Fijkold(Rdu%ieJasold+1:ne,:)
           Rdu%Gi(Rdu%ieJasold) = Rdu%Giold
           Rdu%Gki(:,Rdu%ieJasold) = Rdu%Gkiold(:)
           raibar(1:Rdu%endim1,1:Rdu%endim2,Rdu%ieJasold) = Rdu%raibarOld(:,:)
           Rdu%U = Rdu%U0
           Rdu%Uk = Rdu%Uk0
         endif
         Rdu%Fijold(1:ie-1) = Rdu%Fij(ie,1:ie-1)
         Rdu%Fijold(ie+1:ne) = Rdu%Fij(ie+1:ne,ie)
         Rdu%Fijkold(1:ie-1,:) = Rdu%Fijk(ie,1:ie-1,:)
         Rdu%Fijkold(ie+1:ne,:) = Rdu%Fijk(ie+1:ne,ie,:)
         Rdu%Giold = Rdu%Gi(ie)
         Rdu%Gkiold(:) = Rdu%Gki(:,ie)
         Rdu%raibarOld(:,:) = raibar(1:Rdu%endim1,1:Rdu%endim2,ie)
         Rdu%ieJasold = ie
      ! else same electron as before, simply update ie
      endif

      !!!write(iul,'(ai3,g20.10)') 'DBG:jassmUpdate:ie,Uold:',ie,Rdu%U

      ! subtract current=previous values the three Jastrow terms
      Fsum = 0.d0
      do i=1,ie-1
         Fsum = Fsum + Rdu%Fij(ie,i)
      enddo
      do j=ie+1,ne
         Fsum = Fsum + Rdu%Fij(j,ie)
      enddo
C       Fksum = 0.d0
C       do i=1,ie-1
C          Fksum(:) = Fksum(:) + Rdu%Fijk(ie,i,:)
C       enddo
C       do j=ie+1,ne
C          Fksum(:) = Fksum(:) + Rdu%Fijk(j,ie,:)
C       enddo

      Rdu%U  = Rdu%U  - Rdu%Gi(ie) - Fsum
C      Rdu%Uk = Rdu%Uk - Rdu%Gki(1:k0-1,ie) - Fksum

      !!!write(iul,'(a)') 'DBG:Fsumold,Giold,Uold-Giold-Fsumold'
      !!!write(iul,'(5g20.10)') Fsum,Rdu%Gi(ie),Rdu%U


c---------Precalculation of Terms-------------------------------

      !uu = 0
      t    = ne-1

      i = ie ! only ie !

      do a = 1,nclast
         tmp1          = apg(atoms(a)%sa)*Rdu%rai(a,i)/(1.d0+apg(atoms(a)%sa)*Rdu%rai(a,i))
         raibar(1,a,i) = tmp1
         do s = 2,lmax
            raibar(s,a,i)  = raibar(s-1,a,i)*tmp1
         enddo
      enddo
      do j = 1,i-1
         tmp2          = bpg*Rdu%rij(j,i)/(1.d0+bpg*Rdu%rij(j,i))
         rijbar(1,i,j) = tmp2
         do v = 2,mmax
            rijbar(v,i,j)  = rijbar(v-1,i,j)*tmp2
         enddo
      enddo
      do j = i+1,ne
         tmp2          = bpg*Rdu%rij(i,j)/(1.d0+bpg*Rdu%rij(i,j))
         rijbar(1,j,i) = tmp2
         do v = 2,mmax
            rijbar(v,j,i)  = rijbar(v-1,j,i)*tmp2
         enddo
      enddo

c---------main part-----------------------------------------------------

      ! recalculate the ie terms   
      Rdu%Fij(i,1:i-1) = 0.d0
      Rdu%Fij(i+1:ne,i) = 0.d0
      !!!Rdu%Fijk(i,1:i-1,:) = 0.d0
      !!!Rdu%Fijk(i+1:ne,i,:) = 0.d0
      Rdu%Gi(i) = 0.d0
      !!!Rdu%Gki(:,i) = 0.d0

      k1 = 0
      k2 = 0

      do k = 1,k0

         km1 = k - 1

         ! Electron-Electron-Correlation Terms------------------------
         if (tcontrol(k,1)==1) then
            v = tcontrol(k,5)
            if (v==1) then   ! cusp term, no parameter
               do j = 1,i-1
                  Rdu%Fij(i,j)    = Rdu%Fij(i,j) + cjas(k)*rijbar(v,i,j)
                  !!!write(iul,'(a,2i3,g20.10)') 'DBG:---->',j,i,rijbar(v,i,j)               
               enddo
               do j = i+1,ne
                  Rdu%Fij(j,i)    = Rdu%Fij(j,i) + cjas(k)*rijbar(v,j,i)
               enddo
            else
               k2 = k2 + 1
               do j = 1,i-1
                  Rdu%Fij(i,j)     = Rdu%Fij(i,j) + cjas(k)*rijbar(v,i,j)
                  Rdu%Uk(km1) = Rdu%Uk(km1) - Rdu%Fijk(i,j,k2) + rijbar(v,i,j)
                  Rdu%Fijk(i,j,k2) = rijbar(v,i,j)
                  !!!write(iul,'(a,2i3,g20.10)') 'DBG:---->',j,i,rijbar(v,i,j)               
               enddo
               do j = i+1,ne
                  Rdu%Fij(j,i)     = Rdu%Fij(j,i) + cjas(k)*rijbar(v,j,i)
                  Rdu%Uk(km1) = Rdu%Uk(km1) - Rdu%Fijk(j,i,k2) + rijbar(v,j,i)
                  Rdu%Fijk(j,i,k2) = rijbar(v,j,i)
               enddo
            endif

         endif

         ! satisfy cusp for same and different spin
         if (tcontrol(k,1)==7) then
               !!!write(iul,'(a)') 'DBG: OH OH OH!'               
            v = 1                         ! only linear r_ij term!
            if (i <= nalpha) then
               do j = 1,i-1
                  !uu(k) =  uu(k) + 0.5d0*rijbar(v,i,j)
                  Rdu%Fij(i,j)    = Rdu%Fij(i,j) + cjas(k)*0.5d0*rijbar(v,i,j)
               enddo
               do j = i+1,nalpha
                  !uu(k) =  uu(k) + 0.5d0*rijbar(v,j,i)
                  Rdu%Fij(j,i)    = Rdu%Fij(j,i) + cjas(k)*0.5d0*rijbar(v,j,i)
               enddo
               do j = nalpha+1,ne
                  !uu(k) =  uu(k) + rijbar(v,j,i)
                  Rdu%Fij(j,i)    = Rdu%Fij(j,i) + cjas(k)*rijbar(v,j,i)
               enddo
            else
               do j = i,nalpha
                  !uu(k) =  uu(k) + 0.5d0*rijbar(v,j,i)
                  Rdu%Fij(i,j)    = Rdu%Fij(i,j) + cjas(k)*rijbar(v,i,j)
               enddo
               do j = nalpha+1,i-1
                  !uu(k) =  uu(k) + 0.5d0*rijbar(v,j,i)
                  Rdu%Fij(i,j)    = Rdu%Fij(i,j) + cjas(k)*0.5d0*rijbar(v,i,j)
               enddo
               do j = i+1,ne
                  !uu(k) =  uu(k) + 0.5d0*rijbar(v,j,i)
                  Rdu%Fij(j,i)    = Rdu%Fij(j,i) + cjas(k)*0.5d0*rijbar(v,j,i)
               enddo
            endif
         endif

         ! Electron-Nucleus-Correlation-Terms----------------------------------

         if (tcontrol(k,1)==2) then
            s = tcontrol(k,6)
            k1 = k1 + 1
            Rdu%Uk(km1) = Rdu%Uk(km1) - Rdu%Gki(k1,i)
            Rdu%Gki(k1,i) = 0 
            do a = tcontrol(k,3),tcontrol(k,4)
               tmp = t*raibar(s,a,i)
               Rdu%Gi(i)    = Rdu%Gi(i) + cjas(k)*tmp
               Rdu%Gki(k1,i) = Rdu%Gki(k1,i) + tmp
            enddo
            Rdu%Uk(km1) = Rdu%Uk(km1) + Rdu%Gki(k1,i)
         endif

         ! Electron-Electron-Nucleus-Correlation-Terms----------------------

         if (tcontrol(k,1)==3) then
            s  = tcontrol(k,6)
            k2 = k2 + 1
            do j = 1,i-1
               Rdu%Uk(km1) = Rdu%Uk(km1) - Rdu%Fijk(i,j,k2)
               Rdu%Fijk(i,j,k2) = 0 
               do a = tcontrol(k,3),tcontrol(k,4)
                  !uu(k) =  uu(k) + raibar(s,a,i)*raibar(s,a,j)
                  tmp = raibar(s,a,i)*raibar(s,a,j)
                  Rdu%Fij(i,j)    = Rdu%Fij(i,j) + cjas(k)*tmp
                  Rdu%Fijk(i,j,k2) = Rdu%Fijk(i,j,k2) + tmp
                  !!!write(iul,'(a,3i3,2g20.10)') 'DBG1:---->',j,i,a,raibar(s,a,i),raibar(s,a,j)          
               enddo
               Rdu%Uk(km1) = Rdu%Uk(km1) + Rdu%Fijk(i,j,k2)
            enddo
            do j = i+1,ne
               Rdu%Uk(km1) = Rdu%Uk(km1) - Rdu%Fijk(j,i,k2)
               Rdu%Fijk(j,i,k2) = 0 
               do a = tcontrol(k,3),tcontrol(k,4)
                  !uu(k) =  uu(k) + raibar(s,a,i)*raibar(s,a,j)
                  tmp = raibar(s,a,i)*raibar(s,a,j)
                  Rdu%Fij(j,i)    = Rdu%Fij(j,i) + cjas(k)*tmp
                  Rdu%Fijk(j,i,k2) = Rdu%Fijk(j,i,k2) + tmp
              enddo
               Rdu%Uk(km1) = Rdu%Uk(km1) + Rdu%Fijk(j,i,k2)
            enddo
         endif
   
         if (tcontrol(k,1)==4) then
            v  = tcontrol(k,5)
            s  = tcontrol(k,6)
            k2 = k2 + 1
            do j = 1,i-1
               Rdu%Uk(km1) = Rdu%Uk(km1) - Rdu%Fijk(i,j,k2)
               Rdu%Fijk(i,j,k2) = 0 
               do a = tcontrol(k,3),tcontrol(k,4)
                  tmp1 = rijbar(v,i,j)*(raibar(s,a,i) + raibar(s,a,j))
                  Rdu%Fij(i,j)    = Rdu%Fij(i,j) + cjas(k)*tmp1
                  Rdu%Fijk(i,j,k2) = Rdu%Fijk(i,j,k2) + tmp1
                  !!!write(iul,'(a,3i3,3g20.10)') 'DBG2:---->',j,i,a,raibar(s,a,i),raibar(s,a,j),rijbar(v,i,j)
               enddo
               Rdu%Uk(km1) = Rdu%Uk(km1) + Rdu%Fijk(i,j,k2)
            enddo
            do j = i+1,ne
               Rdu%Uk(km1) = Rdu%Uk(km1) - Rdu%Fijk(j,i,k2)
               Rdu%Fijk(j,i,k2) = 0 
               do a = tcontrol(k,3),tcontrol(k,4)
                  tmp1 = rijbar(v,j,i)*(raibar(s,a,i) + raibar(s,a,j))
                  Rdu%Fij(j,i)    = Rdu%Fij(j,i) + cjas(k)*tmp1
                  Rdu%Fijk(j,i,k2) = Rdu%Fijk(j,i,k2) + tmp1
               enddo
               Rdu%Uk(km1) = Rdu%Uk(km1) + Rdu%Fijk(j,i,k2)
            enddo
         endif
   
      enddo ! K-Schleife


      call assert(k1 == size(Rdu%Gki,1) .and. k2 == size(Rdu%Fijk,3) .and. k1+k2==k0-1, 
     &    "jassmUpdateWithUk: k1,k2 illegal values!")


      ! update of U
      Fsum = 0.d0
      do i=1,ie-1
         Fsum = Fsum + Rdu%Fij(ie,i)
      enddo
      do j=ie+1,ne
         Fsum = Fsum + Rdu%Fij(j,ie)
      enddo
C       Fksum = 0.d0
C       do i=1,ie-1
C          Fksum(:) = Fksum(:) + Rdu%Fijk(ie,i,1:k0-1)
C       enddo
C       do j=ie+1,ne
C          Fksum(:) = Fksum(:) + Rdu%Fijk(j,ie,1:k0-1)
C       enddo

      Rdu%U = Rdu%U + Rdu%Gi(i) + Fsum
C      Rdu%Uk(:) = Rdu%Uk(:) + Rdu%Gki(1:k0-1,ie) + Fksum(:)


C       write(iul,'(a)') 'DBG:JassmUpdate:NewFij'
C       do i=1,ne
C         write(iul,'(i3)') i
C         write(iul,'(5g20.10)') (Rdu%Fij(j,i),j=1,ne)
C       enddo
C       write(iul,'(a)') 'DBG:JassmUpdate:New:Fsum,Ginew,Unew'
      !!!write(iul,'(a,g20.10)') 'DBG:jassmupdate:U,Uk:',Rdu%U
      !!!write(iul,'(5g20.10)') Rdu%Uk(:)
C       write(iul,'(a)') 'Fijk:'
C       do k=1,k0-1
C         do i=1,ne
C           write(iul,'(2i4,4x,g20.10)') k,i,Rdu%Gki(k,i)
C           do j=i+1,ne
C             write(iul,'(3i4,g20.10)') k,i,j,Rdu%Fijk(j,i,k)
C           enddo
C         enddo
C       enddo


      call Rdu%markJastrowValid()

      end subroutine jassmUpdateWithUk


c=======================================================================


c-----------------------------------------------------------------------
      subroutine jassmp(init,ie,rai,rij,ju)
c-----------------------------------------------------------------------
c
c This is essentially the subroutine jassmone without the gradient and laplacian parts.
c It is needed for ECP localisation.

c rai,rij         : el-nucl(a) and el-el distances
c ncenter,nscenter: (params) # nuclei, and # nuclei not symmetry equiv.
c sa(a)           : symmetry equiv. nuclei
c bp(s),ap(s,a)   : Jastrow coefficients in r/(1+br) for rij and rai, resp.
c cjas(k)    : lin. Koeff. c_k im Jastrow-Exponent
c ijuu(i,j)  : Uij's
c uu1(i)     :[sum_j\Uij +sum_j\Uji +sum_(j,a)\Uija +sum_(j,a)\Ujia]
c uu2(i)     : sum_a Uai
c
c rijbar:  (rij/(1+brij))**2
c raibar:  (rai/(1+arai))**2

      logical, intent(in) ::  init   ! .true. for initialization
      integer, intent(in) ::  ie     ! electron with new position
      real*8, intent(in)  ::  rai(:,:),rij(:,:)   ! current distances
      real*8, intent(out) ::  ju     ! returns U rather than exp(U)

      integer i,j,a,k,line
      integer s,v
      real*8  tmp,tmp1,tmp2,t,jexp
      real*8  huu
      real*8  kuu

      save t,jexp

      if(.not. init) call assert(ie > 0,'jassmp: positive ie required')

      cjas(1) = 0.5d0/bpg

c--------Initialisation------------------------------------------------
      if (init) then
c-----------------------------------------------------------------------

      jexp = 0d0
      uu1(1:ne) = 0d0
      uu2(1:ne) = 0d0

c--------precalculation of terms-------------------------------------

      t    = ne-1

ccc necessary since rijbar is defined differently in jassmall
ccc i and j are interchanged. --> This subroutine should be rewritten
ccc with interchanged indices such that it can share raibar aand rijbar
ccc with jassmall

      do i = 1,ne
         do a = 1,nclast
            tmp1     = 1/(1+apg(atoms(a)%sa)*rai(a,i))
            tmp2     = tmp1*apg(atoms(a)%sa)*rai(a,i)
            raibar(1,a,i) = tmp2
            do s = 2,lmax
              raibar(s,a,i) = tmp2*raibar(s-1,a,i)
            enddo
         enddo
         do j = i+1,ne
            tmp1     = 1/(1+bpg*rij(i,j))
            tmp2     = tmp1*bpg*rij(i,j)
            rijbar(1,i,j) = tmp2
            do v = 2,mmax
              rijbar(v,i,j) = tmp2*rijbar(v-1,i,j)
            enddo
         enddo
      enddo

      do i=1,ne
        do j=i+1,ne
          ijuu(i,j) = 0d0
        enddo
      enddo

c---------main part-----------------------------------------------------
      do k = 1,k0
      do i = 1,ne

c---------Electron-Nucleus-Correlation-Terms----------------------------------

c      if (tcontrol(k,1).eq.3) then
      if (tcontrol(k,1).eq.2) then
      s = tcontrol(k,6)
      do a = tcontrol(k,3),tcontrol(k,4)
        uu2(i) = uu2(i)+cjas(k)*t*raibar(s,a,i) !\sum_i^{ne}\sum_{i<j}^{ne}{r_i + r_j} =
                                                !(ne-1)\sum_i^{ne}{r_i}
      enddo
      endif

      do j = i+1,ne

c        ijuu(i,j)    = 0d0

c---------Electron-Electron-Correlation-Terms---------------------------
c        if (tcontrol(k,1).eq.2) then
        if (tcontrol(k,1).eq.1) then
          v = tcontrol(k,5)
          ijuu(i,j) = ijuu(i,j) + cjas(k)*rijbar(v,i,j)
        endif

c---------Electron-Electron-Nucleus-Correlation-Terms----------------------
c        if (tcontrol(k,1).eq.4) then
        if (tcontrol(k,1).eq.3) then
        s  = tcontrol(k,6)
          do a = tcontrol(k,3),tcontrol(k,4)
c            ijuu(i,j) = ijuu(i,j)+2*cjas(k)*raibar(s,a,i)*raibar(s,a,j)
            ijuu(i,j) = ijuu(i,j)+cjas(k)*raibar(s,a,i)*raibar(s,a,j)
          enddo
        endif
c-----------------------------------------------------------------------

c        if (tcontrol(k,1).eq.5) then
        if (tcontrol(k,1).eq.4) then
        v = tcontrol(k,5)
        s = tcontrol(k,6)
          do a = tcontrol(k,3),tcontrol(k,4)
            tmp1 = raibar(s,a,i)+raibar(s,a,j)
            ijuu(i,j) = ijuu(i,j)+cjas(k)*tmp1*rijbar(v,i,j)
          enddo
        endif
c-----------------------------------------------------------------------

      enddo ! j-Schleife
      enddo ! i-Schleife
      enddo ! K-Schleife

c---------Calculation of the Exponent of the Correlation-Term-----------
c---------this is necessary because jexp not equal sum_i\[uu1(i)+uu2(i)]--

      do i = 1,ne
        jexp = jexp + uu2(i)
        do j = i+1,ne
          jexp = jexp + ijuu(i,j)
        enddo
      enddo

c---------Calculation of Terms for the single Electron------------------

      do line = 1,ne
        j = line
        do i = 1,line-1
            uu1(j) = uu1(j) + ijuu(i,j)
        enddo
        i = line
        do j = line+1,ne
            uu1(i) = uu1(i) + ijuu(i,j)
        enddo
      enddo

      ju = jexp

C       write(iul,'(a,g20.10)') 'DBG:jassmp:init',ju
C       write(iul,'(a)') 'DBG:uu1:ee/een'
C       do i=1,ne
C         write(iul,'(i5)') i
C         write(iul,'(5g20.10)') (ijuu(i,j),j=1,ne)
C       enddo
C       write(iul,'(a)') 'DBG:uu2:en'
C       write(iul,'(5g20.10)') (uu2(i),i=1,ne)

c-------------one-electron-move-----------------------------------------
      else  ! init = .false.
c-----------------------------------------------------------------------

      oju = jexp
      do i = 1,ne
        ouu1(i) = uu1(i)
      enddo

      ouu2(ie) = uu2(ie)
      do a = 1,nclast
         do s = 1,lmax
            oraibar(s,a,ie) = raibar(s,a,ie)
         enddo
      enddo

      do i = 1,ie-1
         oijuu(i,ie) = ijuu(i,ie)
         do v = 1,mmax
            orijbar(v,i,ie) = rijbar(v,i,ie)
         enddo
      enddo
      do i = ie+1,ne
         oijuu(ie,i) = ijuu(ie,i)
         do v = 1,mmax
            orijbar(v,ie,i) = rijbar(v,ie,i)
         enddo
      enddo

c--------Calculation of Correlation-Term--------------------------------

      !!!write(iul,'(a,i3,g20.10)') 'DBG:jassmp:update',ie,jexp

      jexp = jexp - uu1(ie) - uu2(ie) !substract contributions of all terms
                                      !arising from the old electron position

      !!!write(iul,'(a)') 'DBG:uu1:old,uu2:old,jexp-uu1-uu2'
      !!!write(iul,'(5g20.10)') uu1(ie),uu2(ie),jexp


      uu1(ie) = 0d0
      uu2(ie) = 0d0

c---------precalculations of terms--------------------------------------
      i = ie
      do a = 1,nclast
         tmp1     = 1/(1+apg(atoms(a)%sa)*rai(a,i))
         tmp2     = tmp1*apg(atoms(a)%sa)*rai(a,i)
         raibar(1,a,i) = tmp2
         do s = 2,lmax
           raibar(s,a,i) = tmp2*raibar(s-1,a,i)
         enddo
      enddo

      j = ie
      do i = 1,ie-1
         tmp1     = 1/(1+bpg*rij(i,j))
         tmp2     = tmp1*bpg*rij(i,j)
         rijbar(1,i,j) = tmp2
         do v = 2,mmax
           rijbar(v,i,j) = tmp2*rijbar(v-1,i,j)
         enddo
      enddo
      i = ie
      do j = ie+1,ne
         tmp1     = 1/(1+bpg*rij(i,j))
         tmp2     = tmp1*bpg*rij(i,j)
         rijbar(1,i,j) = tmp2
         do v = 2,mmax
           rijbar(v,i,j) = tmp2*rijbar(v-1,i,j)
         enddo
      enddo


c---------main part-----------------------------------------------------

      huu    = 0

c-----change of ...(i,ie) terms for one electron move-------------------

      j = ie

      do i = 1,ie-1
         !
         uu1(i) = uu1(i) - ijuu(i,j)
         ijuu(i,j) = 0d0

         do k = 1,k0
           kuu      = 0d0

c---------Electron-Electron-Correlation-Terms---------------------------
c           if (tcontrol(k,1).eq.2) then
           if (tcontrol(k,1).eq.1) then
             v = tcontrol(k,5)
             kuu = rijbar(v,i,j)
             !!!write(iul,'(a,2i3,g20.10)') 'DBG--->',i,ie,rijbar(v,i,j)
           endif


c---------Electron-Electron-Nucleus-Correlation-Terms--------------------
c           if (tcontrol(k,1).eq.4) then
           if (tcontrol(k,1).eq.3) then
             s  = tcontrol(k,6)
             do a = tcontrol(k,3),tcontrol(k,4)
c               kuu = kuu + 2*raibar(s,a,i)*raibar(s,a,j)
               kuu = kuu + raibar(s,a,i)*raibar(s,a,j)
               !!!write(iul,'(a,3i3,2g20.10)') 'DBG1--->',i,ie,a,raibar(s,a,i),raibar(s,a,j)
             enddo
           endif
c-----------------------------------------------------------------------

c           if (tcontrol(k,1).eq.5) then
           if (tcontrol(k,1).eq.4) then
             v  = tcontrol(k,5)
             s  = tcontrol(k,6)
             do a = tcontrol(k,3),tcontrol(k,4)
               tmp1 = raibar(s,a,i)+raibar(s,a,j)
               kuu  = kuu + tmp1*rijbar(v,i,j)
               !!!write(iul,'(a,3i3,3g20.10)') 'DBG2--->',i,ie,a,raibar(s,a,i),raibar(s,a,j),rijbar(v,i,j)
             enddo
           endif
c-----------------------------------------------------------------------

c-----------------------------------------------------------------------

           ijuu(i,j) = ijuu(i,j) + cjas(k)*kuu
         enddo ! K-loop
         uu1(i) = uu1(i) + ijuu(i,j)
         huu    = huu    + ijuu(i,j)
      enddo  ! I-loop


c-------change of ...(ie,j) terms for one electron move----------------

      i = ie

c---------Electron-Nucleus-Correlation-Terms---------------------------
      do k = 1,k0
c        if (tcontrol(k,1).eq.3) then
        if (tcontrol(k,1).eq.2) then
          s   = tcontrol(k,6)
          tmp1 = 0d0
          do a = tcontrol(k,3),tcontrol(k,4)
            tmp1 = tmp1 +raibar(s,a,i)
          enddo
          uu2(i) = uu2(i)+cjas(k)*t*tmp1
        endif
      enddo

      do j = ie+1,ne
         uu1(j)    = uu1(j) - ijuu(i,j)
         ijuu(i,j) = 0d0

         do k = 1,k0
           kuu      = 0d0

c---------Electron-Electron-Correlation-Terms------------------------------
c           if (tcontrol(k,1).eq.2) then
           if (tcontrol(k,1).eq.1) then
             v = tcontrol(k,5)
             kuu = rijbar(v,i,j)
           endif

c---------Electron-Electron-Nucleus-Correlation-Terms----------------------
c           if (tcontrol(k,1).eq.4) then
           if (tcontrol(k,1).eq.3) then
             s = tcontrol(k,6)
             do a = tcontrol(k,3),tcontrol(k,4)
c               kuu = kuu + 2*raibar(s,a,i)*raibar(s,a,j)
               kuu = kuu + raibar(s,a,i)*raibar(s,a,j)
             enddo
           endif
c-----------------------------------------------------------------------

c           if (tcontrol(k,1).eq.5) then
           if (tcontrol(k,1).eq.4) then
             v = tcontrol(k,5)
             s = tcontrol(k,6)
             do a = tcontrol(k,3),tcontrol(k,4)
               tmp1 = raibar(s,a,i)+raibar(s,a,j)
               kuu = kuu + tmp1*rijbar(v,i,j)
             enddo
           endif

c-----------------------------------------------------------------------

           ijuu(i,j) = ijuu(i,j) + cjas(k)*kuu

         enddo ! K-loop

         uu1(j) = uu1(j) + ijuu(i,j)
         huu    = huu    + ijuu(i,j)


      enddo ! J-loop

      uu1(ie) = huu


c--------construction of the correlation-term --------------------------

      jexp = jexp + uu1(ie) + uu2(ie)

C       write(iul,'(a)') 'DBG:newijuu'
C       do i=1,ne
C         write(iul,'(i5)') i
C         write(iul,'(5g20.10)') (ijuu(i,j),j=1,ne)
C       enddo
C       write(iul,'(a)') 'DBG:uu1:new,uu2:new,jexp+uu1+uu2'
C       write(iul,'(a,5g20.10)') 'DBG:jassmp:',jexp

      ju = jexp

ccc restore rijbar and riabar matrices (this differs from jassmone because
ccc all 'pseudo-OEM' in connection with the ECP localisation get rejected.)
ccc --> it makes sense to restore the matrices directly in this subroutine

      jexp = oju

      do i = 1,ne
        uu1(i) = ouu1(i)
      enddo

      uu2(ie) = ouu2(ie)
      do a = 1,nclast
         do s = 1,lmax
            raibar(s,a,ie) = oraibar(s,a,ie)
         enddo
      enddo

      do i = 1,ie-1
         ijuu(i,ie) = oijuu(i,ie)
         do v = 1,mmax
            rijbar(v,i,ie) = orijbar(v,i,ie)
         enddo
      enddo
      do i = ie+1,ne
         ijuu(ie,i) = oijuu(ie,i)
         do v = 1,mmax
            rijbar(v,ie,i) = orijbar(v,ie,i)
         enddo
      enddo

c-----------------------------------------------------------------------
      endif
c-----------------------------------------------------------------------

      end subroutine jassmp



c=======================================================================

c     ----------------------------------------------------------------
      subroutine jassmallOnlyUk(x,y,z,rai,rij,optType,ju,uuk)
c     ----------------------------------------------------------------

c rijbar:  (rij/(1+brij))**m
c raibar:  (rai/(1+arai))**l
c derivi:  term used for calculation of gradient of (rij/(1+brij))**m upper triangle,
c          term used for calculation of laplacian of (rij/(1+brij))**m lower triangle
c fderiva: term used for calculation of gradient of (rai/(1+arai))**l ,
c sdiriva: term used for calculation of laplacian of (rai/(1+arai))**l


      real*8, intent(in)            ::  x(:),y(:),z(:)     ! cartesian coords
      real*8, intent(in)            ::  rai(:,:),rij(:,:)  ! distances r_ai, r_ij
      character(len=*)              ::  optType            ! parameter optimization
      real*8, intent(out)           ::  ju                 ! in exp(U)
      real*8, optional, intent(out) ::  uuk(:)
      integer i,j,a,k, k1
      integer s,v
      real*8  uu(kmax)
      real*8  xij(nmax,nmax),yij(nmax,nmax),zij(nmax,nmax)
      real*8  xai(amax,nmax),yai(amax,nmax),zai(amax,nmax)
      real*8  tmp,tmp1,tmp2,t

      real*8 Fijk(ne,ne,0:k0-1), Gki(0:k0-1,ne)

      Fijk = 0.d0
      Gki = 0.d0


c---------Enforce cusp condition for unlike electrons--------------
c         ((here to allow modification of bpg for optimization))
      if (JPD_enforceECusp) then
         cjas(1) = 0.5d0/bpg
      endif

c---------Calculation of Jastrow factor and its deriviatives-------
      do k = 1,k0
         uu(k) = 0
      enddo

c---------Precalculation of Terms-------------------------------
      t    = ne-1

      do i = 1,ne
         do a = 1,nclast
            xai(a,i) = x(i)-atoms(a)%cx
            yai(a,i) = y(i)-atoms(a)%cy
            zai(a,i) = z(i)-atoms(a)%cz
            tmp1           = 1/(1+apg(atoms(a)%sa)*rai(a,i))
            tmp2           = tmp1*apg(atoms(a)%sa)*rai(a,i)
            raibar(1,a,i)  = tmp2
            do s = 2,lmax
              raibar(s,a,i)  = tmp2*raibar(s-1,a,i)
            enddo
         enddo
         do j = i+1,ne
            xij(j,i) = x(i)-x(j)
            yij(j,i) = y(i)-y(j)
            zij(j,i) = z(i)-z(j)
            tmp1          = 1/(1+bpg*rij(i,j))
            tmp2          = tmp1*bpg*rij(i,j)
            rijbar(1,j,i) = tmp2
            do v = 2,mmax
              rijbar(v,j,i)  = tmp2*rijbar(v-1,j,i)
            enddo
         enddo
      enddo

c---------main part-----------------------------------------------------

      do k = 1,k0

      k1 = k-1

c------------Electron-Electron-Correlation Terms------------------------
      if (tcontrol(k,1)==1) then
      v = tcontrol(k,5)
      do i = 1,ne
         do j = i+1,ne
            uu(k)          =  uu(k)          + rijbar(v,j,i)
            Fijk(j,i,k1) = Fijk(j,i,k1) + rijbar(v,j,i)
            !!write(iul,'(a,i5,g20.10)') 'DBG:k,uuk:',k,uu(k)            
         enddo
      enddo
      endif


      ! satisfy cusp for same and different spin
      if (tcontrol(k,1)==7) then
      v = 1                         ! only linear r_ij term!
      do i = 1,nalpha
         do j = i+1,nalpha
            uu(k)          =  uu(k)          + 0.5d0*rijbar(v,j,i)
         enddo
         do j = nalpha+1,ne
            uu(k)          =  uu(k)          + rijbar(v,j,i)
         enddo
      enddo
      do i = nalpha+1,ne
         do j = i+1,ne
            uu(k)          =  uu(k)          + 0.5d0*rijbar(v,j,i)
         enddo
      enddo
      endif


c---------Electron-Nucleus-Correlation-Terms----------------------------------

      if (tcontrol(k,1)==2) then
      s   = tcontrol(k,6)
      do i = 1,ne
         do a = tcontrol(k,3),tcontrol(k,4)
            uu(k)          =  uu(k)          + t*raibar(s,a,i)
            Gki(k1,i) = Gki(k1,i) + t*raibar(s,a,i)
            !!write(iul,'(a,i5,g20.10)') 'DBG:k,uuk:',k,uu(k)            
       enddo
      enddo
      endif

c---------Electron-Electron-Nucleus-Correlation-Terms----------------------

      if (tcontrol(k,1)==3) then
      s  = tcontrol(k,6)
      do i = 1,ne
        do j = i+1,ne
          do a = tcontrol(k,3),tcontrol(k,4)
            uu(k)          =  uu(k)       + raibar(s,a,i)*raibar(s,a,j)
            Fijk(j,i,k1) = Fijk(j,i,k1) + raibar(s,a,i)*raibar(s,a,j)
            !!write(iul,'(a,i5,g20.10)') 'DBG:k,uuk:',k,uu(k)            
          enddo
        enddo
      enddo
      endif
c-----------------------------------------------------------------------

      if (tcontrol(k,1)==4) then
      v  = tcontrol(k,5)
      s  = tcontrol(k,6)
      do i = 1,ne
        do j = i+1,ne
          do a = tcontrol(k,3),tcontrol(k,4)
            tmp1    = raibar(s,a,i)+raibar(s,a,j)
            uu(k)         =uu(k)         +tmp1*rijbar(v,j,i)
            Fijk(j,i,k1) = Fijk(j,i,k1) + tmp1*rijbar(v,j,i)
            !!write(iul,'(a,i5,g20.10)') 'DBG:k,uuk:',k,uu(k)            
          enddo
        enddo
      enddo
      endif
c-----------------------------------------------------------------------

      enddo ! K-Schleife
c--------Calculation of U and its Deriviatives--------------------------

      ju = 0d0

      do k=1,k0
        ju = ju + cjas(k)*uu(k)
      enddo


      if (optType == 'jastrow' .or. optType == 'jas+ci' .or. optType == 'jas+mo' .or. optType == 'jas+mo+ci') then
         call assert(allocated(uk),
     .               'jassmall: jastrowParamData not allocated')
         call assert(uk_params == k0-1,'jassmall: mismatch of # params')
         uk = uu(2:k0)
      endif
      if (present(uuk)) then
        do i=1,k0-1
           uuk(i)=uu(1+i)
        enddo

C         write(iul,'(a,g20.10)') 'DBG:jassmInitwpd:U,Uk',ju
C         write(iul,'(5g20.10)') uuk(1:k0-1)
C         write(iul,'(a)') 'Fijk:'
C         do k=1,k0-1
C           do i=1,ne
C             write(iul,'(2i4,4x,g20.10)') k,i,Gki(k,i)
C             do j=i+1,ne
C               write(iul,'(3i4,g20.10)') k,i,j,Fijk(j,i,k)
C             enddo
C           enddo
C         enddo

      endif
      end subroutine jassmallOnlyUk


c=======================================================================

c     -----------------------------------------------------
      subroutine getVectorLenSM(optMode, npJ1, npJ2, npJnl)
c     -----------------------------------------------------
      integer, intent(in)    :: optMode
      integer, intent(inout) :: npJ1     ! one-electron linear
      integer, intent(inout) :: npJ2     ! two-electron linear
      integer, intent(inout) :: npJnl    ! nonlinear
      integer k

      npJ1 = 0
      npJ2 = -1
      do k=1,k0
         if (tcontrol(k,1)==2) then
            npJ1 = npJ1 + 1
         else 
            npJ2 = npJ2 + 1
         endif
      enddo

      select case (optMode)
      case (1)   ! linear jastrow parameters
         npJnl = 0
      case (2)   ! linear and nonlinear Jastrow parameters
         npJnl = atoms(nclast)%sa + 1
      case default
         call abortp("getVectorLenSM: optMode not implemented")
      end select

      end subroutine getVectorLenSM


c=======================================================================

c     ---------------------------------
      subroutine getvectorSM(optMode,p)
c     ---------------------------------

      integer, intent(in) ::  optMode       ! optimization mode
      real*8, intent(out) ::  p(:)          ! parameter vector

      select case (optMode)
      case (1)   ! linear jastrow parameters
         if (size(p) /= k0-1) call abortp('getVectorSM: size mismatch')
         p(:) = cjas(2:k0)
      case (2)   ! linear and nonlinear Jastrow parameters
         if (size(p) /= k0+atoms(nclast)%sa)
     .      call abortp('getVectorSM: size mismatch')
         p(1:k0-1) = cjas(2:k0)
         p(k0) = bpg
         p(k0+1:k0+atoms(nclast)%sa) = apg(1:atoms(nclast)%sa)
      case default
         call abortp("getVectorSM: optMode not implemented")
      end select
      end subroutine getvectorSM


c=======================================================================

c     ---------------------------------
      subroutine putvectorSM(optMode,p)
c     ---------------------------------

      integer, intent(in) ::  optMode       ! optimization mode
      real*8, intent(in)  ::  p(:)          ! parameter vector

      select case (optMode)
      case (1)   ! linear jastrow parameters
         if (size(p) /= k0-1) call abortp('getVectorSM: size mismatch')
         cjas(2:k0) = p(:)
      case (2)   ! linear and nonlinear Jastrow parameters
         if (size(p) /= k0+atoms(nclast)%sa)
     .      call abortp('getVectorSM: size mismatch')
         cjas(2:k0) = p(1:k0-1)
         bpg = p(k0)
         apg(1:atoms(nclast)%sa) = p(k0+1:k0+atoms(nclast)%sa)
      case default
         call abortp("getVectorSM: optMode not implemented")
      end select
      end subroutine putvectorSM

c=======================================================================

c     ---------------------------------
      subroutine jas_diffeecusp_sm(pdiff_ee_cusp)
c     ---------------------------------
      logical, intent(in) :: pdiff_ee_cusp

      if (pdiff_ee_cusp .eqv. .TRUE.) then
         diffeecusp = .TRUE.
         if (tcontrol(1,1) == 0) then
            call abortp('jas_diffeecusp_sm: jastrow params not initialized')
         endif
         !change term type for first parameter (cusp)
         tcontrol(1,1) = 7
      else
         diffeecusp = .FALSE.
         if (tcontrol(1,1) == 0) then
            call abortp('jas_diffeecusp_sm: jastrow params not initialized')
         endif
         !change term type for first parameter (cusp)
         tcontrol(1,1) = 1
      endif

      end subroutine jas_diffeecusp_sm

      END MODULE jastrowSM
