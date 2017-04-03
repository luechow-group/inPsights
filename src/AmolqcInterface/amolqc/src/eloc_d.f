
c elocaldata contains values pertaining to local energy
c that are accessed throughout the code
c All "elXXX" data contain just the "current" values
c whatever that means.
c The "elXXX" values are correct *after* a call to eloc(x,y,z)
c for the positions x,y,z.

c Currently: global data
c TODO: access functions get/set
c


c     -----------------
      MODULE elocaldata
c     -----------------

      use wfdata
      use ecpModule
      implicit none


c------- new wf data type (to replace global data) ---------------

      type, public :: WFType
         !!! type(mdetType) :: mDet
         !!! type(jastrowType) :: jas
         type(EcpType) :: ecp
      end type WFType

      type, public :: WFMainResult
         real*8              :: Phi         ! determinantal part
         real*8              :: U           ! Jastrow U
         real*8              :: ELocal      ! E_local
         real*8              :: VPot        ! potential energy without nucl-nucl-repulsion
         real*8, allocatable :: xDrift(:)   ! nabla Psi / Psi,  Psi = Phi * exp(U)
         real*8, allocatable :: yDrift(:)
         real*8, allocatable :: zDrift(:)
      end type WFMainResult


c main return data from eloc, for several elec configurations
      real*8, allocatable :: elPhi(:)        ! Phi_G (determinantal part)
      real*8, allocatable :: elU(:)          ! Jastrow exponent U
      real*8, allocatable :: elEloc(:)       ! E_local(Psi_G)
      real*8, allocatable :: elVPot(:)       ! potential energy (without nucl-nucl-repulsion)
      real*8, allocatable :: elxDrift(:,:)   ! nabla(Psi_G) / Psi_G
      real*8, allocatable :: elyDrift(:,:)   ! nabla(Psi_G) / Psi_G
      real*8, allocatable :: elzDrift(:,:)   ! nabla(Psi_G) / Psi_G

c components of Eloc and Psi (only for one configuration)
      real*8, allocatable :: elEkini(:,:)  ! -0.5*nabla_i**2 Psi_G / Psi_G
      real*8, allocatable :: elUgrad(:,:)  ! grad(U)
      real*8, allocatable :: elULapl(:)    ! laplacian(U)
      real*8, allocatable :: elFgrad(:,:)  ! grad(Phi_G)
      real*8, allocatable :: elFLapl(:)    ! laplacian(Phi_G)
      real*8, allocatable :: elFLapli(:,:) ! laplacian(Phi_G(i,n))
      logical             :: elSingularity ! singular determinant

c additional stuff
      real*8              :: elECPPot,elECPPotNl ! ECP Potential,Non-local part

c E_local cutoff (Rothstein/Vrbik)
      real*8 :: mTauCutOff=0
      real*8 :: mCutOffFactor=0
      logical :: mElocCutOff=.false.

c max and actual number of simultaneous electron configurations
      integer :: mMaxElecConfigs = 0
      integer :: mElecConfigs = 0
c Arrays for energy partitioning
      real*8,allocatable   :: elEkin(:,:)
      real*8,allocatable   :: elVne(:,:,:)
      real*8,allocatable   :: elVee(:,:,:)

c stop when singularity in det is hit?
      logical :: mStopAtSingularity = .true.

      integer*8 :: mElocCut, mElocCutCount
      integer*8 :: mDriftCut, mDriftCutCount


      CONTAINS


c     -------------------------------
      subroutine eloc_initialize(nec)
c     -------------------------------

      integer, intent(in), optional :: nec
      integer nn,alstat
      if (present(nec)) then
         mMaxElecConfigs = nec
      else
         mMaxElecConfigs = 1
      endif
      nn = mMaxElecConfigs

      call assert(ne>0 .and. ncenter>0,
     .     'eloc_initialize: electron and nucleus number not set')
      if (allocated(elEloc) .and. size(elEloc) /= nn) then
         call eloc_deallocate()
      endif

      if (.not. allocated(elEloc)) then
            allocate(elPhi(nn),elU(nn),elEloc(nn),elVpot(nn),
     .             elxDrift(ne,nn),elyDrift(ne,nn),elzDrift(ne,nn),
     .             elEkini(ne,nn),elUgrad(3*ne,nn),elULapl(nn),
     .             elFgrad(3*ne,nn),elFLapl(nn),elFLapli(ne,nn),
     .             stat=alstat)
            if(do_epart) allocate(elEkin(ne,nn),elVne(ncenter,ne,nn)
     .                          ,elVee(ne,ne,nn))
            call assert(alstat==0,'eloc_initialize: allocation3 failed')
            elPhi = 0; elU = 0; elEloc = 0; elVpot = 0; elxDrift = 0
            elyDrift = 0; elzDrift = 0; elEkini = 0; elUgrad = 0
            elFgrad = 0; elFLapli = 0;
      endif
      end subroutine eloc_initialize

      logical function  eloc_getStopAtSingularity()
         eloc_getStopAtSingularity = mStopAtSingularity
      end function  eloc_getStopAtSingularity

      subroutine eloc_setStopAtSingularity(l)
         logical, intent(in) :: l
         mStopAtSingularity = l
      end subroutine eloc_setStopAtSingularity

      logical function eloc_isAtSingularity()
         eloc_isAtSingularity = elSingularity
      end function eloc_isAtSingularity

c     -----------------------------
      subroutine eloc_deallocate()
c     -----------------------------

      call assert(allocated(elEloc),
     .            'eloc: deallocation before allocation')
      deallocate(elPhi,elU,elEloc,elVpot,elxDrift,elyDrift,elzDrift,
     .              elEkini,elUgrad,elULapl,elFgrad,elFLapl,elFLapli)
      if(do_epart) deallocate(elEkin,elVne,elVee)
      end subroutine eloc_deallocate

c     -------------------------------------------------------
      subroutine eloc_getCurrentElocData(nec,phi,u,eloc,vpot,
     .           xdrift,ydrift,zdrift)
c     -------------------------------------------------------

      ! CHANGE to pointer!
      integer, intent(inout):: nec
      real*8, intent(inout) :: phi(:)
      real*8, intent(inout) :: u(:)
      real*8, intent(inout) :: eloc(:)
      real*8, intent(inout) :: vpot(:)
      real*8, intent(inout) :: xdrift(:,:),ydrift(:,:),zdrift(:,:)
      nec = mElecConfigs
      phi = elPhi
      u = elU
      eloc = elEloc
      vpot = elVpot
      xdrift = elxDrift; ydrift = elyDrift; zdrift = elzDrift
      end subroutine eloc_getCurrentElocData

c     -------------------------------------------------------
      subroutine eloc_getCurrentElocData1(phi,u,eloc,vpot,
     .           xdrift,ydrift,zdrift)
c     -------------------------------------------------------

      ! this returns non-array data containing the first elec config
      ! CHANGE to returning pointers to avoid copying
      real*8, intent(inout) :: phi
      real*8, intent(inout) :: u
      real*8, intent(inout) :: eloc
      real*8, intent(inout) :: vpot
      real*8, intent(inout) :: xdrift(:),ydrift(:),zdrift(:)
      phi = elPhi(1)
      u = elU(1)
      eloc = elEloc(1)
      vpot = elVpot(1)
      xdrift = elxDrift(:,1)
      ydrift = elyDrift(:,1)
      zdrift = elzDrift(:,1)
      end subroutine eloc_getCurrentElocData1

      subroutine eloc_getCurrentElocEpart(nec,kin,vee,vne)
      real*8, intent(inout) :: kin(:)
      real*8, intent(inout) :: vee(:,:)
      real*8, intent(inout) :: vne(:,:)
      integer,intent(inout)    :: nec

      kin = elEkin(nec,:)
      vee = elVee(nec,:,:)
      vne = elVne(nec,:,:)

      end subroutine eloc_getCurrentElocEpart

      subroutine eloc_getCurrentElocEpart1(kin,vee,vne)
      real*8, intent(inout) :: kin(:)
      real*8, intent(inout) :: vee(:,:)
      real*8, intent(inout) :: vne(:,:)

      kin = elEkin(1,:)
      vee = elVee(1,:,:)
      vne = elVne(1,:,:)

      end subroutine eloc_getCurrentElocEpart1

      subroutine eloc_getECPContribs1(EECPlocal,EECPnonlocal)
         real*8, intent(inout) :: EECPlocal,EECPnonlocal
         EECPlocal = elECPPot
         EECPnonlocal = elECPPotNl
      end subroutine


c     ---------------------------------------
      subroutine setElocCutOff(cutoff,tau,cf)
c     ---------------------------------------

      logical, intent(in) :: cutoff
      real*8, intent(in) :: tau
      real*8, intent(in) :: cf

      mElocCutOff=cutoff
      mTauCutOff=tau
      mCutOffFactor=cf
      end subroutine setElocCutOff


c     ---------------------------------------
      subroutine getElocCutOff(cutoff,tau,cf)
c     ---------------------------------------

      logical, intent(out) :: cutoff
      real*8, intent(out), optional :: tau
      real*8, intent(out), optional :: cf

      cutoff=mElocCutOff
      if (present(tau)) tau=mTauCutOff
      if (present(cf)) cf=mCutOffFactor
      end subroutine getElocCutOff


      END MODULE elocaldata


