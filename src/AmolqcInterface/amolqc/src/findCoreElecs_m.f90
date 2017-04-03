!
module findCoreModule
!
!  note that there is also Rene's findcore module (findcore_m.f90)
!
!  this module uses simple array data structure to identify
!  and possibly swap electrons close to the nuclei
!
   use global
   use wfdata
   implicit none

   private
   public :: initCoreLists, findCoreElecs, findHCoreElecs, findNucElecs, findAllCoreElecs


contains


   subroutine initCoreLists(naCore,aCoreList,nbCore,bCoreList)
      integer, intent(out)  :: naCore, nbCore    ! alpha, beta elecs at nucleus (not varied)
      integer, intent(inout):: aCoreList(:), bCoreList(:)
      integer a,j
      j = 0
      do a=1,ncenter
         if ( atoms(a)%za > 2 .and. (.not.atoms(a)%ecp) ) then
            j = j + 1
            aCoreList(j) = a
            bCoreList(j) = a
         endif
      enddo
      naCore = j
      nbCore = j
   end subroutine initCoreLists

   subroutine findCoreElecs(naCore,aCoreList,nbCore,bCoreList,x,y,z,idx)
      ! identifies K shell core electrons by comparing distances
      ! swap electrons such the core electrons are the first in alpha and
      ! beta list and put core electrons at nucleus position
      integer, intent(in)   :: naCore, nbCore    ! alpha, beta elecs at nucleus (not varied)
      integer, intent(in)   :: aCoreList(:), bCoreList(:)
      real*8, intent(inout) :: x(:),y(:),z(:)  ! coords
      integer, intent(inout), optional :: idx(:)  ! j=idx(i) means x(i) was before x(j). idx is updated
                                                  ! meaning x(i) was x(j) prior to all swaps
      integer i,ia,ib,a,iamin,ibmin,ica,icb,na,itmp
      real*8 dmin
      real*8 rai(getNNuc(),getNElec())
      real*8, parameter :: EPS=0.d0   ! put electrons slightly outside nucleus to avoid singularities!

      call calcNucElecDists(x,y,z,rai)

      na = getNAlpha()
      do ia=1,naCore
         a = aCoreList(ia)
         iamin = minloc(rai(a,ia:na),dim=1) + ia - 1
         ! copy electron ia to iamin (deleting iamin) and set core electron ia to nucleus pos (slightly outside by EPS)
         ! copy also rai entry
         ica = ia
         x(iamin) = x(ica); y(iamin) = y(ica); z(iamin) = z(ica)
         rai(:,iamin) = rai(:,ica)
         x(ica) = atoms(a)%cx; y(ica) = atoms(a)%cy; z(ica) = atoms(a)%cz+EPS
         if (present(idx)) then
            itmp = idx(ica); idx(ica) = idx(iamin); idx(iamin) = itmp
         end if
      end do
      do ib=1,nbCore
         a = bCoreList(ib)
         ibmin = minloc(rai(a,na+ib:getNElec()),dim=1) + na + ib - 1
         icb = na + ib
         x(ibmin) = x(icb); y(ibmin) = y(icb); z(ibmin) = z(icb)
         rai(:,ibmin) = rai(:,icb)
         x(icb) = atoms(a)%cx; y(icb) = atoms(a)%cy; z(icb) = atoms(a)%cz-EPS
         if (present(idx)) then
            itmp = idx(icb); idx(icb) = idx(ibmin); idx(ibmin) = itmp
         end if
      end do
   end subroutine findCoreElecs

   subroutine findHCoreElecs(naCore,aCoreList,nbCore,bCoreList,thresh,x,y,z,found,idx)
      ! identifies elecs close to H (or He) nucleus.
      ! swap these electrons such that they are the first (after core) in alpha and
      ! beta list and put these electrons at the nucleus position
      integer, intent(inout):: naCore, nbCore    ! alpha, beta elecs at nucleus
      integer, intent(inout):: aCoreList(:), bCoreList(:)
      real*8, intent(in)    :: thresh            ! threshold to identify H core electrons
      real*8, intent(inout) :: x(:),y(:),z(:)  ! coords
      logical, intent(out)   :: found
      integer, intent(inout), optional :: idx(:)
      integer i,ia,ib,a,iamin,ibmin,ica,icb,na,itmp
      real*8 dmin
      real*8 rai(getNNuc(),getNElec())
      real*8, parameter :: EPS=0.d0   ! put electrons slightly outside nucleus to avoid singularities!

      found = .false.
      call calcNucElecDists(x,y,z,rai)   ! not all distances are necessary!

      na = getNAlpha()
      ! alpha electrons
      ia = naCore + 1
      do a=1,getNNuc()
         if (atoms(a)%za > 2) cycle
         if (any(aCoreList(1:naCore)==a)) cycle    ! never two alpha elecs at one nuc
         dmin = minval(rai(a,ia:na))
         iamin = minloc(rai(a,ia:na),dim=1) + ia - 1
         if (dmin < thresh) then
            ! copy electron ia to iamin (deleting iamin) and set core electron ia to nucleus pos (slightly outside by EPS)
            ! copy also rai entry
            ica = ia
            x(iamin) = x(ica); y(iamin) = y(ica); z(iamin) = z(ica)
            rai(:,iamin) = rai(:,ica)
            x(ica) = atoms(a)%cx; y(ica) = atoms(a)%cy; z(ica) = atoms(a)%cz+EPS
            naCore = naCore + 1
            aCoreList(naCore) = a
            ia = ia + 1
            found = .true.
            if (present(idx)) then
               itmp = idx(ica); idx(ica) = idx(iamin); idx(iamin) = itmp
            end if
        end if
      end do

      ! beta electrons
      ib = nbCore + 1
      do a=1,getNNuc()
         if (atoms(a)%za > 2) cycle
         if (any(bCoreList(1:nbCore)==a)) cycle    ! never two alpha elecs at one nuc
         dmin = minval(rai(a,na+ib:getNElec()))
         ibmin = minloc(rai(a,na+ib:getNElec()),dim=1) + na + ib - 1
         if (dmin < thresh) then
            icb = na+ib
            x(ibmin) = x(icb); y(ibmin) = y(icb); z(ibmin) = z(icb)
            rai(:,ibmin) = rai(:,icb)
            x(icb) = atoms(a)%cx; y(icb) = atoms(a)%cy; z(icb) = atoms(a)%cz-EPS
            nbCore = nbCore + 1
            bCoreList(nbCore) = a
            ib = ib + 1
            found = .true.
            if (present(idx)) then
               itmp = idx(icb); idx(icb) = idx(ibmin); idx(ibmin) = itmp
            end if
         end if
      enddo

   end subroutine findHCoreElecs

   subroutine findAllCoreElecs(naCore,aCoreList,nbCore,bCoreList,thresh,x,y,z,found)
      ! identifies elecs close to any nucleus (dist < thresh/Z)
      ! move electrons at nucleus to front and add to the list
      integer, intent(inout):: naCore, nbCore    ! alpha, beta elecs at nucleus
      integer, intent(inout):: aCoreList(:), bCoreList(:)
      real*8, intent(in)    :: thresh            ! threshold to identify H core electrons
      real*8, intent(inout) :: x(:),y(:),z(:)  ! coords
      logical, intent(out)   :: found
      integer i,ia,ib,a,iamin,ibmin,ica,icb,na
      real*8 dmin
      real*8 rai(getNNuc(),getNElec())
      real*8, parameter :: EPS=0.d0   ! put electrons slightly outside nucleus to avoid singularities!

      found = .false.
      call calcNucElecDists(x,y,z,rai)   ! not all distances are necessary!

      na = getNAlpha()
      ! alpha electrons
      ia = naCore + 1
      do a=1,getNNuc()
         if (any(aCoreList(1:naCore)==a)) cycle    ! never two alpha elecs at one nuc
         dmin = minval(rai(a,ia:na))
         iamin = minloc(rai(a,ia:na),dim=1) + ia - 1
         if (dmin < thresh/atoms(a)%za) then
            ! swap electrons and put core electron at nucleus (slightly outside by EPS)
            ica = ia
            x(iamin) = x(ica); y(iamin) = y(ica); z(iamin) = z(ica)
            x(ica) = atoms(a)%cx; y(ica) = atoms(a)%cy; z(ica) = atoms(a)%cz+EPS
            naCore = naCore + 1
            aCoreList(naCore) = a
            ia = ia + 1
            found = .true.
        end if
      end do

      ! beta electrons
      ib = nbCore + 1
      do a=1,getNNuc()
         if (any(bCoreList(1:nbCore)==a)) cycle    ! never two alpha elecs at one nuc
         dmin = minval(rai(a,na+ib:getNElec()))
         ibmin = minloc(rai(a,na+ib:getNElec()),dim=1) + na + ib - 1
         if (dmin < thresh/atoms(a)%za) then
            icb = na+ib
            x(ibmin) = x(icb); y(ibmin) = y(icb); z(ibmin) = z(icb)
            x(icb) = atoms(a)%cx; y(icb) = atoms(a)%cy; z(icb) = atoms(a)%cz-EPS
            nbCore = nbCore + 1
            bCoreList(nbCore) = a
            ib = ib + 1
            found = .true.
         end if
      enddo
   end subroutine findAllCoreElecs

   subroutine findNucElecs(thresh,x,y,z,naNuc,aNucElec,nbNuc,bNucElec)
      ! identifies elecs close to a nucleus core (if dist < thresh)
      ! build a/bNucList(na/bNuc)
      ! do NOT swap these electrons (use findCoreElecs for this)
      ! do NOT check if electrons are sorted
      real*8, intent(in)    :: thresh          ! threshold to identify H core electrons
      real*8, intent(in)    :: x(:),y(:),z(:)  ! coords
      integer, intent(out):: naNuc, nbNuc      ! alpha, beta elecs at nucleus (not varied)
      integer, intent(out):: aNucElec(:), bNucElec(:) ! nuc a -> elec i (at nuc)
      integer i,ia,ib,a,iamin,ibmin,ica,icb,na
      real*8 dmin,tmp
      real*8 rai(getNNuc(),getNElec())
      logical swap

      call assert(size(x)==getNElec() .and. size(aNucElec)==getNNuc(), &
                 '(findNucElecs): illegal sizes on entry')

      call calcNucElecDists(x,y,z,rai)
      na = getNAlpha()

      ! alpha electrons
      naNuc = 0
      aNucElec = 0
      do a=1,getNNuc()
         dmin = minval(rai(a,1:na))
         iamin = minloc(rai(a,1:na),dim=1)
         swap = .false.
         if (dmin < thresh) then
            ! swap electrons
            !if (iamin /= ia) then
            !   swap = .true.
            !   ica = ia
            !   tmp = x(iamin); x(iamin) = x(ica); x(ica) = tmp;
            !   tmp = y(iamin); y(iamin) = y(ica); y(ica) = tmp;
            !   tmp = z(iamin); z(iamin) = z(ica); z(ica) = tmp;
            !end if
            aNucElec(a) = iamin
            naNuc = naNuc + 1
         end if
         !!print'(i3,2f10.5,l4,10i3)',a,dmin,thresh,swap,iamin,naNuc,aNucElec
      end do

      ! beta electrons
      nbNuc = 0
      bNucElec = 0
      do a=1,getNNuc()
         dmin = minval(rai(a,na+1:getNElec()))
         ibmin = minloc(rai(a,na+1:getNElec()),dim=1) + na
         swap = .false.
         if (dmin < thresh) then
            ! swap electrons and put core electron at nucleus (slightly outside by EPS)
            !if (ibmin /= na + ib) then
            !   swap = .true.
            !   icb = na + ib
            !   tmp = x(ibmin); x(ibmin) = x(icb); x(icb) = tmp;
            !   tmp = y(ibmin); y(ibmin) = y(icb); y(icb) = tmp;
            !   tmp = z(ibmin); z(ibmin) = z(icb); z(icb) = tmp;
            !end if
            bNucElec(a) = ibmin
            nbNuc = nbNuc + 1
        end if
        !!print'(i3,2f10.5,l4,10i3)',a,dmin,thresh,swap,ibmin,nbNuc,bNucElec
      end do
   end subroutine findNucElecs

end module findCoreModule
