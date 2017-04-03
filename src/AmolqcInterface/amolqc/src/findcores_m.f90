module findcores
  use wfdata
  use global
  implicit none

   interface findcores_sort
    module procedure findcores_sort_1D_Array
    module procedure findcores_sort_2D_Array
    module procedure findcores_sort_coreel
  end interface


   private
    public ::findcores_find,findcores_sort_Vne,findcores_sort,set_core_rule
    integer                 :: mCore = 1 ! 1= only [He] Core seperated, 2=for Z=3-10,[Ne] = Z=11... so far
    integer                 :: mNumber_of_core_el = 0
    integer                 :: mNumber_of_cores = 0
    integer,allocatable     :: mNumber_of_el_and_coreidx(:,:)
    integer,allocatable     :: mElectron_assign_idx(:)

  contains

 !----------------------------!
  subroutine set_core_rule(n)
!----------------------------
  integer,intent(in) :: n
   mCore = n
  end subroutine set_core_rule


 !----------------------------
  subroutine findcores_find(nocel)
 !---------------------------
  integer,intent(out)  :: nocel
  integer              :: tmp_el_core(ncenter)
  integer              :: i,j

  allocate(mElectron_assign_idx(ne))


  mNumber_of_cores = 0

! we look for the Z of the Atom and define the number of coreelectron*0.5, because we have
! to handle alpha and beta electrons seperately (etc.open shell)

  if(mCore == 1) then
    do i=1,ncenter
      if(atoms(i)%Elemidx <= 2) then
        tmp_el_core(i) = 0
      elseif(atoms(i)%Elemidx >2 ) then
        tmp_el_core(i) = 1
        mNumber_of_cores = mNumber_of_cores + 1
      else
       call abortp('findcores::mCore 1, elemidx assign. not possible')
      endif
    enddo
  elseif(mCore == 2) then
    do i=1,ncenter
      if(atoms(i)%Elemidx <= 2) then
        tmp_el_core(i) = 0
      elseif(atoms(i)%Elemidx >2 .and. atoms(i)%Elemidx <=10 ) then
        tmp_el_core(i) = 1
         mNumber_of_cores = mNumber_of_cores + 1
      elseif(atoms(i)%Elemidx >10) then
        tmp_el_core(i) = 5
         mNumber_of_cores = mNumber_of_cores + 1
      else
        call abortp('findcores::mCore 1, elemidx assign. not possible')
      endif
    enddo
  elseif (mCore == 3) then
        tmp_el_core = 1
        mNumber_of_cores = ncenter
  else
     call abortp('findcores :: mCore undefined')
  endif

  nocel = sum(tmp_el_core(:))
  mNumber_of_core_el = nocel

  allocate(mNumber_of_el_and_coreidx(mNumber_of_cores,2))

! write working array for seperation of corelectron
! doing this before should reduce number of operations in the seperation

  j = 1
  do i=1,ncenter
    if(tmp_el_core(i) /= 0) then
      mNumber_of_el_and_coreidx(j,1) = i
      mNumber_of_el_and_coreidx(j,2) = tmp_el_core(i)
      j=j+1
    endif
  enddo

  end subroutine findcores_find

 !--------------------------------------!
  subroutine findcores_sort_coreel(nall)
 !--------------------------------------
  ! push the alpha/beta core eletrons to the end of the
  ! alpha/beta electrons and write idxtable, for sorting
  ! othere arrays
   type(coord),intent(inout) :: nall
   type(coord)            :: ae,be
   real*8                 :: dist_en_a(nalpha,mNumber_of_cores)
   real*8                 :: dist_en_b(nbeta,mNumber_of_cores)
   real*8                 :: xx,yy,zz
   integer                :: min_val(2)
   integer                :: tmp_idx_a(mNumber_of_cores,ne)
   integer                :: tmp_idx_b(mNumber_of_cores,ne)
   integer                :: el_idx_a(nalpha),el_idx_b(nbeta)
   integer                :: i,tmp

    allocate(ae%x(nalpha),ae%y(nalpha),ae%z(nalpha))
    allocate(be%x(nbeta),be%y(nbeta),be%z(nbeta))

    ae%x = nall%x(1:nalpha)
    ae%y = nall%y(1:nalpha)
    ae%z = nall%z(1:nalpha)

    be%x(1:nbeta) = nall%x(nalpha+1:ne)
    be%y(1:nbeta) = nall%y(nalpha+1:ne)
    be%z(1:nbeta) = nall%z(nalpha+1:ne)

  ! Calculate the distance matrix el - core
   call findcores_get_dist_en(ae,dist_en_a)
   call findcores_get_dist_en(be,dist_en_b)

  ! Get the assignment of core electrons to Atoms
   call findcores_assign(dist_en_a,tmp_idx_a)
   call findcores_assign(dist_en_b,tmp_idx_b)

  ! Sort ae and be and get the index table
  ! index table means (el_idx) with assigns the original ae to the resorted one
  ! we need this for the sorting of arrays related to energy partitioning
   call findcores_sort_array(ae,tmp_idx_a,el_idx_a)
   call findcores_sort_array(be,tmp_idx_b,el_idx_b)

  ! now write the final index table for a complete array(1:ne)

    nall%x(1:nalpha) = ae%x
    nall%y(1:nalpha) = ae%y
    nall%z(1:nalpha) = ae%z

   nall%x(nalpha+1:ne) = be%x(1:nbeta)
   nall%y(nalpha+1:ne) = be%y(1:nbeta)
   nall%z(nalpha+1:ne) = be%z(1:nbeta)

   do i =1,nalpha
    mElectron_assign_idx(i) = el_idx_a(i)
   enddo

   do i=1,nbeta
    tmp =  el_idx_b(i) + nalpha
    mElectron_assign_idx(nalpha+i) = tmp
   enddo

   deallocate(ae%x,ae%y,ae%z)
   deallocate(be%x,be%y,be%z)

  end subroutine findcores_sort_coreel

 !------------------------------- ---!
  subroutine findcores_sort_1D_Array(a)
 !-----------------------------------!
   real*8,intent(inout)    :: a(ne)
   real*8                  :: tmp(ne)
   integer                 :: i

   do i=1,ne
    tmp(i) = a(mElectron_assign_idx(i))
   enddo

    a = tmp

  end subroutine findcores_sort_1D_Array

 !-------------------------------------
  subroutine findcores_sort_2D_Array(a)
 !-------------------------------------
   real*8,intent(inout)    :: a(ne,ne)
   real*8                  :: tmp(ne,ne)
   integer                 :: i,j


    do i=1,ne
     do j=1,ne
      tmp(i,j) = a(mElectron_assign_idx(i),mElectron_assign_idx(j))
     enddo
    enddo

     a = tmp

  end subroutine findcores_sort_2D_Array

 !-------------------------------------
  subroutine findcores_sort_Vne(a)
 !-------------------------------------
   real*8,intent(inout)   :: a(ncenter,ne)
   real*8                 :: tmp(ncenter,ne)
   integer                :: i

   do i=1,ne
     tmp(:,i) = a(:,mElectron_assign_idx(i))
   enddo

     a = tmp

  end subroutine findcores_sort_Vne


!----------------------------------------------------------------------------------------------------
! Routines for sorting the alpha and betha coordinate sets
! Seems to be the best way so far, since we need to know how the whole coordinate set is sorted so far
! without introducing the different subroutines, its pretty messy

 !--------------------------!
  subroutine findcores_get_dist_en(a,dist)
 !--------------------------!
   type(coord),intent(in) :: a
   real*8, intent(out)    :: dist(:,:)
   real*8  :: xx,yy,zz
   integer :: i,j

   do i=1,size(a%x)
    do j=1,mNumber_of_cores
       xx = (a%x(i) - atoms(mNumber_of_el_and_coreidx(j,1))%cx)**2
       yy = (a%y(i) - atoms(mNumber_of_el_and_coreidx(j,1))%cy)**2
       zz = (a%z(i) - atoms(mNumber_of_el_and_coreidx(j,1))%cz)**2
       dist(i,j) = sqrt(xx + yy + zz)
    enddo
   enddo

  end subroutine findcores_get_dist_en

  !--------------------------!
   subroutine findcores_assign(dist,tmp_idx)
  !--------------------------!
    real*8,intent(inout)      :: dist(:,:)
    integer,intent(out)    :: tmp_idx(:,:)
    integer                :: vec_diff(mNumber_of_cores)
    integer                :: tmp_count(mNumber_of_cores)
    integer                :: min_val(2)

    tmp_count = 0
    tmp_idx = 0
    vec_diff = 1
    ! Pretty confusing and messy, but so far the easiest way to assign core electrons to
    ! different atoms
    ! In Principle this is a kind of asymmetric assignment problem,possibly there is a more effective way
    do while(sum(vec_diff) /= 0)
      min_val = minloc(dist, MASK=dist .ge. 0)
     if(tmp_count(min_val(2))/= mNumber_of_el_and_coreidx(min_val(2),2))then
      tmp_count(min_val(2)) = tmp_count(min_val(2)) + 1
      tmp_idx(min_val(2),tmp_count(min_val(2))) = min_val(1)
      dist(min_val(1),:) = -1
      if(tmp_count(min_val(2)) == mNumber_of_el_and_coreidx(min_val(2),2)) dist(:,min_val(2)) = -1
     endif
     vec_diff = mNumber_of_el_and_coreidx(:,2) - tmp_count
    enddo

  end subroutine findcores_assign

 !---------------------------------
  subroutine findcores_sort_array(el,tmp_idx,el_idx_table)
 !---------------------------------
  type(coord),intent(inout) :: el
  integer,intent(in)        :: tmp_idx(:,:)
  integer,intent(inout)     :: el_idx_table(:)
  real*8,allocatable        :: x_tmp(:),y_tmp(:),z_tmp(:)
  integer                   :: idx1,sum_core_el
  integer                   :: c1,i,j
  logical,allocatable       :: assign_table(:) ! to check wether a electron was assigned yet or not
                                            ! this is more secure than the old way, but produces overhead

   allocate(x_tmp(size(el%x)),y_tmp(size(el%x)),z_tmp(size(el%x)))
   allocate(assign_table(size(el%x)))
  !First write core electrons to the end of the arrays

   assign_table = .false.

   sum_core_el = sum(mNumber_of_el_and_coreidx(:,2))
   idx1 = size(el%x) - sum_core_el
   c1 = 0
   el_idx_table = 0

   do i=1,mNumber_of_cores
    do j=1, mNumber_of_el_and_coreidx(i,2)
       c1 = c1 +1
       x_tmp(idx1+c1) = el%x(tmp_idx(i,j))
       y_tmp(idx1+c1) = el%y(tmp_idx(i,j))
       z_tmp(idx1+c1) = el%z(tmp_idx(i,j))
       el_idx_table(idx1+c1) = tmp_idx(i,j)
       assign_table(tmp_idx(i,j)) = .true.
    enddo
   enddo

  ! no write the remaining unassigned electrons vom original el to tmp
     c1 =0
   do i=1,size(el%x)
    if(assign_table(i) .eqv. .false.)then
     c1 = c1 + 1
     x_tmp(c1) = el%x(i)
     y_tmp(c1) = el%y(i)
     z_tmp(c1) = el%z(i)
     el_idx_table(c1) = i
    endif
   enddo

  ! el is now the resorted coordinate and will be given back to main routine
    el%x = x_tmp
    el%y = y_tmp
    el%z = z_tmp

  deallocate(x_tmp,y_tmp,z_tmp,assign_table)
  end subroutine findcores_sort_array

end module findcores
