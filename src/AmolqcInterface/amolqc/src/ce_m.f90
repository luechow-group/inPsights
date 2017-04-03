!rewrite this to make it more flexible (only Eww, Eint can be calculated on the fly so far, add E_int_ne, E_int_ee...)
module chemical_entity
    use global
    use wfdata, only: Vnni,atoms,ncenter
    use statistics
    implicit none

    type ce
        character(5)      ::  ce_type
        integer,pointer   ::  cores(:)
        integer,pointer   ::  electrons(:)
    end type ce

    private
    public :: ce_init,ce_onthefly,ce_blockstat,ce_write_ce_tabletat,ce_get_nce

    real*8,allocatable   :: cx(:),cy(:),cz(:)
    real*8,allocatable   :: rx(:),ry(:),rz(:)
    real*8,allocatable   :: Ekin(:),Vee(:,:),Vne(:,:),Vnn(:,:)
    real*8,allocatable   :: std_Ekin(:),std_Vee(:,:),std_Vne(:,:),std_Vnn(:,:)
    real*8,allocatable   :: E_intern_list(:),std_E_intern_list(:)
    real*8,allocatable   :: E_ww_list(:,:),std_E_ww_list(:,:)
    integer              :: mCores,mElectrons
    integer              :: io_ce = 25
    logical              :: print_full = .false.
    type(ce),allocatable ::  ce_table(:)
    logical              :: mDebug = .false.
    logical              :: mRead=.false.
    integer              :: nce = 0
! Statistics for on the fly ce calculation
   !StepStat
    type(weightStat),allocatable  :: stepEint(:)
    type(weightStat),allocatable  :: stepEww(:,:)
    !BlockStat
    type(simpleStat),allocatable  :: blockEint(:)
    type(simpleStat),allocatable  :: blockEww(:,:)
    !Totalstat
    type(weightStat),allocatable  :: totalEint(:)
    type(weightStat),allocatable  :: totalEww(:,:)
contains

    !---------------------------!
    subroutine ce_init(lines,nl)
        !---------------------------!
  
        character(len=120), intent(in) :: lines(:)
        integer, intent(in)           :: nl
        logical                       :: finda
        integer                       :: ierr,i,j


            if(finda(lines,nl,'ce_print_full')) print_full =.true.
            if(finda(lines,nl,'debug')) mDebug = .true.
            if(finda(lines,nl,'readce')) then
                mRead=.true.
                call ce_read_energy_arrays()
                call ce_generate_table_from_file()
                call ce_calc_E_ce()
                if(print_full .eqv. .false.) then
                    call ce_write_E_ce(E_intern_list,std_E_intern_list,E_ww_list,std_E_ww_list,ce_table)
                endif
            elseif(finda(lines,nl,'autoce')) then
                write(iul,*) "automatic ce definition N.Y.I"
                stop
            else
                mCores = ncenter
                mElectrons = ne
                call ce_generate_table_from_file()
                allocate(Ekin(ne))
                allocate(Vee(ne,ne))
                allocate(Vne(ncenter,ne))
                allocate(Vnn(ncenter,ncenter))
                Vnn = Vnni
                allocate(stepEint(nce),blockEint(nce),totalEint(nce))
                allocate(stepEww(nce,nce),blockEww(nce,nce),totalEww(nce,nce))
                do i=1,nce
                   call reset(stepEint(i))
                   call reset(blockEint(i))
                   call reset(totalEint(i))
                   do j=1,nce
                      call reset(stepEww(i,j))
                      call reset(blockEww(i,j))
                      call reset(totalEww(i,j))
                   enddo
                enddo
            endif

    end subroutine ce_init

    subroutine ce_get_nce(n)
        integer,intent(out) :: n
        n = nce
    end subroutine ce_get_nce

    !--------------------------------------------
    subroutine ce_read_energy_arrays()
        !--------------------------------------------
        integer :: i,j,n,m,dummy

        open(26,file=trim(baseName)//'.ce')
        read(26,*) mElectrons,mCores
        n  = mElectrons
        m  = mCores

        allocate(cx(m),cy(m),cz(m))
        allocate(rx(n),ry(n),rz(n))
        allocate(Ekin(n))
        allocate(Vee(n,n))
        allocate(Vne(m,n))
        allocate(Vnn(m,m))
 
        allocate(std_Ekin(n))
        allocate(std_Vee(n,n))
        allocate(std_Vne(m,n))
        allocate(std_Vnn(m,m))

        Ekin = 0;Vee = 0; Vne = 0;Vnn= 0

        read(26,*)
        do i=1,mCores
            read(26,*) cx(i),cy(i),cz(i)
        enddo
        read(26,*)
        read(26,*)
        do i=1,mElectrons
            read(26,*) rx(i),ry(i),rz(i)
        enddo
        read(26,*)
        read(26,*)
        do i=1,mElectrons
            read(26,*) dummy,Ekin(i),std_Ekin(i)
        enddo
        read(26,*)
        read(26,*)
        do i=1,mElectrons
            do j=1,mElectrons
                read(26,*) dummy,dummy,Vee(i,j),std_Vee(i,j)
            enddo
        enddo

  
        read(26,*)
        read(26,*)
        do i=1,mElectrons
            do j=1,mCores
                read(26,*) dummy,dummy,Vne(j,i),std_Vne(j,i)
            enddo
        enddo
        read(26,*)
        read(26,*)
        do i=1,mCores
            do j=1,mCores
                read(26,*) dummy,dummy,Vnn(i,j)
            enddo
        enddo

        call ce_write_energy_arrays()

    end subroutine ce_read_energy_arrays

    subroutine ce_write_energy_arrays()
        integer :: i,j


        write(iul,*) "Geometry"
        do i=1,mCores
            write(iul,*) cx(i)," ",cy(i)," ",cz(i)
        enddo
        write(iul,*)
        write(iul,*) "Referenz"
        do i=1,mElectrons
            write(iul,*) rx(i)," ",ry(i)," ",rz(i)
        enddo
        write(iul,*)
        write(iul,*) "Electrons"," ","Nuclei"
        write(iul,*) mElectrons," ",mCores
        write(iul,*) "Ekin i"
        do i=1,mElectrons
            write(iul,*) Ekin(i)," ",std_Ekin(i)
        enddo
        write(iul,*)
        write(iul,*) "Vee i j"
        do i=1,mElectrons
            do j=i+1,mElectrons
                write(iul,*) Vee(i,j)," ",std_Vee(i,j)
            enddo
        enddo
        write(iul,*)
        write(iul,*) "Vne i j"
        do i=1,mElectrons
            do j=1,mCores
                write(iul,*) Vne(j,i)," ",std_Vne(j,i)
            enddo
        enddo
        write(iul,*)
        write(iul,*) "Vnn i j"
        do i=1,mCores
            do j=i+1,mCores
                write(iul,*) Vnn(i,j)
            enddo
        enddo
        write(iul,*)
        write(iul,*)

    end subroutine ce_write_energy_arrays

    !--------------------------------------!
    subroutine ce_onthefly(ee_in,ne_in,kin_in,w)
    !--------------------------------------!
     real*8,intent(in)  :: ee_in(:,:),ne_in(:,:),kin_in(:),w
     integer            :: i,j

     Ekin = kin_in
     Vee = ee_in
     Vne = ne_in

     call ce_calc_E_ce()
     do i=1,size(ce_table)
        call addData(stepEint(i),E_intern_list(i),w)
        do j=i+1,size(ce_table)
           call addData(stepEww(i,j),E_ww_list(i,j),w)
        enddo
     enddo

    end subroutine ce_onthefly

    !-----------------------!
    subroutine ce_blockstat(dat)
    !-----------------------!
     integer           :: i,j
     logical           :: dat
      do i=1,size(ce_table)
        call addData(blockEint(i),meanAllNodes(stepEint(i)))
        totalEint(i) = totalEint(i) + stepEint(i)
        !call reset(stepEint(i))
        do j=i+1,size(ce_table)
            call addData(blockEww(i,j),meanAllNodes(stepEww(i,j)))
            totalEww(i,j) = totalEww(i,j) + stepEww(i,j)
            !call reset(stepEww(i,j))
        enddo
     enddo
     if(dat) call ce_write_block_mean()

     do i=1,size(ce_table)
        call reset(stepEint(i))
        do j=i+1,size(ce_table)
            call reset(stepEww(i,j))
        enddo
     enddo
     
    end subroutine ce_blockstat
    subroutine ce_write_block_mean()
        integer      :: i,j
        real*8       :: ceint(nce),ceww(nce,nce)
        real*8       :: sdceint(nce),sdceww(nce,nce)

        do i=1,nce
           ceint(i) = meanAllNodes(stepEint(i))
           sdceint(i) = stdDevMeanAllNodes(stepEint(i))
         do j=i+1,nce
           ceww(i,j) = meanAllNodes(stepEww(i,j))
           sdceww(i,j) = stdDevMeanAllNodes(stepEww(i,j))
         enddo
        enddo
        if(MASTER)then
           do i=1,nce
               write(17) ceint(i),sdceint(i)
           enddo
           do i=1,nce
              do j=i+1,nce
                 write(17) ceww(i,j),sdceww(i,j)
              enddo
           enddo
        endif
    end subroutine ce_write_block_mean

    !--------------------------!
    subroutine ce_write_ce_tabletat()
    !--------------------------!
     integer            :: i,j

     if(MASTER)then
        write(iul,*)
        write(iul,*)
        write(iul,*) "          CE Energy Partitioning -- On the Fly --"
        write(iul,*) "********************************************************"
        write(iul,*)
        write(iul,*) 'CE                 E_int     +/-     stdDev'
        write(iul,*) '___________________________________________'
        do i=1,nce
           write(iul,'(I3,A,F13.5,A,F13.5)') i,trim(ce_table(i)%ce_type), mean(blockEint(i)),'  +/-  ',stdDevMean(blockEint(i))
        enddo
        write(iul,*)
        write(iul,*) 'CE1      CE2       E_ww       +/-     stdDev'
        write(iul,*) '____________________________________________'
        do i=1,nce
            do j=i+1,nce
                 write(iul,'(I2,A,I2,A,F13.5,A,F13.5)') i,trim(ce_table(i)%ce_type),j, &
                 trim(ce_table(j)%ce_type),mean(blockEww(i,j)),'  +/-  ',stdDevMean(blockEww(i,j))
            enddo
        enddo

     endif

    end subroutine ce_write_ce_tabletat

    !---------------------------------!
    subroutine ce_calc_E_ce()
        !--------------------------------!
        integer      :: i,j
        real*8     :: E_ww,std_Eww

        if(print_full)then
            write(iul,*) " "
            write(iul,*) " CE Energy Partitioning"
            write(iul,*) " "
        endif
        do i = 1,size(ce_table)
            if(print_full) write(iul,'(I3,A)') i,trim(ce_table(i)%ce_type)
            if(mDebug) write(*,*) "CE",i
            call ce_calc_E_intern(E_intern_list(i),std_E_intern_list(i),ce_table(i))
            do j= i+1,size(ce_table)
                if(print_full) write(iul,'(I3,A,I3,A)') i,trim(ce_table(i)%ce_type),j,trim(ce_table(j)%ce_type)
                if(mDebug) write(*,*) "CE1",i,"CE2",j
                call ce_calc_E_ww(E_ww_list(i,j) ,std_E_ww_list(i,j),ce_table(i),ce_table(j))
            enddo
        enddo
 

    end subroutine ce_calc_E_ce

    !--------------------------------------------
    subroutine  ce_write_E_ce(E_intern,std_E_intern,E_ww,std_E_ww,ce_table)
        !--------------------------------------------
        real*8,intent(in) :: E_intern(:),E_ww(:,:),std_E_intern(:),std_E_ww(:,:)
        type(ce),intent(in)  :: ce_table(:)
        integer           :: i,j

        write(iul,*) " "
        write(iul,*) " CE Energy Partitioning"
        write(iul,*) " "
        do i = 1,size(ce_table)
            write(iul,'(I3,A,X,F13.5,F13.5)') i,trim(ce_table(i)%ce_type),E_intern(i),std_E_intern(i)
        enddo

        do i = 1,size(ce_table)
            do j= i+1,size(ce_table)
                write(iul,'(I3,A,I3,A,F13.5,F13.5)') i,trim(ce_table(i)%ce_type),j,trim(ce_table(j)%ce_type),E_ww(i,j),std_E_ww(i,j)
            enddo
        enddo
        write(iul,*) " "
        write(iul,*) " "


    end subroutine ce_write_E_ce

    !--------------------------if(mDebug---------------------
    subroutine ce_generate_table_from_file()
        !-----------------------------------------------
        integer     :: no_of_ce,no_of_atoms,no_of_el
        integer     :: i,k

        character(5) :: type_name
        integer,allocatable :: idx_tmp(:)
 
        open(io_ce,file="definition.ce")
        read(io_ce,*) no_of_ce
        nce = no_of_ce
        allocate(idx_tmp(max(mCores,mElectrons)))
 

        !Allocate ce-table
        allocate(ce_table(no_of_ce))
        do i=1,no_of_ce
            allocate(ce_table(i)%cores(mCores))
            allocate(ce_table(i)%electrons(mElectrons))
            ce_table(i)%cores = 0
            ce_table(i)%electrons = 0
            ce_table(i)%ce_type = ""
        enddo

        do i=1,no_of_ce
            read(io_ce,*) type_name
            ce_table(i)%ce_type = trim(type_name)
            read(io_ce,*) no_of_atoms

            !Set Core Assignment Vector
            read(io_ce,*) (idx_tmp(k),k=1,no_of_atoms)
            do k=1,no_of_atoms
                ce_table(i)%cores(idx_tmp(k)) = 1
            enddo

            idx_tmp = 0
            ! Set Electron Assigment Vector
            read(io_ce,*) no_of_el

            read(io_ce,*) (idx_tmp(k),k=1,no_of_el)
            do k=1,no_of_el
                ce_table(i)%electrons(idx_tmp(k)) = 1
            enddo

            idx_tmp = 0
        enddo
        allocate(E_intern_list(nce),std_E_intern_list(nce))
        allocate(E_ww_list(nce,nce),std_E_ww_list(nce,nce))
        call ce_write_ce_table()

    end subroutine ce_generate_table_from_file

    !---------------------------------------
    subroutine ce_write_ce_table()
        !---------------------------------------
        integer              ::  i,k

        write(iul,*) " "
        write(iul,*) " CE Definitions for CE-Partitioning"
        write(iul,*) " "

        do i=1,size(ce_table(:))
            write(iul,*) "CE",i,"Type: ",ce_table(i)%ce_type
            write(iul,*) "Index of Nuclei",(ce_table(i)%cores(k),k=1,size(ce_table(i)%cores(:)) )
            write(iul,*) "Index of Electrons",(ce_table(i)%electrons(k),k=1,size(ce_table(i)%electrons(:)) )
        enddo
    end subroutine ce_write_ce_table

    !--------------------------------!
    subroutine ce_calc_E_intern(E_int,std_Eint,ce1)
        !--------------------------------!
        ! E_int =  n x Ekin  + n x m Ven + (N 2) nxn Vee + (M 2) mxm Vnn
        type(ce) :: ce1
        integer               :: i,j
        real*8                :: E_int,std_Eint
        real*8                :: E_kin,V_ne,V_ee,V_nn
        real*8                :: std_E_kin,std_V_ne,std_V_ee

        E_kin=0;V_ne=0;V_ee=0;V_nn=0;E_int=0
        if(mRead) std_E_kin=0;std_V_ne=0;std_V_ee=0;std_Eint=0

        do i = 1,mElectrons
  
            E_kin = E_kin + ce1%electrons(i)*Ekin(i)
            if(mRead) std_E_kin = std_E_kin + ce1%electrons(i)*std_Ekin(i)**2
            if(mDebug .and. ce1%electrons(i)==1) write(*,*) "Ekin:",i,Ekin(i)
            do j = 1,mCores
                V_ne =  V_ne + ce1%electrons(i)*ce1%cores(j)*Vne(j,i)
                if(mDebug .and. ce1%electrons(i)*ce1%cores(j)==1) write(*,*)"Vne:",j,i,Vne(j,i)
                if(mRead)std_V_ne = std_V_ne +  ce1%electrons(i)*ce1%cores(j)*std_Vne(j,i)**2
            enddo
        enddo

        do i = 1,mElectrons
            do j=i+1,mElectrons
                if(mDebug .and. ce1%electrons(i)*ce1%electrons(j)==1) write(*,*) "Vee:",i,j,Vee(i,j)
                V_ee =  V_ee + ce1%electrons(i)*ce1%electrons(j)*Vee(i,j)
                if(mRead)std_V_ee = std_V_ee + ce1%electrons(i)*ce1%electrons(j)*std_Vee(i,j)**2
            enddo
        enddo

        do i = 1,mCores
            do j=i+1,mCores
                if(mDebug .and. ce1%cores(i)*ce1%cores(j)==1) write(*,*) "Vnn:",i,j,Vnn(i,j)
                V_nn =  V_nn +  ce1%cores(i)*ce1%cores(j)*Vnn(i,j)
            enddo
        enddo

        E_int = E_kin + V_ne + V_ee + V_nn
        if(mRead)std_Eint = std_E_kin + std_V_ne + std_V_ee
        if(mRead)std_Eint = sqrt(std_Eint)
        if(print_full) then
            write(iul,'(5(A))') "E_int ","E_kin ","V_ne ","V_ee ","V_nn "
            write(iul,'(5(F13.5))') E_int,E_kin,V_ne,V_ee,V_nn
           if(mRead) write(iul,'(4(A))') "std_Eint","std_E_kin","std_V_ne","std_Vee"
           if(mRead) write(iul,'(4(F13.6))') std_Eint,sqrt(std_E_kin),sqrt(std_V_ne),sqrt(std_V_ee)
        endif
    end subroutine ce_calc_E_intern

    !----------------------------------------!
    subroutine ce_calc_E_ww(E_ww,std_Eww,ce1,ce2)
        !---------------------------------------!
        ! E_ww  =  E_ce1 x E_ce2 =  (n 2) n1 * n2 Vee + n1 * m2 Ven  + n2 * m1 Ven + (m 2)  m1 x m2 Vnn
        type(ce),intent(in) :: ce1,ce2
        integer               :: i,j
        real*8,intent(out)    :: E_ww,std_Eww
        real*8                :: V_ne,V_ee,V_nn
        real*8                :: std_V_ne,std_V_ee

        V_ne=0;V_ee=0;V_nn=0 ;E_ww=0
        if(mRead) std_V_ne=0;std_V_ee=0;std_Eww=0
        ! Ven
 
        do i = 1,mElectrons
            do j = 1,mCores
                if(ce1%electrons(i)*ce2%cores(j) == 1) then
                    V_ne =  V_ne + Vne(j,i)
                    if(mDebug) write(*,*) "Vne1:",j,i,Vne(j,i)
                    if(mRead) std_V_ne = std_V_ne + std_Vne(j,i)**2
                endif
            enddo
        enddo


        do i = 1,mElectrons
            do j = 1,mCores
                if(ce2%electrons(i)*ce1%cores(j) == 1) then
                    V_ne =  V_ne +  Vne(j,i)
                    if(mDebug) write(*,*) "Vne2:",j,i,Vne(j,i)
                    if(mRead) std_V_ne = std_V_ne + std_Vne(j,i)**2
                endif
            enddo
        enddo

        ! Vee
        do i = 1,mElectrons
            do j = 1,mElectrons
                if( ce1%electrons(i)*ce2%electrons(j)== 1) then
                    if(mDebug) write(*,*) "Vee:",i,j,Vee(i,j)
                    V_ee =  V_ee +  Vee(i,j)
                    if(mRead) std_V_ee = std_V_ee + std_Vee(i,j)**2
                endif
            enddo
        enddo
 
        ! Vnn

        do i = 1,mCores
            do j = 1,mCores
                if(ce1%cores(i)*ce2%cores(j) == 1) then
                    V_nn =  V_nn +  Vnn(i,j)
                    if(mDebug) write(*,*) "Vnn:",i,j,Vnn(i,j)
                endif
            enddo
        enddo


        E_ww = V_ne +  V_ee + V_nn
        if(mRead) std_Eww = std_V_ne + std_V_ee
        if(mRead) std_Eww = sqrt(std_Eww)
        !write(*,*) E_ww
        if(print_full) then
            write(iul,'(4(A))') "E_ww ","V_ne ","V_ee ","V_nn "
            write(iul,'(4(F13.5,X))') E_ww,V_ne,V_ee,V_nn
            if(mRead) write(iul,'(3(A))') "std_Eww","std_V_ne","std_Vee"
            if(mRead) write(iul,'(4(F13.6,X))') std_Eww,sqrt(std_V_ne),sqrt(std_V_ee)
        endif
    end subroutine ce_calc_E_ww

end module chemical_entity
