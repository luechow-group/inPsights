module ce_new_m
    use wfdata, only: Vnni,atoms,ncenter
    use ceestat
    use global
    use primitive_cee
    implicit none
    private
    public :: ce_init,ce_blockstat,ce_onthefly,ce_finalize,ce_get_nce,ce_conv_check

    type ce
        character(5)      ::  name
        integer,pointer   ::  nu(:)
        integer,pointer   ::  el(:)
    end type ce
    logical                        :: mConvCheck = .false.
    real*8                         :: mSigmaTol = 1D-3
    logical                        :: mRead=.false.
    integer                        :: mNCe                       !Number of CEs
    integer                        :: mCores,mElectrons
    integer                        :: io_ce=27
    type(ce),allocatable           :: ces(:)
    type(w_cee),allocatable        :: sEint(:),sEww(:,:)         !Step Statistic
    type(s_cee),allocatable        :: bEint(:),bEww(:,:)         !Block Statistic
    type(cee),allocatable        :: pEint(:),pEww(:,:)           !primitive Cee struct for tmp storage
    type(cee),allocatable        :: pstdDevEint(:),pstdDevEww(:,:)
    real*8,allocatable :: Ekin(:),Vee(:,:),Vne(:,:),Vnn(:,:)  !Energy Arrays
    real*8,allocatable :: stdDevEkin(:),stdDevVee(:,:),stdDevVne(:,:)
    real*8,allocatable :: cx(:),cy(:),cz(:)
    real*8,allocatable :: rx(:),ry(:),rz(:)
contains

    subroutine ce_get_nce(n)
        integer,intent(out) :: n
        n = mNCe
    end subroutine


    !---------------------------!
    subroutine ce_init(lines,nl)
        !---------------------------!
        character(len=120), intent(in) :: lines(:)
        integer, intent(in)           :: nl
        logical                       :: finda
        integer                       :: ierr,i,j,iflag

        if(finda(lines,nl,'readce')) then
            mRead=.true.
            call ce_read_energy_arrays()
            call ce_read_ces_from_file()
            allocate(pEint(mNCe),pEww(mNCe,mNCe))
            allocate(pstdDevEint(mNCe),pstdDevEww(mNCe,mNCe))
            allocate(stdDevEkin(ne),stdDevVee(ne,ne),stdDevVne(ncenter,ne))
            call ce_calc_E_ce()
            call ce_write_ce_energies()
        else
            mCores = ncenter
            mElectrons = ne
            call ce_read_ces_from_file()
            allocate(Ekin(ne),Vee(ne,ne),Vne(ncenter,ne),Vnn(ncenter,ncenter))
            allocate(sEint(mNCe),sEww(mNCe,mNCe),bEint(mNce),bEww(mNCe,mNCe))
            allocate(pEint(mNCe),pEww(mNCe,mNCe))
            allocate(pstdDevEint(mNCe),pstdDevEww(mNCe,mNCe))
            Vnn=Vnni
            call getdbla(lines,nl,'cesigma=',mSigmaTol,iflag)
            if(iflag == 0) mConvCheck = .true.
            do i = 1,mNCe
                call reset(sEint(i))
                call reset(bEint(i))
                do j = 1,mNCe
                    call reset(sEww(i,j))
                    call reset(bEww(i,j))
                enddo
            enddo
        endif
    end subroutine ce_init

    !--------------------------------------!
    subroutine ce_onthefly(ee_in,ne_in,kin_in,w)
        !--------------------------------------!
        real*8,intent(in)  :: ee_in(:,:),ne_in(:,:),kin_in(:),w
        integer            :: i,j

        Ekin = kin_in;Vee = ee_in;Vne = ne_in
        call ce_calc_E_ce()
        do i=1,mNCe
            call addData(sEint(i),pEint(i),w)
            do j=i+1,mNCe
                call addData(sEww(i,j),pEww(i,j),w)
            enddo
        enddo

    end subroutine ce_onthefly

    !-----------------------!
    subroutine ce_blockstat(dat)
        !-----------------------!
        integer           :: i,j
        logical           :: dat
        do i=1,mNCe
            call addData(bEint(i),meanAllNodes(sEint(i)))
            do j=i+1,mNCe
                call addData(bEww(i,j),meanAllNodes(sEww(i,j)))
            enddo
        enddo
        if(dat) call ce_write_block_mean()

        do i=1,mNCe
            call reset(sEint(i))
            do j=i+1,mNCe
                call reset(sEww(i,j))
            enddo
        enddo

    end subroutine ce_blockstat


    !---------------------------------!
    subroutine ce_calc_E_ce()
        !--------------------------------!
        integer      :: i,j

        do i = 1,mNCe
            call ce_E_int(pEint(i),ces(i))
            do j= i+1,mNCe
                call ce_E_ww(pEww(i,j),ces(i),ces(j))
            enddo
        enddo
        if(mRead)then
            !calc approximate stdDev
        endif
    end subroutine ce_calc_E_ce


    !--------------------------------!
    subroutine ce_E_int(e,c)
        !--------------------------------!
        ! E_int =  n x Ekin  + n x m Ven + (N 2) nxn Vee + (M 2) mxm Vnn
        type(ce) :: c
        type(cee) :: e
        integer               :: i,j

        e%kin=0;e%ne=0;e%nn=0;e%ee=0;
        do i = 1,mElectrons
            e%kin = e%kin + c%el(i)*Ekin(i)
            do j = 1,mCores
                e%ne =  e%ne + c%el(i)*c%nu(j)*Vne(j,i)
            enddo
        enddo

        do i = 1,mElectrons
            do j=i+1,mElectrons
                e%ee =  e%ee + c%el(i)*c%el(j)*Vee(i,j)
            enddo
        enddo

        do i = 1,mCores
            do j=i+1,mCores
                e%nn =  e%nn +  c%nu(i)*c%nu(j)*Vnn(i,j)
            enddo
        enddo
    end subroutine ce_E_int

    !----------------------------------------!
    subroutine ce_E_ww(e,ce1,ce2)
        !---------------------------------------!
        ! E_ww  =  E_ce1 x E_ce2 =  (n 2) n1 * n2 Vee + n1 * m2 Ven  + n2 * m1 Ven + (m 2)  m1 x m2 Vnn
        type(ce),intent(in)       :: ce1,ce2
        type(cee),intent(inout)   :: e
        integer                   :: i,j

        e%kin=0;e%ne=0;e%nn=0;e%ee=0
        ! Ven
        do i = 1,mElectrons
            do j = 1,mCores
                if(ce1%el(i)*ce2%nu(j) == 1) then
                    e%ne =  e%ne + Vne(j,i)
                endif
            enddo
        enddo
        do i = 1,mElectrons
            do j = 1,mCores
                if(ce2%el(i)*ce1%nu(j) == 1) then
                    e%ne =  e%ne +  Vne(j,i)
                endif
            enddo
        enddo
        ! Vee
        do i = 1,mElectrons
            do j = 1,mElectrons
                if( ce1%el(i)*ce2%el(j)== 1) then
                    e%ee =  e%ee +  Vee(i,j)
                endif
            enddo
        enddo
        ! Vnn
        do i = 1,mCores
            do j = 1,mCores
                if(ce1%nu(i)*ce2%nu(j) == 1) then
                    e%nn =  e%nn +  Vnn(i,j)
                endif
            enddo
        enddo

    end subroutine ce_E_ww
    !---------------------------------------!
    subroutine ce_read_ces_from_file()
        !---------------------------------------!
        integer     :: no_of_atoms,no_of_el
        integer     :: i,k

        character(5) :: type_name
        integer,allocatable :: idx_tmp(:)

        open(io_ce,file="definition.ce")
        read(io_ce,*) mNCe

        allocate(idx_tmp(max(mCores,mElectrons)))


        !Allocate ce-table
        allocate(ces(mNCe))
        do i=1,mNCe
            allocate(ces(i)%nu(mCores))
            allocate(ces(i)%el(mElectrons))
            ces(i)%nu = 0
            ces(i)%el = 0
            ces(i)%name = ""
        enddo

        do i=1,mNCe
            read(io_ce,*) type_name
            ces(i)%name = trim(type_name)
            read(io_ce,*) no_of_atoms

            !Set Core Assignment Vector
            read(io_ce,*) (idx_tmp(k),k=1,no_of_atoms)
            do k=1,no_of_atoms
                ces(i)%nu(idx_tmp(k)) = 1
            enddo

            idx_tmp = 0
            ! Set Electron Assigment Vector
            read(io_ce,*) no_of_el

            read(io_ce,*) (idx_tmp(k),k=1,no_of_el)
            do k=1,no_of_el
                ces(i)%el(idx_tmp(k)) = 1
            enddo

            idx_tmp = 0
        enddo

    end subroutine ce_read_ces_from_file

    subroutine ce_write_block_mean()
        integer      :: i,j

        do i=1,mNCe
            pEint(i) = meanAllNodes(sEint(i))
            pstdDevEint(i) = stdDevMeanAllNodes(sEint(i))
            do j=i+1,mNCe
                pEww(i,j) = meanAllNodes(sEww(i,j))
                pstdDevEww(i,j) = stdDevMeanAllNodes(sEww(i,j))
            enddo
        enddo

        if(MASTER)then
            do i=1,mNCe
                write(17) Esum(pEint(i)),Esigma(pstdDevEint(i))
                write(17) pEint(i)%kin, pstdDevEint(i)%kin
                write(17) pEint(i)%nn, pstdDevEint(i)%nn
                write(17) pEint(i)%ne, pstdDevEint(i)%ne
                write(17) pEint(i)%ee, pstdDevEint(i)%ee
            enddo
            do i=1,mNCe
                do j=i+1,mNCe
                    write(17) Esum(pEww(i,j)),Esigma(pstdDevEww(i,j))
                    write(17) pEww(i,j)%nn, pstdDevEww(i,j)%nn
                    write(17) pEww(i,j)%ne, pstdDevEww(i,j)%ne
                    write(17) pEww(i,j)%ee, pstdDevEww(i,j)%ee
                enddo
            enddo
        endif
    end subroutine ce_write_block_mean

    !-----------------------!
    subroutine ce_finalize()
    !-----------------------!
        integer :: i,j
        do i=1,mNCe
            pEint(i) = mean(bEint(i))
            pstdDevEint(i) = stdDevMean(bEint(i))
            do j=i+1,mNCe
                pEww(i,j) = mean(bEww(i,j))
                pstdDevEww(i,j) = stdDevMean(bEww(i,j))
            enddo
        enddo
        call ce_write_ce_energies()
    end subroutine ce_finalize
    !-------------------------!
    subroutine ce_conv_check(t)
        !--------------------------!
        type(cee) :: tmp
        logical   :: t
        integer   :: i,j
        real*8    :: r
        real*8    :: ei(mNCe),ew(mNCe,mNCe)

        if(mConvCheck) then
            r=-1
            t = .true.
            do i=1,mNCe
                tmp = stdDevMean(bEint(i))
                ei(i) = Esigma(tmp)
                if(Esigma(tmp) >= mSigmaTol) then
                    t = .false.
                    if(r < Esigma(tmp)) r = Esigma(tmp)
                endif
                do j=i+1,mNCe
                    tmp = stdDevMean(bEww(i,j))
                    ew(i,j) = Esigma(tmp)
                    if(Esigma(tmp) >= mSigmaTol)then
                        t = .false.
                        if(r < Esigma(tmp)) r = Esigma(tmp)
                    endif
                enddo
            enddo
            if(logmode>=3)then
                 write(iul,'(A,F13.5)') "CE Convergence: ",r
                 write(iul,*) "-----------------------------"
                 write(iul,*) "Standard Deviations for CE"
                 write(iul,*)
                 do i=1,mNCe
                   write(iul,'(I2,X,F13.5)') i,ei(i)
                   do j=i+1,mNCe
                    write(iul,'(I2,X,I2,X,F13.5)') i,j,ew(i,j)
                   enddo
                 enddo
            endif
        else
            t = .false.
        endif
    end subroutine ce_conv_check

        !--------------------------!
    subroutine ce_write_ce_energies()
        !--------------------------!
        integer            :: i,j
        real*8             :: etotal,sigmatotal
        etotal = 0
        sigmatotal = 0
        if(MASTER)then
            write(iul,*)
            write(iul,*)
            write(iul,*) "          CE Energy Partitioning -- On the Fly --"
            write(iul,*) "********************************************************"
            write(iul,*)
            write(iul,*) "Internal Energies (Eint):"
            write(iul,*)
            write(iul,*) 'CE                 E        +/-      sigma'
            write(iul,*) '___________________________________________'
            do i=1,mNCe
                write(iul,'(A)') trim(ces(i)%name)
                write(iul,'(A,F13.5,A,F13.5)') "Eint: ",Esum(pEint(i)),'  +/-  ',Esigma(pstdDevEint(i))
                etotal = etotal + Esum(pEint(i))
                sigmatotal = sigmatotal + Esigma(pstdDevEint(i))**2
                if(logmode >= 2)then
                    write(iul,'(A,F13.5,A,F13.5)') "Ekin: ",pEint(i)%kin,'  +/-  ', pstdDevEint(i)%kin
                    write(iul,'(A,F13.5,A,F13.5)') "Vee : ",pEint(i)%ee,'  +/-  ', pstdDevEint(i)%ee
                    write(iul,'(A,F13.5,A,F13.5)') "Vne : ",pEint(i)%ne,'  +/-  ', pstdDevEint(i)%ne
                    write(iul,'(A,F13.5,A,F13.5)') "Vnn : ",pEint(i)%nn,'  +/-  ', pstdDevEint(i)%nn
                endif
            enddo
            write(iul,*)
            write(iul,*) "Interaction Energies (Eww):"
            write(iul,*)
            write(iul,*) 'CE1      CE2       E       +/-     sigma'
            write(iul,*) '____________________________________________'
            do i=1,mNCe
                do j=i+1,mNCe
                    write(iul,'(A,A,A)') trim(ces(i)%name)," ",trim(ces(j)%name)
                    write(iul,'(A,F13.5,A,F13.5)') "Eww:  ",Esum(pEww(i,j)),'  +/-  ',Esigma(pstdDevEww(i,j))
                    etotal = etotal + Esum(pEww(i,j))
                    sigmatotal = sigmatotal + Esigma(pstdDevEww(i,j))**2
                    if(logmode >= 2)then
                        write(iul,'(A,F13.5,A,F13.5)') "Ekin: ",pEww(i,j)%kin,'  +/-  ', pstdDevEww(i,j)%kin
                        write(iul,'(A,F13.5,A,F13.5)') "Vee : ",pEww(i,j)%ee,'  +/-  ', pstdDevEww(i,j)%ee
                        write(iul,'(A,F13.5,A,F13.5)') "Vne : ",pEww(i,j)%ne,'  +/-  ', pstdDevEww(i,j)%ne
                        write(iul,'(A,F13.5,A,F13.5)') "Vnn : ",pEww(i,j)%nn,'  +/-  ', pstdDevEww(i,j)%nn
                    endif
                enddo
            enddo
            write(iul,'(A,F13.4,A,F13.5)') "Total Energy: ",etotal," +/- ",sqrt(sigmatotal)
        endif

    end subroutine ce_write_ce_energies

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

        allocate(stdDevEkin(n))
        allocate(stdDevVee(n,n))
        allocate(stdDevVne(m,n))


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
            read(26,*) dummy,Ekin(i),stdDevEkin(i)
        enddo
        read(26,*)
        read(26,*)
        do i=1,mElectrons
            do j=1,mElectrons
                read(26,*) dummy,dummy,Vee(i,j),stdDevVee(i,j)
            enddo
        enddo


        read(26,*)
        read(26,*)
        do i=1,mElectrons
            do j=1,mCores
                read(26,*) dummy,dummy,Vne(j,i),stdDevVne(j,i)
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
            write(iul,*) Ekin(i)," ",stdDevEkin(i)
        enddo
        write(iul,*)
        write(iul,*) "Vee i j"
        do i=1,mElectrons
            do j=i+1,mElectrons
                write(iul,*) Vee(i,j)," ",stdDevVee(i,j)
            enddo
        enddo
        write(iul,*)
        write(iul,*) "Vne i j"
        do i=1,mElectrons
            do j=1,mCores
                write(iul,*) Vne(j,i)," ",stdDevVne(j,i)
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

end module ce_new_m
