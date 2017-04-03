module epart
    use RandomWalkerModule
    use statistics
    use global
    use findcores
    use assignment
    use ce_new_m
    use wfdata, only: Vnni,atoms,do_epart
    implicit none

    private
    public :: epart_init,epart_write,epart_add,epart_block,epart_get_energies,epart_reset

    integer,allocatable :: mAssign(:)

    !Energystat
    logical            ::  mEpart = .false.
    logical            ::  mCe = .false.
    !StepStat
    type(weightStat),allocatable  :: Ekin(:)
    type(weightStat),allocatable  :: Vne(:,:)
    type(weightStat),allocatable  :: Vee(:,:)
    !BlockStat
    type(simpleStat),allocatable  :: blockEkin(:)
    type(simpleStat),allocatable  :: blockVne(:,:)
    type(simpleStat),allocatable  :: blockVee(:,:)
    !Totalstat
    type(weightStat),allocatable  :: totalEkin(:)
    type(weightStat),allocatable  :: totalVne(:,:)
    type(weightStat),allocatable  :: totalVee(:,:)
    !Blocks
    integer                       ::  bl=1
    logical                       ::  mWriteDat = .false.

contains

    !----------------------------------
    subroutine epart_init(lines,nl)
        !----------------------------------
        integer, intent(in)           :: nl
        character(len=120), intent(in) :: lines(:)
        integer                       :: i,j,k,iflag,alstat,nces
        logical                       :: finda,e


        call assert(do_epart,"epart::init:: epart not set in wf file")
        allocate(mAssign(ne+1),stat=alstat)
        call assert(alstat==0,"epart::init:: allocation 0 failed")
        if (finda(lines,nl,'epart')) then
            if(finda(lines,nl,'cestat'))then
              call ce_init(lines,nl)
              mCE = .true.
            endif
            mEpart = .true.
            if(finda(lines,nl,'write_stat_file')) then
                mWriteDat = .true.
                open(17,file=trim(baseName)//'.dat',form='unformatted')
                call ce_get_nce(nces)
                write(17) ne,ncenter,nces
            endif
            allocate(Ekin(ne),Vne(ncenter,ne),Vee(ne,ne),stat=alstat)
            call assert(alstat==0,"epart::init:: allocation 1 failed")
            allocate(blockEkin(ne),blockVne(ncenter,ne),blockVee(ne,ne),stat=alstat)
            call assert(alstat==0,"epart::init:: allocation 2 failed")
            allocate(totalEkin(ne),totalVne(ncenter,ne),totalVee(ne,ne),stat=alstat)
            call assert(alstat==0,"epart::init:: allocation 3 failed")
            do i=1,ne
                call reset(Ekin(i))
                call reset(blockEkin(i))
                call reset(totalEkin(i))
                do j=1,ne
                    call reset(Vee(i,j))
                    call reset(blockVee(i,j))
                    call reset(totalVee(i,j))
                enddo
                do k=1,ncenter
                    call reset(Vne(k,i))
                    call reset(blockVne(k,i))
                    call reset(totalVne(k,i))
                enddo
            enddo
        endif


    end subroutine epart_init


    !------------------------------!
    subroutine epart_write()
        !-------------------------------!
        !Write Blockstatistic to output- and ce-file
        integer, parameter :: iu=12
        integer :: i,j,k,ierr
        real*8  :: x(ne),y(ne),z(ne)
        real*8  :: etot

        call assignment_get_reference(1,x,y,z)
        etot =0
        if (MASTER) then
            write(iul,*)
            write(iul,*) "Energy partitioning results:"
            write(iul,*)
            if (logmode>2) then
                write(iul,*) ""
                write(iul,*) ""
                write(iul,'(A)') "I    Block_E_kin(I)  +/-  stdDev "
                write(iul,'(A)') "---------------------------------"
                do i=1,ne
                    write(iul,'(I3,2F13.5)') i,mean(blockEkin(i)),stdDevMean(blockEkin(i))
                    etot = etot + mean(blockEkin(i))
                enddo

                write(iul,*) ""
                write(iul,*) ""
                write(iul,'(A)') "I J   Block_E_Vee(I,J)  +/-  stdDev  "
                write(iul,'(A)') "-------------------------------------"
                do i=1,ne
                    do j=i+1,ne
                        write(iul,'(2I3,2F13.5)') i,j,mean(blockVee(i,j)),stdDevMean(blockVee(i,j))
                        etot = etot + mean(blockVee(i,j))
                    enddo
                enddo

                write(iul,*) ""
                write(iul,*) ""
                write(iul,'(A)') "I J   Block_E_Ven(I,J)  +/-  stdDev "
                write(iul,'(A)') "------------------------------------"
                do i=1,ne
                    do j=1,ncenter
                        write(iul,'(2I3,2F13.5)') i,j,mean(blockVne(j,i)),stdDevMean(blockVne(j,i))
                        etot = etot + mean(blockVne(j,i))
                    enddo
                enddo

                write(iul,*) ""
                write(iul,*) ""
                write(iul,'(A)') "I J   E_Vnn(I,J)             "
                write(iul,'(A)') "-----------------------------"
                do i=1,ncenter
                    do j=i+1,ncenter
                        write(iul,'(2I3,F13.5)') i,j,Vnni(i,j)
                        etot = etot + Vnni(i,j)
                    enddo
                enddo
            endif
            open(iu,file=trim(baseName)//'.see',status='unknown')
            write(iu,*) ne," ",ncenter
            write(iu,*) "Geometry"
            do i=1,ncenter
                write(iu,'(i4,x,a2,x,3f12.5)') i,atoms(i)%elem,atoms(i)%cx,atoms(i)%cy,atoms(i)%cz
            enddo
            write(iu,*)
            write(iu,*) "Reference:"
            do i=1,ne
                write(iu,'(3f12.5)') x(i),y(i),z(i)
            enddo
            write(iu,*)
            write(iu,*) "Ekin(i):"
            do i=1,ne
                write(iu,'(i4,2f15.5)') i, mean(blockEkin(i)),stdDevMean(blockEkin(i))
            enddo
            write(iu,*)
            write(iu,*) "Vee(i,j):"
            do i=1,ne
                do j=i+1,ne
                    write(iu,'(2i4,2f15.5)') i,j,mean(blockVee(i,j)),stdDevMean(blockVee(i,j))
                enddo
            enddo
            write(iu,*)
            write(iu,*) "Ven(i,j)"
            do i=1,ne
                do j=1,ncenter
                    write(iu,'(2i4,2f15.5)') i,j,mean(blockVne(j,i)),stdDevMean(blockVne(j,i))
                enddo
            enddo
            write(iu,*)
            write(iu,*) "Vnn(i,j)"
            do i=1,ncenter
                do j=i+1,ncenter
                    write(iu,'(2i4,f15.5)') i,j,Vnni(i,j)
                enddo
            enddo
            close(iu)
            if (logmode>=2) then
               write(iul,*) "Energy partitioning results written to ",trim(baseName)//".see"
            end if
        endif
        if(mCe) call ce_finalize()
        call myMPIBarrier(ierr)
    end subroutine epart_write

    subroutine epart_reset()
        integer :: i,j,k

        do i=1,ne
            call reset(Ekin(i))
            call reset(blockEkin(i))
            call reset(totalEkin(i))
            do j=1,ne
                call reset(Vee(i,j))
                call reset(blockVee(i,j))
                call reset(totalVee(i,j))
            enddo
            do k=1,ncenter
                call reset(Vne(k,i))
                call reset(blockVne(k,i))
                call reset(totalVne(k,i))
            enddo
        enddo

    end subroutine epart_reset

    !--------------------------------------------
    subroutine epart_add(rwp,asgn)
        !-------------------------------------------
        !Get Energy contributions from RW and usw index mapping for adding
        type(RandomWalker),pointer  :: rwp(:)
        integer           :: asgn(:,:)
        real*8            :: Ekin_i(ne)
        real*8            :: Vee_i(ne,ne)
        real*8            :: Vne_i(ncenter,ne)
        real*8			  :: ee_tmp(ne,ne),ne_tmp(ncenter,ne),kin_tmp(ne)
        integer           :: i,j,k,w

        call assert(associated(rwp),'epart::add:: rw pointer not associated')
        do w=1,size(rwp)
            call assert(size(asgn(:,w))==ne+1,'epart::add:asgn size wrong')
            !Get arrays from random walker
            call Ekini(rwp(w),Ekin_i)
            call EVee(rwp(w),Vee_i)
            call EVne(rwp(w),Vne_i)

            !Fill fill up array, this is nec. due to assignment
            do i=1,ne
                do j=i+1,ne
                    Vee_i(j,i)=Vee_i(i,j)
                enddo
                Vee_i(i,i) = 0d0
            enddo

            if(asgn(size(asgn(:,w)),w) /=0)then
                call findcores_sort(Ekin_i)
                call findcores_sort(Vee_i)
                call findcores_sort_Vne(Vne_i)
            endif
            ! Do statistics
            do i=1,ne
                call addData(Ekin(i),Ekin_i(asgn(i,w)),wgt(rwp(w)))
                kin_tmp(i) = Ekin_i(asgn(i,w))
                do j=1,ne
                    call addData(Vee(i,j),Vee_i(asgn(i,w),asgn(j,w)),wgt(rwp(w)))
                    ee_tmp(i,j) = Vee_i(asgn(i,w),asgn(j,w))
                enddo
                do k=1,ncenter
                    call addData(Vne(k,i),Vne_i(k,asgn(i,w)),wgt(rwp(w)))
                    ne_tmp(k,i) = Vne_i(k,asgn(i,w))
                enddo
            enddo
            if(mCe) call ce_onthefly(ee_tmp,ne_tmp,kin_tmp,wgt(rwp(w)))
        enddo
    end subroutine epart_add
    !-----------------------
    subroutine epart_block()
        !-----------------------
        !Do block/total statistics and reset stepstat
        integer :: i,j,k
        real*8 :: vekin(ne),vvee(ne,ne),vvne(ncenter,ne)
        integer, parameter :: iub = 777

        do i=1,ne
            vekin(i) = meanAllNodes(Ekin(i))
            call addData(blockEkin(i),vekin(i))
            totalEkin(i) = totalEkin(i) + Ekin(i)
            do j=1,ne
                vvee(j,i) = meanAllNodes(Vee(i,j))
                call addData(blockVee(i,j),vvee(j,i))
                totalVee(i,j) = totalVee(i,j) + Vee(i,j)
            enddo
            do k=1,ncenter
                vvne(k,i) = meanAllNodes(Vne(k,i))
                call addData(blockVne(k,i),vvne(k,i))
                totalVne(k,i) = totalVne(k,i) + Vne(k,i)
            enddo
        enddo

        if(mWriteDat) call epart_write_stat_to_file()

        do i=1,ne
            call reset(Ekin(i))
            do j=1,ne
                call reset(Vee(i,j))
            enddo
            do k=1,ncenter
                call reset(Vne(k,i))
            enddo
        enddo

        !!if (MASTER) then
        !!   write(iub,*) 'new block:'
        !!   write(iub,'(10e18.10)') vekin
        !!   do i=1,ne
        !!      write(iub,'(10e18.10)') vvee(i+1:ne,i)
        !!   end do
        !!   do i=1,ne
        !!      write(iub,'(10e18.10)') vvne(:,i)
        !!   end do
        !!   flush(iub)
        !!end if

        !!if (mCe) call ce_blockstat()
        if (mCe) call ce_blockstat(mWriteDat)

    end subroutine epart_block

    subroutine epart_write_stat_to_file()
        integer :: i,j,k
        real*8  :: kin(ne),v_ee(ne,ne),v_ne(ncenter,ne)
        real*8  :: stdkin(ne),stdee(ne,ne),stdne(ncenter,ne)

        do i=1,ne
            kin(i) = meanAllNodes(Ekin(i))
            stdkin(i) = stdDevMeanAllNodes(Ekin(i))
            do j=1,ne
                v_ee(i,j) = meanAllNodes(Vee(i,j))
                stdee(i,j) = stdDevMeanAllNodes(Vee(i,j))
            enddo
            do k=1,ncenter
                v_ne(k,i) = meanAllNodes(Vne(k,i))
                stdne(k,i) = stdDevMeanAllNodes(Vne(k,i))
            enddo
        enddo

        if(MASTER)then
            write(17) bl
            do i=1,ne
               write(17) kin(i),stdkin(i)
            enddo
            do i=1,ne
                do j=1,ne
                    write(17) v_ee(i,j),stdee(i,j)
                enddo
            enddo
            do i=1,ncenter
                do k=1,ne
                    write(17) v_ne(i,k),stdne(k,i)
                enddo
            enddo
        endif
        bl = bl +1
    end subroutine epart_write_stat_to_file

    subroutine epart_get_energies(k,ee,vne)
        real*8,intent(out) :: k(ne)
        real*8,intent(out) :: ee(ne,ne)
        real*8,intent(out) :: vne(ncenter,ne)
        integer            :: i,j,a

        do i=1,ne
            k(i) = meanAllNodes(totalEkin(i))
            do j=1,ne
                ee(i,j) = meanAllNodes(totalVee(i,j))
            enddo
            do a=1,ncenter
                vne(a,i) = meanAllNodes(totalVne(a,i))
            enddo
        enddo

    end subroutine epart_get_energies


end module epart
