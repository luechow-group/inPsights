! Module with routines to control the iteration of coc calculations
module coc_control
    use qmc
    use coc
    use global
    use rwsamplemodule
    use assign
    use findCoreModule, only: initCoreLists, findCoreElecs
    use statistics
    implicit none

    private
    public :: center_of_charge_run

contains

    !----------------------------------------------!
    subroutine center_of_charge_run(lines,nl,sample)
    !----------------------------------------------!

#ifdef PARALLEL
        include 'mpif.h'
#endif
        character(len=120), intent(in) :: lines(:)
        integer, intent(in)           :: nl
        type(RWSample), intent(inout) :: sample

        integer                       :: iflag,i,j,ierr
        logical                       :: finda,changed,coresep
        integer                       :: iters,discard,blocklen
        integer                       :: naCore, nbCore    ! alpha, beta elecs at nucleus (not varied)
        integer                       :: aCoreList(ncenter), bCoreList(ncenter)
        logical                       :: convSinglePoint,cocAverage
        real*8                        :: maxSigma,meanSigma,macrotol,microtol
        real*8                        :: coc_it(ne,3)    ! COCs of the current step
        real*8                        :: coc_save(ne,3)  ! COCs of the last step
        real*8, allocatable           :: coc_list(:,:,:) ! all iterations
        real*8                        :: maxError,meanError

        macrotol = 1.d-3
        call getdbla(lines,nl,'coc_tol=',macrotol,iflag)
        microtol = macrotol / 10d0
        call getdbla(lines,nl,'coc_mtol=',microtol,iflag)
        iters = 3
        call getinta(lines,nl,'iters=',iters,iflag)
        blocklen = 100
        call getinta(lines,nl,'steps=',blocklen,iflag)
        coresep = .true.
        if (finda(lines,nl,'nocoresep')) coresep = .false.

        allocate(coc_list(ne,3,iters))
        coc_list = 0

        call initCoreLists(naCore,aCoreList,nbCore,bCoreList)

        call coc_create(blocklen,microtol)

        call qmc_init(lines,nl,sample)

        coc_save = 0
        coc_it   = 0
        do i=1,iters
            call qmc_run(sample)
            convSinglePoint = qmc_cocIsConverged()
            
            ! careful: only MASTER gets coc results!
            call coc_getCOC(coc_it(:,1),coc_it(:,2),coc_it(:,3))
            call coc_getErrors(maxSigma,meanSigma)
            call coc_reset()

            maxError  = maxval(abs(coc_it - coc_save))
            meanError = sum(abs(coc_it - coc_save))/(3*ne)

            if (MASTER .and. logmode>=2) then
               if (convSinglePoint) then
                  write(iul,'(a,i4,a,f12.5)') " --- center of charge iter=",i," is converged with maxSigma=",maxSigma
               else
                  write(iul,'(a,i4,a,f12.5)') " --- center of charge iter=",i," is NOT converged with maxSigma=",maxSigma
               end if
               write(iul,'(2(A,F12.5)/)') "   iteration convergence: max error = ",maxError,"  mean error = ",meanError
               do j=1,ne
                  write(iul,'(3f12.5)') coc_it(j,1),coc_it(j,2),coc_it(j,3)
               end do
            end if

            if (MASTER .and. coresep) then
               ! find core elecs and put them at nuc position
               call findCoreElecs(naCore,aCoreList,nbCore,bCoreList,coc_it(:,1),coc_it(:,2),coc_it(:,3))
            end if
#ifdef PARALLEL
            call MPI_BCAST(coc_it(1,1),3*ne,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
#endif
            call assign_setRef(coc_it(:,1),coc_it(:,2),coc_it(:,3),changed)
            coc_list(:,:,i) = coc_it
            if (MASTER.and.changed.and.logmode>=2) then
               write(iul,*) "Warning: core electrons in COC reference have been permuted"
            end if
            coc_save = coc_it
            !!!! this requires BCAST of maxError:  if (maxError<=macroTol) exit
        enddo

        call coc_writeToFile(coc_list,iters)

        if (MASTER .and. logmode>=2) then
           if (maxError<=macroTol) then
               write(iul,'(/a/3(a,f15.5))') " --- center of charge result: converged", &
                   " maxError = ",maxError," meanError = ",meanError," macroTol=",macroTol
           else
               write(iul,'(/a/3(a,f15.5))') " --- center of charge result: NOT converged", &
                   " maxError = ",maxError," meanError = ",meanError," macroTol=",macroTol
           end if
        end if

        deallocate(coc_list)

    end subroutine center_of_charge_run


    !--------------------------------------------!
    subroutine coc_write_center_of_charge(x,y,z)
        !---------------------------------------------!
        real*8,intent(in)     :: x(ne),y(ne),z(ne)
        integer               :: i,j


        write(iul,'(A)')
        write(iul,'(A)')
        write(iul,'(A)')  "Center of Charge for SEDs"
        write(iul,'(A)')  "-------------------------"
        write(iul,'(A)')
        do i=1,ne
            write(iul,'(3(F15.8))') x(i),y(i),z(i)
        enddo
        write(iul,*)


    end subroutine coc_write_center_of_charge


    !--------------------------------!
    subroutine coc_writeToFile(list,n)
    !--------------------------------!
        ! overwrite the *.ref File
        ! write list backwards!
        real*8,intent(in)     :: list(:,:,:)
        integer,intent(in)    :: n     ! last entry
        integer :: i,j

        call assert(size(list,3)>=n,"coc_writeToFile: illegal size")
        !open(12,file=trim(baseName)//'.ref',status='replace')
        open(12,file=trim(baseName)//'.ref')
        write(12,'(i5,a)') ncenter, " nuclei:"
        do i=1,getNNuc()
           write(12,'(i4,x,a2,x,3f14.7)') i,atoms(i)%elem,atoms(i)%cx,atoms(i)%cy,atoms(i)%cz
        end do
        write(12,'(i3,a)') size(list,3),"  COC"
        do i=n,1,-1
           write(12,'(a,i3,a)') "COC:  ",n-i+1,"  1  F(COC):  0.0  found: 0  0 0 0 0"
           write(12,'(I5)') ne
           do j=1,ne
              write(12,'(3F12.5)') list(j,1,i),list(j,2,i),list(j,3,i)
           end do
        end do
        close(12)

    end subroutine coc_writeToFile

    subroutine coc_averaging(d,maxe,rms)
      integer,intent(in) :: d
      real*8,intent(out) :: rms,maxe
      integer :: i,j
      real*8  :: check(ne,3),dist
      do i=1,ne
         do j=1,3
            !!call addData(coc_stat(i,j),coc_it(i,j))
         enddo
      enddo

      !if(d-mDiscard >2) then
        write(iul,*)
        write(iul,'(A)')  "Step averraged Center of Charge for SEDs"
        write(iul,'(A)')  "----------------------------------------"
        write(iul,'(A)')
        rms = 0
        do i=1,ne
            !!write(iul,'(3(F15.8))') mean(coc_stat(i,1)),mean(coc_stat(i,2)),mean(coc_stat(i,3))
            !!if(d-mDiscard > 3) then
                !!check(i,1) = stddevmean(coc_stat(i,1));check(i,2) = stddevmean(coc_stat(i,2));check(i,3) = stddevmean(coc_stat(i,3))
                !!rms = rms + check(i,j)**2
            !!endif
        enddo
        rms = sqrt(rms/(3*ne))
        maxe = maxval(check)
        write(iul,*)
        !!write(iul,'(A,X,F13.5,X,A,X,F13.5)') "Max displacement",maxe,"Tolerance:",mMacroTol
        !!write(iul,'(A,X,F13.5,X,A,X,F13.5)') "Max displacement",rms,"Tolerance:",mMacroTol
      !endif
    end subroutine coc_averaging

end module coc_control
