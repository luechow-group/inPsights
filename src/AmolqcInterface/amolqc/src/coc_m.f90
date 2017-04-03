! module for calculating the center of charges for all electrons

module coc
    use newstatistics
    use RandomWalkerModule
    use wfdata
    implicit none

    private
    public :: coc_create,coc_reset,coc_destroy,coc_add,coc_isConverged,coc_getErrors,coc_getCOC

    type(blockvectorstat) :: mCOC
    real*8                :: mMicroTol = 1d-4                ! convergence criterium for max(sigma)

contains
    !--------------------------------
    subroutine coc_create(blockLen,t)
    !--------------------------------
        integer, intent(in) :: blockLen    ! blockLen for statistics
        real*8, optional    :: t           ! convergence criterium (max(sigma)<t)

        if (present(t)) mMicroTol = t
        call assert(mMicroTol > 0,"coc_create: negative tolerance")
        call mCOC%create(3*ne,blockLen)
    end subroutine coc_create

    !---------------------!
    subroutine coc_reset(t)
    !---------------------!
        real*8, optional  :: t    ! convergence criterium (max(sigma)<t)
        integer :: i,j

        if (present(t)) mMicroTol = t
        call assert(mMicroTol > 0,"coc_reset: illegal negative tolerance")
        call assert(mCOC%iscreated(),"coc_reset: call coc_create first")
        call mCOC%reset()
    end subroutine coc_reset

    !----------------------!
    subroutine coc_destroy()
    !----------------------!
        if (mCOC%iscreated()) call mCOC%destroy()
    end subroutine coc_destroy


    !--------------------------!
    subroutine coc_add(rwp,asgn)
    !--------------------------!
        type(RandomWalker),pointer,intent(in)  :: rwp(:)
        integer,intent(in)  :: asgn(:,:)
        integer             :: i,w
        integer             :: idx
        real*8 x(ne),y(ne),z(ne),v(3*ne)
        call assert(associated(rwp),'(coc_add): rw pointer not associated')
        do w=1,size(rwp)
            call pos(rwp(w),x,y,z)
            forall (i=1:ne)
                v(i)      = x(asgn(i,w))
                v(i+ne)   = y(asgn(i,w))
                v(i+2*ne) = z(asgn(i,w))
            end forall
            call mCOC%add(v)
        enddo

    end subroutine coc_add


    !--------------------------------!
    logical function coc_isConverged()
    !--------------------------------!

    ! CAREFUL: only MASTER returns convergence

        real*8       :: sigma(3*ne),maxsigma

        coc_isConverged = .false.
        sigma = mCOC%stddev()
        if (MASTER) then
            coc_isConverged = (maxval(sigma) < mMicroTol)
        end if
    end function coc_isConverged

    !------------------------------------------!
    subroutine coc_getErrors(maxError,meanError)
    !------------------------------------------!

        ! CAREFUL: only MASTER returns result
        real*8, intent(out) :: maxError
        real*8, intent(out) :: meanError
        real*8       :: sigma(3*ne)

        maxError = 0
        meanError = 0
        sigma = mCOC%stddev()
        if (MASTER) then
           maxError = maxval(sigma)
           meanError = sum(sigma) / size(sigma)
        end if
    end subroutine

    !--------------------------!
    subroutine coc_getCOC(x,y,z)
    !--------------------------!

        ! MASTER return current mean values 
        real*8, intent(inout) :: x(:),y(:),z(:)
        real*8 mean(3*ne)
        call assert(size(x)==ne,"(coc_getCOC): arguments not properly allocated")
        mean = 0
        mean = mCOC%mean()
        x = mean(1:ne)
        y = mean(ne+1:2*ne)
        z = mean(2*ne+1:3*ne)
    end subroutine coc_getCOC

end module coc
