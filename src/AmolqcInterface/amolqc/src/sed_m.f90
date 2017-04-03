! R.P.: rwbinstat should only contain, init,add
!
module sed

    use bins3d
    use global
    use findcores
    use RandomWalkerModule

    implicit none
    private
    public:: sed_add,sed_init,sed_write,sed_write_init

    !Simple Binstat
    type(bin3d),allocatable           ::    bin(:)            ! grids for the seds
    real*8                            ::    mVal=1d0          ! currently this is the weight of the rw
    real*8                            ::    grid_length       ! Half(!) the length of the "cube" (i.e. cuboid)
    real*8                            ::    ax,ay,az,bx,by,bz ! coordinates of the cuboid
    integer                           ::    mNbin =80         ! bins per grid dimension
    type(coord)                       ::    nall              ! Electron configuration
    logical                           ::    mWriteXYZ = .false.   ! write xyz files with permuted walkers
    integer                           ::    mLCount = 0
  
    ! Bin Manipulation
    logical                           ::    mNormBin = .false.   !Flag for sed norm. calc.
    logical                           ::    mBinOverlap = .false.!Flags  for sed overlap calc.

contains

    !------------------------------------
    subroutine sed_init(lines,nl)
    !-----------------------------------
        integer, intent(in)           :: nl
        character(len=120), intent(in) :: lines(:)
        integer                       :: i,iflag,j
        logical                       :: finda
        real*8                        :: cy(ncenter),cx(ncenter),cz(ncenter)
        integer                       :: ElemIdx(ncenter)

        !Reading everything for the grid initialization
 
        mWriteXYZ = finda(lines,nl,'write_xyz')
        if (mWriteXYZ) then 
            write(500+mytid,'(i5,a)') ncenter,' nuclei:'
            do i=1,ncenter
                write(500+mytid,'(i4,x,a2,x,3f14.7)') i,atoms(i)%elem,atoms(i)%cx,atoms(i)%cy,atoms(i)%cz
            enddo
            write(500+mytid,'(a)') ' UNKNOWN_SIZE  xyz'
            mLCount = 0
        endif

        call getdbla(lines,nl,'grid=',grid_length,iflag)
        if (iflag /= 0) then
            call getdbla(lines,nl,'ax=',ax,iflag)
            if (iflag /= 0) call abortp("(sed_init): either grid or ax,bx,.. required for sed calculation")
            call getdbla(lines,nl,'bx=',bx,iflag)
            call getdbla(lines,nl,'ay=',ay,iflag)
            call getdbla(lines,nl,'by=',by,iflag)
            call getdbla(lines,nl,'az=',az,iflag)
            call getdbla(lines,nl,'bz=',bz,iflag)
            if (ax > bx .or. ay > by .or. az > bz) call abortp('(sed_init): cube definition requires a > b')
        else
            ax = -1d0*grid_length; ay = -1d0*grid_length; az = -1d0*grid_length
            bx = grid_length; by = grid_length; bz = grid_length
        endif
        !!! NOTE: atoms_getBox is available to calculate optimal box

        call getinta(lines,nl,'nbin=',mNbin,iflag)

        allocate(bin(ne))

        do i=1,ne
            call binsAllocate(bin(i),mNbin)             ! create bins
            call binsInit(bin(i),ax,bx,ay,by,az,bz) ! init. bins
        enddo

        ! Set Informaton for Export to Cube
        do i=1,ncenter
          cx(i) = atoms(i)%cx;cy(i) = atoms(i)%cy;cz(i) = atoms(i)%cz
          ElemIdx(i)=atoms(i)%elemIdx
        enddo
        call setncenter(ncenter)
        call setcubevals(cx,cy,cz,elemIdx)
        call setbasename(baseName)

        !Bin Manipulation
        if (finda(lines,nl,'sed_norm')) mNormBin =.true.
        if (finda(lines,nl,'sed_overlap')) mBinOverlap =.true.

        !Allocate
        allocate(nall%x(ne),nall%y(ne),nall%z(ne))
    end subroutine sed_init

    !----------------------------------
    subroutine sed_write_init(iu)
    !----------------------------------
        integer, intent(in) :: iu
        integer :: ierr

        write(iu,'(A/)') "    sed parameters:"
        write(iu,'(1X,A21,F10.4,5X,A21,I4)') "half box length =",grid_length,"# bins =",mNbin
        write(iu,*)
    end subroutine sed_write_init

    !---------------------------------
    subroutine sed_write()
        !---------------------------------
        integer :: i,k

        if(nproc > 1) then
            do i=1,ne
                call binallnodes(bin(i))
            enddo
        endif
        if(mytid == 0) then
            if(mNormBin) call sed_norm()
            if(mBinOverlap) call sed_overlap()
            do i=1,ne
                call writetofile(bin(i),i)
            enddo
        endif

    end subroutine sed_write


    !Routines for adding the configuration to bins

    !------------------------------------
    subroutine sed_add(rwp,asgn)
        !------------------------------------
        type(RandomWalker),pointer,intent(in)  :: rwp(:)
        integer,intent(in)  :: asgn(:,:)
        integer             :: kx,ky,kz
        integer             :: i,w
        integer             :: idx

        call assert(associated(rwp),'sed::add:: rw pointer not associated')
        do w=1,size(rwp)
            call pos(rwp(w),nall%x,nall%y,nall%z)
            mVal = wgt(rwp(w))
            !Add the x,y,z from random walker to the bins by index mapping
            do i=1,ne
                idx = asgn(i,w)
                call binsAddData(bin(i),nall%x(idx),nall%y(idx),nall%z(idx),mVal,kx,ky,kz)
            enddo

            if (mWriteXYZ) then
                mLCount = mLCount + 1
                write(500+mytid,'(a,i9)') 'xyz:   1 ',mLCount
                write(500+mytid,'(i5)') ne
                do i=1,ne
                    idx = asgn(i,w)
                    write(500+mytid,'(3f15.6)') nall%x(idx),nall%y(idx),nall%z(idx)
                enddo
            endif

        enddo


    end subroutine sed_add

    !---------------------------
    subroutine sed_norm()
        !------------------------------
        integer      :: n,i,j,k,b
        real*8       :: norm,val

        ! since sed : Norm**2 INT psi**2 dt -> N = 1/sqrt(INT psi**2 dt)
        ! this should be the correct way
        do b=1,ne
            n = getnobins(bin(b))
            norm=0d0
            do i=0,n-1
                do j=0,n-1
                    do k=0,n-1
                        norm = norm + binsGetData(bin(b),i,j,k)*binsGetData(bin(b),i,j,k)
                    enddo
                enddo
            enddo

            norm = sqrt(norm)

            do i=0,n-1
                do j=0,n-1
                    do k=0,n-1
                        val = binsGetData(bin(b),i,j,k)/norm
                        call binsSetData(bin(b),i,j,k,val)
                    enddo
                enddo
            enddo
        enddo
    end subroutine sed_norm

    !---------------------------------
    subroutine sed_overlap()
        !---------------------------------
        integer      :: n1,n2
        integer      :: b1,b2
        integer      :: i,j,k
        integer      :: i2,j2,k2
        real*8       :: overlap
        real*8       :: x,y,z,x2,y2,z2
        real*8       :: r12

        write(iul,'(A)') ""
        write(iul,'(A)') ""
        write(iul,'(A)') "      Calculating Overlap        "
        write(iul,'(A)') "sed(i) sed(j)     O(i,j)         "
        write(iul,'(A)') "---------------------------------"
        do b1=1,ne
            do b2=b1+1,ne
                n1=getnobins(bin(b1))
                n2=getnobins(bin(b2))
                if(n1 /= n2) call abortp("sedstat_bin_num_int::overlap:: bins must have same size")
                overlap = 0d0
                do i=0,n1-1
                    do j=0,n1-1
                        do k=0,n1-1
                            overlap = overlap + (binsGetData(bin(b1),i,j,k)) * (binsGetData(bin(b2),i,j,k))
                        enddo
                    enddo
                enddo
                write(iul,'(I3,I3,F13.5)') b1,b2,overlap
            enddo
        enddo
    end subroutine sed_overlap


end module sed
