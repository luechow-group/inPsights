!========================================!
!                                        !
!   Modul zum 3-dimensionalen binning    !
!                                        !
!========================================!

! Diese Version definiert eine fast echte
! Klasse, alle Methoden haben "Zeiger" auf Objekt


MODULE bins3D
    use statistics

    implicit none

    type Bin3D
        real*8, pointer  :: bins(:,:,:)        ! actual bins
        real*8           :: ax,bx,ay,by,az,bz  ! lower/upper borders
        real*8           :: hx,hy,hz           ! step size
        integer          :: n                  ! # of bins per dimension
    end type Bin3D

    integer              :: ncenter2
    real*8,allocatable   :: cx2(:),cy2(:),cz2(:)
    integer,allocatable  :: zidx(:)
    character*20         :: filename=''

contains

    subroutine setbasename(n)
        character*20 :: n
        filename = n
    end subroutine setbasename
    !==============================================================
    subroutine setncenter(n)
        integer,intent(in) :: n
        ncenter2 = n
        allocate(cx2(n),cy2(n),cz2(n),zidx(n))
    end subroutine setncenter
    !=============================================================
    subroutine setcubevals(x,y,z,idxz)
        real*8,intent(in)  :: x(:),y(:),z(:)
        integer,intent(in) :: idxz(:)
        cx2=x(1:ncenter2)
        cy2=y(1:ncenter2)
        cz2=z(1:ncenter2)
        zidx=idxz(1:ncenter2)
    end subroutine setcubevals
    !=============================================================
    integer function getnobins(bin)
        type(Bin3d), intent(in) :: bin
        getnobins = bin%n
    end function getnobins
    !=============================================================!
    subroutine binsAllocate(bin,ndim)

        type(Bin3D),intent(inout)  :: bin
        integer ie,ndim

        allocate(bin%bins(0:ndim-1,0:ndim-1,0:ndim-1),stat=ie)
        if (ie /= 0) then
            write(*,*) ' ERROR in binsAllocate: allocation error'
            call exit(1)
        endif
        bin%n = ndim

    end subroutine binsAllocate

    !===========================================================

    subroutine binsInit(bin,ax,bx,ay,by,az,bz)

        type(Bin3D),intent(inout)  :: bin
        real*8,intent(in)       :: ax,bx,ay,by,az,bz
   
    
        bin%ax = ax; bin%bx = bx
        bin%ay = ay; bin%by = by
        bin%az = az; bin%bz = bz
        bin%hx = (bx-ax)/bin%n
        bin%hy = (by-ay)/bin%n
        bin%hz = (bz-az)/bin%n
       
    
        bin%bins = 0.d0

    end subroutine binsInit

!    !===========================================================!
!
!    subroutine binsAddData(bin,x,y,z,value,kx,ky,kz)
!
!        type(Bin3D), intent(inout)  :: bin
!       real*8 , intent(in)         :: x,y,z,value
!        integer, intent(out)        :: kx,ky,kz


!        kx = max(0,int((x - bin%ax)/bin%hx))
!        kx = min(kx,bin%n-1)
!        ky = max(0,int((y - bin%ay)/bin%hy))
!        ky = min(ky,bin%n-1)
!        kz = max(0,int((z - bin%az)/bin%hz))
!        kz = min(kz,bin%n-1)

!        bin%bins(kx,ky,kz) = bin%bins(kx,ky,kz) + value

!    !   write(*,*) kx,ky,kz,value,bin%bins(kx,ky,kz)

!    end subroutine binsAddData

    !===========================================================

    subroutine binsAddData(bin,x,y,z,value,kx,ky,kz)

        type(Bin3D), intent(inout)  :: bin
        real*8 , intent(in)         :: x,y,z,value
        integer, intent(out)        :: kx,ky,kz
    

        kx = int((x - bin%ax)/bin%hx)
        kx = min(kx,bin%n-1)
        ky = max(0,int((y - bin%ay)/bin%hy))
        ky = min(ky,bin%n-1)
        kz = max(0,int((z - bin%az)/bin%hz))
        kz = min(kz,bin%n-1)
   
        if (kx>=0 .and. kx<=bin%n-1 .and. ky>=0 .and. ky<=bin%n-1 .and. kz>=0 .and. kz<=bin%n-1) &
           bin%bins(kx,ky,kz) = bin%bins(kx,ky,kz) + value
        !!! else
        !!!    bin%count_outside_box += 1
        !!! endif
    !   write(*,*) kx,ky,kz,value,bin%bins(kx,ky,kz)

    end subroutine binsAddData

    !=========================================================

    real*8 function binsGetData(bin,ix,iy,iz)

        type(Bin3D)  :: bin
        integer      :: ix,iy,iz,n

        n = bin%n - 1
        if (ix<0.or.ix>n.or.iy<0.or.iy>n.or.iz<0.or.iz>n) then
            write(*,*) ' ERROR in binsGetData: overflow error'
            call exit(1)
        endif
        binsGetData = bin%bins(ix,iy,iz)

    end function binsGetData

    !===========================================================
    subroutine binsSetData(bin,ix,iy,iz,val)

        type(Bin3D)  :: bin
        integer      :: ix,iy,iz,n
        real*8       :: val

        n = bin%n - 1
        if (ix<0.or.ix>n.or.iy<0.or.iy>n.or.iz<0.or.iz>n) then
            write(*,*) ' ERROR in binsGetData: overflow error'
            call exit(1)
        endif
        bin%bins(ix,iy,iz) = val

    end subroutine binsSetData
     !===========================================================

    subroutine binsGetCenter(bin,kx,ky,kz,cx,cy,cz)

        ! returns center of bin (ix,iy,iz)
        type(Bin3D)     :: bin
        real*8, intent(out)      :: cx,cy,cz
        integer, intent(in)      :: kx,ky,kz
        cx = bin%ax + (kx+0.5)*bin%hx
        cy = bin%ay + (ky+0.5)*bin%hy
        cz = bin%az + (kz+0.5)*bin%hz

    end subroutine binsGetCenter

    !===============================================================
    subroutine binallnodes(bin)
        type(Bin3D),intent(inout) :: bin
        real*8,allocatable :: sendbuf(:),recvbuf(:)
        integer            :: c,i,j,k

        allocate(sendbuf(bin%n**3),recvbuf(bin%n**3))
        c = 1
        do i =0,bin%n-1
            do j=0,bin%n-1
                do k=0,bin%n-1
                    sendbuf(c) = bin%bins(i,j,k)
                    c = c + 1
                enddo
            enddo
        enddo
        recvbuf = 0
        call myMPIAllReduceSumDouble(sendbuf,recvbuf,(bin%n)**3)
        c = 1
        do i =0,bin%n-1
            do j=0,bin%n-1
                do k=0,bin%n-1
                    bin%bins(i,j,k) = recvbuf(c)
                    c = c + 1
                enddo
            enddo
        enddo
    end subroutine binallnodes

    !===============================================================!
    subroutine writetofile(bin,idx)
        type(Bin3D),intent(in) :: bin
        integer       :: i,j,k,idx,n_all
        integer       :: i1,i2,i3
        real*8        :: vx,vy,vz
        character(len=3) :: idxc
  
        write(idxc,'(I3)') idx
        idxc = adjustl(idxc)
        open(10,file=trim(filename)//trim(idxc)//'.cube',status='unknown')

        write(10,*) "cube file generated by bin_m"
        write(10,*) "outer loop x, middle y, inner z"
        write(10,100) ncenter2,bin%ax+bin%hx/2,bin%ay+bin%hy/2,bin%az+bin%hz/2
        write(10,100) bin%n,dble(bin%hx),dble(0),dble(0)
        write(10,100) bin%n,dble(0),dble(bin%hy),dble(0)
        write(10,100) bin%n,dble(0),dble(0),dble(bin%hz)

        do i=1,ncenter2
            write(10,101) zidx(i),dble(0),cx2(i),cy2(i),cz2(i)
        end do

        do i1 =0,bin%n-1
            do i2=0,bin%n-1
                write(10,102) (bin%bins(i1,i2,i3),i3=0,bin%n-1)
            enddo
        enddo

        close(10)

100     format((I5,3F12.6))
101     format(((I5,4F12.6)))
102     format((6E13.5))


    end subroutine writetofile

end module bins3D


