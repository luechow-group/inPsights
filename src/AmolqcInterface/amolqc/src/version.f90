
MODULE versionModule

!use global,only: getEnvironmentVariableName
use MachModule, only: mygetEnv
! this file is edited after each commit
implicit none
private


public :: getVersion

character(len=60)           :: versionString

contains

subroutine getVersion(ver)

 character(len=60),intent(out)     :: ver
 character(len=60)                 :: versionString
 character(len=180)                :: basispath
 character(len=190)                :: filePath
 integer                           :: io

 10 format(A60)
  call mygetEnv('AMOLQC',basispath)
  call assert(len(trim(basispath)) < 179,"ecpinputex: amolqc path length exceeds definition")
  filePath = trim(basispath)//"/version"
  open(24,file=filePath,status='old',iostat=io)

      if (io /= 0) then
         call error('version could not be opened')
      end if
      read(24,10) versionString
      close(24)
      ver=versionString
end subroutine getVersion

end MODULE versionModule
