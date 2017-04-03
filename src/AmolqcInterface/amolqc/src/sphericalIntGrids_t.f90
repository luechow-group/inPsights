
program test
   use SphericalIntegrationGrids
   implicit none

   integer n,r,idx0,nGridPoints
   real*8 I
   real*8 v(3),w,sumw

   do n=1,size(RuleIdx)-1
      idx0 = RuleIdx(n) 
      nGridPoints = RuleIdx(n+1) - idx0
      print*, n,nGridPoints,' rule points'
      I = 0.d0
      sumw = 0.d0
      do r=1,nGridPoints
         v = AllGridPoints(:,idx0+r)
         w = AllWeights(idx0+r)
         I = I + w*testfunc(v)
         sumw = sumw + w
         !!!print*,w,v(1)**2 + v(2)**2 + v(3)**2
      enddo
      print*, 'rule sumw int ',n,sumw,I
   enddo

contains

   real*8 function testfunc(v)
      real*8, intent(in) :: v(3)
      testfunc = -v(1)**2 + v(2)**2 - v(3)**2
   end function testfunc

end program test