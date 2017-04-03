program test

  use statistics

  ! this tests generic interface to the many modules in statistics_m.f90

  type(simpleStat)      :: stat
  type(weightStat)      :: wstat
  type(IntHistogram)    :: ihist
  type(DoubleHistogram) :: dhist
  real*8                :: xMin=0.d0,xMax=100.d0
  integer               :: iMin=0,iMax=100,nBins=10
  integer               :: i

  call reset(stat)
  call reset(wstat)

  do i=1,5
     call addData(stat,dble(i))
     call addData(wstat,dble(i),1.d0)
  enddo

  write(*,*) mean(stat),mean(wstat)
  write(*,*) dataCount(stat), dataCount(wstat)

  call createHistogram(ihist,iMin,iMax,nBins)
  write(*,*) 'ihist created with: ',iMin,iMax,nBins
  call createHistogram(dhist,xMin,xMax,nBins)
  write(*,*) 'dhist created with: ',xMin,xMax,nBins

  do i=0,100
     call addData(ihist,i)
     call addData(dhist,dble(i))
  enddo

  write(*,*) mean(ihist),mean(dhist)
  write(*,*) dataCount(ihist), dataCount(dhist)

  write(*,*) 'ihist:'
  call writeHistogram(ihist,6)
  write(*,*) 'dhist:'
  call writeHistogram(dhist,6)
  write(*,*) 'dhist:'
  call writeMeanHistogram(dhist,6)

end program test
