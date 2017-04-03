!
!               ----------------------
!                  A  M  O  L  Q  C
!               ----------------------
!
! Atomic and Molecular Quantum Monte Carlo Calculations
!
! initial version:
!  Arne Luechow, Penn State University, Feb. 1996
!
! main author:
!  Arne Luechow, RWTH Aachen University
!
! with contributions from:
! James B. Anderson
! Sebastian Manten, Christian Diedrich, Annika Bande, Tony Scott
! Rene Petz, Raphael Berner, Alexander Sturm, Kaveh Haghighi Mood
!
!
! main program
!

program amolqc

   use global
   use mainloopModule
   use InitModule, only: initAmolqc, finalizeAmolqc

   implicit none


   call initAmolqc()

   call mainloop()

   call finalizeAmolqc()

end program amolqc


