module OptimizeParamsETRNRModule

! optimize the energy for a fixed sample using the 
! trust-region Newton-Raphson alg. 
! uses nmtr subroutine by More´ and Sorensen
! see: SIAM Journal on Scientific and Statistical Computing, 4 (1983), pages 553-572

   use global
   use subloopModule, only: subloop
   use RWSampleModule
   use ElocAndPsiTermsENRModule
   use WFParamsModule
   use utilsmodule, only: dnmtr
   implicit none

   private
   public :: eminTRNR_optimizeSample


contains



   subroutine eminTRNR_optimizeSample(lines,nl,WFP,sample,converged)
   !---------------------------------------------------------------!

   integer, intent(in)                  :: nl
   character(len=120), intent(in)       :: lines(nl)
   type(WFParamDef), pointer              :: WFP
   type(RWSample), intent(inout)        :: sample  ! (fixed) sample for optimization
   logical, intent(out)                 :: converged
   integer                              :: nParams, np
   real*8, allocatable                  :: p(:)               ! parameter vector
   real*8, allocatable                  :: g(:),g1(:),H(:,:)  ! gradient and Hessian
   real*8                               :: e0,var,pe0,pvar
   real*8, allocatable                  :: fi(:),ELi(:),fiEL(:)
   real*8, allocatable                  :: fifj(:,:),fifjEL(:,:),fiELj(:,:),fij(:,:),fijEL(:,:)
   real*8, allocatable                  :: A(:,:),B(:,:),D(:,:)
   real*8, allocatable                  :: eval(:),evec(:,:)
   integer i,j,info,n,ierr,nSize, gmode, NRMode, iflag
   real*8                               :: maxVar, lambda(6), lambdaOpt, eRef
   type(ElocAndPsiTermsENR)             :: EPsiTENR
   character(len=40)                    :: subName
   ! nmtr variables
   character(len=60)                    :: task
   real*8, allocatable                  :: wa(:)
   real*8                               :: dsave(7),delta, frtol, fatol, fmin, gnorm, dnrm2
   integer                              :: isave(6), nfev, nfev1, ngev, maxfev
   logical scale,forceNewSample

   converged = .true.

   call internal_readInput()        ! internal subroutine after contains

   eRef = 0.d0
   call ElocAndPsiTermsENR_create(EPsiTENR,eRef,wfp)

   if (logmode>=2) then
      write(iul,'(/A/)') '   - -  energy minimization using trusted-region Newton-Raphson: initialization  - -'
      write(iul,'(a,i3)') ' with parameters: method = ',NRMode
   endif

   np = ElocAndPsiTermsENR_nParams(EPsiTENR)
   WFP => ElocAndPsiTermsENR_getWFP(EPsiTENR)
   call assert(np>0,'eminNR_optimizeSample: no parameters')
   allocate(p(np),g(np),g1(np),H(np,np))
   allocate(fi(np),ELi(np),fiEL(np))
   allocate(fifj(np,np),fifjEL(np,np),fiELj(np,np),fij(np,np),fijEL(np,np))
   allocate(A(np,np),B(np,np),D(np,np))
   allocate(eval(np),evec(np,np))
   allocate(wa(7*np))

   p = 0
   wa = 0
   nfev = 0; nfev1 = 0; ngev = 0
   fmin = -1.d30

   if (logmode >= 2) write(iul,*) ' starting wf parameter optimization with optType=',WFP%optType


   ! Start of search.
   task = 'START'
   p = wfparams_get(WFP)
   forceNewSample = .false.

   do
      if (logmode >= 2) write(iul,'(/2a/a,g15.5)') '***** trnr: task=',trim(task),' delta=',delta
      if (task=='NEWX') then
         !if (logmode >= 2) write(iul,*) ' ** NEWX approx e0:',e0
         !!if (logmode >= 2) call internal_writeCurrentVectorToLogFile()
         ! equilibrate sample only before function calls
         ! call "subroutine" subName in .in or macro subName.cmd that should contain
         ! code for equilibrating the sample with the new wave function 
         !call subloop(subName,sample)
         !forceNewSample = .false.
         !nfev = nfev + 1
      end if
      if (task=='F' .or. task=='START') then
         if (task=='F') then
            call subloop(subName,sample)
         end if
         call ElocAndPsiTermsENR_reset(EPsiTENR)
         call internal_calcEPsiTerms()

         e0 = ElocAndPsiTermsENR_EmeanALL(EPsiTENR)
         var = ElocAndPsiTermsENR_varALL(EPsiTENR)
         nSize = getSampleSizeAllNodes(sample)
         call ElocAndPsiTermsENR_resultALL(EPsiTENR,fi,ELi,fiEL,fij,fifj,fijEL,fifjEL,fiELj)
         if (logmode >= 2) then
            write(iul,'(a,f15.5,a,f12.5,a,f12.3,a,i10)') &
               ' calc Emean=',e0,' +/- ',sqrt(var/nSize),' var=',var,' size=',nSize
         end if
         nfev = nfev + 1
       end if

      if (task=='GH' .or. task=='START') then
         call ElocAndPsiTermsENR_resultALL(EPsiTENR,fi,ELi,fiEL,fij,fifj,fijEL,fifjEL,fiELj)
         if (logmode >= 2) then
            write(iul,*) ' calculating grad and Hessian without sample update'
         end if
         if (MASTER) then
            call internal_calcGradAndHessian(g,H)
            call internal_writeDataToLogFile
         end if
         ngev = ngev + 1
      end if

      ! Initialize the trust region bound.

      if (task=='START' .and. MASTER) then
         gnorm = dnrm2(np,g,1)        ! 2-norm from BLAS
         delta = 0.d0
      end if

      if (task=='START') then
         ! dnmtr creates new parameters. Equilibrate sample at next start of loop
         forceNewSample = .true.
      end if 

      if (MASTER) then
         call dnmtr(np,p,e0,g,H,np,frtol,fatol,fmin,task,delta,   &
               wa(5*np+1),scale,isave,dsave,wa(1),wa(np+1),wa(2*np+1), &
               wa(3*np+1),wa(4*np+1))
         if (logmode >= 2) write(iul,*) ' call to dnmtr resulted in task=',trim(task)
         call internal_writeCurrentVectorToLogFile()
      end if
      call myMPIBcastDouble(p,np)
      call wfparams_set(WFP,p)
      call myMPIBcastString(task,60)

      if (nfev > maxfev) then
         task = 'ERROR: NFEV > MAXFEV'
         exit
      end if
      if (task(1:4) == 'CONV') exit 
      if (task(1:4) == 'WARN' .or. task(1:5) == 'ERROR') exit
   end do

   if (logmode >= 2) write (iul,*) 'final task:',task
   if (task(1:4) /= 'CONV') converged = .false.

   call setCurrentResult(e0,0.d0,var)

   deallocate(p,g,g1,H)
   deallocate(fi,ELi,fiEL)
   deallocate(fifj,fifjEL,fiELj,fij,fijEL)
   deallocate(A,B,D)
   deallocate(eval,evec)
   deallocate(wa)

   call ElocAndPsiTermsENR_destroy(EPsiTENR)

   contains

      subroutine internal_readInput()
      !-----------------------------!
         gmode = 1
         call getinta(lines,nl,'gmode=',gmode,iflag)
         NRMode = 1; scale = .true.
         call getinta(lines,nl,'nrmode=',NRMode,iflag)
         if (NRMode==2) scale = .false.
         frtol = 1.d-14
         fatol = 2.d-3
         call getdbla(lines,nl,'tol=',fatol,iflag)
         delta = -1.d0
         call getdbla(lines,nl,'delta=',delta,iflag)
         subName = 'equilibrate'
         call getstra(lines,nl,'eq_call=',subName,iflag)
         maxfev = 5
         call getinta(lines,nl,'eq_iter=',maxfev,iflag)
      end subroutine internal_readInput


      subroutine internal_calcEPsiTerms()
      !---------------------------------!
         ! calculate sample average for E_loc and Psi terms required for linear method
         real*8 :: x(ne),y(ne),z(ne)
         type(RandomWalker), pointer :: rwp
         type(WFParamDerivTerms) :: wfpDT
         type(eConfigArray)  :: ec

         call eConfigArray_new(ec,ne,1)

         rwp => getFirst(sample)
         do
            call pos(rwp,x,y,z)
            call eConfigArray_set(ec,1,x,y,z)
            call eloc(0,ec,WFP%optType,WFP,wfpDT)
            !!!!call resetTo(rwp,x,y,z,WFP%optType)  ! recalculate and reset walker
            call ElocAndPsiTermsENR_add(EPsiTENR,wfpDT)
         if (.not.isNext(sample)) exit
            rwp => getNext(sample)
         enddo

      end subroutine internal_calcEPsiTerms



      subroutine internal_calcGradAndHessian(g,H)
         real*8 :: g(:)
         real*8 :: H(:,:)
         select case (gmode)
         case (1)
            g = 2*( fiEL - e0*fi )
         case (2)
            g = 2*( fiEL - e0*fi ) + ELi
         case (3)
            g = 2*( fiEL - e0*fi ) + 2*ELi
         end select
         g1 = 2*( fiEL - e0*fi )

         A = 2*( fijEL - fij*e0 - fifjEL + fifj*e0 )
         do j=1,np
            B(:,j) = -2*( fi(:)*g1(j) + fi(j)*g1(:) )
            D(:,j) = -fi(:)*ELi(j) - fi(j)*ELi(:)
         enddo
         B = B + 4*( fijEL - fij*e0 )
         D = D + fiELj + transpose(fiELj)

         H = A + B + D
      end subroutine internal_calcGradAndHessian

      subroutine internal_writeDataToLogFile()
         real*8 maxgrad,meangrad
         integer i

         maxgrad = maxval(abs(g))
         meangrad = sum(abs(g))/np
         if (logmode >= 2) then
            write(iul,'(2(A,G12.4))') ' gradient with abs mean = ',meangrad,' and abs max =',maxgrad
            write(iul,'(10g12.4)') (g(i),i=1,np)
         endif
         if (logmode >= 3) then
            write(iul, *)
            write(iul,'(A,G13.5,7G12.4)') 'ONE:',e0,2*(fiEL(1)-e0*fi(1)),2*(fiEL(1)-e0*fi(1))+ELi(1), &
              2*(fiEL(1)-e0*fi(1))+2*ELi(1),H(1,1),A(1,1),B(1,1),D(1,1)
            write(iul,*) 'Hessian:'
            do i=1,np
               write(iul,'(15G10.3)') H(i,:)
            enddo
            if (logmode >= 4) then
               write(iul,*) 'fi:'
               write(iul,'(10G10.3)') fi(:)
               write(iul,*) 'ELi:'
               write(iul,'(10G10.3)') ELi(:)
               write(iul,*) 'fiEL:'
               write(iul,'(10G10.3)') fiEL(:)
               write(iul,*) 'fij:'
               do i=1,np
                  write(iul,'(15G10.3)') fij(i,:)
               enddo
               write(iul,*) 'fijEL:'
               do i=1,np
                  write(iul,'(15G10.3)') fijEL(i,:)
               enddo
               write(iul,*) 'fifj:'
               do i=1,np
                  write(iul,'(15G10.3)') fifj(i,:)
               enddo
               write(iul,*) 'fifjEL:'
               do i=1,np
                  write(iul,'(15G10.3)') fifjEL(i,:)
               enddo
               write(iul,*) 'fiELj:'
               do i=1,np
                  write(iul,'(15G10.3)') fiELj(i,:)
               enddo
               write(iul,*) 'A:'
               do i=1,np
                  write(iul,'(15G10.3)') A(i,:)
               enddo
               write(iul,*) 'B:'
               do i=1,np
                  write(iul,'(15G10.3)') B(i,:)
               enddo
               write(iul,*) 'D:'
               do i=1,np
                  write(iul,'(15G10.3)') D(i,:)
               enddo
            endif
         endif
      end subroutine internal_writeDataToLogFile

      subroutine internal_writeCurrentVectorToLogFile()
         if (logmode >= 2) then
            write(iul,*) ' current parameter vector:'
            write(iul,'(10g12.4)') p 
         endif
      end subroutine internal_writeCurrentVectorToLogFile

   end subroutine eminTRNR_optimizeSample


end module OptimizeParamsETRNRModule
