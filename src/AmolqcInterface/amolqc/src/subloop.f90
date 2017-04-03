module subloopModule

   use global
   use utilsmodule, only: assert, getNextBlock, getinta, finda
   use Utils, only: getToken,getTimeString,readFileParallel,expandMacro
   use InitModule, only: initAmolqc,readInputFile,initGen,finalizeAmolqc
   use RWSampleModule, only: RWSample, writeSampleCommand, recalculateSample
   use RandomWalkerModule, only: setNumberOfElectrons,setNumberOfCenters,setEpart
   use QmcSample, only: sample, initInitialWalker
   use qmc, only: qmc_run, qmc_init, initWalkerStat, initTrajectory
   use properties
   use wfmodule
   use eloctest
   use optimderivstest
   use jastrow, only: jasChangeType
   use wfParamsModule, only: wfparams_change
   use ce_new_m
   use coc_control, only: center_of_charge_run
   use referenceModule, only: references_check, references_check1, references_optimize
   use psimax, only: psimax_init
   use maxanalysisModule, only: maxana_init
   use maxbasinsModule, only: maxbas_init
   implicit none

contains


subroutine subloop_callSubroutine(lines,nl,smpl)
   character(len=120), intent(in) :: lines(:)
   integer, intent(in)            :: nl
   type(RWSample), intent(inout)  :: smpl
   character(len=120) :: subName
   integer iflag

   call getstra(lines,nl,'name=',subName,iflag)
   call subloop(subName,smpl)
end subroutine subloop_callSubroutine



subroutine subloop(subname,smpl)
!------------------------------!
   ! subloop enters the event script loop under the name "subname" which is either
   ! a subroutine (i.e. globally stored series of blocks) or a macro command
   ! (i.e. blocks stored in a file)
   character(len=*), intent(in)  :: subname
   type(RWSample), intent(inout)  :: smpl

   integer, parameter :: MAXLINES=100

   character(len=120)          :: token
   character(len=MAXLEN)       :: inLines(MAXLINES)=''
   character(len=MAXLEN)       :: blockLines(MAXLINES)=''
   character(len=MAXLEN)       :: macroLines(MAXLINES)=''
   character(len=15)           :: si = ' =======>      '
   character(len=14)           :: sf = '      <======='
   character(len=120)          :: macropath,macrofile
   integer                     :: idx,nbl,nil,i,iuf,io,nnew,iflag,fileExistsInt
   integer                     :: loopIdx,loopIter,currentLoopIter,subIdx
   integer                     :: mLines
   real*8                      :: start,myMPIWallTime,startCPU,endCPU
   real*8                      :: tstart,tstartCPU,tendCPU,sendbuf(1),recvbuf(1),t
   logical                     :: wfRead,found,exitLoop,fileExists,wout

   call getAmolqcPath(macropath)
   call assert(len(trim(macropath)) < 116,"amolqc: amolqc path exceeds definition")
   macropath = trim(macropath)//"/cmds"

   do subIdx=1,MAXSUBS 
      if (subNames(subIdx)==subname) then
         if (logmode >= 2) write(iul,'(/2a)') ' ================>     calling subroutine ',trim(subname)
         exit
      end if
   end do
   if (subIdx>MAXSUBS) then ! not found in subroutine list, look for .cmd file
      macrofile = trim(macropath)//'/'//trim(subname)//'.cmd'
      fileExistsInt = 1
      if (MASTER) then
         inquire(file=macrofile,exist=fileExists)
         if (.not.fileExists) then
            fileExistsInt = 0
            write(iul,'(//2a/)') "ERROR in .in file: unknown subroutine ",trim(subname)
         else
            if (logmode >= 2) write(iul,'(/2a)') ' * * *  calling macro ',trim(subname)
         end if
      end if
      call myMPIBcastInteger(fileExistsInt,1)
      if (fileExistsInt==0) call abortp("unkown command in infile")
      call readFileParallel(mytid,macrofile,macroLines,mLines)
      nil = 0
      idx = 1
      call expandMacro(inLines,nil,macroLines,mLines,idx)
      if (logmode >= 2) then
         write(iul,'(//3a)') '  ---  Expanding macro cmd ',trim(subname),' to:'
         do i=1,mLines
            write(iul,'(A)') trim(inlines(idx+i-1))
         enddo
      endif
   else ! found in subroutine list
      do i=1,subLen(subIdx)
         inLines(i) = subLines(i,subIdx)
      end do
      nil = subLen(subIdx)
   end if


   idx = 1 
   wout = .false.   ! write output
   if (MASTER .and. logmode > 1) wout = .true. 
   do
      call getNextBlock(inLines,nil,idx,'$',')','!',MAXLINES,blockLines,nbl)

      if (nbl == 0) exit

      token = getToken(blockLines(1),'$','(')

      start = myMPIWallTime()
      call cpu_time(startCPU)

      if (token=='change_jastrow') then
         if (wout) write(iul,'(/a/)') si//'$change_jastrow - changing Jastrow terms'//sf
         call jasChangeType(blocklines,nbl)
         call recalculateSample(smpl)
      else if (token=='change_parameters') then
         if (wout) write(iul,'(/a/)') si//'$change_parameters'//sf
         call wfparams_change(blocklines,nbl)
      else if (token=='eloctest') then
         if (wout) write(iul,'(/a/)') si//'$eloctest - testing local energies'//sf
         call runEloctest(blockLines,nbl,smpl)
      else if (token=='init_max_analysis') then
         if (wout) write(iul,'(/a/)') si//'$init_max_analysis - initializing maximum analysis'//sf
         call maxana_init(blocklines,nbl)
      else if (token=='init_basin_analysis') then
         if (wout) write(iul,'(/a/)') si//'$init_basin_analysis - initializing basin analysis'//sf
         call maxbas_init(blocklines,nbl)
      else if (token=='init_max_search') then
         if (wout) write(iul,'(/a/)') si//'$init_max_search - initializing maxima search'//sf
         call psimax_init(blocklines,nbl)
      else if (token=='init_walker') then
         if (wout) write(iul,'(/a/)') si//'$init_walker - setting an initial walker'//sf
         call initInitialWalker(blockLines,nbl)
      else if (token=='iterate_coc') then
         if (wout) write(iul,'(/a/)') si//'$iterate_coc - center of charge run'//sf
         call center_of_charge_run(blocklines,nbl,smpl)
      else if (token=='optimize_refs') then
         if (wout) write(iul,'(/a/)') si//'$optimize_refs - optimizing references'//sf
         call references_optimize(blockLines,nbl)
      else if (token=='optimtest') then
         if (wout) write(iul,'(/a/)') si//'$optimtest - testing the optimizers'//sf
         call optimizeTest(blockLines,nbl,smpl)
      !!!else if (token=='params_numderivs') then
      !!!   if (wout) write(iul,'(/a/)') si//'$params_numderivs - ???'//sf
      !!!   call wfparams_derivatives(blocklines,nbl,smpl)
      else if (token=='print_results') then
         if (wout) write(iul,'(/a/)') si//'$print_results - printing stored results'//sf
         call global_printSavedResults(blockLines,nbl)
      else if (token=='props') then
         if (wout) write(iul,'(/a/)') si//'$props - calculting properties'//sf
         call propInit(blocklines,nbl)
      else if (token=='qmc') then
         if (wout) write(iul,'(/a/)') si//'$qmc - running a qmc calculation'//sf
         call qmc_init(blockLines,nbl,smpl)
         call qmc_run(smpl)
      else if (token=='sample') then
         if (wout) write(iul,'(/a/)') si//'$sample - creating or modifying the walker sample'//sf
         call sample(blockLines,nbl,smpl)
      else if (token=='save_result') then
         if (wout) write(iul,'(/a/)') si//'$save_results - storing current results'//sf
         call global_saveResult(blockLines,nbl)
      else if (token=='sed') then
         if (wout) write(iul,'(/a/)') si//'$sed - running a sed calculation'//sf
         call qmc_init(blockLines,nbl,smpl)
         call qmc_run(smpl)
      else if (token=='stop_if') then
         exitLoop = global_exitIf(blocklines,nbl)
         if (exitLoop) then
            if (MASTER .and. logmode>=2) then
               write(iul,'(/a/)') '   ---->   stop on condition   <--------'
            end if
            exit
         end if
      else if (token=='test_balance') then
         call testBalance(smpl)
      else if (token=='write_sample') then
         if (wout) write(iul,'(/a/)') si//'$write_sample - writing the walker sample'//sf
         call writeSampleCommand(blockLines,nbl,smpl)
      ! loop commands
      else if (token=='begin_loop') then
         call assert(nbl==1,"begin_loop calls may consist of only one line")
         loopIdx = idx
         call getinta(blocklines,nbl,'count=',loopIter,iflag)
         currentLoopIter = 1
         call setCurrentLoopIdx(currentLoopIter)
         call assert(iflag==0,'begin_loop requires count value')
      else if (token=='end_loop') then
         if (loopIter > 1) then
            idx = loopIdx
            loopIter = loopIter - 1
            currentLoopIter = currentLoopIter + 1
            call setCurrentLoopIdx(currentLoopIter)
         else
            call setCurrentLoopIdx(1)
         endif
      else if (token=='exit_if') then
         exitLoop = global_exitIf(blocklines,nbl)
         if (exitLoop) then
            if (finda(blocklines,nbl,'stop')) then
               if (MASTER .and. logmode>=2) then
                  write(iul,'(/a/)') '   ---->   exit from loop: stop   <--------'
               end if
               exit
            else 
               found = .false.
               if (MASTER .and. logmode>=2) then
                  write(iul,'(/a/)') '   ---->   exit from loop: continue   <--------'
               end if
               do i=idx,nil
                  if (index(inLines(i),'$end_loop') > 0) then
                     idx = i+1
                     found = .true.
                     exit
                  end if
               end do
               if (.not.found) call abortp("$exit_if only allowed between $begin_loop and $end_loop")
            end if
         end if
      ! macro expansion
      else
         macrofile = trim(macropath)//'/'//trim(token)//'.cmd'
         fileExistsInt = 1
         if (MASTER) then
            inquire(file=macrofile,exist=fileExists)
            if (.not.fileExists) then
               fileExistsInt = 0
               write(iul,'(//2a/)') "ERROR in .in file: unknown command ",trim(token)
            end if
         end if
         call myMPIBcastInteger(fileExistsInt,1)
         if (fileExistsInt==0) call abortp("unkown command in infile")
         call readFileParallel(mytid,macrofile,macroLines,mLines)
         call assert(nbl==1,' currently only one-line macros allowed')
         call expandMacro(inLines,nil,macroLines,mLines,idx)
         if (logmode >= 2) then
            write(iul,'(//3a)') '  ---  Expanding macro cmd ',trim(token),' to:'
            do i=1,mLines
               write(iul,'(A)') trim(inlines(idx+i-1))
            enddo
         endif
      endif

      if (logmode >= 2) then
         if (token=='qmc' .or. token=='sample') then
            write(iul,'(/3A,F18.2,A)') ' wall clock time for   ',trim(token),' : ', &
              myMPIWallTime()-start,' s'
            call cpu_time(endCPU)
            write(iul,'(3A,F18.2,A//)') ' cpu time (master) for ',trim(token),' : ', &
            endCPU-startCPU,' s'
         end if
         call flush(iul)
      end if

   end do

   if (logmode >= 2) then
      if (subIdx>MAXSUBS) then
         write(iul,'(/2a/)') ' =============>    end macro ',trim(subname)
      else
         write(iul,'(/2a/)') ' =============>    end subroutine ',trim(subname)
      end if
   end if      


end subroutine subloop

end module subloopModule
