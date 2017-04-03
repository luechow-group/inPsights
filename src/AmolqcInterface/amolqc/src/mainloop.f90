module mainloopModule

contains 

subroutine mainloop

   use global
   use utilsmodule, only: assert, getNextBlock, getinta, finda
   use Utils, only: getToken,getTimeString,readFileParallel,expandMacro
   use InitModule, only: initAmolqc,readInputFile,initGen,finalizeAmolqc
   use RWSampleModule, only: RWSample, writeSampleCommand, recalculateSample
   use RandomWalkerModule, only: setNumberOfElectrons,setNumberOfCenters,setEpart
   use QmcSample, only: sample, initInitialWalker
   use qmc, only: qmc_run, qmc_init, initWalkerStat, initTrajectory
   use OptimizeParamsModule
   use properties
   use elocal, only: wf_init, wf_ecp_init
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
   use subloopModule
   implicit none

   integer, parameter :: MAXLINES=1000

   character(len=120)          :: token
   character(len=MAXLEN)       :: inLines(MAXLINES)=''
   character(len=MAXLEN)       :: blockLines(MAXLINES)=''
   character(len=MAXLEN)       :: macroLines(MAXLINES)=''
   character(len=120)          :: macropath,macrofile,subName
   character(len=15)           :: si = ' =======>      '
   character(len=14)           :: sf = '      <======='
   type(RWSample)              :: smpl
   integer                     :: idx,nbl,nil,i,iuf,io,nnew,iflag,fileExistsInt
   integer                     :: loopIdx,loopIter,currentLoopIter,subIdx,subLine
   integer                     :: mLines
   real*8                      :: start,myMPIWallTime,startCPU,endCPU
   real*8                      :: tstart,tstartCPU,tendCPU,sendbuf(1),recvbuf(1),t
   logical                     :: wfRead,found,exitLoop,fileExists,subMode,converged,wout

   call readInputFile(inLines,nil)
   call getAmolqcPath(macropath)
   call assert(len(trim(macropath)) < 116,"amolqc: amolqc path exceeds definition")
   macropath = trim(macropath)//"/cmds"

   tstart = myMPIWallTime()
   call cpu_time(tstartCPU)

   idx = 1
   wfRead = .false.
   converged = .true.
   subIdx = 0
   wout = .false.   ! write output
   do

      if (.not.converged) then
         if (MASTER) write(iul,'(/a/)') '   ---->   not converged: terminating run   <--------'
         exit
      end if

      call getNextBlock(inLines,nil,idx,'$',')','!',MAXLINES,blockLines,nbl)

      if (nbl == 0) exit

      token = getToken(blockLines(1),'$','(')

      start = myMPIWallTime()
      call cpu_time(startCPU)

      if (token=='gen') then
         if (wfRead) call abortp('$gen block must precede $wf in .in file')
         call initGen(blockLines,nbl)
         if (MASTER .and. logmode>1) wout = .true.
      !!else if (token=='sparse') then
      !!   if (wfRead) call abortp('$sparse block must precede $wf in .in file')
      !!   call initSparse(blockLines,nbl)
      else if (token=='wf') then
         if (wout) write(iul,'(/a/)') si//'$wf - wave function'//sf
         call wf_init(blockLines,nbl)
         if (.not.wfRead) then
            call setNumberOfElectrons(ne)
            call setNumberOfCenters(ncenter)
            call setEpart(do_epart)
            wfRead = .true.
         endif
      else if (token=='ecp') then
         if (.not.wfRead) call abortp('$ecp block must follow $wf in .in file')
         if (wout) write(iul,'(/a/)') si//'$ecp - effective core potential settings'//sf
         call wf_ecp_init(blockLines,nbl)
      else if (token=='analyze_refs') then
         if (wout) write(iul,'(/a/)') si//'$analyze_refs - analyze references'//sf
         if (finda(blockLines,nbl,'max_mode=')) then
           call references_check1(blockLines,nbl)
         else
           call references_check(blockLines,nbl)
         end if
      else if (token=='call_subroutine') then
         call subloop_callSubroutine(blocklines,nbl,smpl)
      else if (token=='change_jastrow') then
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
      else if (token=='optimize_parameters') then
         if (wout) write(iul,'(/a/)') si//'$optimize_parameters - optimizing wave function parameters'//sf
         call optimizeParameters(blockLines,nbl,smpl,converged)
      else if (token=='optimize_refs') then
         if (wout) write(iul,'(/a/)') si//'$optimize_refs - optimizing references'//sf
         call references_optimize(blockLines,nbl)
      else if (token=='wf_param_deriv_test') then
         if (wout) write(iul,'(/a/)') si//'$wf_param_deriv_test - testing the parameter derivatives'//sf
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
            if (MASTER) then
               write(iul,'(/a/)') '   ---->   stop on condition   <--------'
            end if
            exit
         end if
      else if (token=='test_balance') then
         call testBalance(smpl)
      else if (token=='trajectory') then
         if (wout) write(iul,'(/a/)') si//'$trajectory - initializing a trajectory run'//sf
         call initTrajectory(blockLines,nbl)
      else if (token=='walker_stat') then
         if (wout) write(iul,'(/a/)') si//'$walker_stat - initializing walker statistics'//sf
         call initWalkerStat(blockLines,nbl)
      else if (token=='write_sample') then
         if (wout) write(iul,'(/a/)') si//'$write_sample - writing the walker sample'//sf
         call writeSampleCommand(blockLines,nbl,smpl)

      ! subroutine definition   
      else if (token == 'begin_subroutine') then
         call getstra(blockLines,nbl,'name=',subName,iflag)
         if (iflag /= 0) call abortp('begin_subroutine requires name argument')
         subIdx = subIdx + 1
         subNames(subIdx) = subName
         subLine = 0
         do 
            if (subLine > MAXSUBLINES) call abortp('subroutine has too many lines')
            token = getToken(inLines(idx+subLine),'$','(')
            if (token == 'end_subroutine') then
               idx = idx + subLine + 1
               subLen(subIdx) = subLine
               exit
            end if
            subLine = subLine + 1
            subLines(subLine,subIdx) = inLines(idx+subLine-1)
         end do
         if (wout) write(iul,'(/a)') si//' $subroutine: storing subroutine '//trim(subName)//sf
      
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
               if (MASTER) then
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
         if (token=='qmc' .or. token=='sample' .or. token=='optimize_parameters') then
            write(iul,'(/3A,F18.2,A)') ' wall clock time for   ',trim(token),' : ', &
              myMPIWallTime()-start,' s'
            call cpu_time(endCPU)
            write(iul,'(3A,F18.2,A//)') ' cpu time (master) for ',trim(token),' : ', &
            endCPU-startCPU,' s'
         end if
         call flush(iul)
      end if

   end do

   call cpu_time(tendCPU)
   sendbuf(1) = tendCPU - tstartCPU
   call myMPIReduceSumDouble(sendbuf,recvbuf,1)
   t = myMPIWallTime()-tstart

   if (MASTER .and. logmode >= 2) then
      write(iul,'(//a,a17)') ' wall clock time for run         : ',getTimeString(t)
      write(iul,'(a,f17.4)') ' total cpu time for run (core-h) : ',recvbuf(1)/3600.d0
      write(iul,'(a,f17.4)') ' cpu time per mpi process (h)    : ',recvbuf(1)/3600.d0/nproc
   end if

end subroutine mainloop

end module mainloopModule

