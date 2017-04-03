program test
	use mrg
	use rannum
	use mt19937
	implicit none

	character(len=5) :: rng
	character(len=10) :: tmp
	integer :: i
	real*8 :: cpuStart, cpuEnd, wallStart, wallEnd, s
	integer :: seed = 1633837925
	integer,parameter :: ops = 10000000
	real*8 :: init_mrgran_c, mrg_ran_c
	rng = "mrg"

	if(iargc() .gt. 0) then
		call getarg(1, rng)
		if(iargc() .gt. 1) then
			call getarg(2, tmp)
			read(tmp, "(I10)") seed
			write(0, *) "Using seed", seed
		endif
	endif

	if(rng == "ran2") then
		write(0, *) "Testing ran2"
		call test_ran2()
	elseif(rng == "mrg") then
		write(0, *) "Testing MRG8 (integer)"
		call test_mrg_int()
	elseif(rng == "mrgr") then
		write(0, *) "Testing MRG8 (real)"
		call test_mrg_real()
	elseif(rng == "mt") then
		write(0, *) "Testing MT19937"
		call test_mt()
	elseif(rng == "test") then
		write(0, *) "MRG8 rands"
		write(*,*) init_mrgran(seed)
		!write(*,*) mrg_intran(), mrg_intran(), mrg_intran(), mrg_intran(), mrg_intran()
		write(*,*) mrg_ran(), mrg_ran(), mrg_ran(), mrg_ran(), mrg_ran(), mrg_ran(), &
		           mrg_ran(), mrg_ran(), mrg_ran(), mrg_ran(), mrg_ran(), mrg_ran()
	elseif(rng == "speed") then
		write(*,*) "Testing MRG8 speed"
		s = init_mrgran(seed)
		call startTime()
		do i = 1, ops
			s = s + mrg_ran()
		enddo
		call endTime()
		write(*,*) "Result:", s

		write(*,*) "Testing MRG8 C speed"
		s = init_mrgran_c(seed)
		call startTime()
		do i = 1, ops
			s = s + mrg_ran_c()
		enddo
		call endTime()
		write(*,*) "Result:", s

		write(*,*) "Testing ran2 speed"
		s = init_ran2(seed)
		call startTime()
		do i = 1, ops
			s = s + ran2_ran()
		enddo
		call endTime()
		write(*,*) "Result:", s

		write(*,*) "Testing MT19937 speed"
		s = init_mtran(seed)
		call startTime()
		do i = 1, ops
			s = s + mt_ran()
		enddo
		call endTime()
		write(*,*) "Result:", s
	endif


contains
	subroutine startTime()
		integer*8 :: cnt, rate
		call cpu_time(cpuStart)

		call system_clock(cnt,rate)
		wallStart = dble(cnt)/dble(rate)
	end subroutine

	subroutine endTime()
		integer*8 :: cnt, rate
		call cpu_time(cpuEnd)

		call system_clock(cnt,rate)
		wallEnd = dble(cnt)/dble(rate)
		write(*,*) "CPU time:", (cpuEnd - cpuStart)
		write(*,*) "Walltime:", (wallEnd - wallStart)
		write(*,*) "Ops/sec:", ops/(wallEnd - wallStart)

	end subroutine

	subroutine test_ran2()
		real*8 :: dummy
		dummy = init_ran2(seed)
		do
			call print_reals(ran2_ran(), ran2_ran())
		enddo
	end subroutine

	subroutine test_mrg_int()
		real*8 :: dummy
		integer(kind=8) :: a, b, c, d, rand
		integer, parameter :: m = 1073741823 ! 2**30 - 1
		dummy = init_mrgran(seed)
		do
			a = and(mrg_intran(), m)
			b = and(mrg_intran(), m)
			rand = lshift(a, 30) + b
			call print_int(rand, 7)
			c = and(mrg_intran(), m)
			d = and(mrg_intran(), m)
			rand = lshift(c, 30) + d
			call print_int(rand, 7)
			rand = lshift(rshift(c, 26), 4) + rshift(a, 26)
			call print_int(rand, 1)
		enddo
	end subroutine

	subroutine test_mrg_real()
		real*8 :: dummy
		dummy = init_mrgran(seed)
		do
			call print_reals(mrg_ran(), mrg_ran())
		enddo
	end subroutine

	subroutine test_mt()
		real*8 :: dummy
		dummy = init_mtran(seed)
		do
			call print_reals(mt_ran(), mt_ran())
		enddo
	end subroutine

	subroutine print_int(foo, bytes)
		integer(kind=8), intent(in) :: foo
		integer, intent(in) :: bytes
		integer :: i
		do i = 0, bytes - 1
			write(*,"(A)",advance="no") char(and(rshift(foo, i*8), 255))
		enddo

	end subroutine

	subroutine print_reals(m, n)
		real*8, intent(in) :: m, n
		integer(kind=8) :: a, b, rand
		integer(kind=8), parameter :: k = 4503599627370495_8 ! 2**52 - 1
		a = transfer(m, a)
		a = and(a, k)
		b = transfer(n, b)
		b = and(b, k)
		call print_int(a, 6)
		call print_int(b, 6)
		rand = lshift(rshift(a, 48), 4) + rshift(b, 48)
		call print_int(rand, 1)
	end subroutine

end program

