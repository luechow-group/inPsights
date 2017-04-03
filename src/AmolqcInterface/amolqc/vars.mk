# compilers
ifeq "$(origin FC)" "default"
  ifeq ($(COMPILER),intel)
    FC = ifort
  else ifeq ($(COMPILER),gnu)
    FC = gfortran
  else
    FC = mp95
  endif
endif

ifeq ($(MPIFC),)
  ifeq ($(COMPILER),intel)
    #MPIFC = mpiifort
    MPIFC = mpif90
  else ifeq ($(COMPILER),gnu)
    MPIFC = mpif90
  else
    MPIFC = mpf95
  endif
endif

ifeq ($(PARALLEL),yes)
  FC = $(MPIFC)
endif

# compile cache
ifdef CCACHE
  FC := $(CCACHE) $(FC)
endif

# compiler flags
ifeq ($(COMPILER),intel)
  FPP := $(FC) -fpp
  AR = xiar
  RANLIB = xiar -s
  ifeq ($(OPTLEVEL),no)
    FFLAGS += -132 -O0
    FNOFLAGS += -132 -O0
    F90FLAGS += -O0
  else ifeq ($(OPTLEVEL),debug)
    FFLAGS += -132 -g -check bounds -traceback
    FNOFLAGS += -132 -g -check bounds -traceback
    F90FLAGS += -g -check bounds -traceback
  else ifeq ($(OPTLEVEL),high)
    FFLAGS += -132 -O3
    FNOFLAGS += -132 -O0
    F90FLAGS += -O3
  else ifeq ($(OPTLEVEL),fast)
    ifneq ($(FLAGS_INTEL_FAST_NO_FPOPT),)
      FFLAGS += -132 $(FLAGS_INTEL_FAST_NO_FPOPT)
      FNOFLAGS += -132 -O0
      F90FLAGS += $(FLAGS_INTEL_FAST_NO_FPOPT)
    else
      FFLAGS += -132 -fast -fp-model source
      FNOFLAGS += -132 -O0
      F90FLAGS += -fast -fp-model source
    endif
  else
    FFLAGS += -132
    FNOFLAGS += -132 -O0
    F90FLAGS +=
  endif
  FFLAGS += -module mod -warn all -warn nounused -warn noerrors -warn nointerfaces -WB -diag-disable 8290,8291
  FNOFLAGS += -module mod -warn all -warn nounused -warn noerrors -warn nointerfaces -WB -diag-disable 8290,8291
  F90FLAGS += -module mod -warn all -warn nounused -warn noerrors -warn nointerfaces -WB -diag-disable 8290,8291
else ifeq ($(COMPILER),gnu)
  FPP := $(FC) -cpp
  AR = ar
  RANLIB = ranlib
  ifeq ($(OPTLEVEL),no)
    FFLAGS += -O0
    FNOFLAGS += -O0
    F90FLAGS += -O0
  else ifeq ($(OPTLEVEL),debug)
    FFLAGS += -g -fbounds-check
    FNOFLAGS += -g -fbounds-check
    F90FLAGS += -g -fbounds-check
  else ifeq ($(OPTLEVEL),std)
    FFLAGS += -O2
    FNOFLAGS += -O0
    F90FLAGS += -O2
  else ifeq ($(OPTLEVEL),high)
    FFLAGS += -O3
    FNOFLAGS += -O0
    F90FLAGS += -O3
  endif
  FFLAGS += -ffixed-form -ffixed-line-length-132 -Jmod -Wall -Wtabs -Wno-unused -Wno-maybe-uninitialized -Wno-unused-dummy-argument -Wno-conversion
  FNOFLAGS += -ffixed-form -ffixed-line-length-132 -Jmod -Wall -Wtabs -Wno-unused -Wno-maybe-uninitialized -Wno-unused-dummy-argument -Wno-conversion
  F90FLAGS += -Jmod -Wall -Wtabs -Wno-unused -Wno-maybe-uninitialized -Wno-unused-dummy-argument -Wno-conversion
else
  AR = ar
  RANLIB = ranlib
  FPP := $(FC) -cpp
endif


# link commands
LDFLAGS += -fopenmp
ifeq ($(OPTLEVEL),debug)
  LDFLAGS += -g
endif

ifneq ($(RNG),)
FFLAGS += -D$(RNG)
FNOFLAGS += -D$(RNG)
F90FLAGS += -D$(RNG)
endif

ifeq ($(MLIBS),)
  ifeq ($(LAPACK),mkl)
    MLIBS = -L$(MKLPATH) $(MKLLIBS)
  else ifeq ($(LAPACK),compiled)
    MLIBS = -L$(MATHPATH) -llapack -lblas
  endif
endif
