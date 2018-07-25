if(EXISTS "$ENV{MKLROOT}")
    message("using MKL")
    set(MKL_USE_STATIC_LIBS ON)

    find_library(MKL_INTERFACE_LIBRARY
            NAMES mkl_intel_lp64
            PATHS $ENV{MKLROOT}/lib
            $ENV{MKLROOT}/lib/intel64
            $ENV{INTEL}/mkl/lib/intel64
            NO_DEFAULT_PATH
            )

    find_library(MKL_SEQUENTIAL_LAYER_LIBRARY
            NAMES mkl_sequential
            PATHS $ENV{MKLROOT}/lib
            $ENV{MKLROOT}/lib/intel64
            $ENV{INTEL}/mkl/lib/intel64
            NO_DEFAULT_PATH)

    find_library(MKL_CORE_LIBRARY
            NAMES mkl_core
            PATHS $ENV{MKLROOT}/lib
            $ENV{MKLROOT}/lib/intel64
            $ENV{INTEL}/mkl/lib/intel64
            NO_DEFAULT_PATH)
endif()