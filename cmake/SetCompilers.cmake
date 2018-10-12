
# Determine compiler from environment variables
if(EXISTS $ENV{CC})
    set(CMAKE_C_COMPILER $ENV{CC})
    message("Using specified C Compiler         ${CMAKE_C_COMPILER}")
else()
    message("Using default C Compiler           ${CMAKE_C_COMPILER}")
endif()

if(EXISTS $ENV{CXX})
    set(CMAKE_CXX_COMPILER $ENV{CXX})
    message("Using specified C++ Compiler       ${CMAKE_CXX_COMPILER}")
else()
    message("Using default CXX Compiler         ${CMAKE_CXX_COMPILER}")
endif()

if(EXISTS $ENV{FC})
    set(CMAKE_Fortran_COMPILER $ENV{FC})
    message("Using specified Fortran Compiler   ${CMAKE_Fortran_COMPILER}")
else()
    message("Using default Fortran Compiler     ${CMAKE_Fortran_COMPILER}")
endif()

# Compiler specific optimization settings
if (${CMAKE_CXX_COMPILER_ID} MATCHES "Clang")

    if(${CMAKE_BUILD_TYPE} MATCHES "Release")
        add_definitions(-O3)
    else()
        add_definitions(-Wall -g)
    endif()

elseif ("${CMAKE_CXX_COMPILER_ID}" MATCHES "GNU")

    if(${CMAKE_BUILD_TYPE} MATCHES "Release")
        add_definitions(-O3 -ftracer -floop-optimize -funroll-loops -mtune=native -mmmx -msse2 -mfpmath=sse)
    else()
        add_definitions(-Wall -g)
    endif()

elseif (${CMAKE_CXX_COMPILER_ID} MATCHES "PGI")

    if(${CMAKE_BUILD_TYPE} MATCHES "Release")
        add_definitions(-O3)
    else()
        add_definitions()
    endif()

elseif (${CMAKE_CXX_COMPILER_ID} MATCHES "Intel")

    if(${CMAKE_BUILD_TYPE} MATCHES "Release")
        add_definitions(-O3-ipo -inline-forceinline)
    else()
        add_definitions(-Wall -g)
    endif()

elseif (${CMAKE_CXX_COMPILER_ID} MATCHES "MSVC")

    # using Visual Studio C++
    message(" ## WARNING ##: Amolqcpp was not tested with the Microsoft Visual Studio compiler.
    Consider switching to the Intel or GNU compiler collection ")

    if(${CMAKE_BUILD_TYPE} MATCHES "Release")
        add_definitions(-O3)
    else()
        add_definitions()
    endif()
endif()
