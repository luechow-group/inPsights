
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


# General settings
if(CMAKE_BUILD_TYPE STREQUAL "Release")
    add_definitions(-O3)
endif()

# Compiler specific optimization settings
if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
    # using Clang
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    # using GCC
    add_definitions(-ftracer -floop-optimize -funroll-loops -mtune=native -mmmx -msse2 -mfpmath=sse -g)
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
    # using Intel C++
    add_definitions(-ipo)
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
    # using Visual Studio C++
    message(" ## WARNING ##: Amolqcpp was not tested with the Microsoft Visual Studio compiler.
    Consider switching to the Intel or GNU compiler collection ")
endif()
