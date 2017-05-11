cmake_minimum_required(VERSION 3.5)

project(AmolqcInterface C CXX Fortran)


set ( PROJECT_FOLDER_NAME ${PROJECT_NAME})
set ( PROJECT_PATH ${PROJECT_NAME})
set ( OUTPUT_NAME "${PROJECT_FOLDER_NAME}" )
message ( "-- ++ CREATING PROJECT OR MAKEFILE ---> ${OUTPUT_NAME}" )

set(CMAKE_CXX_STANDARD 11)

find_package(Eigen3 REQUIRED)
include_directories( ${EIGEN3_INCLUDE_DIR} )

# Require LAPACK Library which requires BLAS library implicitly
find_package(LAPACK REQUIRED)

include_directories(
        include
        source
        )

set(HEADER_FILES
        include/ElectronicWaveFunction.h
        include/ElectronicWaveFunctionProblem.h
        )

set(SOURCE_FILES
        source/ElectronicWaveFunction.cpp
        source/ElectronicWaveFunctionProblem.cpp
        )

get_filename_component(PATH_TO_LIBAMOLQC amolqc/src/libamolqc.a REALPATH)
get_filename_component(PATH_TO_LIBUTILS amolqc/utils/libutils.a REALPATH)

add_library(${PROJECT_NAME} STATIC ${SOURCE_FILES} ${HEADER_FILES} dummy.f90)
target_link_libraries(${PROJECT_NAME}
        #eigen
        cppoptlib
        ${PATH_TO_LIBAMOLQC}
        ${PATH_TO_LIBUTILS}
        ${LAPACK_LIBRARIES}
        ${BLAS_LIBRARIES}
        ${LAPACK_LINKER_FLAGS}
        ${BLAS_LINKER_FLAGS}
        )
target_include_directories(${PROJECT_NAME} PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        )

add_executable(${PROJECT_NAME}_main dummy.f90 main.cpp)
target_link_libraries(${PROJECT_NAME}_main
        ${PROJECT_NAME}
        ${LAPACK_LINKER_FLAGS}
        ${BLAS_LINKER_FLAGS}
        )
