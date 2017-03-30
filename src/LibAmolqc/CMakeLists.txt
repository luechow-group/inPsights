cmake_minimum_required(VERSION 3.5)


#============= PROJECT LibAmolqc ===================
project(LibAmolqc)

set ( PROJECT_FOLDER_NAME ${PROJECT_NAME})

set ( PROJECT_PATH ${PROJECT_NAME})
set ( OUTPUT_NAME "${PROJECT_FOLDER_NAME}" )
message ( "-- ++ CREATING PROJECT OR MAKEFILE ---> ${OUTPUT_NAME}" )

find_package(LAPACK REQUIRED)
#find_package(BLAS REQUIRED)
message(STATUS ${BLAS_FOUND})

include_directories(
        include
        source
        )

set(HEADER_FILES
        include/ElectronicWaveFunction.h
        )

set(SOURCE_FILES
        source/ElectronicWaveFunction.cpp
        )

get_filename_component(PATH_TO_LIBAMOLQC amolqc/src/libamolqc.a REALPATH)
get_filename_component(PATH_TO_LIBUTILS amolqc/utils/libutils.a REALPATH)

#set(CMAKE_CXX_FLAGS "-Wall -fopenmp")

add_library(${PROJECT_NAME}_interface STATIC ${SOURCE_FILES} ${HEADER_FILES})
target_link_libraries(${PROJECT_NAME}_interface
        eigen
        )

add_library(libamolqc STATIC IMPORTED GLOBAL)
set_property(TARGET libamolqc PROPERTY IMPORTED_LOCATION ${PATH_TO_LIBAMOLQC})
add_library(libutils STATIC IMPORTED GLOBAL)
set_property(TARGET libutils PROPERTY IMPORTED_LOCATION ${PATH_TO_LIBUTILS})


add_library(${PROJECT_NAME} STATIC dummy.cpp)
target_link_libraries(${PROJECT_NAME}
        ${PROJECT_NAME}_interface
        Externals::cppoptlib
        libamolqc
        libutils
        "-llapack -lblas -lgfortran"
        )
target_include_directories(${PROJECT_NAME} PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        )

add_executable("${PROJECT_NAME}_main" main.cpp)
target_link_libraries("${PROJECT_NAME}_main"
        ${PROJECT_NAME}
        )


