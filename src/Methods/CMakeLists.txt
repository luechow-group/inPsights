cmake_minimum_required(VERSION 3.5)

project(Methods C CXX Fortran)


set ( PROJECT_FOLDER_NAME ${PROJECT_NAME} )
set ( PROJECT_PATH ${PROJECT_NAME})
set ( OUTPUT_NAME "${PROJECT_FOLDER_NAME}" )
message ( "-- ++ CREATING PROJECT OR MAKEFILE ---> ${OUTPUT_NAME}" )

set(CMAKE_CXX_STANDARD 11)

include_directories(
        include
        source
        )

set(HEADER_FILES
        include/QuaternionFit.h
        include/StringMethodProblem.h
        #include/StringMethod.h
        )

set(SOURCE_FILES
        source/QuaternionFit.cpp
        source/StringMethodProblem.cpp
        )

add_library(${PROJECT_NAME} STATIC ${SOURCE_FILES} ${HEADER_FILES})
target_include_directories(${PROJECT_NAME} PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        )
target_link_libraries(${PROJECT_NAME}
        BSplines
        AmolqcInterface # already includes eigen and cppoptlib
        )

add_executable(${PROJECT_NAME}_main main.cpp )
#target_link_libraries(${PROJECT_NAME}_main ${PROJECT_NAME})
target_link_libraries(${PROJECT_NAME}_main
        ${PROJECT_NAME}
        ${LAPACK_LINKER_FLAGS}
        ${BLAS_LINKER_FLAGS}
        )

add_subdirectory(tests)
