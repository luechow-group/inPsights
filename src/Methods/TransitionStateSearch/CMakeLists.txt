cmake_minimum_required(VERSION 3.5)

project(TransitionStateSearch C CXX Fortran)

set ( PROJECT_FOLDER_NAME ${PROJECT_NAME} )
set ( PROJECT_PATH ${PROJECT_NAME})
set ( OUTPUT_NAME "${PROJECT_FOLDER_NAME}" )
message ( "-- ++ CREATING PROJECT OR MAKEFILE ---> ${OUTPUT_NAME}" )

set(HEADER_FILES
        include/StringOptimizationProblem.h
        include/StringMethod.h
        include/ChainOfStates.h
        include/LocalNewtonSearch.h
        )

set(SOURCE_FILES
        source/StringOptimizationProblem.cpp
        source/StringMethod.cpp
        source/ChainOfStates.cpp
        )

add_library(${PROJECT_NAME} STATIC ${SOURCE_FILES} ${HEADER_FILES})

target_include_directories(${PROJECT_NAME}
        PUBLIC include
        PRIVATE source
        )

target_link_libraries(${PROJECT_NAME}
        BSplines
        AmolqcInterface
        cppoptlib
        )

add_executable(${PROJECT_NAME}_main main.cpp)
target_link_libraries(${PROJECT_NAME}_main
        ${PROJECT_NAME}
        ${LAPACK_LINKER_FLAGS}
        ${BLAS_LINKER_FLAGS}
        )

#add_subdirectory(tests) # currently empty


