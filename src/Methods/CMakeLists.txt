cmake_minimum_required(VERSION 3.5)


#============= PROJECT Methods ===================
project(Methods)

set ( PROJECT_FOLDER_NAME ${PROJECT_NAME} )
set ( PROJECT_PATH ${PROJECT_NAME})
set ( OUTPUT_NAME "${PROJECT_FOLDER_NAME}" )
message ( "-- ++ CREATING PROJECT OR MAKEFILE ---> ${OUTPUT_NAME}" )

include_directories(
        include
        source
        )

set(HEADER_FILES
        include/QuaternionFit.h
        #include/Optimization/IOptimizer.h
        #include/Optimization/OptimizerFactory.h
        #include/Optimization/NewtonTypeOptimizer.h
        #include/Optimization/LBFGSOptimizer.h
        #include/Optimization/NewtonRaphsonOptimizer.h
        )

set(SOURCE_FILES
        source/QuaternionFit.cpp
        )

add_library(${PROJECT_NAME} STATIC ${SOURCE_FILES} ${HEADER_FILES})
target_include_directories(${PROJECT_NAME} PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        )
target_link_libraries(${PROJECT_NAME}
        AmolqcGUI::LibAmolqc
        AmolqcGUI::BSplines
        Externals::cppoptlib
        eigen
        )
add_library(AmolqcGUI::${PROJECT_NAME} ALIAS ${PROJECT_NAME})


add_executable("${PROJECT_NAME}_main" main.cpp )
target_link_libraries("${PROJECT_NAME}_main" ${PROJECT_NAME})

add_subdirectory(tests)
