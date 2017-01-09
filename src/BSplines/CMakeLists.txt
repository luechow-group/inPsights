cmake_minimum_required(VERSION 3.6)


#============= PROJECT BSplines ===================
project(BSplines)

set ( PROJECT_FOLDER_NAME ${PROJECT_NAME})

set ( PROJECT_PATH ${PROJECT_NAME})
set ( OUTPUT_NAME "${PROJECT_FOLDER_NAME}" )
message ( "-- ++ CREATING PROJECT OR MAKEFILE ---> ${OUTPUT_NAME}" )

set(CMAKE_CXX_STANDARD 11) # TODO needed?
set(CMAKE_CXX_FLAGS "-g -Wall") # TODO needed?

# Find includes in corresponding build directories
#set(CMAKE_INCLUDE_CURRENT_DIR ON) # TODO needed?
# Instruct CMake to run moc automatically when needed.
#set(CMAKE_AUTOMOC ON) # TODO needed?

# Find Eigen lib
find_package(Eigen3 REQUIRED)

include_directories(
        include
        source
        ${GTEST_INCLUDE_DIRS}
        #${Gmock_INCLUDE_DIR} ${Gtest_INCLUDE_DIR}
        ${EIGEN3_INCLUDE_DIR})

set(HEADER_FILES
        include/ArcLengthParametrizedBSpline.h
        include/BSpline.h
        include/BSplineBasis.h
        include/BSplineFromControlPolygon.h
        #include/BSplineFromPenalizedLeastSquaresFitWithFixedEnds.h
        include/BSplineFromPenalizedLeastSquaresFitWithLooseEnds.h
        include/BSplineFromPointInterpolation.h
        include/BSplineGenerator.h
        include/BSplineKnotInserter.h
        include/BSplineMerger.h
        include/BSplineRootFinder.h
        include/BSplineSegmenter.h
        include/BSplineSplitter.h
        include/BSplineStationaryPointFinder.h
        include/BSplineTools.h
        include/IRCTools.h
        include/LowePathSimplifier.h
        include/LoweRDPCombinedPathSimplifier.h
        include/NumericalIntegration.h
        include/RamerDouglasPeuckerPathSimplifier.h
        include/ContainerConverter.h
        )

set(SOURCE_FILES
        source/ArcLengthParametrizedBSpline.cpp
        source/BSpline.cpp
        source/BSplineBasis.cpp
        source/BSplineFromControlPolygon.cpp
        #source/BSplineFromPenalizedLeastSquaresFitWithFixedEnds.cpp
        source/BSplineFromPenalizedLeastSquaresFitWithLooseEnds.cpp
        source/BSplineFromPointInterpolation.cpp
        source/BSplineGenerator.cpp
        source/BSplineKnotInserter.cpp
        source/BSplineMerger.cpp
        source/BSplineRootFinder.cpp
        source/BSplineSegmenter.cpp
        source/BSplineSplitter.cpp
        source/BSplineStationaryPointFinder.cpp
        source/BSplineTools.cpp
        source/LowePathSimplifier.cpp
        source/LoweRDPCombinedPathSimplifier.cpp
        source/RamerDouglasPeuckerPathSimplifier.cpp
        source/ContainerConverter.cpp
        )

#add_to_source_group(include include INCLUDE_FILES)
#add_to_source_group(source source SOURCE_FILES)

add_library(${PROJECT_NAME} STATIC ${SOURCE_FILES} ${HEADER_FILES})
#add_executable(BSpline_Main main.cpp)
#target_link_libraries(BSpline_Main BSpline)


target_include_directories(${PROJECT_NAME} PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        )

add_library(AmolqcGUI::${PROJECT_NAME} ALIAS ${PROJECT_NAME})

add_subdirectory(tests)
