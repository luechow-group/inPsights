cmake_minimum_required(VERSION 3.5)

project(BSplines CXX)


set ( PROJECT_FOLDER_NAME ${PROJECT_NAME})
set ( PROJECT_PATH ${PROJECT_NAME})
set ( OUTPUT_NAME "${PROJECT_FOLDER_NAME}" )
message ( "-- ++ CREATING PROJECT OR MAKEFILE ---> ${OUTPUT_NAME}" )

set(CMAKE_CXX_STANDARD 11)

include_directories(
        include
        source
)

set(HEADER_FILES
        include/ArcLengthParametrizedBSpline.h
        include/BSpline.h
        include/BSplineBasis.h
        include/BSplineFromControlPolygon.h
        include/BSplineFromPenalizedLeastSquaresFitWithFixedEnds.h
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
        source/BSplineFromPenalizedLeastSquaresFitWithFixedEnds.cpp
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

add_library(${PROJECT_NAME} STATIC ${SOURCE_FILES} ${HEADER_FILES})
target_link_libraries(${PROJECT_NAME}
        eigen
        )
target_include_directories(${PROJECT_NAME} PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        )

add_subdirectory(tests)
