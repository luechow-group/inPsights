cmake_minimum_required(VERSION 3.5)

project(Executables C CXX Fortran)

set ( PROJECT_FOLDER_NAME ${PROJECT_NAME} )

set ( PROJECT_PATH ${PROJECT_NAME})
set ( OUTPUT_NAME "${PROJECT_FOLDER_NAME}" )
message ( "-- ++ CREATING PROJECT OR MAKEFILE ---> ${OUTPUT_NAME}" )

file(GLOB_RECURSE RESOURCE_FILES resources/*.*)
file(COPY ${RESOURCE_FILES} DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

add_executable(StringMethodGUI.exe StringMethodGUI.cpp)
add_executable(PathOptimization.exe PathOptimization.cpp)
add_executable(WavefunctionTransitionStates.exe WavefunctionTransitionStates.cpp)
add_executable(LocalPotentialEnergy.exe LocalPotentialEnergy.cpp)
add_executable(ShowMaxima.exe ShowMaxima.cpp)
add_executable(H6TSsearch.exe H6TSsearch.cpp)
add_executable(StructuralSimilarity.exe StructuralSimilarityExec.cpp)

target_link_libraries(StringMethodGUI.exe
        GUI
        BaseLib
        TransitionStateSearch
        )

target_link_libraries(PathOptimization.exe
        GUI
        BaseLib
        TransitionStateSearch
        )

target_link_libraries(WavefunctionTransitionStates.exe
        GUI
        BaseLib
        TransitionStateSearch
        )

target_link_libraries(WavefunctionTransitionStates.exe
        GUI
        BaseLib
        TransitionStateSearch
        )

target_link_libraries(LocalPotentialEnergy.exe
        GUI
        LocalPotentialEnergy
        AmolqcInterface
        )

target_link_libraries(ShowMaxima.exe
        GUI
        AmolqcInterface
        )

target_link_libraries(H6TSsearch.exe
        GUI
        BaseLib
        TransitionStateSearch
        )

target_link_libraries(StructuralSimilarity.exe
        BaseLib
        AmolqcInterface
        StructuralSimilarity
        )
