cmake_minimum_required(VERSION 3.10)

project(Executables C CXX Fortran)

set ( PROJECT_FOLDER_NAME ${PROJECT_NAME} )

set ( PROJECT_PATH ${PROJECT_NAME})
set ( OUTPUT_NAME "${PROJECT_FOLDER_NAME}" )
message ( "-- ++ CREATING PROJECT OR MAKEFILE ---> ${OUTPUT_NAME}" )

add_executable(BaseLib.exe BaseLib.cpp)
target_link_libraries(BaseLib.exe
        BaseLib
        )

add_executable(inPsights.exe inPsights.cpp)
target_link_libraries(inPsights.exe
        GUI
        AmolqcInterface
        EnergyPartitioning
        )

add_executable(EnergyPartitioning.exe EnergyPartitioning.cpp)
target_link_libraries(EnergyPartitioning.exe
        GUI
        AmolqcInterface
        BaseLib
        EnergyPartitioning
        )

add_executable(StructuralSimilarity.exe StructuralSimilarityExec.cpp)
target_link_libraries(StructuralSimilarity.exe
        AmolqcInterface
        BaseLib
        StructuralSimilarity
        )

if (DEFINED BUILD_GUI AND ON)
    add_executable(StringMethodGUI.exe StringMethodGUI.cpp)
    target_link_libraries(StringMethodGUI.exe
            GUI
            BaseLib
            TransitionStateSearch
            )

    add_executable(PathOptimization.exe PathOptimization.cpp)
    target_link_libraries(PathOptimization.exe
            GUI
            BaseLib
            TransitionStateSearch
            )

    add_executable(WavefunctionTransitionStates.exe WavefunctionTransitionStates.cpp)
    target_link_libraries(WavefunctionTransitionStates.exe
            GUI
            BaseLib
            TransitionStateSearch
            )

    add_executable(LocalPotentialEnergy.exe LocalPotentialEnergy.cpp)
    target_link_libraries(LocalPotentialEnergy.exe
            GUI
            LocalPotentialEnergy
            AmolqcInterface
            )

    add_executable(ShowMaxima.exe ShowMaxima.cpp)
    target_link_libraries(ShowMaxima.exe
            GUI
            AmolqcInterface
            )

    add_executable(ShowPath.exe ShowPath.cpp)
    target_link_libraries(ShowPath.exe
            GUI
            AmolqcInterface
            )

    add_executable(H6TSsearch.exe H6TSsearch.cpp)
    target_link_libraries(H6TSsearch.exe
            GUI
            BaseLib
            TransitionStateSearch
            )
endif()
