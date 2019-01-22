cmake_minimum_required(VERSION 3.10)

project(Executables C CXX Fortran)

set ( PROJECT_FOLDER_NAME ${PROJECT_NAME} )

set ( PROJECT_PATH ${PROJECT_NAME})
set ( OUTPUT_NAME "${PROJECT_FOLDER_NAME}" )
message ( "-- ++ CREATING PROJECT OR MAKEFILE ---> ${OUTPUT_NAME}" )


add_executable(MaximaProcessing.exe MaximaProcessing.cpp)
target_link_libraries(MaximaProcessing.exe
        BaseLib
        MaximaProcessing
        )

add_executable(StructuralSimilarity.exe StructuralSimilarityExec.cpp)
target_link_libraries(StructuralSimilarity.exe
        BaseLib
        StructuralSimilarity
        AmolqcFileImport
        )

if (DEFINED BUILD_GUI AND ON)
    add_executable(CustomMesh CustomMesh.cpp)
    target_link_libraries(CustomMesh
            VoxelProcessing # contains dualmc
            GUI
            )

    add_executable(inPsights inPsights.cpp)
    target_link_libraries(inPsights
            GUI
            MaximaProcessing
            )

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
            AmolqcElocInterface
            )

    add_executable(ShowMaxima.exe ShowMaxima.cpp)
    target_link_libraries(ShowMaxima.exe
            GUI
            AmolqcFileImport
            )

    add_executable(ShowPath.exe ShowPath.cpp)
    target_link_libraries(ShowPath.exe
            GUI
            AmolqcFileImport
            )

endif()
