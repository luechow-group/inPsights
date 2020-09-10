cmake_minimum_required(VERSION 3.11)

project(Executables CXX)

set ( PROJECT_FOLDER_NAME ${PROJECT_NAME} )

set ( PROJECT_PATH ${PROJECT_NAME})
set ( OUTPUT_NAME "${PROJECT_FOLDER_NAME}" )
message ( "-- ++ CREATING PROJECT OR MAKEFILE ---> ${OUTPUT_NAME}" )


add_executable(ProcessMaxima ProcessMaxima.cpp)
target_link_libraries(ProcessMaxima
        BaseLib
        MaximaProcessing
        )

add_executable(FetchOverlaps fetchOverlaps.cpp)
target_link_libraries(FetchOverlaps
        BaseLib
        MaximaProcessing
        )

add_executable(CalculateSpinStatisticsOverMultipleMolecules CalculateSpinStatisticsOverMultipleMolecules.cpp)
target_link_libraries(CalculateSpinStatisticsOverMultipleMolecules
        BaseLib
        MaximaProcessing
        )

if (${BUILD_GUI})
    add_executable(inPsights inPsights.cpp)
    target_link_libraries(inPsights
            GUI
            MaximaProcessing
            )

endif()
