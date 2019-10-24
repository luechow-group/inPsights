cmake_minimum_required(VERSION 3.11)

project(Executables C CXX Fortran)

set ( PROJECT_FOLDER_NAME ${PROJECT_NAME} )

set ( PROJECT_PATH ${PROJECT_NAME})
set ( OUTPUT_NAME "${PROJECT_FOLDER_NAME}" )
message ( "-- ++ CREATING PROJECT OR MAKEFILE ---> ${OUTPUT_NAME}" )


add_executable(ProcessMaxima ProcessMaxima.cpp)
target_link_libraries(ProcessMaxima
        BaseLib
        MaximaProcessing
        )

if (${BUILD_GUI})
    add_executable(inPsights inPsights.cpp)
    target_link_libraries(inPsights
            GUI
            MaximaProcessing
            )

    add_executable(ShowMaximum ShowMaximum.cpp)
    target_link_libraries(ShowMaximum
            GUI
            AmolqcFileImport
            )

    add_executable(ShowPath ShowPath.cpp)
    target_link_libraries(ShowPath
            GUI
            AmolqcFileImport
            )

endif()
