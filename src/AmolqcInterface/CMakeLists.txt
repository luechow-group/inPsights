cmake_minimum_required(VERSION 3.11)

add_subdirectory(AmolqcFileImport)

if (${BUILD_AMOLQC})
    message("BUILD_AMOLQC=${BUILD_AMOLQC} -> Building Amolqc and the eloc interface.")

    option(PFUNIT "" OFF)
    option(CONFIG_FILE "" OFF)
    option(WARNINGS "" OFF)

    add_subdirectory(Amolqc)
    add_subdirectory(AmolqcElocInterface)
else()
    message("BUILD_AMOLQC=${BUILD_AMOLQC} -> Amolqc and the eloc interface are not built.")
endif()
