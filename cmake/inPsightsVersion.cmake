# creates the version flag
execute_process(
        COMMAND git describe --abbrev=6 --dirty --always --tags
        WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
        OUTPUT_VARIABLE INPSIGHTS_VERSION
        OUTPUT_STRIP_TRAILING_WHITESPACE
)

if ("${CMAKE_CXX_COMPILER_ID}" MATCHES "Intel")
    set(INPSIGHTS_VERSION_FLAG "-DINPSIGHTS_VERSION=\\\"${INPSIGHTS_VERSION}\\\"")
else()
    set(INPSIGHTS_VERSION_FLAG "-DINPSIGHTS_VERSION=${INPSIGHTS_VERSION}")
endif()

message("Building inPsights: ${INPSIGHTS_VERSION_FLAG}")


# creates the CXX compiler flag
execute_process(
        COMMAND bash "-c" "${CMAKE_CXX_COMPILER} --version | head -n 1"
        WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
        OUTPUT_VARIABLE INPSIGHTS_CXX_COMPILER
        OUTPUT_STRIP_TRAILING_WHITESPACE
)

if ("${CMAKE_CXX_COMPILER_ID}" MATCHES "Intel")
    set(INPSIGHTS_CXX_COMPILER_FLAG
            "-DINPSIGHTS_CXX_COMPILER=\\\"${INPSIGHTS_CXX_COMPILER}; flags:${CMAKE_CXX_FLAGS}\\\"")
else()
    set(INPSIGHTS_CXX_COMPILER_FLAG
            "-DINPSIGHTS_CXX_COMPILER='${INPSIGHTS_CXX_COMPILER}; flags:${CMAKE_CXX_FLAGS}'")
endif()

message("C++ compiler: ${INPSIGHTS_CXX_COMPILER}; flags:${CMAKE_CXX_FLAGS}")
