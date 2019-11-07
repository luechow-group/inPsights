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