# Download and unpack yaml-cpp Library at configure time

message("Configuring yaml-cpp...")

set(YAML_CPP_BUILD_TESTS OFF CACHE BOOL "Don't build the YAML tests.")
set(YAML_CPP_BUILD_CONTRIB OFF CACHE BOOL "Don't build YAML contributions.")

configure_file(cmake/downloadYamlCpp.in
        yaml-cpp-download/CMakeLists.txt)
execute_process(COMMAND ${CMAKE_COMMAND} -G "${CMAKE_GENERATOR}" .
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/yaml-cpp-download )
execute_process(COMMAND ${CMAKE_COMMAND} --build .
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/yaml-cpp-download )

# Add yaml-cpp library directly to our build. This adds
# the following targets: yaml-cpp

add_subdirectory(
        ${CMAKE_BINARY_DIR}/yaml-cpp-src
        ${CMAKE_BINARY_DIR}/yaml-cpp-build
        )