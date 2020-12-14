# Download and unpack Spdlog Library at configure time

message("Configuring spdlog...")

configure_file(cmake/downloadSpdlog.in
        spdlog-download/CMakeLists.txt)
execute_process(COMMAND ${CMAKE_COMMAND} -G "${CMAKE_GENERATOR}" .
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/spdlog-download )
execute_process(COMMAND ${CMAKE_COMMAND} --build .
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/spdlog-download )

# Add Spdlog library directly to our build. This adds
# the following targets: spdlog

add_subdirectory(
        ${CMAKE_BINARY_DIR}/spdlog-src
        ${CMAKE_BINARY_DIR}/spdlog-build
        )
