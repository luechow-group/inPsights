# Download and unpack Eigen Library at configure time
configure_file(cmake/DownloadEigen.in
        eigen-download/CMakeLists.txt)
execute_process(COMMAND ${CMAKE_COMMAND} -G "${CMAKE_GENERATOR}" .
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/eigen-download )
execute_process(COMMAND ${CMAKE_COMMAND} --build .
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/eigen-download )

# Add Eigen library directly to our build. This adds
# the following targets: eigen

add_subdirectory(
        ${CMAKE_BINARY_DIR}/eigen-src
        ${CMAKE_BINARY_DIR}/eigen-build
        )
