add_subdirectory(BaseLib)
add_subdirectory(Externals)
add_subdirectory(Methods)

if (${BUILD_GUI})
    message("BUILD_GUI=${BUILD_GUI} -> Trying to build the GUI. Make sure, Qt and Qt3D is available.")
    add_subdirectory(GUI)
else()
    message("BUILD_GUI=${BUILD_GUI} -> Don't build the GUI and dependent executables.")
endif()

add_subdirectory(Executables)
