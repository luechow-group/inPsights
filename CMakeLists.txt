cmake_minimum_required(VERSION 3.7)
project(LocalSpinMultiplicity)

set(CMAKE_CXX_STANDARD 11)

include_directories(include)

set(SOURCE_FILES src/main.cpp include/pse.h src/Particle.cpp include/Particle.h src/Electron.cpp src/Electron.h src/Core.cpp src/Core.h)
add_executable(LocalSpinMultiplicity ${SOURCE_FILES})