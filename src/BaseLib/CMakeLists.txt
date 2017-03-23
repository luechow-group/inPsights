cmake_minimum_required(VERSION 3.5)


#============= PROJECT Analysis ===================
project(BaseLib)

set ( PROJECT_FOLDER_NAME ${PROJECT_NAME} )

set ( PROJECT_PATH ${PROJECT_NAME})
set ( OUTPUT_NAME "${PROJECT_FOLDER_NAME}" )
message ( "-- ++ CREATING PROJECT OR MAKEFILE ---> ${OUTPUT_NAME}" )

#set(CMAKE_AUTOMOC ON)

include_directories(
        include
        source
        )

find_package(Qt5 REQUIRED COMPONENTS
        Core
        )

set(HEADER_FILES
        include/ChemicalSystem.h
        include/MaxRefParser.h
        )

set(SOURCE_FILES
        source/ChemicalSystem.cpp
        source/MaxRefParser.cpp
        )

add_library(${PROJECT_NAME} STATIC ${SOURCE_FILES} ${HEADER_FILES})

target_link_libraries(${PROJECT_NAME}
        Qt5::Core
        eigen
        )

target_include_directories(${PROJECT_NAME} PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        )

add_executable(${PROJECT_NAME}_main main.cpp)

add_library(AmolqcGUI::${PROJECT_NAME} ALIAS ${PROJECT_NAME})

add_subdirectory(tests)
