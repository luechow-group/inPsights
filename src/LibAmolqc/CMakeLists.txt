cmake_minimum_required(VERSION 3.5)


#============= PROJECT LibAmolqc ===================
project(LibAmolqc)

set(CMAKE_CXX_FLAGS "-g -Wall -fopenmp -llapack -lblas -lgfortran")
set(CMAKE_Fortran_FLAGS "-g -fopenmp -llapack -lblas -lstdc++")


set ( PROJECT_FOLDER_NAME ${PROJECT_NAME})

set ( PROJECT_PATH ${PROJECT_NAME})
set ( OUTPUT_NAME "${PROJECT_FOLDER_NAME}" )
message ( "-- ++ CREATING PROJECT OR MAKEFILE ---> ${OUTPUT_NAME}" )

include_directories(
        include
        source
        )

set(HEADER_FILES
        include/ElectronicWaveFunction.h
        )

set(SOURCE_FILES
        source/ElectronicWaveFunction.cpp
        )

get_filename_component(PATH_TO_LIBAMOLQC amolqc/src/libamolqc.a REALPATH)
get_filename_component(PATH_TO_LIBUTILS amolqc/utils/libutils.a REALPATH)


add_library(libamolqc STATIC IMPORTED GLOBAL)
set_property(TARGET libamolqc PROPERTY IMPORTED_LOCATION ${PATH_TO_LIBAMOLQC})

add_library(libutils STATIC IMPORTED GLOBAL)
set_property(TARGET libutils PROPERTY IMPORTED_LOCATION ${PATH_TO_LIBUTILS})

add_library(${PROJECT_NAME} STATIC ${SOURCE_FILES} ${HEADER_FILES})
target_include_directories(${PROJECT_NAME} PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        )
target_link_libraries(${PROJECT_NAME}
        libamolqc
        libutils
        eigen
        )
add_library(AmolqcGUI::${PROJECT_NAME} ALIAS ${PROJECT_NAME})

add_executable("${PROJECT_NAME}_main" main.cpp)
target_link_libraries("${PROJECT_NAME}_main" ${PROJECT_NAME}
        Externals::cppoptlib
        )

#add_custom_target(exportAmolqcPath
#        DEPENDS wrapper.sh
#        COMMAND /bin/sh ${CMAKE_CURRENT_SOURCE_DIR}/wrapper.sh #<real_command> args...
#        #COMMAND ${PROJECT_NAME}_main
#        COMMENT "exporting amolqc path"
#        )


#add_custom_target(LibAmolqcWithExport
#        COMMAND ${CMAKE_COMMAND} -E env "AMOLQC=/Users/michaelheuer/amolqcGUI/src/LibAmolqc/amolqc;$ENV{AMOLQC}"
#        export
#        )

#set($ENV{AMOLQC} "/Users/michaelheuer/amolqcGUI/src/LibAmolqc/amolqc")


#add_custom_target(LibAmolqcWithExport ALL -E env
#        "AMOLQC=/Users/michaelheuer/amolqcGUI/src/LibAmolqc/amolqc;$ENV{AMOLQC}" #export
#        )

#add_subdirectory(tests)
