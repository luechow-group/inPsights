cmake_minimum_required(VERSION 3.11)

project(BaseLib CXX)
set ( PROJECT_FOLDER_NAME ${PROJECT_NAME} )
set ( PROJECT_PATH ${PROJECT_NAME})
set ( OUTPUT_NAME "${PROJECT_FOLDER_NAME}" )
message ( "-- ++ CREATING PROJECT OR MAKEFILE ---> ${OUTPUT_NAME}" )

set_source_files_properties(source/BuildInfo.cpp PROPERTIES COMPILE_FLAGS
        "${INPSIGHTS_VERSION_INFO_FLAG} ${INPSIGHTS_CXX_COMPILER_INFO_FLAG}"
        )

add_library(inPsightsInfo STATIC include/BuildInfo.h source/BuildInfo.cpp)
target_include_directories(inPsightsInfo
        PUBLIC include
        PRIVATE source
        )

set(HEADER_FILES
        resources/TestMolecules.h

        include/ElementColor.h
        include/ElementInfo.h
        include/pElementInfo.h
        include/ElementType.h
        include/NaturalConstants.h
        include/Particle.h
        include/ParticlesVector.h
        include/ParticlesVectorCollection.h
        include/SpinType.h
        include/Importer.h
        include/Serialization.h
        include/AbstractVector.h
        include/PositionsVector.h
        include/PositionsVectorCollection.h
        include/ToString.h
        include/PositionsVectorTransformer.h
        include/TypesVector.h
        include/Metrics.h
        include/MolecularGeometry.h
        include/EnumeratedType.h
        include/EigenYamlConversion.h
        include/SpecialMathFunctions/BoostSphericalHarmonics.h
        include/SpecialMathFunctions/ModifiedSphericalBessel1stKind.h
        include/SpecialMathFunctions/Gaussian.h
        include/CoulombPotential.h
        include/BinaryFileReader.h
        include/Interval.h
        include/DataVector.h
        include/InsertableVector.h
        include/Statistics.h
        include/SpinCorrelation.h
        include/ISettings.h
        include/Property.h
        include/Varname.h
        include/Enumerate.h
        include/Combinatorics.h
        include/Logo.h
        include/ElectronSelection.h
        include/SpinPairClassification.h
        include/CameraSettings.h
        include/ErrorHandling.h
        include/MapUtils.h
        include/ParticleIndices.h
        include/MolecularSelection.h
        include/PermutationHandling.h
        include/GraphAnalysis.h
        include/Formatting.h
        )


set(SOURCE_FILES
        source/pElementInfo.cpp
        source/ElementInfo.cpp
        source/Importer.cpp
        source/Serialization.cpp
        source/SpinType.cpp
        source/ElementType.cpp
        source/PositionsVector.cpp
        source/AbstractVector.cpp
        source/PositionsVectorCollection.cpp
        source/ToString.cpp
        source/PositionsVectorTransformer.cpp
        source/TypesVector.cpp
        source/MolecularGeometry.cpp
        source/Particle.cpp
        source/EigenYamlConversion.cpp
        source/SpecialMathFunctions/BoostSphericalHarmonics.cpp
        source/SpecialMathFunctions/ModifiedSphericalBessel1stKind.cpp
        source/SpecialMathFunctions/Gaussian.cpp
        source/BinaryFileReader.cpp
        source/ISettings.cpp
        source/Combinatorics.cpp
        source/ElectronSelection.cpp
        source/SpinPairClassification.cpp
        source/CameraSettings.cpp
        source/ErrorHandling.cpp
        source/Metrics.cpp
        source/ParticleIndices.cpp
        source/MolecularSelection.cpp
        source/PermutationHandling.cpp
        source/GraphAnalysis.cpp
        )

add_library(${PROJECT_NAME} STATIC ${SOURCE_FILES} ${HEADER_FILES})

target_include_directories(${PROJECT_NAME}
        PUBLIC include
        PRIVATE source
)

target_link_libraries(${PROJECT_NAME}
        inPsightsInfo
        yaml-cpp
        spdlog
        )

if (${BUILD_EIGEN})
        message("BUILD_EIGEN build eigen is ON.")
        target_link_libraries(${PROJECT_NAME} eigen)
endif()

target_include_directories(${PROJECT_NAME}
        PUBLIC
                include
                resources
        PRIVATE
                source
        )

add_gtest_subdirectory_with_dependency(tests "gmock;pthread")
