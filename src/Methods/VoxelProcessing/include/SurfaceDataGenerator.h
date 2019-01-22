//
// Created by Michael Heuer on 2019-01-07.
//

#ifndef INPSIGHTS_SURFACEDATAGENERATOR_H
#define INPSIGHTS_SURFACEDATAGENERATOR_H

#include <DualMC.h>
#include "VoxelCube.h"
#include "SurfaceData.h"

class SurfaceDataGenerator {
public:

    explicit SurfaceDataGenerator(VoxelCube volume);

    SurfaceData computeSurfaceData(double volumeThreshold = 0.5);

    VoxelCube::VolumeDataType getIsoValue(double volumeThreshold = 0.5, unsigned maxSteps = 10, double eps = 1e-2);

    VoxelCube cube_;
    std::vector<dualmc::Vertex> dualMcVertices_;
    std::vector<dualmc::Quad> dualMcQuads_;
};

#endif //INPSIGHTS_SURFACEDATAGENERATOR_H
