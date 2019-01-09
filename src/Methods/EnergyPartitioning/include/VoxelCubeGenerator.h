//
// Created by Michael Heuer on 2019-01-08.
//

#ifndef INPSIGHTS_VOXELCUBEGENERATOR_H
#define INPSIGHTS_VOXELCUBEGENERATOR_H

#include <VoxelCube.h>

class Sample;
class SimilarReferences;

namespace VoxelCubeGeneration{

    VoxelCube<uint16_t> fromCluster(
            const std::vector<SimilarReferences> &cluster,
            const std::vector<Sample> &samples);
};

#endif //INPSIGHTS_VOXELCUBEGENERATOR_H
