//
// Created by Michael Heuer on 2019-01-08.
//

#ifndef INPSIGHTS_VOXELCUBEGENERATION_H
#define INPSIGHTS_VOXELCUBEGENERATION_H

#include <VoxelCube.h>

class Sample;
class SimilarReferences;

namespace VoxelCubeGeneration{

    VoxelCube<uint16_t> fromCluster(
            const std::vector<SimilarReferences> &cluster,
            const std::vector<Sample> &samples);
};

#endif //INPSIGHTS_VOXELCUBEGENERATION_H
