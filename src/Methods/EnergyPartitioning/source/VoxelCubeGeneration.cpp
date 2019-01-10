//
// Created by Michael Heuer on 2019-01-09.
//

#include <VoxelCubeGeneration.h>
#include <Sample.h>
#include <SimilarReferences.h>

//YAML_SETTINGS_DEFINITION(Settings::VoxelCubeGeneration)*/

std::vector<VoxelCube<uint16_t>> VoxelCubeGeneration::fromCluster(const std::vector<SimilarReferences> &cluster,
                                                     const std::vector<Sample> &samples) {
    std::vector<VoxelCube<uint16_t>> voxels;

    auto firstMax = cluster[0].representativeReference().maximum();//cluster[0].similarReferencesIterators()[0].base()->maximum();//TODO use averaged point

    for (long i = 0; i < firstMax.numberOfEntities(); ++i) {
        VoxelCube<uint16_t> voxel(16, 2.0*ConversionFactors::angstrom2bohr, firstMax[i].position().cast<VoxelCube<uint16_t>::VertexComponentsType>()); //TODO use averaged point
        for (const auto &simRef : cluster) {
            for (const auto &ref : simRef.similarReferencesIterators()) {

                for (auto sample : ref.base()->sampleIds()) {
                    auto electrons = samples[sample].sample_;
                    auto p = electrons[i].position();
                    voxel.add(p);
                }
            }
        }
        voxels.emplace_back(voxel);
    }

    return voxels;
}
