//
// Created by Michael Heuer on 2019-01-09.
//

#include <VoxelCubeGeneration.h>
#include <Sample.h>
#include <SimilarReferences.h>

VoxelCube<uint16_t> VoxelCubeGeneration::fromCluster(const std::vector<SimilarReferences> &cluster,
                                                     const std::vector<Sample> &samples) {
    VoxelCube<uint16_t> volume;

    auto nElectrons = cluster[0].similarReferencesIterators()[0].base()->maximum().numberOfEntities();

    long i = 0;
    //for (long i = 0; i < nElectrons; ++i) {
    for (const auto &simRef : cluster)
        for (const auto &ref : simRef.similarReferencesIterators())
            for (auto sample : ref.base()->sampleIds()) {
                auto electrons = samples[sample].sample_;
                auto p = electrons[i].position();
                volume.add(p);
            }
    //}

    return volume;
}
