//
// Created by Michael Heuer on 2019-01-08.
//

#ifndef INPSIGHTS_VOXELDATAGENERATOR_H
#define INPSIGHTS_VOXELDATAGENERATOR_H

#include <Volume.h>


class VoxelDataGenerator{
public:
    VoxelDataGenerator(const std::vector<SimilarReferences>& cluster, const std::vector<Sample>& samples)
    : cluster_(cluster) {

        Volume<uint16_t> volume;

        auto nElectrons = cluster_[0].similarReferencesIterators()[0].base()->maximum().numberOfEntities();

        long i = 0;
        //for (long i = 0; i < nElectrons; ++i) {
            for(const auto& simRef : cluster_)
                for(const auto& ref : simRef.similarReferencesIterators())
                    for(auto sample : ref.base()->sampleIds()) {
                        auto electrons = samples[sample].sample_;
                        auto p = electrons[i].position();
                        volume.add(p);
                    }
        //}
    }


    const std::vector<SimilarReferences>& cluster_;
};

#endif //INPSIGHTS_VOXELDATAGENERATOR_H
