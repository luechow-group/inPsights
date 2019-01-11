//
// Created by Michael Heuer on 28.08.18.
//

#ifndef INPSIGHTS_SAMPLE_H
#define INPSIGHTS_SAMPLE_H

#include <utility>
#include <ParticlesVector.h>

class Sample{
public:
    Sample(ElectronsVector sample, Eigen::VectorXd kineticEnergies)
    : sample_(std::move(sample)), kineticEnergies_(std::move(kineticEnergies)) {
        assert(sample.numberOfEntities() == kineticEnergies_.size());
    }

    void permute(const Eigen::PermutationMatrix<Eigen::Dynamic>& perm) {
        assert(perm.indices().size() == sample_.numberOfEntities()
        && "Permutation vector size must match with the number of entities");
        sample_.permute(perm);
        kineticEnergies_ = perm*kineticEnergies_;
    }

public:
    ElectronsVector sample_;
    Eigen::VectorXd kineticEnergies_;
};

#endif //INPSIGHTS_SAMPLE_H
