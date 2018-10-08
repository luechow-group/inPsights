//
// Created by Michael Heuer on 28.08.18.
//

#ifndef AMOLQCPP_SAMPLE_H
#define AMOLQCPP_SAMPLE_H

#include <utility>
#include <ParticlesVector.h>

class Sample{
public:
    Sample(ElectronsVector sample, Eigen::VectorXd kineticEnergies)
    : sample_(std::move(sample)), kineticEnergies_(std::move(kineticEnergies)) {
        assert(sample.numberOfEntities() == kineticEnergies_.size());
    }

public:
    ElectronsVector sample_;
    Eigen::VectorXd kineticEnergies_;
};

#endif //AMOLQCPP_SAMPLE_H