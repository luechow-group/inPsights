//
// Created by Michael Heuer on 28.08.18.
//

#ifndef INPSIGHTS_SAMPLE_H
#define INPSIGHTS_SAMPLE_H

#include <utility>
#include <ParticlesVector.h>

class Sample{
public:
    Sample(ElectronsVector sample, Eigen::VectorXd kineticEnergies);

    void permute(const Eigen::PermutationMatrix<Eigen::Dynamic>& perm);

public:
    ElectronsVector sample_;
    Eigen::VectorXd kineticEnergies_;
};

#endif //INPSIGHTS_SAMPLE_H
