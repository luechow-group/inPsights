// Copyright (C) 2018-2019 Michael Heuer.
// SPDX-License-Identifier: GPL-3.0-or-later

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
