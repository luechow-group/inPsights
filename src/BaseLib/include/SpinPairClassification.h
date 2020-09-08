// Copyright (C) 2018-2019 Michael Heuer.
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef INPSIGHTS_SPINPAIRCLASSIFICATION_H
#define INPSIGHTS_SPINPAIRCLASSIFICATION_H

#include <ParticlesVector.h>

namespace SpinPairClassification{
    enum class PairType { atSamePosition, closeBy, distant };

    using PairTypeMap = std::map<std::pair<Eigen::Index, Eigen::Index>, SpinPairClassification::PairType>;

    PairTypeMap classify(const ElectronsVector& electronsVector,
            double maxDistance = std::numeric_limits<double>::max(),
            double identicalThreshold = 0.01);

    bool atSamePositionQ(const std::pair<Eigen::Index, Eigen::Index>& pair,
            const ElectronsVector& electrons, double identicalThreshold = 0.01);

    bool closeByQ(const std::pair<Eigen::Index, Eigen::Index>& pair,
                         const ElectronsVector& electronsVector, double maxDistance = 1.6);

    bool isAtSamePositionQ(const PairTypeMap& pairTypes , Eigen::Index i);
}

#endif //INPSIGHTS_SPINPAIRCLASSIFICATION_H
