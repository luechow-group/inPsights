// Copyright (C) 2018-2019 Michael Heuer.
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef INPSIGHTS_CUTOFF_H
#define INPSIGHTS_CUTOFF_H

#include <Eigen/Core>

namespace SOAP {
    namespace Cutoff {
        bool withinCutoffRadiusQ(double distance);

        double getWeight(double distanceFromExpansionCenter);

        double getWeight(const Eigen::Vector3d &position,
                         const Eigen::Vector3d &expansionCenter = Eigen::Vector3d::Zero());

        Eigen::Vector3d getWeightGradient(const Eigen::Vector3d &position);

        double getCenterWeight();

        double distance(const Eigen::Vector3d &position,
                        const Eigen::Vector3d &expansionCenter = Eigen::Vector3d::Zero());
    }
}

#endif //INPSIGHTS_CUTOFF_H
