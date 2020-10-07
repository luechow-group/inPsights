// Copyright (C) 2018-2019 Michael Heuer.
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef INPSIGHTS_CUTOFF_H
#define INPSIGHTS_CUTOFF_H

#include <Eigen/Core>

namespace SOAP {
    namespace Cutoff {
        bool withinCutoffRadiusQ(double distance);

        double value(double distanceFromExpansionCenter);

        double value(const Eigen::Vector3d &position,
                     const Eigen::Vector3d &expansionCenter = Eigen::Vector3d::Zero());

        Eigen::Vector3d gradient(const Eigen::Vector3d &position);

        Eigen::Vector3d gradient(const Eigen::Vector3d &position,
                                 const Eigen::Vector3d &expansionCenter);
    }
}

#endif //INPSIGHTS_CUTOFF_H
