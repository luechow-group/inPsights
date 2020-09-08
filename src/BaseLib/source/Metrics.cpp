// Copyright (C) 2020 Michael Heuer.
// SPDX-License-Identifier: GPL-3.0-or-later

#include <Metrics.h>

Eigen::Vector3d Metrics::averagedPosition(const PositionsVector &positions){
    Eigen::Vector3d sum = Eigen::Vector3d::Zero();

    for (long i = 0; i < positions.numberOfEntities(); ++i)
        sum += positions[i];

    return sum / positions.numberOfEntities();;
}
