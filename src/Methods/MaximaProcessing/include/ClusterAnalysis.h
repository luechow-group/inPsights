// Copyright (C) 2020 Michael Heuer.
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef INPSIGHTS_CLUSTERANALYSIS_H
#define INPSIGHTS_CLUSTERANALYSIS_H

#include "Cluster.h"

namespace ClusterAnalysis {
    Eigen::MatrixXd calculateBestMatchDistanceMatrix(const Cluster &cluster);
}

#endif //INPSIGHTS_CLUSTERANALYSIS_H
