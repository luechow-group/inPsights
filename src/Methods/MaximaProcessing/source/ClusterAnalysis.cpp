// Copyright (C) 2020 Michael Heuer.
// SPDX-License-Identifier: GPL-3.0-or-later

#include "ClusterAnalysis.h"
#include <DistanceBasedMetric.h>
#include <Maximum.h>

// constructs an adjacency matrix with indices being determined from the order in the cluster
Eigen::MatrixXd  ClusterAnalysis::calculateBestMatchDistanceMatrix(const Cluster& cluster) {
    assert(!cluster.empty() && "The cluster cannot be empty.");

    // construct matrix
    Eigen::MatrixXd mat = Eigen::MatrixXd::Zero(cluster.size(),cluster.size());
    Eigen::Index nClusters = cluster.size();
    for (Eigen::Index i = 0; i < nClusters-1; ++i) {
        for (Eigen::Index j = i+1; j < nClusters; ++j) {
            auto [norm, perm] = Metrics::Similarity::DistanceBased::compare< Eigen::Infinity, 2>(
                    cluster[i].representative()->maximum().positionsVector(),
                    cluster[j].representative()->maximum().positionsVector());
            mat(i,j) = norm;
        }
    }
    // symmetrization
    return mat.selfadjointView<Eigen::Upper>();
}