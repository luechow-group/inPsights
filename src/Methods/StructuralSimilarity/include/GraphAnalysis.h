// Copyright (C) 2019-2020 Michael Heuer.
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef INPSIGHTS_GRAPHANALYSIS_H
#define INPSIGHTS_GRAPHANALYSIS_H

#include <Eigen/Core>
#include <queue>
#include <set>
#include <algorithm>
#include <map>

namespace Eigen {
    using MatrixXb = Matrix<bool, Eigen::Dynamic, Eigen::Dynamic>;
    using VectorXb = Matrix<bool, Eigen::Dynamic, 1>;
}
namespace GraphAnalysis {
    // filter to obtain an adjacency matrix
    Eigen::MatrixXb filter(const Eigen::MatrixXd &matrix, double threshold = 1.0);

    Eigen::MatrixXb lowerOrEqualFilter(const Eigen::MatrixXd &matrix, double threshold = 1.0);

    std::set<Eigen::Index> findVerticesOfOutgoingEdges(const Eigen::MatrixXb &adjacencyMatrix, Eigen::Index vertes);

    std::set<Eigen::Index> findVerticesOfIncomingEdges(const Eigen::MatrixXb &adjacencyMatrix, Eigen::Index vertex);

    std::set<Eigen::Index> findConnectedVertices(const Eigen::MatrixXb &adjacencyMatrix, Eigen::Index startVertex);

    std::vector<std::set<Eigen::Index>> findGraphClusters(const Eigen::MatrixXb &adjacencyMatrix);


    std::map<std::size_t, std::size_t> findMergeMap(
            const std::vector<std::set<Eigen::Index>> &subsets,
            const std::vector<std::set<Eigen::Index>> &referenceSets);
}

#endif //INPSIGHTS_GRAPHANALYSIS_H
