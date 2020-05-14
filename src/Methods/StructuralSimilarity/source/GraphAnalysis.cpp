/* Copyright (C) 2019 Michael Heuer.
 *
 * This file is part of inPsights.
 * inPsights is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * inPsights is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with inPsights. If not, see <https://www.gnu.org/licenses/>.
 */

#include <GraphAnalysis.h>
#include <numeric>
#include <Enumerate.h>
#include <algorithm>

#include <iostream>
namespace GraphAnalysis {
    Eigen::MatrixXb filter(const Eigen::MatrixXd &matrix, double threshold) {
        assert((matrix.array() >= 0.0).all());
        assert(threshold >= 0.0);

        return matrix.unaryExpr([&](const double x) { return (x >= threshold) ? 1.0 : 0.0; }).cast<bool>();
    }

    Eigen::MatrixXb lowerOrEqualFilter(const Eigen::MatrixXd &matrix, double threshold) {
        assert((matrix.array() >= 0.0).all());
        assert(threshold >= 0.0);

        return matrix.unaryExpr([&](const double x) { return (x <= threshold) ? 1.0 : 0.0; }).cast<bool>();
    }

    std::set<Eigen::Index> findVerticesOfIncomingEdges(const Eigen::MatrixXb &adjacencyMatrix, Eigen::Index vertex) {
        std::set<Eigen::Index> incomingVertexIndices;

        for (Eigen::Index i = 0; i < adjacencyMatrix.rows(); ++i)
            if (adjacencyMatrix(i, vertex))
                incomingVertexIndices.emplace(i);

        return incomingVertexIndices;
    }

    std::set<Eigen::Index> findVerticesOfOutgoingEdges(const Eigen::MatrixXb &adjacencyMatrix, Eigen::Index vertex) {
        std::set<Eigen::Index> outgoingVertexIndices;

        for (Eigen::Index j = 0; j < adjacencyMatrix.cols(); ++j)
            if (adjacencyMatrix(j, vertex))
                outgoingVertexIndices.emplace(j);

        return outgoingVertexIndices;
    }

    std::set<Eigen::Index> findConnectedVertices(const Eigen::MatrixXb &adjacencyMatrix, Eigen::Index startVertex) {
        assert(adjacencyMatrix.rows() == adjacencyMatrix.cols());

        Eigen::Index vertexCount = adjacencyMatrix.rows();
        assert(vertexCount > 0);
        assert(startVertex < vertexCount);

        Eigen::VectorXb marks = Eigen::VectorXb::Constant(vertexCount, false);
        marks(startVertex) = true;


        std::queue<Eigen::Index> queue;
        queue.push(startVertex);

        // breadth first search
        while (!queue.empty()) {
            for (Eigen::Index i = 0; i < vertexCount; ++i)
                if (adjacencyMatrix(queue.front(), i) && !marks(i)) {
                    marks(i) = true;
                    queue.push(i);
                }
            queue.pop();
        }

        std::set<Eigen::Index> connectedVertices;

        for (Eigen::Index i = 0; i < vertexCount; ++i)
            if (marks(i)) connectedVertices.emplace(i);

        return connectedVertices;
    }

    // returns a vector of sets containing connected components
    std::vector<std::set<Eigen::Index>> findGraphClusters(const Eigen::MatrixXb &adjacencyMatrix) {
        assert(adjacencyMatrix.rows() == adjacencyMatrix.cols());

        Eigen::Index vertexCount = adjacencyMatrix.rows();
        assert(vertexCount > 0);

        // create
        std::set<Eigen::Index> remainingVertices;
        for (Eigen::Index i = 0; i < vertexCount; ++i) {
            remainingVertices.emplace(i);
        }

        std::vector<std::set<Eigen::Index>> clusters;

        while (!remainingVertices.empty()) {
            // breadth-first search of connected vertices in the adjacency matrix
            // starting at the first of the remaining vertices
            auto connectedVertices = findConnectedVertices(adjacencyMatrix, *std::begin(remainingVertices));
            clusters.push_back(connectedVertices);

            std::set<Eigen::Index> difference;

            // determine remaining vertices from the difference to the newly found connected vertices
            std::set_difference(
                    remainingVertices.begin(), remainingVertices.end(),
                    connectedVertices.begin(), connectedVertices.end(),
                    std::inserter(difference, std::end(difference)));
            remainingVertices = difference;
        };

        return clusters;
    }
}

std::map<std::size_t, std::size_t> GraphAnalysis::findMergeMap(
        const std::vector<std::set<Eigen::Index>> &subsets,
        const std::vector<std::set<Eigen::Index>> &referenceSets) {
    // identify, which sets are subsets of the previous ones
    std::map<std::size_t, std::size_t> map;

    std::vector<bool> foundQ(subsets.size(), false);
    for (const auto &[referenceSetIndex, referenceSet]  : enumerate(referenceSets)) {

        bool matchedQ = false;
        for (const auto &[subsetIndex, subset] : enumerate(subsets)) {
            auto isSubsetQ = std::includes(referenceSet.begin(), referenceSet.end(), subset.begin(), subset.end());
            if (isSubsetQ) {
                map[subsetIndex] = referenceSetIndex;
                foundQ[subsetIndex] = true;
                matchedQ = true;
            }
        }
        assert(matchedQ && "A set of the reference sets has no matching subset.");
    }
    assert(std::all_of(foundQ.begin(), foundQ.end(), [](bool foundQ) {
        return foundQ == true;
    }) && "All subsets should appear in the reference set.");

    return map;
};