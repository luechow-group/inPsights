/* Copyright (C) 2020 Michael Heuer.
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

#include "GroupAnalysis.h"
#include <BestMatchDistance.h>
#include <Reference.h>

// constructs an adjacency matrix with indices being determined from the order in the group
Eigen::MatrixXd  GroupAnalysis::calculateAdjacencyMatrix(const Group& group) {
    assert(!group.empty() && "The group cannot be empty.");

    // construct matrix
    Eigen::MatrixXd mat = Eigen::MatrixXd::Zero(group.size(),group.size());
    Eigen::Index nClusters = group.size();
    for (Eigen::Index i = 0; i < nClusters-1; ++i) {
        for (Eigen::Index j = i+1; j < nClusters; ++j) {

            auto [norm, perm] = BestMatch::Distance::compare<Spin, Eigen::Infinity, 2>(
                    group[i].representative()->maximum(),
                    group[j].representative()->maximum());
            mat(i,j) = norm;
        }
    }
    // symmetrization
    return mat.selfadjointView<Eigen::Upper>();
}