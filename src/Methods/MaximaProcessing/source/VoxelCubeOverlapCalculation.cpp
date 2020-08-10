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

#include "VoxelCubeOverlapCalculation.h"
#include <VoxelCubeGeneration.h>

namespace Settings {
    VoxelCubeOverlapCalculation::VoxelCubeOverlapCalculation()
            : ISettings(VARNAME(VoxelCubeOverlapCalculation::)) {}

    VoxelCubeOverlapCalculation::VoxelCubeOverlapCalculation(const YAML::Node &node)
            : VoxelCubeOverlapCalculation() {
        boolProperty::decode(node[className], calculateOverlapQ);
        unsignedShortProperty::decode(node[className], dimension);
        YAML::convert<Property<VoxelCube::VertexComponentsType>>::decode(node[className], length);
    };

    void VoxelCubeOverlapCalculation::appendToNode(YAML::Node &node) const {
        node[className][calculateOverlapQ.name()] = calculateOverlapQ.get();
        node[className][dimension.name()] = dimension.get();
        node[className][length.name()] = length.get();
    };
}
YAML_SETTINGS_DEFINITION(Settings::VoxelCubeOverlapCalculation)


Eigen::MatrixXd VoxelCubeOverlapCalculation::fromCluster(const Cluster &maxima, const std::vector<Sample> &samples){
    auto dimension = settings.dimension.get();
    auto length = settings.length.get();

     auto voxels = VoxelCubeGeneration::getVoxels(maxima, samples, dimension, length,
             false, false,0);

    return calculateOverlaps(voxels);
}

Eigen::MatrixXd VoxelCubeOverlapCalculation::calculateOverlaps(const std::vector<VoxelCube> &voxels) {
    assert(!voxels.empty() && "The voxel vector cannot be empty.");
    auto dimension = voxels.front().dimension_;
    Eigen::MatrixXd overlap = Eigen::MatrixXd::Identity(voxels.size(), voxels.size());

    // iterate over upper triangular matrix
    for (Eigen::Index a = 0; a < overlap.rows() - 1; ++a) {
        for (Eigen::Index b = a + 1; b < overlap.cols(); ++b) {

            assert(voxels[a].origin_.isApprox(voxels[b].origin_) && "Voxels must have the same center.");

            for (VoxelCube::IndexType i = 0; i < dimension; ++i) {
                for (VoxelCube::IndexType j = 0; j < dimension; ++j) {
                    for (VoxelCube::IndexType k = 0; k < dimension; ++k) {
                        overlap(a, b) += std::sqrt(double(voxels[a].getData(i, j, k)))
                                * std::sqrt(double(voxels[b].getData(i, j, k)));
                    }
                }
            }
            overlap(a, b) /= std::sqrt(voxels[a].totalWeight_ * voxels[b].totalWeight_);
        }
    }
    // symmetrization
    overlap = overlap.selfadjointView<Eigen::Upper>();

    return overlap;
}