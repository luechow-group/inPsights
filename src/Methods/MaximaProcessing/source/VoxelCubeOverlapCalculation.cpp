// Copyright (C) 2019 Michael Heuer.
// SPDX-License-Identifier: GPL-3.0-or-later

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

     auto voxels = VoxelCubeGeneration::getVoxels(maxima, samples, {dimension,dimension,dimension},
                                                  {length, length, length},
                                                  false, false,0,
                                                  {-length/2,-length/2,-length/2});

    return calculateOverlaps(voxels);
}

Eigen::MatrixXd VoxelCubeOverlapCalculation::calculateOverlaps(const std::vector<VoxelCube> &voxels) {
    assert(!voxels.empty() && "The voxel vector cannot be empty.");
    auto dimension = voxels.front().dimensions_;
    Eigen::MatrixXd overlap = Eigen::MatrixXd::Identity(voxels.size(), voxels.size());

    // iterate over upper triangular matrix
    for (Eigen::Index a = 0; a < overlap.rows() - 1; ++a) {
        for (Eigen::Index b = a + 1; b < overlap.cols(); ++b) {

            assert(voxels[a].center_.isApprox(voxels[b].center_) && "Voxels must have the same center.");

            for (VoxelCube::IndexType i = 0; i < dimension[0]; ++i) {
                for (VoxelCube::IndexType j = 0; j < dimension[1]; ++j) {
                    for (VoxelCube::IndexType k = 0; k < dimension[2]; ++k) {
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