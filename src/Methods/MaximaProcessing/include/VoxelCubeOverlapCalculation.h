// Copyright (C) 2019 Michael Heuer.
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef INPSIGHTS_VOXELCUBEOVERLAPCALCULATION_H
#define INPSIGHTS_VOXELCUBEOVERLAPCALCULATION_H

#include <VoxelCube.h>
#include <ISettings.h>
#include <Property.h>
#include <Cluster.h>
#include <Sample.h>

namespace Settings {
    class VoxelCubeOverlapCalculation : public ISettings {
        inline static const std::string className = {VARNAME(VoxelCubeOverlapCalculation)};
    public:
        Property<bool> calculateOverlapQ = {true, VARNAME(calculateOverlapQ)};
        Property<uint16_t> dimension = {100, VARNAME(dimension)};
        Property<VoxelCube::VertexComponentsType > length = {20, VARNAME(length)};


        VoxelCubeOverlapCalculation();
        explicit VoxelCubeOverlapCalculation(const YAML::Node &node);
        void appendToNode(YAML::Node &node) const override;
    };
}
YAML_SETTINGS_DECLARATION(Settings::VoxelCubeOverlapCalculation)

namespace VoxelCubeOverlapCalculation {
    inline Settings::VoxelCubeOverlapCalculation settings {};

    Eigen::MatrixXd fromCluster(const Cluster &maxima, const std::vector<Sample> &samples);
    
    Eigen::MatrixXd calculateOverlaps(const std::vector<VoxelCube>& voxels);
};

#endif //INPSIGHTS_VOXELCUBEOVERLAPCALCULATION_H
