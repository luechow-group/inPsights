// Copyright (C) 2019 Michael Heuer.
// Copyright (C) 2021 Leonard Reuter.
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef INPSIGHTS_VOXELCUBEGENERATION_H
#define INPSIGHTS_VOXELCUBEGENERATION_H

#include <VoxelCube.h>

class Sample;
class Cluster;

#include <ISettings.h>
#include <Property.h>
#include <Formatting.h>

namespace Settings {
    class VoxelCubeGeneration : public ISettings {
        inline static const std::string className = {VARNAME(VoxelCubeGeneration)};
    public:
        Property<bool> generateVoxelCubesQ = {false, VARNAME(generateVoxelCubesQ)};
        Property<bool> centerCubesAtElectronsQ = {true, VARNAME(centerCubesAtElectronsQ)};
        Property<Eigen::Matrix<VoxelCube::IndexType , 3, 1>> dimensions = {{16, 16, 16}, VARNAME(dimensions)};
        Property<bool> smoothingQ = {false, VARNAME(smoothingQ)};
        Property<uint16_t> smoothingNeighbors = {0, VARNAME(smoothingNeighbors)};
        Property<Eigen::Matrix<VoxelCube::VertexComponentsType , 3, 1>> lengths = {{4, 4, 4}, VARNAME(lengths)};
        Property<Eigen::Matrix<VoxelCube::VertexComponentsType , 3, 1>> center = {{0, 0, 0}, VARNAME(center)};


        VoxelCubeGeneration();
        explicit VoxelCubeGeneration(const YAML::Node &node);
        void appendToNode(YAML::Node &node) const override;
    };
}
YAML_SETTINGS_DECLARATION(Settings::VoxelCubeGeneration)


namespace VoxelCubeGeneration{
    inline Settings::VoxelCubeGeneration settings {};

    std::vector<VoxelCube> fromCluster(const Cluster &maxima, const std::vector<Sample> &samples);

    std::vector<VoxelCube> getVoxels(const Cluster &maxima, const std::vector<Sample> &samples,
                                     const Eigen::Matrix<VoxelCube::IndexType, 3, 1>& dimensions,
                                     const Eigen::Matrix<VoxelCube::VertexComponentsType , 3, 1>& lengths,
                                     bool centerCubesAtElectronsQ, bool smoothingQ,
                                     VoxelCube::IndexType smoothingNeighbors,
                                     const Eigen::Matrix<VoxelCube::VertexComponentsType , 3, 1>& center);
};

#endif //INPSIGHTS_VOXELCUBEGENERATION_H
