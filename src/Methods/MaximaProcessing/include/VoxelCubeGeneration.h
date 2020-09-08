// Copyright (C) 2019 Michael Heuer.
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef INPSIGHTS_VOXELCUBEGENERATION_H
#define INPSIGHTS_VOXELCUBEGENERATION_H

#include <VoxelCube.h>

class Sample;
class Cluster;

#include <ISettings.h>
#include <Property.h>

namespace Settings {
    class VoxelCubeGeneration : public ISettings {
        inline static const std::string className = {VARNAME(VoxelCubeGeneration)};
    public:
        Property<bool> generateVoxelCubesQ = {false, VARNAME(generateVoxelCubesQ)};
        Property<bool> centerCubesAtElectronsQ = {true, VARNAME(centerCubesAtElectronsQ)};
        Property<uint16_t> dimension = {16, VARNAME(dimension)};
        Property<bool> smoothingQ = {false, VARNAME(smoothingQ)};
        Property<uint16_t> smoothingNeighbors = {0, VARNAME(smoothingNeighbors)};
        Property<VoxelCube::VertexComponentsType > length = {4, VARNAME(length)};


        VoxelCubeGeneration();
        explicit VoxelCubeGeneration(const YAML::Node &node);
        void appendToNode(YAML::Node &node) const override;
    };
}
YAML_SETTINGS_DECLARATION(Settings::VoxelCubeGeneration)


namespace VoxelCubeGeneration{
    inline Settings::VoxelCubeGeneration settings {};

    std::vector<VoxelCube> fromCluster(const Cluster &maxima, const std::vector<Sample> &samples);

    std::vector<VoxelCube> getVoxels(const Cluster &maxima, const std::vector<Sample> &samples, uint16_t dimension,
                                     VoxelCube::VertexComponentsType length, bool centerCubesAtElectronsQ, bool smoothingQ,
                                     uint16_t smoothingNeighbors);
};

#endif //INPSIGHTS_VOXELCUBEGENERATION_H
