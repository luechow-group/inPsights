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
