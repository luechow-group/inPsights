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

#ifndef INPSIGHTS_VOXELCUBEOVERLAPCALCULATION_H
#define INPSIGHTS_VOXELCUBEOVERLAPCALCULATION_H

#include <VoxelCube.h>
#include <ISettings.h>
#include <Property.h>
#include <Group.h>
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

    Eigen::MatrixXd fromCluster(const Group &maxima, const std::vector<Sample> &samples);
    
    Eigen::MatrixXd calculateOverlaps(const std::vector<VoxelCube>& voxels);
};

#endif //INPSIGHTS_VOXELCUBEOVERLAPCALCULATION_H
