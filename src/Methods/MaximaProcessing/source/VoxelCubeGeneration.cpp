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

#include <VoxelCubeGeneration.h>
#include <Sample.h>
#include <Reference.h>
#include <Group.h>


namespace Settings {
    VoxelCubeGeneration::VoxelCubeGeneration()
    : ISettings(VARNAME(VoxelCubeGeneration::)) {}

    VoxelCubeGeneration::VoxelCubeGeneration(const YAML::Node &node)
    : VoxelCubeGeneration() {
        boolProperty::decode(node[className], generateVoxelCubesQ);
        boolProperty::decode(node[className], centerCubesAtElectronsQ);
        unsignedShortProperty::decode(node[className], dimension);
        YAML::convert<Property<VoxelCube::VertexComponentsType>>::decode(node[className], length);
    };

    void VoxelCubeGeneration::appendToNode(YAML::Node &node) const {
        node[className][generateVoxelCubesQ.name()] = generateVoxelCubesQ.get();
        node[className][centerCubesAtElectronsQ.name()] = centerCubesAtElectronsQ.get();
        node[className][dimension.name()] = dimension.get();
        node[className][length.name()] = length.get();
    };
}
YAML_SETTINGS_DEFINITION(Settings::VoxelCubeGeneration)

std::vector<VoxelCube> VoxelCubeGeneration::fromCluster(const Group &maxima, const std::vector<Sample> &samples) {
    auto dimension = settings.dimension.get();
    auto length = settings.length.get();

    ElectronsVector representativeMax = maxima.representative()->maximum(); //TODO use averaged point
    std::vector<VoxelCube> voxels(static_cast<unsigned long>(representativeMax.numberOfEntities()));

    for (long i = 0; i <representativeMax.numberOfEntities(); ++i) {
        Eigen::Vector3f cubeOrigin = Eigen::Vector3f::Zero();
        if(settings.centerCubesAtElectronsQ.get())
            cubeOrigin = representativeMax.positionsVector()[i].cast<float>();

        VoxelCube voxel(dimension, length, cubeOrigin);

        auto allSampleIds = maxima.allSampleIds();

        for (auto id : allSampleIds)
            voxel.add(samples[id].sample_[i].position());

        voxels[i] = voxel;
    }

    return voxels;
}
