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
#include <Cluster.h>

namespace Settings {
    VoxelCubeGeneration::VoxelCubeGeneration()
    : ISettings(VARNAME(VoxelCubeGeneration::)) {}

    VoxelCubeGeneration::VoxelCubeGeneration(const YAML::Node &node)
    : VoxelCubeGeneration() {
        boolProperty::decode(node[className], generateVoxelCubesQ);
        boolProperty::decode(node[className], centerCubesAtElectronsQ);
        unsignedShortProperty::decode(node[className], dimension);
        boolProperty::decode(node[className], smoothingQ);
        unsignedShortProperty::decode(node[className], smoothingNeighbors);
        YAML::convert<Property<VoxelCube::VertexComponentsType>>::decode(node[className], length);
    };

    void VoxelCubeGeneration::appendToNode(YAML::Node &node) const {
        node[className][generateVoxelCubesQ.name()] = generateVoxelCubesQ.get();
        node[className][centerCubesAtElectronsQ.name()] = centerCubesAtElectronsQ.get();
        node[className][dimension.name()] = dimension.get();
        node[className][smoothingQ.name()] = smoothingQ.get();
        node[className][smoothingNeighbors.name()] = smoothingNeighbors.get();
        node[className][length.name()] = length.get();
    };
}
YAML_SETTINGS_DEFINITION(Settings::VoxelCubeGeneration)

std::vector<VoxelCube> VoxelCubeGeneration::fromCluster(const Cluster &maxima, const std::vector<Sample> &samples) {
    auto dimension = settings.dimension.get();
    auto length = settings.length.get();
    auto smoothingQ = settings.smoothingQ.get();
    auto smoothingNeighbors = settings.smoothingNeighbors.get();
    auto centerCubesAtElectronsQ = settings.centerCubesAtElectronsQ.get();

    return getVoxels(maxima, samples, dimension, length, centerCubesAtElectronsQ, smoothingQ, smoothingNeighbors);
}

std::vector<VoxelCube>
VoxelCubeGeneration::getVoxels(const Cluster &maxima, const std::vector<Sample> &samples, uint16_t dimension,
                               VoxelCube::VertexComponentsType length, bool centerCubesAtElectronsQ,
                               bool smoothingQ, uint16_t smoothingNeighbors) {

    ElectronsVector representativeMax = maxima.representative()->maximum(); //TODO use averaged point
    std::vector<VoxelCube> voxels(static_cast<unsigned long>(representativeMax.numberOfEntities()));

    auto allSampleIds = maxima.allSampleIds();

    for (long i = 0; i <representativeMax.numberOfEntities(); ++i) {
        Eigen::Vector3f cubeOrigin = Eigen::Vector3f::Zero();
        if(centerCubesAtElectronsQ)
            cubeOrigin = representativeMax.positionsVector()[i].cast<float>();

        VoxelCube voxel(dimension, length, cubeOrigin, smoothingQ);

        for (auto id : allSampleIds)
            voxel.add(samples[id].sample_[i].position());

        if(smoothingQ)
            voxel.smooth(smoothingNeighbors);

        voxels[i] = voxel;
    }

    return voxels;
}
