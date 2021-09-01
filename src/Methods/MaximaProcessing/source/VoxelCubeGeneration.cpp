// Copyright (C) 2019 Michael Heuer.
// Copyright (C) 2021 Leonard Reuter.
// SPDX-License-Identifier: GPL-3.0-or-later

#include <VoxelCubeGeneration.h>
#include <Sample.h>
#include <Maximum.h>
#include <Cluster.h>

namespace Settings {
    VoxelCubeGeneration::VoxelCubeGeneration()
    : ISettings(VARNAME(VoxelCubeGeneration::)) {}

    VoxelCubeGeneration::VoxelCubeGeneration(const YAML::Node &node)
    : VoxelCubeGeneration() {
        boolProperty::decode(node[className], generateVoxelCubesQ);
        boolProperty::decode(node[className], centerCubesAtElectronsQ);
        boolProperty::decode(node[className], smoothingQ);
        unsignedShortProperty::decode(node[className], smoothingNeighbors);

        // if statements for backwards compatibility with old keywords dimension and length
        if (node[className]["dimension"]) {
            auto dimension = node[className]["dimension"].as<VoxelCube::IndexType>();
            dimensions = {{dimension, dimension, dimension}, VARNAME(dimensions)};
        }
        else{
            YAML::convert<Property<Eigen::Matrix<VoxelCube::IndexType, 3, 1>>>::decode(node[className], dimensions);
        }

        if (node[className]["length"]) {
            auto length = node[className]["length"].as<VoxelCube::VertexComponentsType>();
            lengths = {{length, length, length}, VARNAME(lengths)};
        }
        else{
            YAML::convert<Property<Eigen::Matrix<VoxelCube::VertexComponentsType, 3, 1>>>::decode(node[className], lengths);
        }

        YAML::convert<Property<Eigen::Matrix<VoxelCube::VertexComponentsType, 3, 1>>>::decode(node[className], center);
    };

    void VoxelCubeGeneration::appendToNode(YAML::Node &node) const {
        node[className][generateVoxelCubesQ.name()] = generateVoxelCubesQ.get();
        node[className][centerCubesAtElectronsQ.name()] = centerCubesAtElectronsQ.get();
        node[className][dimensions.name()] = dimensions.get();
        node[className][smoothingQ.name()] = smoothingQ.get();
        node[className][smoothingNeighbors.name()] = smoothingNeighbors.get();
        node[className][lengths.name()] = lengths.get();
        node[className][center.name()] = center.get();
    };
}
YAML_SETTINGS_DEFINITION(Settings::VoxelCubeGeneration)

std::vector<VoxelCube> VoxelCubeGeneration::fromCluster(const Cluster &maxima, const std::vector<Sample> &samples) {
    auto dimensions = settings.dimensions.get();
    auto lengths = settings.lengths.get();
    auto origin = settings.center.get();
    auto smoothingQ = settings.smoothingQ.get();
    auto smoothingNeighbors = settings.smoothingNeighbors.get();
    auto centerCubesAtElectronsQ = settings.centerCubesAtElectronsQ.get();

    return getVoxels(maxima, samples, dimensions, lengths, centerCubesAtElectronsQ, smoothingQ, smoothingNeighbors, origin);
}

std::vector<VoxelCube>
VoxelCubeGeneration::getVoxels(const Cluster &maxima, const std::vector<Sample> &samples,
                               Eigen::Matrix<VoxelCube::IndexType, 3, 1> dimensions,
                               Eigen::Matrix<VoxelCube::VertexComponentsType, 3, 1> lengths, bool centerCubesAtElectronsQ,
                               bool smoothingQ, VoxelCube::IndexType smoothingNeighbors,
                               Eigen::Matrix<VoxelCube::VertexComponentsType , 3, 1> center) {

    ElectronsVector representativeMax = maxima.representative()->maximum(); //TODO use averaged point
    std::vector<VoxelCube> voxels(static_cast<unsigned long>(representativeMax.numberOfEntities()));

    auto allSampleIds = maxima.allSampleIds();

    for (long i = 0; i <representativeMax.numberOfEntities(); ++i) {
        Eigen::Vector3f voxelCenter = center;
        if(centerCubesAtElectronsQ)
            voxelCenter = representativeMax.positionsVector()[i].cast<float>();

        VoxelCube voxel(dimensions, lengths, voxelCenter, smoothingQ);

        for (auto id : allSampleIds)
            voxel.add(samples[id].sample_[i].position());

        if(smoothingQ)
            voxel.smooth(smoothingNeighbors);

        voxels[i] = voxel;
    }

    return voxels;
}
