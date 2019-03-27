//
// Created by Michael Heuer on 2019-01-09.
//

#include <VoxelCubeGeneration.h>
#include <Sample.h>
#include <SimilarReferences.h>

namespace Settings {
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

std::vector<VoxelCube> VoxelCubeGeneration::fromCluster(const std::vector<SimilarReferences> &cluster,
                                                     const std::vector<Sample> &samples) {
    auto dimension = settings.dimension.get();
    auto length = settings.length.get();

    ElectronsVector firstMax = cluster[0].representativeReference().maximum(); //TODO use averaged point
    std::vector<VoxelCube> voxels(static_cast<unsigned long>(firstMax.numberOfEntities()));

    for (long i = 0; i <firstMax.numberOfEntities(); ++i) {
        Eigen::Vector3f cubeOrigin = Eigen::Vector3f::Zero();
        if(settings.centerCubesAtElectronsQ.get())
            cubeOrigin = firstMax.positionsVector()[i].cast<float>();

        VoxelCube voxel(dimension, length, cubeOrigin);
        for (const auto &simRef : cluster) {
            for (const auto &ref : simRef.similarReferencesIterators()) {
                for (auto sampleId : ref.base()->sampleIds()) {
                    voxel.add(samples[sampleId].sample_[i].position());
                }
            }
        }
        voxels[i] = voxel;
    }

    return voxels;
}
