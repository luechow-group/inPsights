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

    std::vector<VoxelCube> voxels;

    auto firstMax = cluster[0].representativeReference().maximum();//cluster[0].similarReferencesIterators()[0].base()->maximum();//TODO use averaged point

    for (long i = 0; i < firstMax.numberOfEntities(); ++i) {
        Eigen::Matrix<VoxelCube::VertexComponentsType,3,1> cubeOrigin({0,0,0});

        if(settings.centerCubesAtElectronsQ.get())  //TODO use averaged point
            cubeOrigin = firstMax[i].position().cast<VoxelCube::VertexComponentsType>();

        VoxelCube voxel(dimension, length, cubeOrigin);
        for (const auto &simRef : cluster) {
            for (const auto &ref : simRef.similarReferencesIterators()) {

                for (auto sample : ref.base()->sampleIds()) {
                    auto electrons = samples[sample].sample_;
                    auto p = electrons[i].position();
                    voxel.add(p);
                }
            }
        }
        voxels.emplace_back(voxel);
    }

    return voxels;
}
