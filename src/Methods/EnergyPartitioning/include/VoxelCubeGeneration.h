//
// Created by Michael Heuer on 2019-01-08.
//

#ifndef INPSIGHTS_VOXELCUBEGENERATION_H
#define INPSIGHTS_VOXELCUBEGENERATION_H

#include <VoxelCube.h>

class Sample;
class SimilarReferences;

//#include <ISettings.h>
//#include <Property.h>
//namespace Settings {
//    class VoxelCubeGeneration : public ISettings {
//        inline static const std::string className = {VARNAME(VoxelCubeGeneration)};
//    public:
//        Property<uint16_t > dimension = {16, VARNAME(dimension)};
//        Property<double> length = {2, VARNAME(length)};
//        Property<Eigen::Vector3f> origin = {{0,0,0}, VARNAME(origin)};
//
//        VoxelCubeGeneration();
//        explicit VoxelCubeGeneration(const YAML::Node &node);
//        void appendToNode(YAML::Node &node) const override;
//    };
//}
//YAML_SETTINGS_DECLARATION(Settings::VoxelCubeGeneration)
//

namespace VoxelCubeGeneration{

    std::vector<VoxelCube<uint16_t>> fromCluster(
            const std::vector<SimilarReferences> &cluster,
            const std::vector<Sample> &samples);
};

#endif //INPSIGHTS_VOXELCUBEGENERATION_H
