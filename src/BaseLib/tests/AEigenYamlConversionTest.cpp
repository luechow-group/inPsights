//
// Created by Michael Heuer on 2019-01-08.
//

#include <gmock/gmock.h>
#include <yaml-cpp/yaml.h>
#include <Eigen/Core>
#include <iostream>
#include <EigenYamlConversion.h>

TEST(AEigenYamlConversionTest, Vector3d){
    Eigen::Vector3d vec({1,2,3});

    YAML::Emitter out;
    out << vec;
    auto node = YAML::convert<Eigen::Vector3d>::encode(vec);

    Eigen::Vector3d decodedVector;
    YAML::convert<Eigen::Vector3d>::decode(node,decodedVector);

    ASSERT_TRUE(decodedVector.isApprox(vec));
};