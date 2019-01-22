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

    std::cout << node << std::endl;
    std::cout << out.c_str() << std::endl;
};