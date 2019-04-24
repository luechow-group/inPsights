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