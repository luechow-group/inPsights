// Copyright (C) 2019 Michael Heuer.
// SPDX-License-Identifier: GPL-3.0-or-later

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

TEST(AEigenYamlConversionTest, MatrixXd){
    Eigen::MatrixXd mat(3,3);
    mat << 1,2,3,4,5,6,7,8,9;

    YAML::Emitter out;
    out << mat;
    auto node = YAML::convert<Eigen::MatrixXd>::encode(mat);
    
    Eigen::MatrixXd decodedMatrix;
    YAML::convert<Eigen::MatrixXd>::decode(node,decodedMatrix);

    ASSERT_TRUE(decodedMatrix.isApprox(mat));
};

TEST(AEigenYamlConversionTest, TriangularMatrixXd){
    Eigen::MatrixXd mat(3,3);
    mat << 1,2,3,4,5,6,7,8,9;

    YAML::Emitter out;
    out << mat.selfadjointView<Eigen::Upper>();
    auto node = YAML::convert<TriangularMatrixXd>::encode(mat.selfadjointView<Eigen::Upper>());

    Eigen::MatrixXd id = Eigen::MatrixXd::Identity(3,3);
    TriangularMatrixXd decodedMatrix = id.selfadjointView<Eigen::Upper>();
    YAML::convert<TriangularMatrixXd>::decode(node,decodedMatrix);

    Eigen::MatrixXd expected(3,3);
    expected << 0,2,3,0,0,6,0,0,0;

    ASSERT_TRUE(decodedMatrix.nestedExpression().isApprox(expected));
};