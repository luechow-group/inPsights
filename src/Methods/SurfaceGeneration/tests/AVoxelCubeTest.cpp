// Copyright (C) 2021 Leonard Reuter.
// SPDX-License-Identifier: GPL-3.0-or-later

#include <gmock/gmock.h>
#include <spdlog/spdlog.h>
#include <VoxelCube.h>

using namespace testing;

class AVoxelCubeTest : public ::testing::Test {
};

TEST_F(AVoxelCubeTest, Wxmacmolplt) {
    int dimension = 30;
    float length = 3.221078;
    auto origin = -length / 2;
    auto a = VoxelCube(dimension, length);
    float x, y, z, f;

    uint64_t weight;
    for (VoxelCube::IndexType i = 0; i < dimension; ++i) {
        for (VoxelCube::IndexType j = 0; j < dimension; ++j) {
            for (VoxelCube::IndexType k = 0; k < dimension; ++k) {
                x = origin + float(i) * length / float(dimension);
                y = origin + float(j) * length / float(dimension);
                z = origin + float(k) * length / float(dimension);
                // 1s orbital stretched in z-direction
                f = float(1.0 / sqrt(2.0 * 3.14) * exp(-sqrt(pow(x, 2) + pow(y, 2) + pow(z / 2, 2))));
                weight = unsigned(10000 * f);
                a.add(i, j, k, weight);
            }
        }
    }
    a.exportMacmolplt("test.txt", "this is a test");
}

TEST_F(AVoxelCubeTest, Wxmacmolplt2) {
    Eigen::Matrix<VoxelCube::IndexType, 3, 1> dimensions({30, 30, 30 / 2});
    Eigen::Matrix<VoxelCube::VertexComponentsType, 3, 1> lengths({3.221078, 3.221078, 3.221078 / 2});
    auto a = VoxelCube(dimensions, lengths);
    float x, y, z, f;

    auto origin = -lengths / 2;
    uint64_t weight;
    for (VoxelCube::IndexType i = 0; i < dimensions[0]; ++i) {
        for (VoxelCube::IndexType j = 0; j < dimensions[1]; ++j) {
            for (VoxelCube::IndexType k = 0; k < dimensions[2]; ++k) {
                x = origin[0] + float(i) * lengths[0] / float(dimensions[0]);
                y = origin[1] + float(j) * lengths[1] / float(dimensions[1]);
                z = origin[2] + float(k) * lengths[2] / float(dimensions[2]);

                // 1s orbital stretched in z-direction
                f = float(1.0 / sqrt(2.0 * 3.14) * exp(-sqrt(pow(x, 2) + pow(y, 2) + pow(z / 2, 2))));
                weight = unsigned(10000 * f);
                a.add(i, j, k, weight);
            }
        }
    }
    a.exportMacmolplt("test2.txt", "this is another test");
}

TEST_F(AVoxelCubeTest, Wxmacmolplt3) {
    Eigen::Matrix<VoxelCube::IndexType, 3, 1> dimensions({30, 30, 30 / 2});
    Eigen::Matrix<VoxelCube::VertexComponentsType, 3, 1> lengths({3.221078, 3.221078, 3.221078 / 2});
    auto a = VoxelCube(dimensions, lengths);
    float x, y, z, f;
    Eigen::Vector3d pos;

    auto origin = -lengths / 2;
    uint64_t weight;
    for (VoxelCube::IndexType i = 0; i < dimensions[0]; ++i) {
        for (VoxelCube::IndexType j = 0; j < dimensions[1]; ++j) {
            for (VoxelCube::IndexType k = 0; k < dimensions[2]; ++k) {
                x = origin[0] + float(i) * lengths[0] / float(dimensions[0]);
                y = origin[1] + float(j) * lengths[1] / float(dimensions[1]);
                z = origin[2] + float(k) * lengths[2] / float(dimensions[2]);
                pos = Eigen::Vector3d({x, y, z});
                // 1s orbital stretched in z-direction
                f = float(1.0 / sqrt(2.0 * 3.14) * exp(-sqrt(pow(x, 2) + pow(y, 2) + pow(z / 2, 2))));
                weight = unsigned(10000 * f);
                a.add(pos, weight);
            }
        }
    }
    a.exportMacmolplt("test3.txt", "this is another test");
}
