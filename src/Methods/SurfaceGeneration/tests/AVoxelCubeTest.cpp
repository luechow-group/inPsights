// Copyright (C) 2021 Leonard Reuter.
// SPDX-License-Identifier: GPL-3.0-or-later

#include <gmock/gmock.h>
#include <spdlog/spdlog.h>
#include <VoxelCube.h>

using namespace testing;

class AVoxelCubeTest : public ::testing::Test {
public:
    VoxelCube a;
    void SetUp() override {
        int dimension = 30;
        float length = 3.221078;
        float origin = -1.58753;
        a = VoxelCube(dimension,length,{origin + length/2, origin + length/2, origin + length/2});
        float x,y,z,f;

        auto data = a.getData();
        uint64_t weight;
        uint64_t totalWeight = 0;
        for (VoxelCube::IndexType i = 0; i < dimension; ++i) {
            for (VoxelCube::IndexType j = 0; j < dimension; ++j) {
                for (VoxelCube::IndexType k = 0; k < dimension; ++k) {
                        z = origin + k * length / dimension;
                        y = origin + j * length / dimension;
                        x = origin + i * length / dimension;
                        // 1s orbital stretched in z-direction
                        f = 1.0/sqrt(2.0*3.14)*exp(-sqrt(pow(x,2) + pow(y,2)  + pow(z/2,2)));
                        weight = unsigned(10000 * f);
                        // std::cout << x << " " << y << " " << z << " " << f << std::endl;
                        data[a.index(i,j,k)] = weight;
                }
            }
        }
        a.setData(data);
    }
};

TEST_F(AVoxelCubeTest, Wxmacmolplt) {
    a.exportMacmolplt("test.txt", "this is a test");
}
