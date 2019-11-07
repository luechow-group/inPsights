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
#include <spdlog/spdlog.h>
#include <VoxelCubeOverlapCalculation.h>

using namespace testing;

class AVoxelCubeOverlapCalculationTest : public ::testing::Test {
public:
    VoxelCube a, b, c;
    void SetUp() override {
        auto dimension = 10;
        auto length = 10.0;
        a = VoxelCube(dimension,length,{0,0,0});
        b = VoxelCube(dimension,length,{0,0,0});
        c = VoxelCube(dimension,length,{0,0,0});

        for (VoxelCube::IndexType i = 0; i < dimension; ++i) {
            for (VoxelCube::IndexType j = 0; j < dimension; ++j) {
                for (VoxelCube::IndexType k = 0; k < dimension; ++k) {
                    if(i < dimension/2)
                        a.add(i, j, k,10);

                    if(j < dimension/2)
                        b.add(i, j, k, 15);

                    if(k < dimension/2)
                        c.add(i, j, k,5);
                }
            }
        }
    }
};

TEST_F(AVoxelCubeOverlapCalculationTest, OverlappingCubeHalves) {
    auto overlaps = VoxelCubeOverlapCalculation::calculateOverlaps({a,b,c});

    Eigen::MatrixXd expectedOverlaps(3,3);
    expectedOverlaps <<
    1.0, 0.5, 0.5,\
    0.5, 1.0, 0.5,\
    0.5, 0.5, 1.0;

    ASSERT_TRUE(overlaps.isApprox(expectedOverlaps));
}
