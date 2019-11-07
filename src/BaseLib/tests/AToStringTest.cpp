/* Copyright (C) 2018 Leonard Reuter.
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
#include "ToString.h"

using namespace testing;

class AToStringTest : public Test {
public:
    void SetUp() override {
    }
};

TEST_F(AToStringTest, UsignedLong){
    int a = 19;
    ASSERT_EQ(ToString::longToString(a, 2)," 19");
}

TEST_F(AToStringTest, Double){
    double a = 7.845;
    ASSERT_EQ(ToString::doubleToString(a, 5, 0)," 7.84500");
}

TEST_F(AToStringTest, Vector3d){
    Eigen::Vector3d b(0,-2.01238,7.13);
    Eigen::Vector3d c(1e-10,-9.8987,-12.238);
    Eigen::Vector3d d(-99.9,99.9,9.99);
    std::ostringstream os;

    os << ToString::vector3dToString(b)
       << std::endl
       << ToString::vector3dToString(c)
       << std::endl
       << ToString::vector3dToString(d)
       << std::endl;

    ASSERT_EQ(os.str(),"   0.00000  -2.01238   7.13000\n   0.00000  -9.89870 -12.23800\n -99.90000  99.90000   9.99000\n");
}

