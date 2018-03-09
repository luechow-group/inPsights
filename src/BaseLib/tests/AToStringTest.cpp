//
// Created by Leonard Reuter on 09.03.18.
//

#include <gtest/gtest.h>
#include "ToString.h"

using namespace testing;

class AToStringTest : public Test {
public:
    void SetUp() override {
    }
};

TEST_F(AToStringTest, Double){
    double a = 7.845;
    std::ostringstream os;
    os << ToString::double2string(a,5,0);
    ASSERT_EQ(os.str()," 7.84500");
}

TEST_F(AToStringTest, Vector3d){
    Eigen::Vector3d b(0,-2.01238,7.13);
    Eigen::Vector3d c(1e-10,-9.8987,-12.238);
    Eigen::Vector3d d(-99.9,99.9,9.99);
    std::ostringstream os;

    os << ToString::vector3d2string(b)
       << std::endl
       << ToString::vector3d2string(c)
       << std::endl
       << ToString::vector3d2string(d)
       << std::endl;

    ASSERT_EQ(os.str(),"   0.00000  -2.01238   7.13000\n   0.00000  -9.89870 -12.23800\n -99.90000  99.90000   9.99000\n");
}

