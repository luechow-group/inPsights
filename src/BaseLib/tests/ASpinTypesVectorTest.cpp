//
// Created by Michael Heuer on 29.10.17.
//

#include <gtest/gtest.h>
#include "SpinTypesVector.h"

using namespace testing;
using namespace Eigen;

class ASpinTypesVectorTest : public Test {
public:
    void SetUp() override {}
};

TEST_F(ASpinTypesVectorTest, Constructors){
    SpinTypesVector spinTypesVector1(3);
    VectorXi spinTypesNone = VectorXi::Constant(3,0);
    ASSERT_EQ(spinTypesVector1.spinTypesAsEigenVector(),spinTypesNone);

    VectorXi spinTypesAlpha(3);
    spinTypesAlpha << 1,1,1;
    SpinTypesVector spinTypesVector2(spinTypesAlpha);
    ASSERT_EQ(spinTypesVector2.spinTypesAsEigenVector(),spinTypesAlpha);
}

TEST_F(ASpinTypesVectorTest, CopyConstructor){
    VectorXi spinTypesAlpha(3);
    spinTypesAlpha << 1,1,1;
    SpinTypesVector spinTypesVector(spinTypesAlpha);
    SpinTypesVector copySpinTypesVector(spinTypesVector);
}

TEST_F(ASpinTypesVectorTest, IndexOperator){
    VectorXi spinTypes(3);
    spinTypes << 1,-1,0;
    SpinTypesVector spinTypesVector(spinTypes);

    ASSERT_EQ(spinTypesVector[0],Spin::SpinType::alpha);
    ASSERT_EQ(spinTypesVector[1],Spin::SpinType::beta);
    ASSERT_EQ(spinTypesVector[2],Spin::SpinType::none);
}

TEST_F(ASpinTypesVectorTest, IndexOutOfBoundsDeaths){
    VectorXi spinTypes1(3);
    spinTypes1 << 0,1,3;
    EXPECT_DEATH(SpinTypesVector spinTypesVector(spinTypes1),"");

    VectorXi spinTypes2(3);
    spinTypes2 << -1,1,3;
    EXPECT_DEATH(SpinTypesVector spinTypesVector(spinTypes2),"");
}
