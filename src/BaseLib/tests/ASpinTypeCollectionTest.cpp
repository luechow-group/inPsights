//
// Created by Michael Heuer on 29.10.17.
//

#include <gtest/gtest.h>
#include "SpinTypeCollection.h"

using namespace testing;

class ASpinTypeCollectionTest : public Test {
public:
    void SetUp() override {}
};

TEST_F(ASpinTypeCollectionTest, Constructors){
    SpinTypeCollection spinTypeCollection1(3);
    VectorXi spinTypesNone = VectorXi::Constant(3,2);
    ASSERT_EQ(spinTypeCollection1.asVectorXi(),spinTypesNone);

    VectorXi spinTypesAlpha(3);
    spinTypesAlpha << 0,0,0;
    SpinTypeCollection spinTypeCollection2(spinTypesAlpha);
    ASSERT_EQ(spinTypeCollection2.asVectorXi(),spinTypesAlpha);
}

TEST_F(ASpinTypeCollectionTest, CopyConstructor){
    VectorXi spinTypesAlpha(3);
    spinTypesAlpha << 0,0,0;
    SpinTypeCollection spinTypeCollection(spinTypesAlpha);
    SpinTypeCollection copySpinTypeCollection(spinTypeCollection);
}

TEST_F(ASpinTypeCollectionTest, IndexOperator){
    VectorXi spinTypes(3);
    spinTypes << 0,1,2;
    SpinTypeCollection spinTypeCollection(spinTypes);

    ASSERT_EQ(spinTypeCollection.spinType(0),Spin::SpinType::alpha);
    ASSERT_EQ(spinTypeCollection.spinType(1),Spin::SpinType::beta);
    ASSERT_EQ(spinTypeCollection.spinType(2),Spin::SpinType::none);
}

TEST_F(ASpinTypeCollectionTest, IndexOutOfBoundsDeaths){
    VectorXi spinTypes1(3);
    spinTypes1 << 0,1,3;
    EXPECT_DEATH(SpinTypeCollection spinTypeCollection(spinTypes1),"");

    VectorXi spinTypes2(3);
    spinTypes2 << -1,1,3;
    EXPECT_DEATH(SpinTypeCollection spinTypeCollection(spinTypes2),"");
}
