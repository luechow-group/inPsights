//
// Created by Michael Heuer on 29.10.17.
//

#include <gtest/gtest.h>
#include "SpinType.h"

using namespace testing;
using namespace Spins;

class ASpinTypeTest : public Test {
public:
    void SetUp() override {
    }
};

TEST_F(ASpinTypeTest , IntegerValues){
    ASSERT_EQ(int(SpinType::alpha),-2);
    ASSERT_EQ(int(SpinType::beta),-1);
    ASSERT_EQ(int(SpinType::none),0);
}

TEST_F(ASpinTypeTest , QuantumNumbers){
    ASSERT_EQ(quantumNumber(),1/2.0);
    ASSERT_EQ(magneticQuantumNumber(SpinType::alpha),1/2.0);
    ASSERT_EQ(magneticQuantumNumber(SpinType::beta),-1/2.0);
    EXPECT_DEATH(magneticQuantumNumber(SpinType::none),"");
}
