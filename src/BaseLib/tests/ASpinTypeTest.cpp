//
// Created by Michael Heuer on 29.10.17.
//

#include <gtest/gtest.h>
#include "SpinType.h"

using namespace testing;


class ASpinTypeTest : public Test {
public:
    void SetUp() override {
    }
};

TEST_F(ASpinTypeTest , IntegerValues){
    ASSERT_EQ(int(Spin::SpinType::alpha),0);
    ASSERT_EQ(int(Spin::SpinType::beta),1);
    ASSERT_EQ(int(Spin::SpinType::none),2);
}


