//
// Created by dahl on 21.11.17.
//

#include <gtest/gtest.h>

using namespace testing;

class AExampleTest : public Test {};


TEST_F(AExampleTest, SimpleTest) {
    int a = 1;
    double b = 1.0;
    ASSERT_TRUE(a == b);
}
