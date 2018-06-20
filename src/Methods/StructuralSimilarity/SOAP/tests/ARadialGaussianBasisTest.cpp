//
// Created by Michael Heuer on 16.05.18.
//


#include <gtest/gtest.h>
#include <RadialGaussianBasis.h>
#include "ExpansionSettings.h"

using namespace testing;

class ARadialGaussianBasisTest : public ::testing::Test {};


TEST_F(ARadialGaussianBasisTest, EquispacedCenters) {
    //ExpansionSettings::defaults();
    RadialGaussianBasis();
}
