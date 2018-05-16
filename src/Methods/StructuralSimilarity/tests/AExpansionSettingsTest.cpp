//
// Created by Michael Heuer on 15.05.18.
//

#include <gtest/gtest.h>
#include "ExpansionSettings.h"

using namespace testing;

class AExpansionSettingsTest : public ::testing::Test {
public:
};

TEST_F(AExpansionSettingsTest , uninitialized) {
    ASSERT_EQ(ExpansionSettings::Radial::nmax,0);
    ASSERT_EQ(ExpansionSettings::Radial::basisType,RadialGaussianBasisType(0));
    ASSERT_EQ(ExpansionSettings::Radial::sigmaAtom,0);
    ASSERT_EQ(ExpansionSettings::Radial::basisType,RadialGaussianBasisType::equispaced);

    ASSERT_EQ(ExpansionSettings::Angular::lmax,0);

    ASSERT_EQ(ExpansionSettings::Cutoff::cutoffRadius,0);


    ExpansionSettings::Radial::nmax = 4;
    ExpansionSettings::Angular::lmax = 4;

    ASSERT_EQ(ExpansionSettings::Radial::nmax,4);
    ASSERT_EQ(ExpansionSettings::Angular::lmax,4);
}

TEST_F(AExpansionSettingsTest, defaults) {
    ExpansionSettings::defaults();

    ASSERT_EQ(ExpansionSettings::Radial::nmax,2);
    ASSERT_EQ(ExpansionSettings::Radial::basisType,RadialGaussianBasisType::equispaced);
    ASSERT_EQ(ExpansionSettings::Radial::sigmaAtom,0.5);

    ASSERT_EQ(ExpansionSettings::Angular::lmax,1);

    ASSERT_EQ(ExpansionSettings::Cutoff::cutoffRadius,4.0);

    ExpansionSettings::Radial::nmax = 4;
    ExpansionSettings::Angular::lmax = 4;

    ASSERT_EQ(ExpansionSettings::Radial::nmax,4);
    ASSERT_EQ(ExpansionSettings::Angular::lmax,4);
}
