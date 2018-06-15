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
    ASSERT_EQ(ExpansionSettings::Radial::nmax,5);
    ASSERT_EQ(ExpansionSettings::Radial::basisType, ExpansionSettings::Radial::BasisType::equispaced);
    ASSERT_EQ(ExpansionSettings::Radial::sigmaAtom,0.5);

    ASSERT_EQ(ExpansionSettings::Radial::integrationSteps,100);
    ASSERT_EQ(ExpansionSettings::Radial::desiredAbsoluteError,0.0);
    ASSERT_EQ(ExpansionSettings::Radial::desiredRelativeError,1e-6);


    ASSERT_EQ(ExpansionSettings::Angular::lmax,5);

    ASSERT_EQ(ExpansionSettings::Cutoff::radius,8.0);
    ASSERT_EQ(ExpansionSettings::Cutoff::width,2.0);
    ASSERT_EQ(ExpansionSettings::Cutoff::centerWeight,1.0);


    ExpansionSettings::Radial::nmax = 4;
    ExpansionSettings::Angular::lmax = 4;

    ASSERT_EQ(ExpansionSettings::Radial::nmax,4);
    ASSERT_EQ(ExpansionSettings::Angular::lmax,4);
}

TEST_F(AExpansionSettingsTest, defaults) {
    ExpansionSettings::defaults();

    ASSERT_EQ(ExpansionSettings::Radial::nmax,5);
    ASSERT_EQ(ExpansionSettings::Radial::basisType,ExpansionSettings::Radial::BasisType::equispaced);
    ASSERT_EQ(ExpansionSettings::Radial::sigmaAtom,1.0);

    ASSERT_EQ(ExpansionSettings::Angular::lmax,5);

    ASSERT_EQ(ExpansionSettings::Cutoff::radius,8.0);
    ASSERT_EQ(ExpansionSettings::Cutoff::width,2.0);
    ASSERT_EQ(ExpansionSettings::Cutoff::centerWeight,1.0);

    ExpansionSettings::Radial::nmax = 4;
    ExpansionSettings::Angular::lmax = 4;

    ASSERT_EQ(ExpansionSettings::Radial::nmax,4);
    ASSERT_EQ(ExpansionSettings::Angular::lmax,4);
}

TEST_F(AExpansionSettingsTest, toString) {
    ExpansionSettings::defaults();

    auto s = ExpansionSettings::toString();
    std::string ref =
            "General:\n"
            "--------\n"
            "Expansion mode\t\t: chemical\n"
            "Sharpness zeta\t\t: 2\n"
            "Regularization gamma: 0.1\n"
            "\n"
            "Radial:\n"
            "-------\n"
            "BasisType\t\t\t: equispaced\n"
            "n_max\t\t\t\t: 5\n"
            "sigma_atom\t\t\t: 1 angstrom\n"
            "Integration steps\t: 100\n"
            "Desired abs. err\t: 0\n"
            "Desired rel. err\t: 2.22045e-14\n"
            "\n"
            "Angular:\n"
            "--------\n"
            "l_max\t\t\t\t: 5\n"
            "\n"
            "Cutoff:\n"
            "-------\n"
            "Radius\t\t\t\t: 8 angstrom\n"
            "Width\t\t\t\t: 2 angstrom\n"
            "Center weight\t\t: 1\n\n"
            "Alchemical Similarities:\n"
            "------------------------\n"
            "a <-> b: 0.5\n\n";
    ASSERT_EQ(s,ref);
}