//
// Created by Michael Heuer on 15.05.18.
//

#include "ExpansionSettings.h"
#include <NaturalConstants.h>
#include <gtest/gtest.h>

using namespace testing;


TEST(AExpansionSettingsTest, toString) {

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
            "sigma_atom\t\t\t: 0.944863 bohr\n"
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
            "Radius\t\t\t\t: 7.5589 bohr\n"
            "Width\t\t\t\t: 1.88973 bohr\n"
            "Center weight\t\t: 1\n\n"
            "Alchemical Similarities:\n"
            "------------------------\n"
            "a <-> b: 0.5\n\n";
    ASSERT_EQ(s,ref);
}