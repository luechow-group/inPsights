//
// Created by Michael Heuer on 16.03.18.
//

#include <gtest/gtest.h>
#include <LebedevSphericalIntegration/TestFunctions.h>
#include "SphericalHarmonicsRadialBasisExpander.h"


using namespace testing;

class ASphericalHarmonicsRadialBasisExpanderTest : public ::testing::Test {};

TEST_F(ASphericalHarmonicsRadialBasisExpanderTest, CoefficientsVector) {

    UnitSphere f;

    unsigned nmax = 4;
    unsigned lmax = 2;
    double rCutoff = 1.0;

    SphericalHarmonicsRadialBasisExpander expander(nmax,lmax,rCutoff);

    auto coeffs = expander.coefficients(f);

    for (unsigned n = 1; n <= nmax; ++n) {
        for (unsigned l = 0; l < lmax; ++l) {
            std::cout << "\n";
            for (int m = -lmax; m < +lmax; ++m) {
                std::cout << coeffs[n-1][l][m+lmax] << "";
            }
            std::endl;
        }
    }


}