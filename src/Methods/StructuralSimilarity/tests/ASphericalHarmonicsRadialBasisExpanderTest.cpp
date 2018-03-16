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

    int nmax = 4;
    int lmax = 2;
    double rCutoff = 1.0;

    SphericalHarmonicsRadialBasisExpander expander(nmax,lmax,rCutoff);

    auto coeffs = expander.coefficients(f);

    for (int n = 1; n <= nmax; ++n) {
        std::cout << n << " ";
        for (int l = 0; l <= lmax; ++l) {
            std::cout << "  " << l << "\n";
            for (int m = -l; m <= +l; ++m) {
                std::cout << "   " << m << " m+l=" << m+l;
                std::cout << coeffs[n-1][l][m+l] << std::endl;
            }
        }
        std::cout << std::endl;
    }


}