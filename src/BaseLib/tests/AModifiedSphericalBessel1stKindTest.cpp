// Copyright (C) 2020 Michael Heuer.
// SPDX-License-Identifier: GPL-3.0-or-later

#include <gmock/gmock.h>
#include <SpecialMathFunctions/ModifiedSphericalBessel1stKind.h>
#include "NaturalConstants.h"

using namespace testing;


/* Modified Bessel Function of the First Kind
 * Tested with Mathematica 11.1.1.0 as a reference.
 * Reference values for x in {0.0, 1.0, pi} and l in {0,1,2,3} were produced with the following Wolfram code:
 * >    ModifiedSphericalBessel1stKind[n_, x_] := Re[I^(-n)*SphericalBesselJ[n, I*x]];
 * >    Table[N[ModifiedSphericalBessel1stKind[n, x], 16], {x, {0, 1, Pi}}, {n, 0, 3}]
 *  which is an implementation of eq (2) from https://mathworld.wolfram.com/ModifiedBesselFunctionoftheFirstKind.html
 */

class AModifiedSphericalBessel1stKindTest : public ::testing::Test {
public:
    double eps;
    void SetUp() override {
        eps = std::numeric_limits<double>::epsilon();
    }
};

TEST_F(AModifiedSphericalBessel1stKindTest, EvaulateAtX0) {
    unsigned maxDegree = 3;
    std::vector<double> reference = {1.0, 0.0, 0.0, 0.0};
    auto results = ModifiedSphericalBessel1stKind::evaluateToMaxDegree(maxDegree, 0, eps, eps);

    for (unsigned n = 0; n <= maxDegree; ++n)
        ASSERT_NEAR(results[n],reference[n], eps);
}

TEST_F(AModifiedSphericalBessel1stKindTest, EvaulateAtX1) {
    unsigned maxDegree = 3;
    std::vector<double> reference = {1.175201193643801, 0.3678794411714423, 0.07156287012947449, 0.01006509052406986};
    auto results = ModifiedSphericalBessel1stKind::evaluateToMaxDegree(maxDegree, 1.0, eps, eps);

    for (unsigned n = 0; n <= maxDegree; ++n)
        ASSERT_NEAR(results[n],reference[n], eps*10);
}

TEST_F(AModifiedSphericalBessel1stKindTest, EvaulateAtXPi) {
    unsigned maxDegree = 3;
    std::vector<double> reference = {3.676077910374978, 2.519701386524868, 1.269940325689366, 0.4985285838729271};
    auto results = ModifiedSphericalBessel1stKind::evaluateToMaxDegree(maxDegree, Constant::pi, eps, eps);

    for (unsigned n = 0; n <= maxDegree; ++n)
        ASSERT_NEAR(results[n],reference[n], eps*10);
}