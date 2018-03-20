//
// Created by Michael Heuer on 20.03.18.
//

#include <gtest/gtest.h>
#include <BoostSphericalHarmonics.h>
#include <sh/spherical_harmonics.h>

using namespace testing;

class ABoostSphericalHarmonicsTest : public ::testing::Test {};

TEST_F(ABoostSphericalHarmonicsTest, RealSphericalHarmonicY) {

    double reference = 3.0*std::sqrt(3.0/M_PI)/8.0; // from Mathematica

    double theta = M_PI/3.0;
    double phi = theta;
    double boostCalculated = BoostSphericalHarmonics::realSphericalHarmonicY(1,-1,theta,phi);
    double googleCalculated = sh::EvalSH(1,-1,theta,phi);

    ASSERT_NEAR(boostCalculated,reference,std::numeric_limits<double>::epsilon());
    ASSERT_NEAR(-googleCalculated,reference,1e-7); // the google library is imprecise and has a sign change?!
}
