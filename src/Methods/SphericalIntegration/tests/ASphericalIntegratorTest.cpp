//
// Created by Michael Heuer on 12.03.18.
//

#include <gtest/gtest.h>
#include "SphericalIntegrator.h"
#include "LebedevSphericalIntegration/TestFunctions.h"

class ASphericalIntegratorTest : public ::testing::Test {};

TEST_F(ASphericalIntegratorTest, UnitSphereVolume) {

    UnitSphere f;
    SphericalIntegrator integrator(SphericalIntegratorSettings::defaults());

    double calculated = integrator.integrate(f,0,1);
    double reference = 4/3.0*M_PI;

    ASSERT_NEAR(calculated, reference, 1e-15);
}

TEST_F(ASphericalIntegratorTest, GaussVolume) {

    Gauss3d f;
    SphericalIntegrator sphericalIntegrator(SphericalIntegratorSettings::defaults());

    double calculated = sphericalIntegrator.integrate(f,0,10);
    double reference = std::pow(M_PI,3.0/2.0);

    ASSERT_NEAR(calculated, reference, 1e-15);
}
