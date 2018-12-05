//
// Created by Michael Heuer on 20.03.18.
//

#include <gmock/gmock.h>
#include <SpecialMathFunctions/BoostSphericalHarmonics.h>

using namespace testing;

class ABoostSphericalHarmonicsTest : public ::testing::Test {
public:
    double absError, maxR;
    unsigned NpointsTheta, NpointsPhi, NpointsR;
    int logScaleMin;
    int logScaleMax;

    void SetUp() override {
        absError = 1e2*std::numeric_limits<double>::epsilon();
        maxR = 3.0;
        NpointsTheta = 181;
        NpointsPhi = 361;
        NpointsR = 100;
        logScaleMin = -15;
        logScaleMax = 15;
    }
};

TEST_F(ABoostSphericalHarmonicsTest, RealSphericalHarmonicY) {

    double reference = 3.0*std::sqrt(3.0/M_PI)/8.0; // from Mathematica

    double theta = M_PI/3.0;
    double phi = theta;
    double boostCalculated = BoostSphericalHarmonics::realSphericalHarmonicY(1,-1,theta,phi);

    ASSERT_NEAR(boostCalculated,reference,std::numeric_limits<double>::epsilon());
}

TEST_F(ABoostSphericalHarmonicsTest, AngleBoundaries) {

    for (int k = logScaleMin; k <= logScaleMax; ++k) { //ignored r = 0 case here
        double r = pow(1,k);

        for (unsigned i = 1; i <= NpointsTheta; ++i) { //ignored theta = 0 case here
            double theta = double(i) / double(NpointsTheta) * M_PI;

            for (unsigned j = 0; j < NpointsPhi; ++j) {
                double phi = double(j) / double(NpointsPhi - 1) * 2 * M_PI;

                Eigen::Vector3d vec = {r * sin(theta) * cos(phi),
                                       r * sin(theta) * sin(phi),
                                       r * cos(theta)};

                double rCalc, thetaCalc, phiCalc;
                BoostSphericalHarmonics::toSphericalCoordsStandardizedWith2PiShift(vec, rCalc, thetaCalc, phiCalc);

                ASSERT_NEAR(rCalc, r, absError);
                ASSERT_NEAR(thetaCalc, theta, absError);
                ASSERT_NEAR(phiCalc, phi, absError);
            }
        }
    }
}

TEST_F(ABoostSphericalHarmonicsTest, ZAxisCase) {

    double theta = 0.0;

    for (int k = logScaleMin; k <= logScaleMax; ++k) { //ignored r = 0 case here
        double r = pow(1,k);

        for (unsigned j = 0; j < NpointsPhi; ++j) {
            double phi = double(j) / double(NpointsPhi - 1) * 2 * M_PI;

            Eigen::Vector3d vec = {r * sin(theta) * cos(phi),
                                   r * sin(theta) * sin(phi),
                                   r * cos(theta)};

            double rCalc, thetaCalc, phiCalc;
            BoostSphericalHarmonics::toSphericalCoordsStandardizedWith2PiShift(vec, rCalc, thetaCalc, phiCalc);
            //std::cout << phi/double(2 * M_PI)*360 << ",  "  << phiCalc/double(2 * M_PI)*360 << std::endl;
            ASSERT_NEAR(rCalc, r, absError);
            ASSERT_NEAR(thetaCalc, theta, absError);
            ASSERT_NEAR(phiCalc, 0.0, absError);
        }
    }
}

TEST_F(ABoostSphericalHarmonicsTest, OriginCase) {

    double r = 0.0;

    for (unsigned i = 0; i < NpointsTheta; ++i) { //ignored theta = 0 case here
        double theta = double(i) / double(NpointsTheta-1) * M_PI;

        for (unsigned j = 0; j < NpointsPhi; ++j) {
            double phi = double(j) / double(NpointsPhi - 1) * 2 * M_PI;

            Eigen::Vector3d vec = {r * sin(theta) * cos(phi),
                                   r * sin(theta) * sin(phi),
                                   r * cos(theta)};

            double rCalc, thetaCalc, phiCalc;
            BoostSphericalHarmonics::toSphericalCoordsStandardizedWith2PiShift(vec, rCalc, thetaCalc, phiCalc);

            ASSERT_NEAR(rCalc, r, absError);
            ASSERT_NEAR(thetaCalc, 0.0, absError);
            ASSERT_NEAR(phiCalc, 0.0, absError);
        }
    }
}
