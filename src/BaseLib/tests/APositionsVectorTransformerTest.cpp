//
// Created by Michael Heuer on 29.10.17.
//


#include <gtest/gtest.h>
#include <ParticlesVector.h>
#include <PositionsVectorTransformer.h>
#include <NaturalConstants.h>

using namespace testing;
using namespace Eigen;

class APositionsVectorTransformerTest : public Test {
public:
    ParticlesVector<Spins::SpinType> ev;
    void SetUp() override {
        Particle<Spins::SpinType > e1 = {Spins::SpinType::alpha,{0, 0, 0}};
        Particle<Spins::SpinType > e2 = {Spins::SpinType::alpha,{1, 0, 0}};
        Particle<Spins::SpinType > e3 = {Spins::SpinType::beta ,{0, 1, 0}};
        Particle<Spins::SpinType > e4 = {Spins::SpinType::beta ,{0, 0, 1}};
        ev.append(e1);
        ev.append(e2);
        ev.append(e3);
        ev.append(e4);
    }
};

TEST_F(APositionsVectorTransformerTest, counterclockwiseRotation){

    double angle = 120.*ConversionFactors::deg2rad; // 120° counterclockwise rotiation
    Eigen::Vector3d axis = {1,1,1};

    auto rotmat =PositionsVectorTransformer::rotationMatrixFromQuaternion(
            PositionsVectorTransformer::quaternionFromAngleAndAxis(angle,axis));

    Eigen::Matrix3d expectedRotationMatrix;
    expectedRotationMatrix << 0,0,1, 1,0,0, 0,1,0;

    ASSERT_TRUE(rotmat.isApprox(expectedRotationMatrix));

    PositionsVectorTransformer::rotateAroundAxis(ev.positionsVector(), angle, axis);

    Eigen::VectorXd expectedPositions(12);
    expectedPositions << 0,0,0, 0,0,1, 1,0,0, 0,1,0;

    ASSERT_TRUE(ev.positionsVector().positionsAsEigenVector().isApprox(expectedPositions));
}

TEST_F(APositionsVectorTransformerTest, clockwiseRotation){

    double angle = -120.*ConversionFactors::deg2rad; // 120° clockwise rotiation
    Eigen::Vector3d axis = {1,1,1};

    auto rotmat =PositionsVectorTransformer::rotationMatrixFromQuaternion(
            PositionsVectorTransformer::quaternionFromAngleAndAxis(angle,axis));

    Eigen::Matrix3d expectedRotationMatrix;
    expectedRotationMatrix << 0,1,0, 0,0,1, 1,0,0;

    ASSERT_TRUE(rotmat.isApprox(expectedRotationMatrix));

    PositionsVectorTransformer::rotateAroundAxis(ev.positionsVector(), angle, axis);

    Eigen::VectorXd expectedPositions(12);
    expectedPositions << 0,0,0, 0,1,0, 0,0,1, 1,0,0;

    ASSERT_TRUE(ev.positionsVector().positionsAsEigenVector().isApprox(expectedPositions));
}

TEST_F(APositionsVectorTransformerTest, centerOfMassTranslation){

    auto centerOfMass = PositionsVectorTransformer::calculateCenterOfMass(ev.positionsVector());
    Eigen::Vector3d expectedCenterOfMass = {1./4.,1./4.,1./4.};
    ASSERT_TRUE(centerOfMass.isApprox(expectedCenterOfMass));

    PositionsVectorTransformer::translateCenterOfMassToOrigin(ev.positionsVector());

    Eigen::VectorXd expectedPositions(12);
    expectedPositions << \
    -1./4.,-1./4.,-1./4.,\
    3./4.,-1./4.,-1./4.,\
    -1./4.,3./4.,-1./4.,\
    -1./4.,-1./4.,3./4.;

    ASSERT_TRUE(ev.positionsVector().positionsAsEigenVector().isApprox(expectedPositions));

}
