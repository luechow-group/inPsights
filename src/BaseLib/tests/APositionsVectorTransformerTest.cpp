//
// Created by Michael Heuer on 29.10.17.
//

#include <gtest/gtest.h>
#include <ElectronsVector.h>
#include <PositionsVectorTransformer.h>
#include <limits>

using namespace testing;
using namespace Eigen;

class APositionsVectorTransformerTest : public Test {
public:
    ElectronsVector ev;
    void SetUp() override {
        Electron e1 = {{0,0,0},Spin::SpinType::alpha};
        Electron e2 = {{1,0,0},Spin::SpinType::alpha};
        Electron e3 = {{0,1,0},Spin::SpinType::beta};
        Electron e4 = {{0,0,1},Spin::SpinType::beta};
        ev.append(e1);
        ev.append(e2);
        ev.append(e3);
        ev.append(e4);
    }
};



TEST_F(APositionsVectorTransformerTest, counterclockwiseRotation){
    ElectronsVector electronsVector = ev;

    double angle = 2.*M_PI/3.; // 120° counterclockwise rotiation
    Eigen::Vector3d axis = {1,1,1};

    auto rotmat =PositionsVectorTransformer::rotationMatrixFromQuaternion(
            PositionsVectorTransformer::quaternionFromAngleAndAxis(angle,axis));

    Eigen::Matrix3d expectedRotationMatrix;
    expectedRotationMatrix << 0,0,1, 1,0,0, 0,1,0;

    ASSERT_TRUE(rotmat.isApprox(expectedRotationMatrix));

    PositionsVectorTransformer::rotateAroundAxis(electronsVector.positionsVector(), angle, axis);

    Eigen::VectorXd expectedPositions(12);
    expectedPositions << 0,0,0, 0,0,1, 1,0,0, 0,1,0;

    ASSERT_TRUE(electronsVector.positionsVector().positionsAsEigenVector().isApprox(expectedPositions));
}

TEST_F(APositionsVectorTransformerTest, clockwiseRotation){
    ElectronsVector electronsVector = ev;

    double angle = -2.*M_PI/3.; // 120° clockwise rotiation
    Eigen::Vector3d axis = {1,1,1};

    auto rotmat =PositionsVectorTransformer::rotationMatrixFromQuaternion(
            PositionsVectorTransformer::quaternionFromAngleAndAxis(angle,axis));

    Eigen::Matrix3d expectedRotationMatrix;
    expectedRotationMatrix << 0,1,0, 0,0,1, 1,0,0;

    ASSERT_TRUE(rotmat.isApprox(expectedRotationMatrix));

    PositionsVectorTransformer::rotateAroundAxis(electronsVector.positionsVector(), angle, axis);

    Eigen::VectorXd expectedPositions(12);
    expectedPositions << 0,0,0, 0,1,0, 0,0,1, 1,0,0;

    ASSERT_TRUE(electronsVector.positionsVector().positionsAsEigenVector().isApprox(expectedPositions));
}

TEST_F(APositionsVectorTransformerTest, centerOfMassTranslation){
    ElectronsVector electronsVector = ev;

    auto centerOfMass = PositionsVectorTransformer::calculateCenterOfMass(electronsVector.positionsVector());
    Eigen::Vector3d expectedCenterOfMass = {1./4.,1./4.,1./4.};
    ASSERT_TRUE(centerOfMass.isApprox(expectedCenterOfMass));

    PositionsVectorTransformer::translateCenterOfMassToOrigin(electronsVector.positionsVector());

    Eigen::VectorXd expectedPositions(12);
    expectedPositions << \
    -1./4.,-1./4.,-1./4.,\
    3./4.,-1./4.,-1./4.,\
    -1./4.,3./4.,-1./4.,\
    -1./4.,-1./4.,3./4.;

    ASSERT_TRUE(electronsVector.positionsVector().positionsAsEigenVector().isApprox(expectedPositions));

}