// Copyright (C) 2017-2019 Michael Heuer.
// SPDX-License-Identifier: GPL-3.0-or-later


#include <gmock/gmock.h>
#include <ParticlesVector.h>
#include <PositionsVectorTransformer.h>
#include <NaturalConstants.h>

using namespace testing;
using namespace Eigen;

class APositionsVectorTransformerTest : public Test {
public:
    ParticlesVector<Spin> ev;
    void SetUp() override {
        Particle<Spin > e1 = {Spin::alpha,{0, 0, 0}};
        Particle<Spin > e2 = {Spin::alpha,{1, 0, 0}};
        Particle<Spin > e3 = {Spin::beta ,{0, 1, 0}};
        Particle<Spin > e4 = {Spin::beta ,{0, 0, 1}};
        ev.append(e1);
        ev.append(e2);
        ev.append(e3);
        ev.append(e4);
    }
};

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

    ASSERT_TRUE(ev.positionsVector().asEigenVector().isApprox(expectedPositions));
}
