//
// Created by Michael Heuer on 29.10.17.
//

#include <gtest/gtest.h>
#include "ParticleCollection.h"

using namespace testing;

class AParticleCollectionTest : public Test {
public:
    void SetUp() override {}
};

TEST_F(AParticleCollectionTest, IndexOperator){

    Vector3d position0 = Vector3d(1,2,3);
    Vector3d position1 = Vector3d(4,5,6);
    Vector3d position2 = Vector3d(7,8,9);

    Particle particle0(position0);

    VectorXd positions(9);
    positions<< position0, position1, position2;
    ParticleCollection particleCollection(positions);

    ASSERT_EQ(particleCollection[0].position(),position0);
    ASSERT_EQ(particleCollection[1].position(),position1);
    ASSERT_EQ(particleCollection[2].position(),position2);
    ASSERT_EQ(particleCollection[-1].position(),position2);
    ASSERT_EQ(particleCollection[-2].position(),position1);
    ASSERT_EQ(particleCollection[-3].position(),position0);

    EXPECT_DEATH(particleCollection[3],"");
    EXPECT_DEATH(particleCollection[-4],"");
}

TEST_F(AParticleCollectionTest, extractSet){
    VectorXd positions;
    ASSERT_TRUE(false);
}